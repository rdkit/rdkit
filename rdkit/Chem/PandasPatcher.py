import importlib
import logging
import re
from io import StringIO
from xml.dom import minidom
from xml.parsers.expat import ExpatError

from rdkit.Chem import Mol

log = logging.getLogger(__name__)
RDK_MOLS_AS_IMAGE_ATTR = "__rdkitMolAsImage"
InteractiveRenderer = None
PrintAsImageString = None
molJustify = None
pandas_frame = None

pandas_formats = None
for pandas_formats_name in ("pandas.io", "pandas"):
  try:
    pandas_formats = importlib.import_module(f"{pandas_formats_name}.formats")
  except ModuleNotFoundError:
    continue
  break
if pandas_formats is None:
  log.warning("Failed to import pandas formats module")
  raise ModuleNotFoundError
to_html_class = None
for to_html_class_name in ("DataFrameRenderer", "DataFrameFormatter"):
  if hasattr(pandas_formats, "format") and hasattr(pandas_formats.format, to_html_class_name):
    to_html_class = getattr(pandas_formats.format, to_html_class_name)
    if hasattr(to_html_class, "to_html"):
      break
    else:
      to_html_class = None
if to_html_class is None:
  log.warning("Failed to find the pandas to_html method to patch")
  raise AttributeError
dataframeformatter_class = None
orig_get_formatter = None
if (hasattr(pandas_formats.format, "DataFrameFormatter")
    and hasattr(pandas_formats.format.DataFrameFormatter, "_get_formatter")):
  dataframeformatter_class = pandas_formats.format.DataFrameFormatter
  orig_get_formatter = getattr(dataframeformatter_class, "_get_formatter")
if orig_get_formatter is None:
  log.warning("Failed to find the pandas _get_formatter() function to patch")
orig_get_adjustment = None
for get_adjustment_module_name in ("format", "printing"):
  if hasattr(pandas_formats, get_adjustment_module_name):
    get_adjustment_module = getattr(pandas_formats, get_adjustment_module_name)
    for get_adjustment_name in ("get_adjustment", "_get_adjustment"):
      if hasattr(get_adjustment_module, get_adjustment_name):
        orig_get_adjustment = getattr(get_adjustment_module, get_adjustment_name)
        break
    if orig_get_adjustment is not None:
      break
if orig_get_adjustment is None:
  log.warning("Failed to find the pandas get_adjustment() function to patch")
  raise AttributeError
adj = orig_get_adjustment()
if not hasattr(adj, "len"):
  log.warning(f"Failed to find the pandas {adj.__class.__name__}.len() method to patch")
  raise AttributeError
html_formatter_module = None
html_formatter_class = None
for html_formatter_module_name in ("format", "html"):
  try:
    html_formatter_module = importlib.import_module(
      f"{pandas_formats.__name__}.{html_formatter_module_name}")
  except ModuleNotFoundError:
    continue
  if hasattr(html_formatter_module, "HTMLFormatter"):
    html_formatter_class = getattr(html_formatter_module, "HTMLFormatter")
    break
if html_formatter_class is None:
  log.warning("Failed to find the pandas HTMLFormatter class to patch")
  raise AttributeError
orig_write_cell = None
if not hasattr(html_formatter_class, "_write_cell"):
  log.warning("Failed to find the HTMLFormatter._write_cell() method to patch")
  raise AttributeError
orig_write_cell = getattr(html_formatter_class, "_write_cell")
if not (hasattr(pandas_formats, "printing") and hasattr(pandas_formats.printing, "pprint_thing")):
  log.warning("Failed to find the pprint_thing function")
  raise AttributeError
try:
  import pandas as pd
except ImportError:
  log.warning("Failed to import pandas")
  raise

dataframe_applymap = pd.DataFrame.applymap
try:
  if tuple(map(int, (pd.__version__.split(".")))) >= (2, 1, 0):
    dataframe_applymap = pd.DataFrame.map
except:
  pass

orig_to_html = getattr(to_html_class, "to_html")
pprint_thing = pandas_formats.printing.pprint_thing


def is_molecule_image(s):
  result = False
  try:
    # is text valid XML / HTML?
    xml = minidom.parseString(s)
    root_node = xml.firstChild
    # check data-content attribute
    if (root_node.nodeName in ['svg', 'img', 'div']
        and 'data-content' in root_node.attributes.keys()
        and root_node.attributes['data-content'].value == 'rdkit/molecule'):
      result = True
  except ExpatError:
    pass  # parsing xml failed and text is not a molecule image
  return result


styleRegex = re.compile("^(.*style=[\"'][^\"^']*)([\"'].*)$")


class MolFormatter:
  """Format molecules as images"""

  def __init__(self, orig_formatter=None):
    """Store original formatters (if any)"""
    self.orig_formatter = orig_formatter

  @staticmethod
  def default_formatter(x):
    """Default formatter function"""
    return pprint_thing(x, escape_chars=("\t", "\r", "\n"))

  @staticmethod
  def is_mol(x):
    """Return True if x is a Chem.Mol"""
    return isinstance(x, Mol)

  @classmethod
  def get_formatters(cls, df, orig_formatters):
    """Return an instance of MolFormatter for each column that contains Chem.Mol objects"""
    df_subset = df.select_dtypes("object")
    return {
      col: cls(orig_formatters.get(col, None))
      for col in df_subset.columns[dataframe_applymap(df_subset, MolFormatter.is_mol).any()]
    }

  def __call__(self, x):
    """Return x formatted based on its type"""
    if self.is_mol(x):
      return PrintAsImageString(x)
    if callable(self.orig_formatter):
      return self.orig_formatter(x)
    return self.default_formatter(x)


def check_rdk_attr(frame, attr):
  return hasattr(frame, attr) and getattr(frame, attr)


def set_rdk_attr(frame, attr):
  setattr(frame, attr, True)


def patched_to_html(self, *args, **kwargs):
  """A patched version of the to_html method
     that allows rendering molecule images in data frames.
  """
  frame = None
  if self.__class__.__name__ == "DataFrameRenderer":
    fmt = self.fmt
  elif self.__class__.__name__ == "DataFrameFormatter":
    fmt = self
  else:
    raise ValueError(f"patched_to_html: unexpected class {self.__class__.__name__}")
  frame = fmt.frame
  if not check_rdk_attr(frame, RDK_MOLS_AS_IMAGE_ATTR):
    return orig_to_html(self, *args, **kwargs)
  orig_formatters = fmt.formatters
  try:
    formatters = orig_formatters or {}
    if not isinstance(formatters, dict):
      formatters = {col: formatters[i] for i, col in enumerate(self.columns)}
    else:
      formatters = dict(formatters)
    formatters.update(MolFormatter.get_formatters(frame, formatters))
    fmt.formatters = formatters
    res = orig_to_html(self, *args, **kwargs)
    # in pandas 0.25 DataFrameFormatter.to_html() returns None
    if (res is None and not hasattr(html_formatter_class, "get_result") and hasattr(self, "buf")
        and hasattr(self.buf, "getvalue")):
      res = self.buf.getvalue()
    should_inject = res and InteractiveRenderer and InteractiveRenderer.isEnabled()
    if should_inject:
      res = InteractiveRenderer.injectHTMLFooterAfterTable(res)
      # in pandas 0.25 we need to make sure to update buf as return value will be ignored
      if hasattr(self, "buf") and isinstance(self.buf, StringIO):
        self.buf.seek(0)
        self.buf.write(res)
    return res
  finally:
    fmt.formatters = orig_formatters


def patched_get_formatter(self, i, *args, **kwargs):
  if (isinstance(self.formatters, dict) and isinstance(i, int) and i not in self.columns
      and hasattr(self, "tr_col_num") and i >= self.tr_col_num):
    max_cols = 0
    if hasattr(self, "max_cols_fitted"):
      max_cols = self.max_cols_fitted
    elif hasattr(self, "max_cols_adj"):
      max_cols = self.max_cols_adj
    n_trunc_cols = len(self.columns) - max_cols
    if n_trunc_cols > 0:
      i += n_trunc_cols
  return orig_get_formatter(self, i, *args, **kwargs)


def patched_write_cell(self, s, *args, **kwargs):
  """ Disable escaping of HTML in order to render img / svg tags """
  styleTags = f"text-align: {molJustify};"
  def_escape = self.escape
  try:
    if hasattr(self.frame, RDK_MOLS_AS_IMAGE_ATTR) and is_molecule_image(s):
      self.escape = False
      kind = kwargs.get('kind', None)
      if kind == 'td':
        tags = kwargs.get('tags', None) or ''
        match = styleRegex.match(tags)
        if match:
          tags = styleRegex.sub(f'\\1 {styleTags}\\2', tags)
        else:
          if tags:
            tags += ' '
          tags += f'style="{styleTags}"'
        kwargs['tags'] = tags
    return orig_write_cell(self, s, *args, **kwargs)
  finally:
    self.escape = def_escape


def patched_get_adjustment():
  """ Avoid truncation of data frame values that contain HTML content """
  adj = orig_get_adjustment()
  orig_len = adj.len
  adj.len = lambda text, *args, **kwargs: (0 if is_molecule_image(text) else orig_len(
    text, *args, **kwargs))
  return adj


def renderImagesInAllDataFrames(images=True):
  if images:
    set_rdk_attr(pd.core.frame.DataFrame, RDK_MOLS_AS_IMAGE_ATTR)
  elif hasattr(pd.core.frame.DataFrame, RDK_MOLS_AS_IMAGE_ATTR):
    delattr(pd.core.frame.DataFrame, RDK_MOLS_AS_IMAGE_ATTR)


def changeMoleculeRendering(frame, renderer='image'):
  if not renderer.lower().startswith('str'):
    set_rdk_attr(frame, RDK_MOLS_AS_IMAGE_ATTR)
  elif hasattr(frame, RDK_MOLS_AS_IMAGE_ATTR):
    delattr(frame, RDK_MOLS_AS_IMAGE_ATTR)


def patchPandas():
  if getattr(to_html_class, "to_html") != patched_to_html:
    setattr(to_html_class, "to_html", patched_to_html)
  if getattr(html_formatter_class, "_write_cell") != patched_write_cell:
    setattr(html_formatter_class, "_write_cell", patched_write_cell)
  if getattr(get_adjustment_module, get_adjustment_name) != patched_get_adjustment:
    setattr(get_adjustment_module, get_adjustment_name, patched_get_adjustment)
  if (orig_get_formatter
      and getattr(dataframeformatter_class, "_get_formatter") != patched_get_formatter):
    setattr(dataframeformatter_class, "_get_formatter", patched_get_formatter)


def unpatchPandas():
  if orig_to_html:
    setattr(to_html_class, "to_html", orig_to_html)
  if orig_write_cell:
    setattr(html_formatter_class, "_write_cell", orig_write_cell)
  if orig_get_adjustment:
    setattr(pandas_formats.format, get_adjustment_name, orig_get_adjustment)
  if orig_get_formatter:
    setattr(dataframeformatter_class, "_get_formatter", orig_get_formatter)
