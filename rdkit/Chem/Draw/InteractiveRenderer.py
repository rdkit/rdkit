#
#  Copyright (C) 2022 Novartis Institute of BioMedical Research
#  and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Interactive molecule rendering through rdkit-structure-renderer.js """
from xml.dom import minidom
import uuid
import json
import re
from . import rdMolDraw2D
from IPython.display import HTML
from rdkit import Chem

rdkitStructureRendererJsUrl = "https://unpkg.com/rdkit-structure-renderer/dist/rdkit-structure-renderer-module.js"
minimalLibJsUrl = "https://unpkg.com/rdkit-structure-renderer/public/RDKit_minimal.js"
parentNodeQuery = "div[class*=jp-NotebookPanel-notebook]"
_enabled = False
_camelCaseOptToTagRe = re.compile("[A-Z]")
_opts = "__rnrOpts"
_disabled = "_disabled"
_defaultDrawOptions = None
_defaultDrawOptionsDict = None

_enabled = False


def _isAcceptedKeyValue(key, value):
  return not key.startswith("__") and (type(value) in (bool, int, float, str, tuple, list) or
                                       (callable(value) and key.startswith("get")))


def _getValueFromKey(d, key):
  value = getattr(d, key)
  return value() if callable(value) else value


def _molDrawOptionsToDict(molDrawOptions=rdMolDraw2D.MolDrawOptions()):
  return {
    key: _getValueFromKey(molDrawOptions, key)
    for key in dir(molDrawOptions) if _isAcceptedKeyValue(key, getattr(molDrawOptions, key))
  }


def _isDrawOptionEqual(v1, v2):
  if type(v1) != type(v2):
    return False
  if type(v1) in (tuple, list):
    if len(v1) != len(v2):
      return False
    return all(_isDrawOptionEqual(item1, v2[i]) for i, item1 in enumerate(v1))
  if type(v1) == float:
    return abs(v1 - v2) < 1.e-5
  return v1 == v2


def filterDefaultDrawOpts(molDrawOptions):
  global _defaultDrawOptionsDict
  if not isinstance(molDrawOptions, rdMolDraw2D.MolDrawOptions):
    raise ValueError(f"Bad args for {__name__}.filterDefaultDrawOpts("
                     "molDrawOptions: Chem.Draw.rdMolDraw2D.MolDrawOptions)")
  molDrawOptionsAsDict = _molDrawOptionsToDict(molDrawOptions)
  if (_defaultDrawOptionsDict is None):
    _defaultDrawOptionsDict = _molDrawOptionsToDict()
  return {
    key: value
    for key, value in molDrawOptionsAsDict.items()
    if not _isDrawOptionEqual(value, _defaultDrawOptionsDict[key])
  }


def setEnabled():
  """ Enable interactive molecule rendering """
  global _enabled
  if _enabled:
    return
  _enabled = True
  div_uuid = str(uuid.uuid1())
  return HTML(f"""\
<div
class="lm-Widget p-Widget jp-RenderedText jp-mod-trusted jp-OutputArea-output"
id="{div_uuid}"
>Loading rdkit-structure-renderer.js...</div>
<script type="module">
const jsLoader = document.getElementById('{div_uuid}') || {{}};
const setError = e => jsLoader.innerHTML = 'Failed to load rdkit-structure-renderer.js:<br>' +
e.toString() + '<br>' +
'Interactive molecule rendering will not be available in this Jupyter Notebook.'
try {{
import('{rdkitStructureRendererJsUrl}').then(
  ({{ default: Renderer }}) =>
    Renderer.init('{minimalLibJsUrl}').then(
      () => jsLoader.innerHTML = 'Interactive molecule rendering is available in this Jupyter Notebook.'
    ).catch(
      e => setError(e)
    )
  ).catch(
    e => setError(e)
  );
}} catch(e) {{
setError(e);
}}
</script>""")


def isEnabled(mol=None):
  return _enabled and (mol is None or not isNoInteractive(mol))


def getOpts(mol):
  if not isinstance(mol, Chem.Mol):
    raise ValueError(f"Bad args for {__name__}.getOpts(mol: Chem.Mol)")
  opts = {}
  if hasattr(mol, _opts):
    opts = getattr(mol, _opts)
    if not isinstance(opts, dict):
      opts = {}
  return opts


def setOpts(mol, opts):
  if not isinstance(mol, Chem.Mol) or not isinstance(opts, dict):
    raise ValueError(f"Bad args for {__name__}.setOpts(mol: Chem.Mol, opts: dict)")
  if not all(opts.keys()):
    raise ValueError(
      f"{__name__}.setOpts(mol: Chem.Mol, opts: dict): no key in opts should be null")
  if opts:
    setattr(mol, _opts, opts)
  elif hasattr(mol, _opts):
    delattr(mol, _opts)


def setOpt(mol, key, value):
  if not isinstance(mol, Chem.Mol) or not isinstance(key, str) or not key:
    raise ValueError(f"Bad args for {__name__}.setOpt(mol: Chem.Mol, key: str, value: Any)")
  opts = getOpts(mol)
  opts[key] = value
  setOpts(mol, opts)


def clearOpts(mol):
  if not isinstance(mol, Chem.Mol):
    raise ValueError(f"Bad args for {__name__}.clearOpts(mol: Chem.Mol)")
  setOpts(mol, {})


def clearOpt(mol, key):
  if not isinstance(mol, Chem.Mol) or not isinstance(key, str):
    raise ValueError(f"Bad args for {__name__}.clearOpt(mol: Chem.Mol, key: str)")
  opts = getOpts(mol)
  if key in opts:
    opts.pop(key)
  setOpts(mol, opts)


def setNoInteractive(mol, shouldSet=True):
  opts = getOpts(mol)
  if shouldSet:
    opts[_disabled] = True
  elif _disabled in opts:
    opts.pop(_disabled)
  setOpts(mol, opts)


def isNoInteractive(mol):
  opts = getOpts(mol)
  return _disabled in opts


def injectHTMLHeaderBeforeTable(html):
  doc = minidom.parseString(html.replace(" scoped", ""))
  table = doc.getElementsByTagName("table")
  if len(table) != 1:
    return html
  table = table.pop(0)
  tbody = table.getElementsByTagName("tbody")
  if tbody:
    if len(tbody) != 1:
      return html
    tbody = tbody.pop(0)
  else:
    tbody = table
  div_list = tbody.getElementsByTagName("div")
  if any(div.getAttribute("class") == "rdk-str-rnr-mol-container" for div in div_list):
    return generateHTMLHeader(doc, table)
  return html

def generateHTMLHeader(doc, element):
  element_parent = element.parentNode
  script = doc.createElement("script")
  script.setAttribute("type", "module")
  # avoid arrow function as minidom encodes => into HTML (grrr)
  # Also use single quotes to avoid the &quot; encoding
  cmd = doc.createTextNode(f"""\
try {{
import('{rdkitStructureRendererJsUrl}').then(
  function({{ default: Renderer }}) {{
    Renderer.init('{minimalLibJsUrl}').then(
      function(Renderer) {{
        Renderer.updateMolDrawDivs();
      }}
    ).catch()
  }}
).catch();
}} catch {{}}""")
  script.appendChild(cmd)
  element_parent.appendChild(script)
  html = doc.toxml()
  return html


def newlineToXml(molblock):
  return molblock.replace("\n", "<\\n>")


def xmlToNewline(xmlblock):
  return xmlblock.replace("&lt;\\n&gt;", "&#10;")


def toDataMol(mol):
  return (mol.GetNumConformers() and newlineToXml(Chem.MolToMolBlock(mol)) or Chem.MolToSmiles(mol))


def _dashLower(m):
  return "-" + m.group(0).lower()


def camelCaseOptToDataTag(opt):
  tag = _camelCaseOptToTagRe.sub(_dashLower, opt)
  if not tag.startswith("data-"):
    tag = "data-" + tag
  return tag


def generateHTMLBody(useSvg, mol, size, drawOptions=None):
  if drawOptions is None:
    drawOptions = _defaultDrawOptions
  doc = minidom.Document()
  unique_id = str(uuid.uuid1())
  div = doc.createElement("div")
  for key, value in [
    ("style", f"width: {size[0]}px; height: {size[1]}px;"),
    ("class", "rdk-str-rnr-mol-container"),
    ("id", f"rdk-str-rnr-mol-{unique_id}"),
    ("data-mol", toDataMol(mol)),
    ("data-content", "rdkit/molecule"),
    ("data-parent-node", parentNodeQuery),
  ]:
    div.setAttribute(key, value)
  userDrawOpts = filterDefaultDrawOpts(drawOptions)
  molOpts = getOpts(mol)
  molOptsDashed = {}
  for key, value in molOpts.items():
    keyDashed = camelCaseOptToDataTag(key)
    if keyDashed == "data-draw-opts":
      if isinstance(value, str):
        value = json.loads(value)
      if not isinstance(value, dict):
        raise ValueError(f"data-draw-opts: expected dict, found {str(type(value))}")
      userDrawOpts.update(value)
    else:
      molOptsDashed[keyDashed] = value
  if "addAtomIndices" in userDrawOpts:
    addAtomIndices = userDrawOpts["addAtomIndices"]
    if addAtomIndices:
      molOptsDashed["data-atom-idx"] = True
    userDrawOpts.pop('addAtomIndices')
  for key, value in molOptsDashed.items():
    if isinstance(value, Chem.Mol):
      value = toDataMol(value)
    elif not isinstance(value, str):
      value = json.dumps(value, separators=(',', ':'))
    div.setAttribute(key, value)
  if userDrawOpts:
    div.setAttribute("data-draw-opts", json.dumps(userDrawOpts, separators=(',', ':')))
  if useSvg:
    div.setAttribute("data-use-svg", "true")
  doc.appendChild(div)
  return xmlToNewline(div.toxml())
