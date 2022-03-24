#
#  Copyright (C) 2011-2021 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from IPython.display import SVG
from xml.dom import minidom
import uuid
import json
import re
import base64
import copy
import warnings
from io import BytesIO
import IPython
from IPython.display import HTML, SVG
from rdkit import Chem
from rdkit.Chem import Draw, rdchem, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

if IPython.release.version < '0.11':
  raise ImportError('this module requires at least v0.11 of IPython')
try:
  import py3Dmol
  _canUse3D = True
except ImportError:
  _canUse3D = False

from PIL import Image
from PIL.PngImagePlugin import PngInfo

molSize = (450, 150)
highlightSubstructs = True
kekulizeStructures = True
highlightByReactant = False
ipython_useSVG = False
ipython_showProperties = True
ipython_maxProperties = 10
ipython_3d = False
molSize_3d = (400, 400)
drawing_type_3d = 'stick'  # default drawing type for 3d structures
bgcolor_3d = '0xeeeeee'
drawOptions = rdMolDraw2D.MolDrawOptions()

class InteractiveRenderer:
  """ Interactive molecule rendering through rdkit-structure-renderer.js """

  rdkitStructureRendererJsUrl = "https://unpkg.com/rdkit-structure-renderer/dist/rdkit-structure-renderer-module.js"
  minimalLibJsUrl = "https://unpkg.com/rdkit-structure-renderer/public/RDKit_minimal.js"
  parentNodeQuery = "div[class*=jp-NotebookPanel-notebook]"
  _enabled = False
  _camelCaseOptToTagRe = re.compile("[A-Z]")
  _opts = "__rnrOpts"
  _disabled = "_disabled"
  _defaultDrawOptions = None

  @staticmethod
  def _isAcceptedKeyValue(key, value):
    return not key.startswith("__") and (
      type(value) in (bool, int, float, str, tuple, list) or (
        callable(value) and key.startswith("get")
      )
    )

  @staticmethod
  def _getValueFromKey(d, key):
    value = getattr(d, key)
    return value() if callable(value) else value

  @classmethod
  def _molDrawOptionsToDict(cls, molDrawOptions=rdMolDraw2D.MolDrawOptions()):
    return {key: cls._getValueFromKey(molDrawOptions, key) for key in dir(molDrawOptions)
      if cls._isAcceptedKeyValue(key, getattr(molDrawOptions, key))}

  @classmethod
  def _isDrawOptionEqual(cls, v1, v2):
    if type(v1) != type(v2):
      return False
    if type(v1) in (tuple, list):
      if len(v1) != len(v2):
        return False
      return all(cls._isDrawOptionEqual(item1, v2[i]) for i, item1 in enumerate(v1))
    if type(v1) == float:
      return abs(v1 - v2) < 1.e-5
    return v1 == v2

  @classmethod
  def filterDefaultDrawOpts(cls, molDrawOptions):
    if not isinstance(molDrawOptions, rdMolDraw2D.MolDrawOptions):
      raise ValueError(f"Bad args for {cls.__name__}.filterDefaultDrawOpts("
        "molDrawOptions: Chem.Draw.rdMolDraw2D.MolDrawOptions)")
    molDrawOptionsAsDict = cls._molDrawOptionsToDict(molDrawOptions)
    if (cls._defaultDrawOptions is None):
      cls._defaultDrawOptions = cls._molDrawOptionsToDict()
    return {key: value for key, value in molDrawOptionsAsDict.items()
      if not cls._isDrawOptionEqual(value, cls._defaultDrawOptions[key])}

  @classmethod
  def setEnabled(cls):
    """ Enable interactive molecule rendering """
    if cls._enabled:
      return
    cls._enabled = True
    div_uuid = str(uuid.uuid1())
    return HTML(f"""\
<div
  class="lm-Widget p-Widget jp-RenderedText jp-mod-trusted jp-OutputArea-output"
  id="{div_uuid}"
>Loading rdkit-structure-renderer.js</div>
<script type="module">
const jsLoader = document.getElementById('{div_uuid}') || {{}};
const setError = e => jsLoader.innerHTML = 'Failed to load rdkit-structure-renderer.js:<br>' +
  e.toString() + '<br>' +
  'Interactive molecule rendering will not be available in this Jupyter Notebook.<br>'
try {{
  import('{cls.rdkitStructureRendererJsUrl}').then(
    ({{ default: Renderer }}) =>
      Renderer.init('{cls.minimalLibJsUrl}').then(
        () => jsLoader.innerHTML = 'Using rdkit-structure-renderer.js'
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

  @classmethod
  def isEnabled(cls, mol=None):
    return cls._enabled and (mol is None or not cls.isNoInteractive(mol))

  @classmethod
  def getOpts(cls, mol):
    if not isinstance(mol, Chem.Mol):
      raise ValueError(f"Bad args for {cls.__name__}.getOpts(mol: Chem.Mol)")
    opts = {}
    if hasattr(mol, cls._opts):
      opts = getattr(mol, cls._opts)
      if not isinstance(opts, dict):
        opts = {}
    return opts

  @classmethod
  def setOpts(cls, mol, opts):
    if not isinstance(mol, Chem.Mol) or not isinstance(opts, dict):
      raise ValueError(f"Bad args for {cls.__name__}.setOpts(mol: Chem.Mol, opts: dict)")
    if not all(opts.keys()):
      raise ValueError(f"{cls.__name__}.setOpts(mol: Chem.Mol, opts: dict): no key in opts should be null")
    if opts:
      setattr(mol, cls._opts, opts)
    else:
      delattr(mol, cls._opts)

  @classmethod
  def setOpt(cls, mol, key, value):
    if not isinstance(mol, Chem.Mol) or not isinstance(key, str) or not key:
      raise ValueError(f"Bad args for {cls.__name__}.setOpt(mol: Chem.Mol, key: str, value: Any)")
    opts = cls.getOpts(cls, mol)
    opts[key] = value
    cls.setOpts(mol, opts)

  @classmethod
  def clearOpts(cls, mol):
    if not isinstance(mol, Chem.Mol):
      raise ValueError(f"Bad args for {cls.__name__}.clearOpts(mol: Chem.Mol)")
    cls.setOpts(mol, {})

  @classmethod
  def clearOpt(cls, mol, key):
    if not isinstance(mol, Chem.Mol) or not isinstance(key, str):
      raise ValueError(f"Bad args for {cls.__name__}.clearOpt(mol: Chem.Mol, key: str)")
    opts = cls.getOpts(cls, mol)
    if key in opts:
      opts.pop(key)
    cls.setOpts(mol, opts)

  @classmethod
  def setNoInteractive(cls, mol, shouldSet=True):
    opts = cls.getOpts(mol)
    if shouldSet:
      opts[cls._disabled] = True
    elif cls._disabled in opts:
      opts.pop(cls._disabled)
    cls.setOpts(mol, opts)

  @classmethod
  def isNoInteractive(cls, mol):
    opts = cls.getOpts(mol)
    return cls._disabled in opts

  @classmethod
  def injectHTMLHeaderBeforeTable(cls, html):
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
      return cls.generateHTMLHeader(doc, table)
    return html

  @classmethod
  def wrapHTMLIntoTable(cls, html):
    return cls.injectHTMLHeaderBeforeTable(
      f'<div><table><tbody><tr><td style="width: {molSize[0]}px; ' +
      f'height: {molSize[1]}px; text-align: center;">' +
      html.replace(" scoped", "") +
      '</td></tr></tbody></table></div>')

  @classmethod
  def generateHTMLHeader(cls, doc, element):
    element_parent = element.parentNode
    script = doc.createElement("script")
    script.setAttribute("type", "module")
    # avoid arrow function as minidom encodes => into HTML (grrr)
    # Also use single quotes to avoid the &quot; encoding
    cmd = doc.createTextNode(f"""\
try {{
  import('{cls.rdkitStructureRendererJsUrl}').then(
    function({{ default: Renderer }}) {{
      Renderer.init('{cls.minimalLibJsUrl}').then(
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

  @staticmethod
  def newlineToXml(molblock):
    return molblock.replace("\n", "<\\n>")

  @staticmethod
  def xmlToNewline(xmlblock):
    return xmlblock.replace("&lt;\\n&gt;", "&#10;")

  @classmethod
  def toDataMol(cls, mol):
    return (mol.GetNumConformers() and
      cls.newlineToXml(Chem.MolToMolBlock(mol)) or
      Chem.MolToSmiles(mol))

  @staticmethod
  def _dashLower(m):
    return "-" + m.group(0).lower()

  @classmethod
  def camelCaseOptToDataTag(cls, opt):
    tag = cls._camelCaseOptToTagRe.sub(cls._dashLower, opt)
    if not tag.startswith("data-"):
      tag = "data-" + tag
    return tag

  @classmethod
  def generateHTMLBody(cls, useSvg, mol, size):
    doc = minidom.Document()
    unique_id = str(uuid.uuid1())
    div = doc.createElement("div")
    for key, value in [
      ("style", f"width: {size[0]}px; height: {size[1]}px;"),
      ("class", "rdk-str-rnr-mol-container"),
      ("id", f"rdk-str-rnr-mol-{unique_id}"),
      ("data-mol", cls.toDataMol(mol)),
      ("data-content", "rdkit/molecule"),
      ("data-parent-node", cls.parentNodeQuery),
    ]:
      div.setAttribute(key, value)
    userDrawOpts = cls.filterDefaultDrawOpts(drawOptions)
    molOpts = cls.getOpts(mol)
    molOptsDashed = {}
    for key, value in molOpts.items():
      keyDashed = cls.camelCaseOptToDataTag(key)
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
        value = cls.toDataMol(value)
      elif not isinstance(value, str):
        value = json.dumps(value, separators=(',', ':'))
      div.setAttribute(key, value)
    if userDrawOpts:
      div.setAttribute("data-draw-opts", json.dumps(userDrawOpts, separators=(',', ':')))
    if useSvg:
      div.setAttribute("data-use-svg", "true")
    doc.appendChild(div)
    return cls.xmlToNewline(div.toxml())


def addMolToView(mol, view, confId=-1, drawAs=None):
  if mol.GetNumAtoms() >= 999 or drawAs == 'cartoon':
    # py3DMol is happier with TER and MASTER records present
    pdb = Chem.MolToPDBBlock(mol, flavor=0x20 | 0x10)
    view.addModel(pdb, 'pdb')
  else:
    # py3Dmol does not currently support v3k mol files, so
    # we can only provide those with "smaller" molecules
    mb = Chem.MolToMolBlock(mol, confId=confId)
    view.addModel(mb, 'sdf')
  if drawAs is None:
    drawAs = drawing_type_3d
  view.setStyle({drawAs: {}})


def drawMol3D(m, view=None, confId=-1, drawAs=None, bgColor=None, size=None):
  if bgColor is None:
    bgColor = bgcolor_3d
  if size is None:
    size = molSize_3d
  if view is None:
    view = py3Dmol.view(width=size[0], height=size[1])
  view.removeAllModels()

  try:
    ms = iter(m)
    for m in ms:
      addMolToView(m, view, confId, drawAs)
  except TypeError:
    addMolToView(m, view, confId, drawAs)
  view.setBackgroundColor(bgColor)
  view.zoomTo()
  return view.show()


def _toJSON(mol):
  """For IPython notebook, renders 3D webGL objects."""
  if not ipython_3d or not mol.GetNumConformers():
    return None
  conf = mol.GetConformer()
  if not conf.Is3D():
    return None
  res = drawMol3D(mol)
  if hasattr(res, 'data'):
    return res.data
  return ""


def _toHTML(mol):
  useInteractiveRenderer = InteractiveRenderer.isEnabled(mol)
  if _canUse3D and ipython_3d and mol.GetNumConformers():
    return _toJSON(mol)
  props = mol.GetPropsAsDict()
  if not ipython_showProperties or not props:
    if useInteractiveRenderer:
      return InteractiveRenderer.wrapHTMLIntoTable(
        InteractiveRenderer.generateHTMLBody(ipython_useSVG, mol, molSize))
    else:
      return _toSVG(mol)
  if mol.HasProp('_Name'):
    nm = mol.GetProp('_Name')
  else:
    nm = ''

  res = []
  if useInteractiveRenderer:
    content = InteractiveRenderer.generateHTMLBody(ipython_useSVG, mol, molSize)
  else:
    if not ipython_useSVG:
      png = Draw._moltoimg(mol, molSize, [], nm, returnPNG=True, drawOptions=drawOptions)
      png = base64.b64encode(png)
      res.append(f'<tr><td colspan=2 style="text-align:center"><image src="data:image/png;base64,{png.decode()}"></td></tr>')
    else:
      content = Draw._moltoSVG(mol, molSize, [], nm, kekulize=kekulizeStructures, drawOptions=drawOptions)
  res.append(f'<tr><td colspan=2 style="text-align:center">{content}</td></tr>')

  for i,(pn, pv) in enumerate(props.items()):
    if ipython_maxProperties >= 0 and i >= ipython_maxProperties:
      res.append('<tr><td colspan=2 style="text-align:center">Property list truncated.<br />Increase IPythonConsole.ipython_maxProperties (or set it to -1) to see more properties.</td></tr>')
      break
    res.append(
      f'<tr><th style="text-align:right">{pn}</th><td style="text-align:left">{pv}</td></tr>')
  res = '\n'.join(res)
  res = f'<table>{res}</table>'
  if useInteractiveRenderer:
    res = InteractiveRenderer.injectHTMLHeaderBeforeTable(res)
  return res


def _toPNG(mol):
  if hasattr(mol, '__sssAtoms'):
    highlightAtoms = mol.__sssAtoms
  else:
    highlightAtoms = []
  kekulize = kekulizeStructures
  return Draw._moltoimg(mol, molSize, highlightAtoms, "", returnPNG=True, kekulize=kekulize,
                        drawOptions=drawOptions)


def _toSVG(mol):
  if not ipython_useSVG:
    return None
  if hasattr(mol, '__sssAtoms'):
    highlightAtoms = mol.__sssAtoms
  else:
    highlightAtoms = []
  kekulize = kekulizeStructures
  return Draw._moltoSVG(mol, molSize, highlightAtoms, "", kekulize, drawOptions=drawOptions)


def _toReactionPNG(rxn):
  rc = copy.deepcopy(rxn)
  return Draw.ReactionToImage(rc, subImgSize=(int(molSize[0] / 3), molSize[1]),
                              highlightByReactant=highlightByReactant, drawOptions=drawOptions,
                              returnPNG=True)


def _toReactionSVG(rxn):
  if not ipython_useSVG:
    return None
  rc = copy.deepcopy(rxn)
  return Draw.ReactionToImage(rc, subImgSize=(int(molSize[0] / 3), molSize[1]), useSVG=True,
                              highlightByReactant=highlightByReactant, drawOptions=drawOptions)


def _toMolBundlePNG(bundle):
  if _MolsToGridImageSaved is not None:
    fn = _MolsToGridImageSaved
  else:
    fn = Draw.MolsToGridImage
  return fn(bundle, subImgSize=molSize, drawOptions=drawOptions, useSVG=False, returnPNG=True)


def _toMolBundleSVG(bundle):
  if not ipython_useSVG:
    return None
  if _MolsToGridImageSaved is not None:
    fn = _MolsToGridImageSaved
  else:
    fn = Draw.MolsToGridImage
  return fn(bundle, subImgSize=molSize, drawOptions=drawOptions, useSVG=True)


def _GetSubstructMatch(mol, query, *args, **kwargs):
  res = mol.__GetSubstructMatch(query, *args, **kwargs)
  if highlightSubstructs:
    mol.__sssAtoms = list(res)
  else:
    mol.__sssAtoms = []
  return res


_GetSubstructMatch.__doc__ = rdchem.Mol.GetSubstructMatch.__doc__


def _GetSubstructMatches(mol, query, *args, **kwargs):
  res = mol.__GetSubstructMatches(query, *args, **kwargs)
  mol.__sssAtoms = []
  if highlightSubstructs:
    for entry in res:
      mol.__sssAtoms.extend(list(entry))
  return res


_GetSubstructMatches.__doc__ = rdchem.Mol.GetSubstructMatches.__doc__


# code for displaying PIL images directly,
def display_pil_image(img):
  """displayhook function for PIL Images, rendered as PNG"""
  # pull metadata from the image, if there
  metadata = PngInfo()
  for k, v in img.info.items():
    metadata.add_text(k, v)
  bio = BytesIO()
  img.save(bio, format='PNG', pnginfo=metadata)
  return bio.getvalue()


_MolsToGridImageSaved = None

from IPython import display


def ShowMols(mols, maxMols=50, **kwargs):
  global _MolsToGridImageSaved
  if 'useSVG' not in kwargs:
    kwargs['useSVG'] = ipython_useSVG
  if 'returnPNG' not in kwargs:
    kwargs['returnPNG'] = True
  if _MolsToGridImageSaved is not None:
    fn = _MolsToGridImageSaved
  else:
    fn = Draw.MolsToGridImage

  if len(mols) > maxMols:
    warnings.warn(
      "Truncating the list of molecules to be displayed to %d. Change the maxMols value to display more."
      % (maxMols))
    mols = mols[:maxMols]
    for prop in ('legends', 'highlightAtoms', 'highlightBonds'):
      if prop in kwargs:
        kwargs[prop] = kwargs[prop][:maxMols]
  if "drawOptions" not in kwargs:
    kwargs["drawOptions"] = drawOptions

  res = fn(mols, **kwargs)
  if kwargs['useSVG']:
    return SVG(res)
  if kwargs['returnPNG']:
    return display.Image(data=res, format='png')
  return res


ShowMols.__doc__ = Draw.MolsToGridImage.__doc__


def _DrawBit(fn, *args, **kwargs):
  if 'useSVG' not in kwargs:
    kwargs['useSVG'] = ipython_useSVG

  res = fn(*args, **kwargs)
  if kwargs['useSVG']:
    return SVG(res)
  sio = BytesIO(res)
  return Image.open(sio)


def _DrawBits(fn, *args, **kwargs):
  if 'useSVG' not in kwargs:
    kwargs['useSVG'] = ipython_useSVG

  res = fn(*args, **kwargs)
  if kwargs['useSVG']:
    return SVG(res)
  sio = BytesIO(res)
  return Image.open(sio)


_DrawMorganBitSaved = None


def DrawMorganBit(mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs):
  global _DrawMorganBitSaved
  if _DrawMorganBitSaved is not None:
    fn = _DrawMorganBitSaved
  else:
    fn = Draw.DrawMorganBit
  return _DrawBit(fn, mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs)


DrawMorganBit.__doc__ = Draw.DrawMorganBit.__doc__

_DrawMorganBitsSaved = None


def DrawMorganBits(*args, drawOptions=drawOptions, **kwargs):
  global _DrawMorganBitsSaved
  if _DrawMorganBitsSaved is not None:
    fn = _DrawMorganBitsSaved
  else:
    fn = Draw.DrawMorganBits
  return _DrawBit(fn, *args, drawOptions=drawOptions, **kwargs)


DrawMorganBits.__doc__ = Draw.DrawMorganBits.__doc__

_DrawRDKitBitSaved = None


def DrawRDKitBit(mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs):
  global _DrawRDKitBitSaved
  if _DrawRDKitBitSaved is not None:
    fn = _DrawRDKitBitSaved
  else:
    fn = Draw.DrawRDKitBit
  return _DrawBit(fn, mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs)


DrawRDKitBit.__doc__ = Draw.DrawRDKitBit.__doc__

_DrawRDKitBitsSaved = None


def DrawRDKitBits(*args, drawOptions=drawOptions, **kwargs):
  global _DrawRDKitBitsSaved
  if _DrawRDKitBitsSaved is not None:
    fn = _DrawRDKitBitsSaved
  else:
    fn = Draw.DrawRDKitBits
  return _DrawBit(fn, *args, drawOptions=drawOptions, **kwargs)


DrawRDKitBits.__doc__ = Draw.DrawRDKitBits.__doc__

_rendererInstalled = False


def EnableSubstructMatchRendering():
  if not hasattr(rdchem.Mol, '__GetSubstructMatch'):
    rdchem.Mol.__GetSubstructMatch = rdchem.Mol.GetSubstructMatch
  rdchem.Mol.GetSubstructMatch = _GetSubstructMatch
  if not hasattr(rdchem.Mol, '__GetSubstructMatches'):
    rdchem.Mol.__GetSubstructMatches = rdchem.Mol.GetSubstructMatches
  rdchem.Mol.GetSubstructMatches = _GetSubstructMatches


_methodsToDelete = []


def InstallIPythonRenderer():
  global _MolsToGridImageSaved, _DrawRDKitBitSaved, _DrawRDKitBitsSaved, _DrawMorganBitSaved, _DrawMorganBitsSaved
  global _rendererInstalled
  if _rendererInstalled:
    return
  rdchem.Mol._repr_png_ = _toPNG
  rdchem.Mol._repr_svg_ = _toSVG
  _methodsToDelete.append((rdchem.Mol, '_repr_png_'))
  _methodsToDelete.append((rdchem.Mol, '_repr_svg_'))
  rdchem.Mol._repr_html_ = _toHTML
  _methodsToDelete.append((rdchem.Mol, '_repr_html_'))

  rdChemReactions.ChemicalReaction._repr_png_ = _toReactionPNG
  rdChemReactions.ChemicalReaction._repr_svg_ = _toReactionSVG
  _methodsToDelete.append((rdChemReactions.ChemicalReaction, '_repr_png_'))
  _methodsToDelete.append((rdChemReactions.ChemicalReaction, '_repr_svg_'))

  rdchem.MolBundle._repr_png_ = _toMolBundlePNG
  rdchem.MolBundle._repr_svg_ = _toMolBundleSVG
  _methodsToDelete.append((rdchem.MolBundle, '_repr_png_'))
  _methodsToDelete.append((rdchem.MolBundle, '_repr_svg_'))

  EnableSubstructMatchRendering()
  Image.Image._repr_png_ = display_pil_image
  _methodsToDelete.append((Image.Image, '_repr_png_'))
  _MolsToGridImageSaved = Draw.MolsToGridImage
  Draw.MolsToGridImage = ShowMols
  _DrawRDKitBitSaved = Draw.DrawRDKitBit
  Draw.DrawRDKitBit = DrawRDKitBit
  _DrawRDKitBitsSaved = Draw.DrawRDKitBits
  Draw.DrawRDKitBits = DrawRDKitBits
  _DrawMorganBitSaved = Draw.DrawMorganBit
  Draw.DrawMorganBit = DrawMorganBit
  _DrawMorganBitsSaved = Draw.DrawMorganBits
  Draw.DrawMorganBits = DrawMorganBits
  rdchem.Mol.__DebugMol = rdchem.Mol.Debug
  rdchem.Mol.Debug = lambda self, useStdout=False: self.__DebugMol(useStdout=useStdout)
  _rendererInstalled = True


InstallIPythonRenderer()


def DisableSubstructMatchRendering():
  if hasattr(rdchem.Mol, '__GetSubstructMatch'):
    rdchem.Mol.GetSubstructMatch = rdchem.Mol.__GetSubstructMatch
    del rdchem.Mol.__GetSubstructMatch
  if hasattr(rdchem.Mol, '__GetSubstructMatches'):
    rdchem.Mol.GetSubstructMatches = rdchem.Mol.__GetSubstructMatches
    del rdchem.Mol.__GetSubstructMatches


def UninstallIPythonRenderer():
  global _MolsToGridImageSaved, _DrawRDKitBitSaved, _DrawMorganBitSaved, _DrawMorganBitsSaved
  global _rendererInstalled, _methodsToDelete
  if not _rendererInstalled:
    return

  for cls, attr in _methodsToDelete:
    delattr(cls, attr)

  _methodsToDelete = []
  DisableSubstructMatchRendering()
  if _MolsToGridImageSaved is not None:
    Draw.MolsToGridImage = _MolsToGridImageSaved
  if _DrawRDKitBitSaved is not None:
    Draw.DrawRDKitBit = _DrawRDKitBitSaved
  if _DrawRDKitBitsSaved is not None:
    Draw.DrawRDKitBits = _DrawRDKitBitsSaved
  if _DrawMorganBitSaved is not None:
    Draw.DrawMorganBit = _DrawMorganBitSaved
  if _DrawMorganBitsSaved is not None:
    Draw.DrawMorganBits = _DrawMorganBitsSaved
  if hasattr(rdchem.Mol, '__DebugMol'):
    rdchem.Mol.Debug = rdchem.Mol.__DebugMol
    del rdchem.Mol.__DebugMol
  _rendererInstalled = False
