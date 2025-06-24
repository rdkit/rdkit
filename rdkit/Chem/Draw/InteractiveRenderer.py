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
import base64
import json
import logging
import re
import uuid
from xml.dom import minidom

from IPython.display import HTML, display

from rdkit import Chem
from rdkit.Chem import Draw

from . import rdMolDraw2D

log = logging.getLogger(__name__)

rdkitStructureRendererJsUrl = "https://unpkg.com/rdkit-structure-renderer/dist/rdkit-structure-renderer-module.js"
minimalLibJsUrl = "https://unpkg.com/rdkit-structure-renderer/public/RDKit_minimal.js"
parentNodeQuery = "div[class*=jp-NotebookPanel-notebook]"
_enabled_div_uuid = None
_camelCaseOptToTagRe = re.compile("[A-Z]")
_opts = "__rnrOpts"
_disabled = "_disabled"
_defaultDrawOptions = None
_defaultDrawOptionsDict = None


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
    raise ValueError(f"Bad args ({str(type(molDrawOptions))}) for {__name__}.filterDefaultDrawOpts("
                     "molDrawOptions: Chem.Draw.rdMolDraw2D.MolDrawOptions)")
  molDrawOptionsAsDict = _molDrawOptionsToDict(molDrawOptions)
  if (_defaultDrawOptionsDict is None):
    _defaultDrawOptionsDict = _molDrawOptionsToDict()
  return {
    key: value
    for key, value in molDrawOptionsAsDict.items()
    if not _isDrawOptionEqual(value, _defaultDrawOptionsDict[key])
  }


def setEnabled(shouldEnable=True, quiet=False):
  """ Enable interactive molecule rendering """

  def _wrapMsgIntoDiv(uuid, msg, quiet):
    return ('<div '
            'class="lm-Widget p-Widget jp-RenderedText jp-mod-trusted jp-OutputArea-output"'
            f'id="{uuid}">{"" if quiet else msg}</div>')

  global _enabled_div_uuid
  loadingMsg = "Loading rdkit-structure-renderer.js..."
  failedToLoadMsg = "Failed to load rdkit-structure-renderer.js:"
  renderingMsg = "Interactive molecule rendering is {0:s} in this Jupyter Notebook."
  renderingUnavailableMsg = renderingMsg.format("not available")
  renderingEnabledMsg = renderingMsg.format("enabled")
  renderingDisabledMsg = renderingMsg.format("disabled")
  if not shouldEnable:
    _enabled_div_uuid = None
    return display(HTML(_wrapMsgIntoDiv("", renderingDisabledMsg, quiet)))
  if _enabled_div_uuid:
    return display(HTML(_wrapMsgIntoDiv(_enabled_div_uuid, renderingEnabledMsg, quiet)))
  _enabled_div_uuid = str(uuid.uuid1())
  return display(
    HTML(
      _wrapMsgIntoDiv(_enabled_div_uuid, loadingMsg, quiet) + f"""<script type="module">
const jsLoader = document.getElementById('{_enabled_div_uuid}') || {{}};
const setError = (e, resolve) => {{
  jsLoader.innerHTML = (
    '{failedToLoadMsg}<br>' +
    e.toString() + '<br>' +
    '{renderingUnavailableMsg}'
  );
  resolve && resolve();
}};
if (window.rdkStrRnr) {{
  jsLoader.innerHTML = '{renderingEnabledMsg}';
}} else {{
  window.rdkStrRnr = new Promise(resolve => {{
    try {{
      fetch('{rdkitStructureRendererJsUrl}').then(
        r => r.text().then(
          t => import(URL.createObjectURL(new Blob([t], {{ type: 'application/javascript' }}))).then(
            ({{ default: Renderer }}) => {{
              const res = Renderer.init('{minimalLibJsUrl}');
              return res.then(
                Renderer => {{
                  jsLoader.innerHTML = '{renderingEnabledMsg}';
                  resolve(Renderer);
                }}
              ).catch(
                e => setError(e, resolve)
              );
            }}
          ).catch(
              e => setError(e, resolve)
          )
        ).catch(
          e => setError(e, resolve)
        )
      ).catch(
        e => setError(e, resolve)
      );
    }} catch(e) {{
      setError(e, resolve);
    }}
  }});
}}
</script>"""))


def isEnabled(mol=None):
  return _enabled_div_uuid and (mol is None or not isNoInteractive(mol))


def getOpts(mol):
  if not isinstance(mol, Chem.Mol):
    raise ValueError(f"Bad args ({str(type(mol))}) for {__name__}.getOpts(mol: Chem.Mol)")
  opts = {}
  if hasattr(mol, _opts):
    opts = getattr(mol, _opts)
    if not isinstance(opts, dict):
      opts = {}
  return opts


def setOpts(mol, opts):
  if not isinstance(mol, Chem.Mol) or not isinstance(opts, dict):
    raise ValueError(
      f"Bad args ({str(type(mol)), str(type(opts))}) for {__name__}.setOpts(mol: Chem.Mol, opts: dict)"
    )
  if not all(opts.keys()):
    raise ValueError(
      f"{__name__}.setOpts(mol: Chem.Mol, opts: dict): no key in opts should be null")
  if opts:
    setattr(mol, _opts, opts)
  elif hasattr(mol, _opts):
    delattr(mol, _opts)


def setOpt(mol, key, value):
  if not isinstance(mol, Chem.Mol) or not isinstance(key, str) or not key:
    raise ValueError(
      f"Bad args ({str(type(mol))}, {str(type(key))}) for {__name__}.setOpt(mol: Chem.Mol, key: str, value: Any)"
    )
  opts = getOpts(mol)
  opts[key] = value
  setOpts(mol, opts)


def clearOpts(mol):
  if not isinstance(mol, Chem.Mol):
    raise ValueError(f"Bad args ({str(type(mol))}) for {__name__}.clearOpts(mol: Chem.Mol)")
  setOpts(mol, {})


def clearOpt(mol, key):
  if not isinstance(mol, Chem.Mol) or not isinstance(key, str):
    raise ValueError(
      f"Bad args ({str(type(mol))}, {str(type(key))}) for {__name__}.clearOpt(mol: Chem.Mol, key: str)"
    )
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


def injectHTMLFooterAfterTable(html):
  doc = minidom.parseString(html.replace(" scoped", ""))
  tables = doc.getElementsByTagName("table")
  for table in tables:
    tbody = table.getElementsByTagName("tbody")
    if tbody:
      if len(tbody) != 1:
        return html
      tbody = tbody.pop(0)
    else:
      tbody = table
    div_list = tbody.getElementsByTagName("div")
    if any(div.getAttribute("class") == "rdk-str-rnr-mol-container" for div in div_list):
      return generateHTMLFooter(doc, table)
  return html


def generateHTMLFooter(doc, element):
  element_parent = element.parentNode
  if element_parent.nodeName.lower() != "div":
    element_parent = doc.createElement("div")
    element_parent.appendChild(element)
  script = doc.createElement("script")
  script.setAttribute("type", "module")
  # avoid arrow function as minidom encodes => into HTML (grrr)
  # Also use single quotes to avoid the &quot; encoding
  cmd = doc.createTextNode(f"""\
if (window.rdkStrRnr) {{
  window.rdkStrRnr.then(
    function(Renderer) {{
      if (Renderer) {{
        Renderer.updateMolDrawDivs();
      }}
    }}
  ).catch();
}}""")
  script.appendChild(cmd)
  element_parent.appendChild(script)
  html = element_parent.toxml()
  return html


def newlineToXml(molblock):
  return molblock.replace("\n", "<\\n>")


def xmlToNewline(xmlblock):
  return xmlblock.replace("&lt;\\n&gt;", "&#10;")


def toDataMol(mol):
  return "pkl_" + base64.b64encode(
    mol.ToBinary(Chem.PropertyPickleOptions.AllProps
                 ^ Chem.PropertyPickleOptions.ComputedProps)).decode("utf-8")


def _dashLower(m):
  return "-" + m.group(0).lower()


def camelCaseOptToDataTag(opt):
  tag = _camelCaseOptToTagRe.sub(_dashLower, opt)
  if not tag.startswith("data-"):
    tag = "data-" + tag
  return tag


def generateHTMLBody(mol, size, **kwargs):

  def toJson(x):
    return json.dumps(x, separators=(",", ":"))

  drawOptions = kwargs.get("drawOptions", _defaultDrawOptions)
  legend = kwargs.get("legend", None)
  useSVG = kwargs.get("useSVG", False)
  kekulize = Draw.shouldKekulize(mol, kwargs.get("kekulize", True))
  highlightAtoms = kwargs.get("highlightAtoms", []) or []
  highlightBonds = kwargs.get("highlightBonds", []) or []
  if not highlightAtoms and hasattr(mol, '__sssAtoms'):
    highlightAtoms = mol.__sssAtoms
    highlightBonds = [
      b.GetIdx() for b in mol.GetBonds()
      if b.GetBeginAtomIdx() in highlightAtoms and b.GetEndAtomIdx() in highlightAtoms
    ]
  doc = minidom.Document()
  unique_id = str(uuid.uuid1())
  div = doc.createElement("div")
  for key, value in [
    ("style", f"margin: auto;"),
    ("class", "rdk-str-rnr-mol-container"),
    ("id", f"rdk-str-rnr-mol-{unique_id}"),
    ("data-mol", toDataMol(mol)),
    ("data-content", "rdkit/molecule"),
    ("data-parent-node", parentNodeQuery),
    ("data-width", str(size[0])),
    ("data-height", str(size[1])),
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
    userDrawOpts.pop("addAtomIndices")
  if highlightAtoms or highlightBonds:
    molOptsDashed["data-scaffold-highlight"] = True
    userDrawOpts["atoms"] = highlightAtoms
    userDrawOpts["bonds"] = highlightBonds
  if not kekulize:
    userDrawOpts["kekulize"] = False
  for key, value in molOptsDashed.items():
    if isinstance(value, Chem.Mol):
      value = toDataMol(value)
    elif not any(isinstance(value, t) for t in (str, int, float)):
      value = toJson(value)
    div.setAttribute(key, str(value))
  if userDrawOpts:
    div.setAttribute("data-draw-opts", toJson(userDrawOpts))
  if useSVG:
    div.setAttribute("data-use-svg", "true")
  if legend:
    outerTable = doc.createElement("table")
    outerTable.setAttribute("style", f"margin: auto;")
    molTr = doc.createElement("tr")
    molTd = doc.createElement("td")
    molTd.setAttribute("style", f"padding: 0;")
    molTd.appendChild(div)
    molTr.appendChild(molTd)
    nameTr = doc.createElement("tr")
    nameTh = doc.createElement("th")
    legendText = doc.createTextNode(legend)
    nameTh.appendChild(legendText)
    nameTh.setAttribute("style", f"text-align: center; background: white;")
    nameTr.appendChild(nameTh)
    outerTable.appendChild(molTr)
    outerTable.appendChild(nameTr)
    div = outerTable
  return xmlToNewline(div.toxml())


def MolsToHTMLTable(mols, molsPerRow=3, subImgSize=(200, 200), legends=None,
                    highlightAtomLists=None, highlightBondLists=None, useSVG=False,
                    drawOptions=None, **kwargs):
  if legends and len(legends) > len(mols):
    legends = legends[:len(mols)]
  if highlightAtomLists and len(highlightAtomLists) > len(mols):
    highlightAtomLists = highlightAtomLists[:len(mols)]
  if highlightBondLists and len(highlightBondLists) > len(mols):
    highlightBondLists = highlightBondLists[:len(mols)]
  nRows = (len(mols) - 1) // molsPerRow + 1 if mols else 0
  if nRows:
    doc = minidom.Document()
    table = doc.createElement("table")
  res = ""
  i = 0
  for _ in range(nRows):
    tr = doc.createElement("tr")
    for _ in range(molsPerRow):
      td = doc.createElement("td")
      td.setAttribute("style", f"padding: 0; background-color: white;")
      highlights = None
      mol = mols[i]
      legend = legends[i] if legends else None
      highlights = highlightAtomLists[i] if highlightAtomLists else None
      kwargs["highlightBonds"] = highlightBondLists[i] if highlightBondLists else None
      content = None
      if isinstance(mol, Chem.Mol):
        if isEnabled(mol):
          content = generateHTMLBody(mol, subImgSize, drawOptions=drawOptions, legend=legend,
                                     useSVG=useSVG, highlightAtoms=highlights, **kwargs)
        else:
          fn = Draw._moltoSVG if useSVG else Draw._moltoimg
          content = fn(mol, subImgSize, highlights, legend, **kwargs)
      try:
        content = minidom.parseString(content)
        td.appendChild(content.firstChild)
      except Exception as e:
        log.warning("Failed to parse HTML returned by generateHTMLBody()")
      tr.appendChild(td)
      i += 1
      if i == len(mols):
        break
    table.appendChild(tr)
  if nRows:
    doc.appendChild(table)
    res = doc.toxml()
  return injectHTMLFooterAfterTable(res)
