#
#  Copyright (C) 2011-2021 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import base64
import copy
import html
import warnings
from io import BytesIO

import IPython
from IPython.display import HTML, SVG

from rdkit import Chem
from rdkit.Chem import Draw, rdchem, rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D

from . import InteractiveRenderer

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
highlightSubstructs = True      # highlight substructure matches when drawing 2D structures
kekulizeStructures = True       # try to kekulize structures before drawing
highlightByReactant = False     # highlight reactions by reactant
ipython_useSVG = False          # use the SVG renderer for 2D depictions (otherwise PNG is used)
ipython_showProperties = True   # display molecule properties when rendering them in 2D
ipython_maxProperties = 10      # maximum number of properties to display
ipython_3d = False              # try to use 3D rendering for molecules with a 3D conformer
molSize_3d = (400, 400)         # default size of 3D structures
drawing_type_3d = 'stick'       # default drawing type for 3d structures
bgcolor_3d = '0xeeeeee'         # default background color for 3d structures
drawOptions = rdMolDraw2D.MolDrawOptions()   # drawing options for 2D structures
InteractiveRenderer._defaultDrawOptions = drawOptions


def addMolToView(mol, view, confId=-1, drawAs=None):
  ''' adds a single molecule to a py3Dmol view

Parameters
----------
m : RDKit molecule
    the molecule to draw
view : py3Dmol.view
    the view to draw into
confId : int
    the conformer ID to draw. If not provided, the default conformer (-1) will be used
drawAs : str
    the drawing type to use. If not provided, 'stick' will be used for small molecules, 
    'cartoon' for large (more than 999 atoms) molecules
    
Returns
-------
None'''
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
  ''' draws a single molecule in 3D using py3Dmol

Parameters
----------
m : RDKit molecule
    the molecule to draw
view : py3Dmol.view
    the view to draw into. If not provided, one will be created
confId : int
    the conformer ID to draw. If not provided, the default conformer (-1) will be used
drawAs : str
    the drawing type to use. If not provided, 'stick' will be used
bgColor : str
    the background color to use for the drawing. If not provided, the default background color will be used
size : tuple of int
    the size of the drawing. If not provided, molSize_3d will be used

Returns
-------
the py3Dmol.view object containing the drawing'''
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

def drawMols3D(mols, view=None, confIds=None, drawAs=None, bgColor=None, size=None, removeHs=False,
               colors=('cyanCarbon','redCarbon','blueCarbon','magentaCarbon','yellowCarbon','cornflowerblueCarbon')):
  ''' draws a list/tuple of molecules in 3D using py3Dmol

Parameters
----------
mols : sequence (list, tuple, etc.) of RDKit molecules
    the molecules to draw
view : py3Dmol.view
    the view to draw into. If not provided, one will be created
confIds : sequence (list, tuple, etc.) of int
    the conformer ID to draw for each molecule. If not provided, the default conformer (-1) will be used
drawAs : str or sequence of str
    the drawing type to use for each molecule. If not provided, 'stick' will be used
bgColor : str
    the background color to use for the drawing. If not provided, the default background color will be used
size : tuple of int
    the size of the drawing. If not provided, molSize_3d will be used
removeHs : bool
    whether or not to remove Hs from the molecules before drawing, the default is False
colors : sequence of str
    the colors to use for drawing the molecules. 

Returns
-------
the py3Dmol.view object containing the drawing'''
  if bgColor is None:
    bgColor = bgcolor_3d
  if size is None:
    size = molSize_3d
  if view is None:
    view = py3Dmol.view(width=size[0], height=size[1])

  if drawAs is None:
    drawAs = ['stick']*len(mols)
  if confIds is None:
    confIds = [-1] * len(mols)

  for m,confId in zip(mols, confIds):
    if removeHs:
      m = Chem.RemoveHs(m)
    addMolToView(m, view, confId)
  for i in range(len(mols)):
    view.setStyle({'model': i,}, {drawAs[i]: {'colorscheme': colors[i % len(colors)]}})

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


def _wrapHTMLIntoTable(html):
  return InteractiveRenderer.injectHTMLFooterAfterTable(
    f'<div><table><tbody><tr><td style="width: {molSize[0]}px; ' +
    f'height: {molSize[1]}px; text-align: center;">' + html.replace(" scoped", "") +
    '</td></tr></tbody></table></div>')


def _toHTML(mol):
  useInteractiveRenderer = InteractiveRenderer.isEnabled(mol)
  if _canUse3D and ipython_3d and mol.GetNumConformers() and mol.GetConformer().Is3D():
    return _toJSON(mol)
  props = mol.GetPropsAsDict()
  if not ipython_showProperties or not props:
    if useInteractiveRenderer:
      return _wrapHTMLIntoTable(
        InteractiveRenderer.generateHTMLBody(mol, molSize, useSVG=ipython_useSVG))
    else:
      return _toSVG(mol)
  if mol.HasProp('_Name'):
    nm = mol.GetProp('_Name')
  else:
    nm = ''

  res = []

  if useInteractiveRenderer:
    content = InteractiveRenderer.generateHTMLBody(mol, molSize, legend=nm, useSVG=ipython_useSVG)
  else:
    if not ipython_useSVG:
      png = _toPNG(mol)
      png = base64.b64encode(png)
      content = f'<image src="data:image/png;base64,{png.decode()}">'
    else:
      content = _toSVG(mol)
  res.append(f'<tr><td colspan="2" style="text-align: center;">{content}</td></tr>')

  for i, (pn, pv) in enumerate(props.items()):
    if ipython_maxProperties >= 0 and i >= ipython_maxProperties:
      res.append(
        '<tr><td colspan="2" style="text-align: center">Property list truncated.<br />Increase IPythonConsole.ipython_maxProperties (or set it to -1) to see more properties.</td></tr>'
      )
      break
    pv = html.escape(str(pv))
    res.append(
      f'<tr><th style="text-align: right">{pn}</th><td style="text-align: left">{pv}</td></tr>')
  res = '\n'.join(res)
  res = f'<table>{res}</table>'
  if useInteractiveRenderer:
    res = InteractiveRenderer.injectHTMLFooterAfterTable(res)
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
  if Draw._MolsToGridImageSaved is not None:
    fn = Draw._MolsToGridImageSaved
  else:
    fn = Draw.MolsToGridImage
  return fn(bundle, subImgSize=molSize, drawOptions=drawOptions, useSVG=False, returnPNG=True)


def _toMolBundleSVG(bundle):
  if not ipython_useSVG:
    return None
  if Draw._MolsToGridImageSaved is not None:
    fn = Draw._MolsToGridImageSaved
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


from IPython import display


def ShowMols(mols, maxMols=50, **kwargs):
  if 'useSVG' not in kwargs:
    kwargs['useSVG'] = ipython_useSVG
  if 'returnPNG' not in kwargs:
    kwargs['returnPNG'] = True
  if InteractiveRenderer.isEnabled():
    fn = InteractiveRenderer.MolsToHTMLTable
  elif Draw._MolsToGridImageSaved is not None:
    fn = Draw._MolsToGridImageSaved
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
  if InteractiveRenderer.isEnabled():
    return HTML(res)
  elif kwargs['useSVG']:
    return SVG(res)
  elif kwargs['returnPNG']:
    return display.Image(data=res, format='png')
  return res

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

def DrawMorganBit(mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs):
  if Draw._DrawMorganBitSaved is not None:
    fn = Draw._DrawMorganBitSaved
  else:
    fn = Draw.DrawMorganBit
  return _DrawBit(fn, mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs)

def DrawMorganBits(*args, drawOptions=drawOptions, **kwargs):
  if Draw._DrawMorganBitsSaved is not None:
    fn = Draw._DrawMorganBitsSaved
  else:
    fn = Draw.DrawMorganBits
  return _DrawBit(fn, *args, drawOptions=drawOptions, **kwargs)

def DrawRDKitBit(mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs):
  if Draw._DrawRDKitBitSaved is not None:
    fn = Draw._DrawRDKitBitSaved
  else:
    fn = Draw.DrawRDKitBit
  return _DrawBit(fn, mol, bitId, bitInfo, drawOptions=drawOptions, **kwargs)

def DrawRDKitBits(*args, drawOptions=drawOptions, **kwargs):
  if Draw._DrawRDKitBitsSaved is not None:
    fn = Draw._DrawRDKitBitsSaved
  else:
    fn = Draw.DrawRDKitBits
  return _DrawBit(fn, *args, drawOptions=drawOptions, **kwargs)


def EnableSubstructMatchRendering():
  if not hasattr(rdchem.Mol, '__GetSubstructMatch'):
    rdchem.Mol.__GetSubstructMatch = rdchem.Mol.GetSubstructMatch
  rdchem.Mol.GetSubstructMatch = _GetSubstructMatch
  if not hasattr(rdchem.Mol, '__GetSubstructMatches'):
    rdchem.Mol.__GetSubstructMatches = rdchem.Mol.GetSubstructMatches
  rdchem.Mol.GetSubstructMatches = _GetSubstructMatches


_methodsToDelete = []
_rendererInstalled = False

def InstallIPythonRenderer():
  global _rendererInstalled
  if _rendererInstalled:
    return

  # in case of a reload, there may be some cleanup necessary in other modules that we have
  # monkey patched. So go ahead and do an uninstall to cleanup before we start. This was
  # Github #8082
  UninstallIPythonRenderer()

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
  Draw._MolsToGridImageSaved = Draw.MolsToGridImage
  Draw.MolsToGridImage = ShowMols
  Draw.MolsToGridImage.__doc__ = Draw._MolsToGridImageSaved.__doc__
  Draw._DrawRDKitBitSaved = Draw.DrawRDKitBit
  Draw.DrawRDKitBit = DrawRDKitBit
  Draw.DrawRDKitBit.__doc__ = Draw._DrawRDKitBitSaved.__doc__
  Draw._DrawRDKitBitsSaved = Draw.DrawRDKitBits
  Draw.DrawRDKitBits = DrawRDKitBits
  Draw.DrawRDKitBits.__doc__ = Draw._DrawRDKitBitsSaved.__doc__
  Draw._DrawMorganBitSaved = Draw.DrawMorganBit
  Draw.DrawMorganBit = DrawMorganBit
  Draw.DrawMorganBit.__doc__ = Draw.DrawMorganBit.__doc__
  Draw._DrawMorganBitsSaved = Draw.DrawMorganBits
  Draw.DrawMorganBits = DrawMorganBits
  Draw.DrawMorganBits.__doc__ = Draw._DrawMorganBitsSaved.__doc__
  if not hasattr(rdchem.Mol, '__DebugMol'):
    rdchem.Mol.__DebugMol = rdchem.Mol.Debug
    rdchem.Mol.Debug = lambda self, useStdout=False: self.__DebugMol(useStdout=useStdout)
        
  _rendererInstalled = True

def DisableSubstructMatchRendering():
  if hasattr(rdchem.Mol, '__GetSubstructMatch'):
    rdchem.Mol.GetSubstructMatch = rdchem.Mol.__GetSubstructMatch
    del rdchem.Mol.__GetSubstructMatch
  if hasattr(rdchem.Mol, '__GetSubstructMatches'):
    rdchem.Mol.GetSubstructMatches = rdchem.Mol.__GetSubstructMatches
    del rdchem.Mol.__GetSubstructMatches


def UninstallIPythonRenderer():
  global _rendererInstalled, _methodsToDelete

  for cls, attr in _methodsToDelete:
    if hasattr(cls,attr):
      delattr(cls, attr)

  _methodsToDelete = []

  DisableSubstructMatchRendering()

  if hasattr(Draw,'_MolsToGridImageSaved') and Draw._MolsToGridImageSaved is not None:
    Draw.MolsToGridImage = Draw._MolsToGridImageSaved
    del Draw._MolsToGridImageSaved
  if hasattr(Draw,'_DrawRDKitBitSaved') and Draw._DrawRDKitBitSaved is not None:
    Draw.DrawRDKitBit = Draw._DrawRDKitBitSaved
    del Draw._DrawRDKitBitSaved
  if hasattr(Draw,'_DrawRDKitBitsSaved') and Draw._DrawRDKitBitsSaved is not None:
    Draw.DrawRDKitBits = Draw._DrawRDKitBitsSaved
    del Draw._DrawRDKitBitsSaved
  if hasattr(Draw,'_DrawMorganBitSaved') and Draw._DrawMorganBitSaved is not None:
    Draw.DrawMorganBit = Draw._DrawMorganBitSaved
    del Draw._DrawMorganBitSaved
  if hasattr(Draw,'_DrawMorganBitsSaved') and Draw._DrawMorganBitsSaved is not None:
    Draw.DrawMorganBits = Draw._DrawMorganBitsSaved
    del Draw._DrawMorganBitsSaved

  if hasattr(rdchem.Mol, '__DebugMol'):
    rdchem.Mol.Debug = rdchem.Mol.__DebugMol
    del rdchem.Mol.__DebugMol
    
  _rendererInstalled = False

InstallIPythonRenderer()

