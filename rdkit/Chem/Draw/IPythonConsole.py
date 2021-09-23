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
import numpy
import warnings
import uuid
import json
import os
import base64
import copy
from io import BytesIO, StringIO
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Draw
from rdkit.Chem import rdchem, rdChemReactions
from rdkit import Chem
import sys
import IPython

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

# expose RDLogs to Python StdErr so they are shown
#  in the IPythonConsole not the server logs.
Chem.WrapLogs()


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
    iter(m)
  except TypeError:
    addMolToView(m, view, confId, drawAs)
  else:
    ms = m
    for m in ms:
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
  if _canUse3D and ipython_3d and mol.GetNumConformers():
    return _toJSON(mol)
  props = mol.GetPropsAsDict()
  if not ipython_showProperties or not props:
    return _toSVG(mol)
  if mol.HasProp('_Name'):
    nm = mol.GetProp('_Name')
  else:
    nm = ''
  res = []
  if not ipython_useSVG:
    png = Draw._moltoimg(mol, molSize, [], nm, returnPNG=True, drawOptions=drawOptions)
    png = base64.b64encode(png)
    res.append(f'<tr><td colspan=2 style="text-align:center"><image src="data:image/png;base64,{png.decode()}"></td></tr>')
  else:
    svg = Draw._moltoSVG(mol, molSize, [], nm, kekulize=kekulizeStructures, drawOptions=drawOptions)
    res.append(f'<tr><td colspan=2 style="text-align:center">{svg}</td></tr>')

  for i,(pn, pv) in enumerate(props.items()):
    if ipython_maxProperties>=0 and i>= ipython_maxProperties:
      res.append('<tr><td colspan=2 style="text-align:center">Property list truncated.<br />Increase IPythonConsole.ipython_maxProperties (or set it to -1) to see more properties.</td></tr>')
      break
    res.append(
      f'<tr><th style="text-align:right">{pn}</th><td style="text-align:left">{pv}</td></tr>')
  res = "\n".join(res)
  return f'<table>{res}</table>'


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
  if not "drawOptions" in kwargs:
    kwargs["drawOptions"] = drawOptions
  res = fn(mols, **kwargs)
  if kwargs['useSVG']:
    return SVG(res)
  else:
    if kwargs['returnPNG']:
      return display.Image(data=res, format='png')
    else:
      return res


ShowMols.__doc__ = Draw.MolsToGridImage.__doc__


def _DrawBit(fn, *args, **kwargs):
  if 'useSVG' not in kwargs:
    kwargs['useSVG'] = ipython_useSVG
  res = fn(*args, **kwargs)
  if kwargs['useSVG']:
    return SVG(res)
  else:
    sio = BytesIO(res)
    return Image.open(sio)


def _DrawBits(fn, *args, **kwargs):
  if 'useSVG' not in kwargs:
    kwargs['useSVG'] = ipython_useSVG
  res = fn(*args, **kwargs)
  if kwargs['useSVG']:
    return SVG(res)
  else:
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
