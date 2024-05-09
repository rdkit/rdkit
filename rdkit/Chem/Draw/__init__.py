#
# Copyright (C) 2006-2021 Greg Landrum
#  All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os
import warnings
from collections import namedtuple
from importlib.util import find_spec
from io import BytesIO

import numpy
from rdkit import Chem
from rdkit import RDConfig
from rdkit import rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.rdMolDraw2D import *


def _sip_available():
  try:
    from rdkit.Chem.Draw.rdMolDraw2DQt import rdkitQtVersion
  except ImportError:
    return False
  pyqt_pkg = f'PyQt{rdkitQtVersion[0]}'
  if find_spec(pyqt_pkg) and find_spec(f'{pyqt_pkg}.sip'):
    return True
  elif find_spec('sip'):
    return True
  return False


if _sip_available():

  def MolDraw2DFromQPainter(qpainter, width=-1, height=-1, panelWidth=-1, panelHeight=-1):
    from rdkit.Chem.Draw import rdMolDraw2DQt
    if rdMolDraw2DQt.rdkitQtVersion.startswith('6'):
      from PyQt6.QtGui import QPainter
    else:
      from PyQt5.Qt import QPainter
    try:
      # Prefer the PyQt-bundled sip
      if rdMolDraw2DQt.rdkitQtVersion.startswith('6'):
        from PyQt6 import sip
      else:
        from PyQt5 import sip
    except ImportError:
      # No bundled sip, try the standalone package
      import sip

    if not isinstance(qpainter, QPainter):
      raise ValueError("argument must be a QPainter instance")
    if width <= 0:
      width = qpainter.viewport().width()
    if height <= 0:
      height = qpainter.viewport().height()
    ptr = sip.unwrapinstance(qpainter)
    d2d = rdMolDraw2DQt.MolDraw2DFromQPainter_(width, height, ptr, panelWidth, panelWidth)
    # tie the lifetime of the QPainter to this MolDraw2D object
    d2d._qptr = qpainter
    return d2d


def MolToImage(mol, size=(300, 300), kekulize=True, wedgeBonds=True, fitImage=False, options=None,
               **kwargs):
  """Returns a PIL image containing a drawing of the molecule

      ARGUMENTS:

        - kekulize: run kekulization routine on input `mol` (default True)

        - size: final image size, in pixel (default (300,300))

        - wedgeBonds: draw wedge (stereo) bonds (default True)

        - highlightAtoms: list of atoms to highlight (default [])

        - highlightBonds: list of bonds to highlight (default [])

        - highlightColor: RGB color as tuple (default [1, 0, 0])

      NOTE:

            use 'matplotlib.colors.to_rgb()' to convert string and
            HTML color codes into the RGB tuple representation, eg.

              from matplotlib.colors import ColorConverter
              img = Draw.MolToImage(m, highlightAtoms=[1,2], highlightColor=ColorConverter().to_rgb('aqua'))
              img.save("molecule.png")

      RETURNS:

        a PIL Image object
  """
  if not mol:
    raise ValueError('Null molecule provided')
  if not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    raise RuntimeError("MolToImage requires that the RDKit be built with Cairo support")
  return _moltoimg(mol, size, kwargs.get('highlightAtoms', []), kwargs.get('legend', ''),
                   highlightBonds=kwargs.get('highlightBonds',
                                             []), drawOptions=options, kekulize=kekulize,
                   wedgeBonds=wedgeBonds, highlightColor=kwargs.get('highlightColor', None))


def MolToFile(mol, filename, size=(300, 300), kekulize=True, wedgeBonds=True, imageType=None,
              fitImage=False, options=None, **kwargs):
  """ Generates a drawing of a molecule and writes it to a file
  """
  # original contribution from Uwe Hoffmann
  if not filename:
    raise ValueError('no filename provided')
  if not mol:
    raise ValueError('Null molecule provided')

  if imageType is None:
    imageType = os.path.splitext(filename)[1][1:]

  if imageType not in ('svg', 'png'):
    raise ValueError('only PNG and SVG image types are supported')

  if imageType == 'png':
    drawfn = _moltoimg
    mode = 'b'
  elif imageType == 'svg':
    drawfn = _moltoSVG
    mode = 't'
  else:
    raise ValueError("unsupported output format")
  data = drawfn(mol, size, kwargs.get('highlightAtoms', []), kwargs.get('legend', ''),
                highlightBonds=kwargs.get('highlightBonds', []), drawOptions=options,
                kekulize=kekulize, wedgeBonds=wedgeBonds, returnPNG=True)
  with open(filename, 'w+' + mode) as outf:
    outf.write(data)
    outf.close()


def ShowMol(mol, size=(300, 300), kekulize=True, wedgeBonds=True, title='RDKit Molecule',
            stayInFront=True, **kwargs):
  """ Generates a picture of a molecule and displays it in a Tkinter window
  """
  import tkinter

  from PIL import ImageTk

  img = MolToImage(mol, size, kekulize, wedgeBonds, **kwargs)

  tkRoot = tkinter.Tk()
  tkRoot.title(title)
  tkPI = ImageTk.PhotoImage(img)
  tkLabel = tkinter.Label(tkRoot, image=tkPI)
  tkLabel.place(x=0, y=0, width=img.size[0], height=img.size[1])
  tkRoot.geometry('%dx%d' % (img.size))
  tkRoot.lift()
  if stayInFront:
    tkRoot.attributes('-topmost', True)
  tkRoot.mainloop()


def _bivariate_normal(X, Y, sigmax=1.0, sigmay=1.0, mux=0.0, muy=0.0, sigmaxy=0.0):
  """

    This is the implementation from matplotlib:
    https://github.com/matplotlib/matplotlib/blob/81e8154dbba54ac1607b21b22984cabf7a6598fa/lib/matplotlib/mlab.py#L1866
    it was deprecated in v2.2 of matplotlib, so we are including it here.


    Bivariate Gaussian distribution for equal shape *X*, *Y*.
    See `bivariate normal
    <http://mathworld.wolfram.com/BivariateNormalDistribution.html>`_
    at mathworld.
    """
  Xmu = X - mux
  Ymu = Y - muy

  rho = sigmaxy / (sigmax * sigmay)
  z = Xmu**2 / sigmax**2 + Ymu**2 / sigmay**2 - 2 * rho * Xmu * Ymu / (sigmax * sigmay)
  denom = 2 * numpy.pi * sigmax * sigmay * numpy.sqrt(1 - rho**2)
  return numpy.exp(-z / (2 * (1 - rho**2))) / denom


def calcAtomGaussians(mol, a=0.03, step=0.02, weights=None):
  """
useful things to do with these:
fig.axes[0].imshow(z,cmap=cm.gray,interpolation='bilinear',origin='lower',extent=(0,1,0,1))
fig.axes[0].contour(x,y,z,20,colors='k')

fig=Draw.MolToMPL(m);
contribs=Crippen.rdMolDescriptors._CalcCrippenContribs(m)
logps,mrs=zip(*contribs)
x,y,z=Draw.calcAtomGaussians(m,0.03,step=0.01,weights=logps)
fig.axes[0].imshow(z,cmap=cm.jet,interpolation='bilinear',origin='lower',extent=(0,1,0,1))
fig.axes[0].contour(x,y,z,20,colors='k',alpha=0.5)
fig.savefig('coumlogps.colored.png',bbox_inches='tight')


  """
  x = numpy.arange(0, 1, step)
  y = numpy.arange(0, 1, step)
  X, Y = numpy.meshgrid(x, y)
  if weights is None:
    weights = [1.] * mol.GetNumAtoms()
  Z = _bivariate_normal(X, Y, a, a, mol._atomPs[0][0], mol._atomPs[0][1]) * weights[0]
  for i in range(1, mol.GetNumAtoms()):
    Zp = _bivariate_normal(X, Y, a, a, mol._atomPs[i][0], mol._atomPs[i][1])
    Z += Zp * weights[i]
  return X, Y, Z


def MolsToImage(mols, subImgSize=(200, 200), legends=None, **kwargs):
  """
  """
  from PIL import Image
  if legends is None:
    legends = [None] * len(mols)
  res = Image.new("RGBA", (subImgSize[0] * len(mols), subImgSize[1]))
  for i, mol in enumerate(mols):
    res.paste(MolToImage(mol, subImgSize, legend=legends[i], **kwargs), (i * subImgSize[0], 0))
  return res


def _drawerToImage(d2d):
  from PIL import Image
  sio = BytesIO(d2d.GetDrawingText())
  return Image.open(sio)


def shouldKekulize(mol, kekulize):
  if kekulize:
    return not any(bond.GetIsAromatic() and bond.HasQuery() for bond in mol.GetBonds())
  return kekulize


def _moltoimg(mol, sz, highlights, legend, returnPNG=False, drawOptions=None, **kwargs):
  try:
    with rdBase.BlockLogs():
      mol.GetAtomWithIdx(0).GetExplicitValence()
  except RuntimeError:
    mol.UpdatePropertyCache(False)

  kekulize = shouldKekulize(mol, kwargs.get('kekulize', True))
  wedge = kwargs.get('wedgeBonds', True)

  if not drawOptions or drawOptions.prepareMolsBeforeDrawing:
    try:
      with rdBase.BlockLogs():
        mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize, wedgeBonds=wedge)
    except ValueError:  # <- can happen on a kekulization failure
      mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False, wedgeBonds=wedge)
  d2d = rdMolDraw2D.MolDraw2DCairo(sz[0], sz[1])
  if drawOptions is not None:
    d2d.SetDrawOptions(drawOptions)
  if 'highlightColor' in kwargs and kwargs['highlightColor']:
    d2d.drawOptions().setHighlightColour(kwargs['highlightColor'])
  # we already prepared the molecule:
  d2d.drawOptions().prepareMolsBeforeDrawing = False
  bondHighlights = kwargs.get('highlightBonds', None)
  if bondHighlights is not None:
    d2d.DrawMolecule(mol, legend=legend or "", highlightAtoms=highlights or [],
                     highlightBonds=bondHighlights)
  else:
    d2d.DrawMolecule(mol, legend=legend or "", highlightAtoms=highlights or [])
  d2d.FinishDrawing()
  if returnPNG:
    img = d2d.GetDrawingText()
  else:
    img = _drawerToImage(d2d)
  return img


def _moltoSVG(mol, sz, highlights, legend, kekulize, drawOptions=None, **kwargs):
  try:
    with rdBase.BlockLogs():
      mol.GetAtomWithIdx(0).GetExplicitValence()
  except RuntimeError:
    mol.UpdatePropertyCache(False)

  kekulize = shouldKekulize(mol, kekulize)

  if not drawOptions or drawOptions.prepareMolsBeforeDrawing:
    try:
      with rdBase.BlockLogs():
        mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize)
    except ValueError:  # <- can happen on a kekulization failure
      mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
  d2d = rdMolDraw2D.MolDraw2DSVG(sz[0], sz[1])
  if drawOptions is not None:
    d2d.SetDrawOptions(drawOptions)
  # we already prepared the molecule:
  d2d.drawOptions().prepareMolsBeforeDrawing = False
  bondHighlights = kwargs.get('highlightBonds', None)
  if bondHighlights is not None:
    d2d.DrawMolecule(mol, legend=legend or "", highlightAtoms=highlights or [],
                     highlightBonds=bondHighlights)
  else:
    d2d.DrawMolecule(mol, legend=legend or "", highlightAtoms=highlights or [])
  d2d.FinishDrawing()
  svg = d2d.GetDrawingText()
  return svg


def _MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=None,
                     highlightAtomLists=None, highlightBondLists=None, drawOptions=None,
                     returnPNG=False, **kwargs):
  """ returns a PIL Image of the grid
  """
  if legends is None:
    legends = [''] * len(mols)

  nRows = len(mols) // molsPerRow
  if len(mols) % molsPerRow:
    nRows += 1

  if not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    raise RuntimeError("MolsToGridImage requires that the RDKit be built with Cairo support")
  fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
  d2d = rdMolDraw2D.MolDraw2DCairo(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])
  if drawOptions is not None:
    d2d.SetDrawOptions(drawOptions)
  else:
    dops = d2d.drawOptions()
    for k, v in list(kwargs.items()):
      if hasattr(dops, k):
        setattr(dops, k, v)
        del kwargs[k]
  d2d.DrawMolecules(list(mols), legends=legends or None, highlightAtoms=highlightAtomLists,
                    highlightBonds=highlightBondLists, **kwargs)
  d2d.FinishDrawing()
  if not returnPNG:
    res = _drawerToImage(d2d)
  else:
    res = d2d.GetDrawingText()

  return res


def _padList(inputList, lengthShouldBe, padWith=""):
  length = len(inputList)
  paddingCount = lengthShouldBe - length
  paddedList = inputList + [padWith] * paddingCount
  return paddedList


def _padMatrix(inputMatrix, rowLength, padWith=""):
  paddedMatrix = [_padList(row, rowLength, padWith) for row in inputMatrix]
  return paddedMatrix


def _flattenTwoDList(twoDList):
  return [item for sublist in twoDList for item in sublist]


def _MolsNestedToLinear(molsMatrix, legendsMatrix=None, highlightAtomListsMatrix=None,
                        highlightBondListsMatrix=None):
  """Converts a nested data structure (where each data substructure represents a row in mol grid image)
  to a linear one, padding rows as needed so all rows are the length of the longest row
  """
  # Check that other matrices (if provided) are same length,
  #   and each element (sub-iterable) is the same length, as molsMatrix

  nMolsRows = len(molsMatrix)

  if legendsMatrix is not None:
    nLegendsRows = len(legendsMatrix)
    if nLegendsRows != nMolsRows:
      err = f"If legendsMatrix is provided it must be the same length (have the same number "
      err += f"of sub-iterables) as molsMatrix, {nMolsRows}; its length is {nLegendsRows}."
      raise ValueError(err)
    for rowIndex, row in enumerate(legendsMatrix):
      if len(row) != len(molsMatrix[rowIndex]):
        err = f"If legendsMatrix is provided each of its sub-iterables must be the same length "
        err += f"as the corresponding sub-iterable of molsMatrix. For sub-iterable of index "
        err += f"{rowIndex}, its length in molsMatrix is {len(molsMatrix[rowIndex])} "
        err += f"while its length in legendsMatrix is {len(row)}."
        raise ValueError(err)

  if highlightAtomListsMatrix is not None:
    nHighlightAtomListsRows = len(highlightAtomListsMatrix)
    if nHighlightAtomListsRows != nMolsRows:
      err = f"If highlightAtomListsMatrix is provided it must be the same length (have the same number "
      err += f"of sub-iterables) as molsMatrix, {nMolsRows}; its length is {nHighlightAtomListsRows}."
      raise ValueError(err)
    for rowIndex, row in enumerate(highlightAtomListsMatrix):
      if len(row) != len(molsMatrix[rowIndex]):
        err = f"If highlightAtomListsMatrix is provided each of its sub-iterables must be the same length "
        err += f"as the corresponding sub-iterable of molsMatrix. For sub-iterable of index "
        err += f"{rowIndex}, its length in molsMatrix is {len(molsMatrix[rowIndex])} "
        err += f"while its length in highlightAtomListsMatrix is {len(row)}."
        raise ValueError(err)

  if highlightBondListsMatrix is not None:
    nHighlightBondListsRows = len(highlightBondListsMatrix)
    if nHighlightBondListsRows != nMolsRows:
      err = f"If highlightBondListsMatrix is provided it must be the same length (have the same number "
      err += f"of sub-iterables) as molsMatrix, {nMolsRows}; its length is {nHighlightBondListsRows}."
      raise ValueError(err)
    for rowIndex, row in enumerate(highlightBondListsMatrix):
      if len(row) != len(molsMatrix[rowIndex]):
        err = f"If highlightBondListsMatrix is provided each of its sub-iterables must be the same length "
        err += f"as the corresponding sub-iterable of molsMatrix. For sub-iterable of index "
        err += f"{rowIndex}, its length in molsMatrix is {len(molsMatrix[rowIndex])} "
        err += f"while its length in highlightBondListsMatrix is {len(row)}."
        raise ValueError(err)

  molsPerRow = max(len(row) for row in molsMatrix)

  # Pad matrices so they're rectangular (same length for each sublist),
  #   then convert to 1D lists
  # Pad using None for molecule for empty cells
  molsMatrixPadded = _padMatrix(molsMatrix, molsPerRow, None)
  mols = _flattenTwoDList(molsMatrixPadded)

  if legendsMatrix is not None:
    legendsMatrixPadded = _padMatrix(legendsMatrix, molsPerRow, "")
    legends = _flattenTwoDList(legendsMatrixPadded)
  else:
    legends = None

  if highlightAtomListsMatrix is not None:
    highlightAtomListsPadded = _padMatrix(highlightAtomListsMatrix, molsPerRow, [])
    highlightAtomLists = _flattenTwoDList(highlightAtomListsPadded)
  else:
    highlightAtomLists = None

  if highlightBondListsMatrix is not None:
    highlightBondListsPadded = _padMatrix(highlightBondListsMatrix, molsPerRow, [])
    highlightBondLists = _flattenTwoDList(highlightBondListsPadded)
  else:
    highlightBondLists = None

  return mols, molsPerRow, legends, highlightAtomLists, highlightBondLists


def MolsMatrixToGridImage(molsMatrix, subImgSize=(200, 200), legendsMatrix=None,
                          highlightAtomListsMatrix=None, highlightBondListsMatrix=None,
                          useSVG=False, returnPNG=False, **kwargs):
  r"""Creates a mol grid image from a nested data structure (where each data substructure represents a row),
  padding rows as needed so all rows are the length of the longest row
          ARGUMENTS:

        - molsMatrix: A two-deep nested data structure of RDKit molecules to draw,
         iterable of iterables (for example list of lists) of RDKit molecules

        - subImgSize: The size of a cell in the drawing; passed through to MolsToGridImage (default (200, 200))

        - legendsMatrix: A two-deep nested data structure of strings to label molecules with,
         iterable of iterables (for example list of lists) of strings (default None)

        - highlightAtomListsMatrix: A three-deep nested data structure of integers of atoms to highlight,
         iterable of iterables (for example list of lists) of integers (default None)

        - highlightBondListsMatrix: A three-deep nested data structure of integers of bonds to highlight,
         iterable of iterables (for example list of lists) of integers (default None)

        - useSVG: Whether to return an SVG (if true) or PNG (if false);
         passed through to MolsToGridImage (default false)

        - returnPNG: Whether to return PNG data (if true) or a PIL object for a PNG image file (if false);
         has no effect if useSVG is true; passed through to MolsToGridImage (default false)

        - kwargs: Any other keyword arguments are passed to MolsToGridImage

      NOTES:

            To include a blank cell in the middle of a row, supply None for that entry in molsMatrix.
            You do not need to do that for empty cells at the end of a row; 
            this function will automatically pad rows so that all rows are the same length.
            
            This function is useful when each row has some meaning,
            for example the generation in a mass spectrometry fragmentation tree--refer to 
            example at https://en.wikipedia.org/wiki/Fragmentation_(mass_spectrometry).
            If you want to display a set molecules where each row does not have any specific meaning,
            use MolsToGridImage instead.

            This function nests data structures one additional level beyond the analogous function MolsToGridImage
            (in which the molecules and legends are non-nested lists, 
            and the highlight parameters are two-deep nested lists) 

      RETURNS:

        A grid of molecular images in one of these formats:
        
        - useSVG=False and returnPNG=False (default): A PIL object for a PNG image file

        - useSVG=False and returnPNG=True: PNG data

        - useSVG=True: An SVG string

      EXAMPLES:

        from rdkit import Chem
        from rdkit.Chem.Draw import MolsMatrixToGridImage, rdMolDraw2D
        FCl = Chem.MolFromSmiles("FCl")
        molsMatrix = [[FCl, FCl], [FCl, None, FCl]]

        # Minimal example: Only molsMatrix is supplied,
        # result will be a drawing containing (where each row contains molecules):
        # F-Cl    F-Cl
        # F-Cl            F-Cl
        img = MolsMatrixToGridImage(molsMatrix)
        img.save("MolsMatrixToGridImageMinimal.png")
        # img is a PIL object for a PNG image file like:
        # <PIL.PngImagePlugin.PngImageFile image mode=RGB size=600x200 at 0x1648CC390>
        # Drawing will be saved as PNG file MolsMatrixToGridImageMinimal.png

        # Exhaustive example: All parameters are supplied,
        # result will be a drawing containing (where each row of molecules is followed by a row of legends):
        # 1 F-Cl 0              1 F-Cl 0
        # no highlighting       bond highlighted         
        # 1 F-Cl 0                                  1 F-Cl 0
        # sodium highlighted                        chloride and bond highlighted
        legendsMatrix = [["no highlighting", "bond highlighted"], 
        ["F highlighted", "", "Cl and bond highlighted"]]
        highlightAtomListsMatrix = [[[],[]], [[0], None, [1]]]
        highlightBondListsMatrix = [[[],[0]], [[], None, [0]]]

        dopts = rdMolDraw2D.MolDrawOptions()
        dopts.addAtomIndices = True

        img_binary = MolsMatrixToGridImage(molsMatrix=molsMatrix, subImgSize=(300, 400), 
        legendsMatrix=legendsMatrix, highlightAtomListsMatrix=highlightAtomListsMatrix, 
        highlightBondListsMatrix=highlightBondListsMatrix, useSVG=False, returnPNG=True, drawOptions=dopts)
        print(img_binary[:20])
        # Prints a binary string: b'\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x03\x84'
  """
  mols, molsPerRow, legends, highlightAtomLists, highlightBondLists = _MolsNestedToLinear(
    molsMatrix, legendsMatrix, highlightAtomListsMatrix, highlightBondListsMatrix)

  return MolsToGridImage(mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=legends,
                         highlightAtomLists=highlightAtomLists, useSVG=useSVG, returnPNG=returnPNG,
                         highlightBondLists=highlightBondLists, **kwargs)


def _MolsToGridSVG(mols, molsPerRow=3, subImgSize=(200, 200), legends=None, highlightAtomLists=None,
                   highlightBondLists=None, drawOptions=None, **kwargs):
  """ returns an SVG of the grid
  """
  if legends is None:
    legends = [''] * len(mols)

  nRows = len(mols) // molsPerRow
  if len(mols) % molsPerRow:
    nRows += 1

  blocks = [''] * (nRows * molsPerRow)

  fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])

  d2d = rdMolDraw2D.MolDraw2DSVG(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])
  if drawOptions is not None:
    d2d.SetDrawOptions(drawOptions)
  else:
    dops = d2d.drawOptions()
    for k, v in list(kwargs.items()):
      if hasattr(dops, k):
        setattr(dops, k, v)
        del kwargs[k]
  d2d.DrawMolecules(list(mols), legends=legends or None, highlightAtoms=highlightAtomLists or [],
                    highlightBonds=highlightBondLists or [], **kwargs)
  d2d.FinishDrawing()
  res = d2d.GetDrawingText()
  return res


def MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=None,
                    highlightAtomLists=None, highlightBondLists=None, useSVG=False, returnPNG=False,
                    **kwargs):
  if legends and len(legends) > len(mols):
    legends = legends[:len(mols)]
  if highlightAtomLists and len(highlightAtomLists) > len(mols):
    highlightAtomLists = highlightAtomLists[:len(mols)]
  if highlightBondLists and len(highlightBondLists) > len(mols):
    highlightBondLists = highlightBondLists[:len(mols)]

  if useSVG:
    return _MolsToGridSVG(mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=legends,
                          highlightAtomLists=highlightAtomLists,
                          highlightBondLists=highlightBondLists, **kwargs)
  else:
    return _MolsToGridImage(mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=legends,
                            highlightAtomLists=highlightAtomLists,
                            highlightBondLists=highlightBondLists, returnPNG=returnPNG, **kwargs)


def ReactionToImage(rxn, subImgSize=(200, 200), useSVG=False, drawOptions=None, returnPNG=False,
                    **kwargs):
  if not useSVG and not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    raise RuntimeError("PNG output requires that the RDKit be built with Cairo support")
  width = subImgSize[0] * (rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates() + 1)
  if useSVG:
    d = rdMolDraw2D.MolDraw2DSVG(width, subImgSize[1])
  else:
    d = rdMolDraw2D.MolDraw2DCairo(width, subImgSize[1])
  if drawOptions is not None:
    d.SetDrawOptions(drawOptions)
  d.DrawReaction(rxn, **kwargs)
  d.FinishDrawing()
  if useSVG or returnPNG:
    return d.GetDrawingText()
  else:
    return _drawerToImage(d)


def DebugDraw(mol, size=(350, 350), drawer=None, asSVG=True, useBW=True, includeHLabels=True,
              addAtomIndices=True, addBondIndices=False):
  if drawer is None:
    if asSVG:
      drawer = rdMolDraw2D.MolDraw2DSVG(size[0], size[1])
    else:
      drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
  if useBW:
    drawer.drawOptions().useBWAtomPalette()

  drawer.drawOptions().addAtomIndices = addAtomIndices
  drawer.drawOptions().addBondIndices = addBondIndices

  drawer.drawOptions().annotationFontScale = 0.75

  if includeHLabels:
    for atom in mol.GetAtoms():
      if atom.GetTotalNumHs():
        atom.SetProp('atomNote', f'H{atom.GetTotalNumHs()}')

  aromAtoms = [x.GetIdx() for x in mol.GetAtoms() if x.GetIsAromatic()]
  clrs = {x: (.9, .9, .2) for x in aromAtoms}
  aromBonds = [x.GetIdx() for x in mol.GetBonds() if x.GetIsAromatic()]

  rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False, addChiralHs=False)
  drawer.drawOptions().prepareMolsBeforeDrawing = False

  drawer.DrawMolecule(mol, highlightAtoms=aromAtoms, highlightAtomColors=clrs,
                      highlightBonds=aromBonds)

  drawer.FinishDrawing()
  return drawer.GetDrawingText()


def DrawMorganBit(mol, bitId, bitInfo, whichExample=0, **kwargs):
  atomId, radius = bitInfo[bitId][whichExample]
  return DrawMorganEnv(mol, atomId, radius, **kwargs)


def DrawMorganBits(tpls, **kwargs):
  envs = []
  for tpl in tpls:
    if len(tpl) == 4:
      mol, bitId, bitInfo, whichExample = tpl
    else:
      mol, bitId, bitInfo = tpl
      whichExample = 0

    atomId, radius = bitInfo[bitId][whichExample]
    envs.append((mol, atomId, radius))
  return DrawMorganEnvs(envs, **kwargs)


# adapted from the function drawFPBits._drawFPBit() from the CheTo package
# original author Nadine Schneider
FingerprintEnv = namedtuple(
  'FingerprintEnv',
  ('submol', 'highlightAtoms', 'atomColors', 'highlightBonds', 'bondColors', 'highlightRadii'))


def _getMorganEnv(mol, atomId, radius, baseRad, aromaticColor, ringColor, centerColor, extraColor,
                  **kwargs):
  if not mol.GetNumConformers():
    rdDepictor.Compute2DCoords(mol)
  bitPath = Chem.FindAtomEnvironmentOfRadiusN(mol, radius, atomId)

  # get the atoms for highlighting
  atomsToUse = set((atomId, ))
  for b in bitPath:
    atomsToUse.add(mol.GetBondWithIdx(b).GetBeginAtomIdx())
    atomsToUse.add(mol.GetBondWithIdx(b).GetEndAtomIdx())

  #  enlarge the environment by one further bond
  enlargedEnv = set()
  for atom in atomsToUse:
    a = mol.GetAtomWithIdx(atom)
    for b in a.GetBonds():
      bidx = b.GetIdx()
      if bidx not in bitPath:
        enlargedEnv.add(bidx)
  enlargedEnv = list(enlargedEnv)
  enlargedEnv += bitPath

  # set the coordinates of the submol based on the coordinates of the original molecule
  amap = {}
  if enlargedEnv:
    submol = Chem.PathToSubmol(mol, enlargedEnv, atomMap=amap)
  else:
    # generate submol from fragments with no bonds
    submol = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, atomsToUse=atomsToUse))
  Chem.FastFindRings(submol)
  conf = Chem.Conformer(submol.GetNumAtoms())
  confOri = mol.GetConformer(0)
  for i1, i2 in amap.items():
    conf.SetAtomPosition(i2, confOri.GetAtomPosition(i1))
  submol.AddConformer(conf)
  envSubmol = []
  for i1, i2 in amap.items():
    for b in bitPath:
      beginAtom = amap[mol.GetBondWithIdx(b).GetBeginAtomIdx()]
      endAtom = amap[mol.GetBondWithIdx(b).GetEndAtomIdx()]
      envSubmol.append(submol.GetBondBetweenAtoms(beginAtom, endAtom).GetIdx())

  # color all atoms of the submol in gray which are not part of the bit
  # highlight atoms which are in rings
  atomcolors, bondcolors = {}, {}
  highlightAtoms, highlightBonds = [], []
  highlightRadii = {}
  for aidx in amap.keys():
    if aidx in atomsToUse:
      color = None
      if centerColor and aidx == atomId:
        color = centerColor
      elif aromaticColor and mol.GetAtomWithIdx(aidx).GetIsAromatic():
        color = aromaticColor
      elif ringColor and mol.GetAtomWithIdx(aidx).IsInRing():
        color = ringColor
      if color is not None:
        atomcolors[amap[aidx]] = color
        highlightAtoms.append(amap[aidx])
        highlightRadii[amap[aidx]] = baseRad
    else:
      #drawopt.atomLabels[amap[aidx]] = '*'
      submol.GetAtomWithIdx(amap[aidx]).SetAtomicNum(0)
      submol.GetAtomWithIdx(amap[aidx]).UpdatePropertyCache()
  color = extraColor
  for bid in submol.GetBonds():
    bidx = bid.GetIdx()
    if bidx not in envSubmol:
      bondcolors[bidx] = color
      highlightBonds.append(bidx)
  return FingerprintEnv(submol, highlightAtoms, atomcolors, highlightBonds, bondcolors,
                        highlightRadii)


def DrawMorganEnvs(envs, molsPerRow=3, subImgSize=(150, 150), baseRad=0.3, useSVG=True,
                   aromaticColor=(0.9, 0.9, 0.2), ringColor=(0.8, 0.8, 0.8),
                   centerColor=(0.6, 0.6, 0.9), extraColor=(0.9, 0.9, 0.9), legends=None,
                   drawOptions=None, **kwargs):
  submols = []
  highlightAtoms = []
  atomColors = []
  highlightBonds = []
  bondColors = []
  highlightRadii = []
  for mol, atomId, radius in envs:
    menv = _getMorganEnv(mol, atomId, radius, baseRad, aromaticColor, ringColor, centerColor,
                         extraColor, **kwargs)
    submols.append(menv.submol)
    highlightAtoms.append(menv.highlightAtoms)
    atomColors.append(menv.atomColors)
    highlightBonds.append(menv.highlightBonds)
    bondColors.append(menv.bondColors)
    highlightRadii.append(menv.highlightRadii)

  if legends is None:
    legends = [''] * len(envs)

  nRows = len(envs) // molsPerRow
  if len(envs) % molsPerRow:
    nRows += 1

  fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
  # Drawing
  if useSVG:
    drawer = rdMolDraw2D.MolDraw2DSVG(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])
  else:
    drawer = rdMolDraw2D.MolDraw2DCairo(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])

  if drawOptions is None:
    drawOptions = drawer.drawOptions()
  drawOptions.continuousHighlight = False
  drawOptions.includeMetadata = False
  drawer.SetDrawOptions(drawOptions)
  drawer.DrawMolecules(submols, legends=legends, highlightAtoms=highlightAtoms,
                       highlightAtomColors=atomColors, highlightBonds=highlightBonds,
                       highlightBondColors=bondColors, highlightAtomRadii=highlightRadii, **kwargs)
  drawer.FinishDrawing()
  return drawer.GetDrawingText()


def DrawMorganEnv(mol, atomId, radius, molSize=(150, 150), baseRad=0.3, useSVG=True,
                  aromaticColor=(0.9, 0.9, 0.2), ringColor=(0.8, 0.8, 0.8),
                  centerColor=(0.6, 0.6, 0.9), extraColor=(0.9, 0.9, 0.9), drawOptions=None,
                  **kwargs):
  menv = _getMorganEnv(mol, atomId, radius, baseRad, aromaticColor, ringColor, centerColor,
                       extraColor, **kwargs)

  # Drawing
  if useSVG:
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
  else:
    drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])

  if drawOptions is None:
    drawOptions = drawer.drawOptions()
  drawOptions.continuousHighlight = False
  drawOptions.includeMetadata = False
  drawer.SetDrawOptions(drawOptions)
  drawer.DrawMolecule(menv.submol, highlightAtoms=menv.highlightAtoms,
                      highlightAtomColors=menv.atomColors, highlightBonds=menv.highlightBonds,
                      highlightBondColors=menv.bondColors, highlightAtomRadii=menv.highlightRadii,
                      **kwargs)
  drawer.FinishDrawing()
  return drawer.GetDrawingText()


def DrawRDKitBits(tpls, **kwargs):
  envs = []
  for tpl in tpls:
    if len(tpl) == 4:
      mol, bitId, bitInfo, whichExample = tpl
    else:
      mol, bitId, bitInfo = tpl
      whichExample = 0

    bondpath = bitInfo[bitId][whichExample]
    envs.append((mol, bondpath))
  return DrawRDKitEnvs(envs, **kwargs)


def DrawRDKitBit(mol, bitId, bitInfo, whichExample=0, **kwargs):
  bondPath = bitInfo[bitId][whichExample]
  return DrawRDKitEnv(mol, bondPath, **kwargs)


def _getRDKitEnv(mol, bondPath, baseRad, aromaticColor, extraColor, nonAromaticColor, **kwargs):
  if not mol.GetNumConformers():
    rdDepictor.Compute2DCoords(mol)

  # get the atoms for highlighting
  atomsToUse = set()
  for b in bondPath:
    atomsToUse.add(mol.GetBondWithIdx(b).GetBeginAtomIdx())
    atomsToUse.add(mol.GetBondWithIdx(b).GetEndAtomIdx())

  # set the coordinates of the submol based on the coordinates of the original molecule
  amap = {}
  submol = Chem.PathToSubmol(mol, bondPath, atomMap=amap)
  Chem.FastFindRings(submol)
  conf = Chem.Conformer(submol.GetNumAtoms())
  confOri = mol.GetConformer(0)
  for i1, i2 in amap.items():
    conf.SetAtomPosition(i2, confOri.GetAtomPosition(i1))
  submol.AddConformer(conf)
  envSubmol = []
  for i1, i2 in amap.items():
    for b in bondPath:
      beginAtom = amap[mol.GetBondWithIdx(b).GetBeginAtomIdx()]
      endAtom = amap[mol.GetBondWithIdx(b).GetEndAtomIdx()]
      envSubmol.append(submol.GetBondBetweenAtoms(beginAtom, endAtom).GetIdx())

  # color all atoms of the submol in gray which are not part of the bit
  # highlight atoms which are in rings
  atomcolors, bondcolors = {}, {}
  highlightAtoms, highlightBonds = [], []
  highlightRadii = {}
  for aidx in amap.keys():
    if aidx in atomsToUse:
      color = None
      if aromaticColor and mol.GetAtomWithIdx(aidx).GetIsAromatic():
        color = aromaticColor
      elif nonAromaticColor and not mol.GetAtomWithIdx(aidx).GetIsAromatic():
        color = nonAromaticColor
      if color is not None:
        atomcolors[amap[aidx]] = color
        highlightAtoms.append(amap[aidx])
        highlightRadii[amap[aidx]] = baseRad
  color = extraColor
  for bid in submol.GetBonds():
    bidx = bid.GetIdx()
    if bidx not in envSubmol:
      bondcolors[bidx] = color
      highlightBonds.append(bidx)
  return FingerprintEnv(submol, highlightAtoms, atomcolors, highlightBonds, bondcolors,
                        highlightRadii)


def DrawRDKitEnvs(envs, molsPerRow=3, subImgSize=(150, 150), baseRad=0.3, useSVG=True,
                  aromaticColor=(0.9, 0.9, 0.2), extraColor=(0.9, 0.9, 0.9), nonAromaticColor=None,
                  legends=None, drawOptions=None, **kwargs):
  submols = []
  highlightAtoms = []
  atomColors = []
  highlightBonds = []
  bondColors = []
  highlightRadii = []
  for mol, bondpath in envs:
    menv = _getRDKitEnv(mol, bondpath, baseRad, aromaticColor, extraColor, nonAromaticColor,
                        **kwargs)
    submols.append(menv.submol)
    highlightAtoms.append(menv.highlightAtoms)
    atomColors.append(menv.atomColors)
    highlightBonds.append(menv.highlightBonds)
    bondColors.append(menv.bondColors)
    highlightRadii.append(menv.highlightRadii)

  if legends is None:
    legends = [''] * len(envs)

  nRows = len(envs) // molsPerRow
  if len(envs) % molsPerRow:
    nRows += 1

  fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
  # Drawing
  if useSVG:
    drawer = rdMolDraw2D.MolDraw2DSVG(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])
  else:
    drawer = rdMolDraw2D.MolDraw2DCairo(fullSize[0], fullSize[1], subImgSize[0], subImgSize[1])

  if drawOptions is None:
    drawOptions = drawer.drawOptions()
  drawOptions.continuousHighlight = False
  drawOptions.includeMetadata = False
  drawer.SetDrawOptions(drawOptions)
  drawer.DrawMolecules(submols, legends=legends, highlightAtoms=highlightAtoms,
                       highlightAtomColors=atomColors, highlightBonds=highlightBonds,
                       highlightBondColors=bondColors, highlightAtomRadii=highlightRadii, **kwargs)
  drawer.FinishDrawing()
  return drawer.GetDrawingText()


def DrawRDKitEnv(mol, bondPath, molSize=(150, 150), baseRad=0.3, useSVG=True,
                 aromaticColor=(0.9, 0.9, 0.2), extraColor=(0.9, 0.9, 0.9), nonAromaticColor=None,
                 drawOptions=None, **kwargs):
  menv = _getRDKitEnv(mol, bondPath, baseRad, aromaticColor, extraColor, nonAromaticColor, **kwargs)

  # Drawing
  if useSVG:
    drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0], molSize[1])
  else:
    drawer = rdMolDraw2D.MolDraw2DCairo(molSize[0], molSize[1])

  if drawOptions is None:
    drawOptions = drawer.drawOptions()
  drawOptions.continuousHighlight = False
  drawOptions.includeMetadata = False
  drawer.SetDrawOptions(drawOptions)
  drawer.DrawMolecule(menv.submol, highlightAtoms=menv.highlightAtoms,
                      highlightAtomColors=menv.atomColors, highlightBonds=menv.highlightBonds,
                      highlightBondColors=menv.bondColors, highlightAtomRadii=menv.highlightRadii,
                      **kwargs)
  drawer.FinishDrawing()
  return drawer.GetDrawingText()


def SetComicMode(opts):
  opts.fontFile = os.path.join(RDConfig.RDDataDir, "Fonts", "ComicNeue-Regular.ttf")
  opts.comicMode = True
