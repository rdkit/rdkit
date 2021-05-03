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
from importlib.util import find_spec

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem.Draw.rdMolDraw2D import *
from rdkit.Chem import rdDepictor
from rdkit import Chem, rdBase
from rdkit import RDConfig


def _getCanvas():
  useAGG = False
  useCairo = False
  useSping = False
  Canvas = None
  if not os.environ.get('RDKIT_CANVAS', ''):
    try:
      from rdkit.Chem.Draw.cairoCanvas import Canvas
      useCairo = True
    except ImportError:
      try:
        from rdkit.Chem.Draw.aggCanvas import Canvas
        useAGG = True
      except ImportError:
        from rdkit.Chem.Draw.spingCanvas import Canvas
        useSping = True
  else:
    canv = os.environ['RDKIT_CANVAS'].lower()
    if canv == 'cairo':
      from rdkit.Chem.Draw.cairoCanvas import Canvas
      useCairo = True
    elif canv == 'agg':
      from rdkit.Chem.Draw.aggCanvas import Canvas
      useAGG = True
    else:
      from rdkit.Chem.Draw.spingCanvas import Canvas
      useSping = True
  if useSping:
    # <- the sping canvas doesn't support unicode well
    DrawingOptions.radicalSymbol = '.'
  return useAGG, useCairo, Canvas


def _createCanvas(size):
  useAGG, useCairo, Canvas = _getCanvas()
  if useAGG or useCairo:
    from PIL import Image
    img = Image.new("RGBA", size, (0, 0, 0, 0))
    canvas = Canvas(img)
  else:
    from rdkit.Chem.Draw.spingCanvas import Canvas
    canvas = Canvas(size=size, name='MolToImageFile')
    img = canvas._image
  return img, canvas


def _legacyMolToImage(mol, size, kekulize, wedgeBonds, fitImage, options, canvas, **kwargs):
  """Returns a PIL image containing a drawing of the molecule using the legacy drawing code

      ARGUMENTS:

        - kekulize: run kekulization routine on input `mol` (default True)

        - size: final image size, in pixel (default (300,300))

        - wedgeBonds: draw wedge (stereo) bonds (default True)

        - highlightAtoms: list of atoms to highlight (default [])

        - highlightMap: dictionary of (atom, color) pairs (default None)

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
  if canvas is None:
    img, canvas = _createCanvas(size)
  else:
    img = None

  options = options or DrawingOptions()
  if fitImage:
    options.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds = wedgeBonds
  if 'highlightColor' in kwargs:
    color = kwargs.pop('highlightColor', (1, 0, 0))
    options.selectColor = color

  drawer = MolDrawing(canvas=canvas, drawingOptions=options)

  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)

  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)

  if 'legend' in kwargs:
    legend = kwargs['legend']
    del kwargs['legend']
  else:
    legend = ''

  drawer.AddMol(mol, **kwargs)

  if legend:
    from rdkit.Chem.Draw.MolDrawing import Font
    bbox = drawer.boundingBoxes[mol]
    pos = size[0] / 2, int(.94 * size[1]), 0  # the 0.94 is extremely empirical
    # canvas.addCanvasPolygon(((bbox[0],bbox[1]),(bbox[2],bbox[1]),(bbox[2],bbox[3]),(bbox[0],bbox[3])),
    #                         color=(1,0,0),fill=False,stroke=True)
    # canvas.addCanvasPolygon(((0,0),(0,size[1]),(size[0],size[1]),(size[0],0)   ),
    #                         color=(0,0,1),fill=False,stroke=True)
    font = Font(face='sans', size=12)
    canvas.addCanvasText(legend, pos, font)

  if kwargs.get('returnCanvas', False):
    return img, canvas, drawer
  else:
    canvas.flush()
    return img


def _sip_available():
  if find_spec('PyQt5') and find_spec('PyQt5.sip'):
    return True
  elif find_spec('sip'):
    return True
  return False


if find_spec('rdkit.Chem.Draw.rdMolDraw2DQt') and _sip_available():

  def MolDraw2DFromQPainter(qpainter, width=-1, height=-1, panelWidth=-1, panelHeight=-1):
    from PyQt5.Qt import QPainter
    try:
      # Prefer the PyQt5-bundled sip
      from PyQt5 import sip
    except ImportError:
      # No bundled sip, try the standalone package
      import sip
    from rdkit.Chem.Draw import rdMolDraw2DQt

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
               canvas=None, **kwargs):
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
  if canvas is not None or not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    return _legacyMolToImage(mol, size, kekulize, wedgeBonds, fitImage, options, canvas, **kwargs)
  if type(options) == DrawingOptions:
    warnings.warn(
      "legacy DrawingOptions not translated for new drawing code, please update manually",
      DeprecationWarning)
    options = None
  return _moltoimg(mol, size, kwargs.get('highlightAtoms', []), kwargs.get('legend', ''),
                   highlightBonds=kwargs.get('highlightBonds',
                                             []), drawOptions=options, kekulize=kekulize,
                   wedgeBonds=wedgeBonds, highlightColor=kwargs.get('highlightColor', None))


def _legacyMolToFile(mol, fileName, size, kekulize, wedgeBonds, imageType, fitImage, options,
                     **kwargs):
  """ Generates a drawing of a molecule and writes it to a file
  """
  if options is None:
    options = DrawingOptions()
  useAGG, useCairo, Canvas = _getCanvas()
  if fitImage:
    options.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds = wedgeBonds
  if useCairo or useAGG:
    canvas = Canvas(size=size, imageType=imageType, fileName=fileName)
  else:
    options.radicalSymbol = '.'  # <- the sping canvas doesn't support unicode well
    canvas = Canvas(size=size, name=fileName, imageType=imageType)
  drawer = MolDrawing(canvas=canvas, drawingOptions=options)
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)

  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)

  drawer.AddMol(mol, **kwargs)
  if useCairo or useAGG:
    canvas.flush()
  else:
    canvas.save()


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
    _legacyMolToFile(mol, filename, size, kekulize, wedgeBonds, imageType, fitImage, options,
                     **kwargs)

  if type(options) == DrawingOptions:
    warnings.warn(
      "legacy DrawingOptions not translated for new drawing code, please update manually",
      DeprecationWarning)
    options = None
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


def MolToImageFile(mol, filename, size=(300, 300), kekulize=True, wedgeBonds=True, **kwargs):
  """  DEPRECATED:  please use MolToFile instead

  """
  warnings.warn("MolToImageFile is deprecated, please use MolToFile instead", DeprecationWarning)
  img = MolToImage(mol, size=size, kekulize=kekulize, wedgeBonds=wedgeBonds, **kwargs)
  img.save(filename)


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


def MolToMPL(mol, size=(300, 300), kekulize=True, wedgeBonds=True, imageType=None, fitImage=False,
             options=None, **kwargs):
  """ Generates a drawing of a molecule on a matplotlib canvas
  """
  if not mol:
    raise ValueError('Null molecule provided')
  from rdkit.Chem.Draw.mplCanvas import Canvas
  canvas = Canvas(size)
  if options is None:
    options = DrawingOptions()
    options.bgColor = None
  if fitImage:
    options.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds = wedgeBonds
  drawer = MolDrawing(canvas=canvas, drawingOptions=options)
  omol = mol
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)

  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)

  drawer.AddMol(mol, **kwargs)
  omol._atomPs = drawer.atomPs[mol]
  for k, v in omol._atomPs.items():
    omol._atomPs[k] = canvas.rescalePt(v)
  canvas._figure.set_size_inches(float(size[0]) / 100, float(size[1]) / 100)
  return canvas._figure


import numpy


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


from io import BytesIO


def _drawerToImage(d2d):
  from PIL import Image
  sio = BytesIO(d2d.GetDrawingText())
  return Image.open(sio)


def _okToKekulizeMol(mol, kekulize):
  if kekulize:
    for bond in mol.GetBonds():
      if bond.GetIsAromatic() and bond.HasQuery():
        return False
    return True
  return kekulize


def _moltoimg(mol, sz, highlights, legend, returnPNG=False, drawOptions=None, **kwargs):
  try:
    blocker = rdBase.BlockLogs()
    mol.GetAtomWithIdx(0).GetExplicitValence()
  except RuntimeError:
    mol.UpdatePropertyCache(False)

  kekulize = _okToKekulizeMol(mol, kwargs.get('kekulize', True))
  wedge = kwargs.get('wedgeBonds', True)

  try:
    blocker = rdBase.BlockLogs()
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize, wedgeBonds=wedge)
  except ValueError:  # <- can happen on a kekulization failure
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False, wedgeBonds=wedge)
  if not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    img = MolToImage(mc, sz, legend=legend, highlightAtoms=highlights, **kwargs)
    if returnPNG:
      bio = BytesIO()
      img.save(bio, format='PNG')
      img = bio.getvalue()
  else:
    d2d = rdMolDraw2D.MolDraw2DCairo(sz[0], sz[1])
    if drawOptions is not None:
      d2d.SetDrawOptions(drawOptions)
    if 'highlightColor' in kwargs and kwargs['highlightColor']:
      d2d.drawOptions().setHighlightColour(kwargs['highlightColor'])
    # we already prepared the molecule:
    d2d.drawOptions().prepareMolsBeforeDrawing = False
    bondHighlights = kwargs.get('highlightBonds', None)
    if bondHighlights is not None:
      d2d.DrawMolecule(mc, legend=legend or "", highlightAtoms=highlights or [],
                       highlightBonds=bondHighlights)
    else:
      d2d.DrawMolecule(mc, legend=legend or "", highlightAtoms=highlights or [])
    d2d.FinishDrawing()
    if returnPNG:
      img = d2d.GetDrawingText()
    else:
      img = _drawerToImage(d2d)
  return img


def _moltoSVG(mol, sz, highlights, legend, kekulize, drawOptions=None, **kwargs):
  try:
    blocker = rdBase.BlockLogs()
    mol.GetAtomWithIdx(0).GetExplicitValence()
  except RuntimeError:
    mol.UpdatePropertyCache(False)

  kekulize = _okToKekulizeMol(mol, kekulize)

  try:
    blocker = rdBase.BlockLogs()
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize)
  except ValueError:  # <- can happen on a kekulization failure
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
  d2d = rdMolDraw2D.MolDraw2DSVG(sz[0], sz[1])
  if drawOptions is not None:
    d2d.SetDrawOptions(drawOptions)

  bondHighlights = kwargs.get('highlightBonds', None)
  if bondHighlights is not None:
    d2d.DrawMolecule(mc, legend=legend or "", highlightAtoms=highlights or [],
                     highlightBonds=bondHighlights)
  else:
    d2d.DrawMolecule(mc, legend=legend or "", highlightAtoms=highlights or [])
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
    from PIL import Image
    res = Image.new("RGBA", (molsPerRow * subImgSize[0], nRows * subImgSize[1]), (255, 255, 255, 0))
    for i, mol in enumerate(mols):
      row = i // molsPerRow
      col = i % molsPerRow
      highlights = None
      if highlightAtomLists and highlightAtomLists[i]:
        highlights = highlightAtomLists[i]
      if mol is not None:
        img = _moltoimg(mol, subImgSize, highlights, legends[i], **kwargs)
        res.paste(img, (col * subImgSize[0], row * subImgSize[1]))
  else:
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


def _legacyReactionToImage(rxn, subImgSize=(200, 200), **kwargs):
  from PIL import Image

  mols = []
  for i in range(rxn.GetNumReactantTemplates()):
    tmpl = rxn.GetReactantTemplate(i)
    tmpl.UpdatePropertyCache(False)
    mols.append(tmpl)
  mols.append(None)
  for i in range(rxn.GetNumProductTemplates()):
    tmpl = rxn.GetProductTemplate(i)
    tmpl.UpdatePropertyCache(False)
    mols.append(tmpl)

  res = Image.new("RGBA", (subImgSize[0] * len(mols), subImgSize[1]), (255, 255, 255, 0))
  for i, mol in enumerate(mols):
    if mol is not None:
      nimg = MolToImage(mol, subImgSize, kekulize=False, **kwargs)
    else:
      nimg, canvas = _createCanvas(subImgSize)
      p0 = (10, subImgSize[1] // 2)
      p1 = (subImgSize[0] - 10, subImgSize[1] // 2)
      p3 = (subImgSize[0] - 20, subImgSize[1] // 2 - 10)
      p4 = (subImgSize[0] - 20, subImgSize[1] // 2 + 10)
      canvas.addCanvasLine(p0, p1, lineWidth=2, color=(0, 0, 0))
      canvas.addCanvasLine(p3, p1, lineWidth=2, color=(0, 0, 0))
      canvas.addCanvasLine(p4, p1, lineWidth=2, color=(0, 0, 0))
      if hasattr(canvas, 'flush'):
        canvas.flush()
      else:
        canvas.save()
    res.paste(nimg, (i * subImgSize[0], 0))
  return res


def ReactionToImage(rxn, subImgSize=(200, 200), useSVG=False, drawOptions=None, returnPNG=False,
                    **kwargs):
  if not useSVG and not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    return _legacyReactionToImage(rxn, subImgSize=subImgSize, **kwargs)
  else:
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


def MolToQPixmap(mol, size=(300, 300), kekulize=True, wedgeBonds=True, fitImage=False, options=None,
                 **kwargs):
  """ Generates a drawing of a molecule on a Qt QPixmap
    """
  if not mol:
    raise ValueError('Null molecule provided')
  from rdkit.Chem.Draw.qtCanvas import Canvas
  canvas = Canvas(size)
  if options is None:
    options = DrawingOptions()
  options.bgColor = None
  if fitImage:
    options.dotsPerAngstrom = int(min(size) / 10)
  options.wedgeDashedBonds = wedgeBonds
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
  drawer = MolDrawing(canvas=canvas, drawingOptions=options)
  drawer.AddMol(mol, **kwargs)
  canvas.flush()
  return canvas.pixmap


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
from collections import namedtuple
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
  submol = Chem.PathToSubmol(mol, enlargedEnv, atomMap=amap)
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
