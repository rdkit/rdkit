#
# Copyright (C) 2006-2016 Greg Landrum
#  All Rights Reserved
#
import os
import re

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem.Draw.rdMolDraw2D import *
from rdkit.six import iteritems


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
    try:
      import Image
    except ImportError:
      from PIL import Image
    img = Image.new("RGBA", size, (0, 0, 0, 0))
    canvas = Canvas(img)
  else:
    from rdkit.Chem.Draw.spingCanvas import Canvas
    canvas = Canvas(size=size, name='MolToImageFile')
    img = canvas._image
  return img, canvas


def MolToImage(mol, size=(300, 300), kekulize=True, wedgeBonds=True, fitImage=False, options=None,
               canvas=None, **kwargs):
  """Returns a PIL image containing a drawing of the molecule

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


def MolToFile(mol, fileName, size=(300, 300), kekulize=True, wedgeBonds=True, imageType=None,
              fitImage=False, options=None, **kwargs):
  """ Generates a drawing of a molecule and writes it to a file
  """
  # original contribution from Uwe Hoffmann
  if not fileName:
    raise ValueError('no fileName provided')
  if not mol:
    raise ValueError('Null molecule provided')

  if imageType is None:
    imageType = os.path.splitext(fileName)[1][1:]

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


def MolToImageFile(mol, filename, size=(300, 300), kekulize=True, wedgeBonds=True, **kwargs):
  """  DEPRECATED:  please use MolToFile instead

  """
  img = MolToImage(mol, size=size, kekulize=kekulize, wedgeBonds=wedgeBonds, **kwargs)
  img.save(filename)


tkRoot = None
tkLabel = None
tkPI = None


def ShowMol(mol, size=(300, 300), kekulize=True, wedgeBonds=True, title='RDKit Molecule', **kwargs):
  """ Generates a picture of a molecule and displays it in a Tkinter window
  """
  global tkRoot, tkLabel, tkPI
  try:
    import Tkinter
  except ImportError:
    import tkinter as Tkinter
  try:
    import ImageTk
  except ImportError:
    from PIL import ImageTk

  img = MolToImage(mol, size, kekulize, wedgeBonds, **kwargs)

  if not tkRoot:
    tkRoot = Tkinter.Tk()
    tkRoot.title(title)
    tkPI = ImageTk.PhotoImage(img)
    tkLabel = Tkinter.Label(tkRoot, image=tkPI)
    tkLabel.place(x=0, y=0, width=img.size[0], height=img.size[1])
  else:
    tkPI.paste(img)
  tkRoot.geometry('%dx%d' % (img.size))


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
  for k, v in iteritems(omol._atomPs):
    omol._atomPs[k] = canvas.rescalePt(v)
  canvas._figure.set_size_inches(float(size[0]) / 100, float(size[1]) / 100)
  return canvas._figure


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
  import numpy
  from matplotlib import mlab
  x = numpy.arange(0, 1, step)
  y = numpy.arange(0, 1, step)
  X, Y = numpy.meshgrid(x, y)
  if weights is None:
    weights = [1.] * mol.GetNumAtoms()
  Z = mlab.bivariate_normal(X, Y, a, a, mol._atomPs[0][0], mol._atomPs[0][1]) * weights[0]
  for i in range(1, mol.GetNumAtoms()):
    Zp = mlab.bivariate_normal(X, Y, a, a, mol._atomPs[i][0], mol._atomPs[i][1])
    Z += Zp * weights[i]
  return X, Y, Z


def MolsToImage(mols, subImgSize=(200, 200), legends=None, **kwargs):
  """
  """
  try:
    import Image
  except ImportError:
    from PIL import Image
  if legends is None:
    legends = [None] * len(mols)
  res = Image.new("RGBA", (subImgSize[0] * len(mols), subImgSize[1]))
  for i, mol in enumerate(mols):
    res.paste(MolToImage(mol, subImgSize, legend=legends[i], **kwargs), (i * subImgSize[0], 0))
  return res

from io import BytesIO
def _drawerToImage(d2d):
  try:
    import Image
  except ImportError:
    from PIL import Image
  sio = BytesIO(d2d.GetDrawingText())
  return Image.open(sio)

def _okToKekulizeMol(mol,kekulize):
  if kekulize:
    for bond in mol.GetBonds():
      if bond.GetIsAromatic() and bond.HasQuery():
        return False
    return True
  return kekulize

def _moltoimg(mol, sz, highlights, legend, returnPNG=False, **kwargs):
  try:
    mol.GetAtomWithIdx(0).GetExplicitValence()
  except RuntimeError:
    mol.UpdatePropertyCache(False)

  kekulize=_okToKekulizeMol(mol,kwargs.get('kekulize', True))

  try:
    mc = rdMolDraw2D.PrepareMolForDrawing(mol,kekulize=kekulize)
  except ValueError:  # <- can happen on a kekulization failure
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
  if not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    img = MolToImage(mc, sz, legend=legend, highlightAtoms=highlights, **kwargs)
    if returnPNG:
      bio = BytesIO()
      img.save(bio, format='PNG')
      img = bio.getvalue()
  else:
    d2d = rdMolDraw2D.MolDraw2DCairo(sz[0], sz[1])
    d2d.DrawMolecule(mc, legend=legend, highlightAtoms=highlights)
    d2d.FinishDrawing()
    if returnPNG:
        img = d2d.GetDrawingText()
    else:
        img = _drawerToImage(d2d)
  return img

def _moltoSVG(mol, sz, highlights, legend, kekulize, **kwargs):
  try:
    mol.GetAtomWithIdx(0).GetExplicitValence()
  except RuntimeError:
    mol.UpdatePropertyCache(False)

  kekulize=_okToKekulizeMol(mol,kekulize)

  try:
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=kekulize)
  except ValueError:  # <- can happen on a kekulization failure
    mc = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
  d2d = rdMolDraw2D.MolDraw2DSVG(sz[0], sz[1])
  d2d.DrawMolecule(mc, legend=legend, highlightAtoms=highlights)
  d2d.FinishDrawing()
  svg = d2d.GetDrawingText()
  return svg.replace("svg:", "")



def _MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=None,
                     highlightAtomLists=None, highlightBondLists=None, **kwargs):
  """ returns a PIL Image of the grid
  """
  if legends is None:
    legends = [''] * len(mols)

  nRows = len(mols) // molsPerRow
  if len(mols) % molsPerRow:
    nRows += 1

  if not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
    try:
      import Image
    except ImportError:
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
    d2d = rdMolDraw2D.MolDraw2DCairo(fullSize[0],fullSize[1],subImgSize[0], subImgSize[1])
    d2d.DrawMolecules(list(mols),legends=legends,highlightAtoms=highlightAtomLists,
                      highlightBonds=highlightBondLists,**kwargs)
    d2d.FinishDrawing()
    res = _drawerToImage(d2d)

  return res


def _MolsToGridSVG(mols, molsPerRow=3, subImgSize=(200, 200), legends=None,
                   highlightAtomLists=None, highlightBondLists=None,
                   stripSVGNamespace=True, **kwargs):
  """ returns an SVG of the grid
  """
  if legends is None:
    legends = [''] * len(mols)

  nRows = len(mols) // molsPerRow
  if len(mols) % molsPerRow:
    nRows += 1

  blocks = [''] * (nRows * molsPerRow)

  fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])

  d2d = rdMolDraw2D.MolDraw2DSVG(fullSize[0],fullSize[1],subImgSize[0], subImgSize[1])
  d2d.DrawMolecules(mols,legends=legends,highlightAtoms=highlightAtomLists,
                    highlightBonds=highlightBondLists,**kwargs)
  d2d.FinishDrawing()
  res = d2d.GetDrawingText()
  if stripSVGNamespace:
    res = res.replace('svg:', '')
  return res


def MolsToGridImage(mols, molsPerRow=3, subImgSize=(200, 200), legends=None,
                    highlightAtomLists=None, highlightBondLists=None,
                    useSVG=False, **kwargs):
  if legends and len(legends)>len(mols):
    legends = legends[:len(mols)]
  if highlightAtomLists and len(highlightAtomLists)>len(mols):
    highlightAtomLists = highlightAtomLists[:len(mols)]
  if highlightBondLists and len(highlightBondLists)>len(mols):
    highlightBondLists = highlightBondLists[:len(mols)]

  if useSVG:
    return _MolsToGridSVG(mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=legends,
                          highlightAtomLists=highlightAtomLists,
                          highlightBondLists=highlightBondLists,
                          **kwargs)
  else:
    return _MolsToGridImage(mols, molsPerRow=molsPerRow, subImgSize=subImgSize, legends=legends,
                            highlightAtomLists=highlightAtomLists,
                            highlightBondLists=highlightBondLists,
                            **kwargs)


def _legacyReactionToImage(rxn, subImgSize=(200, 200), **kwargs):
  try:
    import Image
  except ImportError:
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


def ReactionToImage(rxn, subImgSize=(200, 200), useSVG=False, **kwargs):
  if not useSVG and not hasattr(rdMolDraw2D, 'MolDraw2DCairo'):
      return _legacyReactionToImage(rxn,subImgSize=subImgSize,**kwargs)
  else:
      width = subImgSize[0] * (rxn.GetNumReactantTemplates() + rxn.GetNumProductTemplates() + 1)
      if useSVG:
          d = rdMolDraw2D.MolDraw2DSVG(width, subImgSize[1])
      else:
          d = rdMolDraw2D.MolDraw2DCairo(width, subImgSize[1])
      d.DrawReaction(rxn,**kwargs)
      d.FinishDrawing()
      if useSVG:
          return d.GetDrawingText().replace('svg:', '')
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
