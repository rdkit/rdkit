# $Id$
#
# Copyright (C) 2006-2011 Greg Landrum
#  All Rights Reserved
#
import os

from MolDrawing import MolDrawing

def _getCanvas():
  useAGG=False
  useCairo=False
  Canvas=None
  if not os.environ.get('RDKIT_CANVAS',''):
    try:
      from cairoCanvas import Canvas
      useCairo=True
    except ImportError:
      try:
        from aggCanvas import Canvas
        useAGG=True
      except ImportError:
        from spingCanvas import Canvas
  else:
    canv=os.environ['RDKIT_CANVAS'].lower()
    if canv =='cairo':
      from cairoCanvas import Canvas
      useCairo=True
    elif canv =='agg':
      from aggCanvas import Canvas
      useAGG=True
    else:
      from spingCanvas import Canvas      
  return useAGG,useCairo,Canvas

def MolToImage(mol, size=(300,300), kekulize=True, wedgeBonds=True,
               canvas=None, **kwargs):
  """ returns a PIL image containing a drawing of the molecule

    Keyword arguments:
    kekulize -- run kekulization routine on input `mol` (default True)
    size -- final image size, in pixel (default (300,300))
    wedgeBonds -- draw wedge (stereo) bonds (default True)
    highlightAtoms -- list of atoms to highlight (default [])
    highlightMap -- dictionary of (atom, color) pairs (default None)
  """
  if not mol:
    raise ValueError,'Null molecule provided'
  if canvas is None:
    useAGG,useCairo,Canvas=_getCanvas()
    if useAGG or useCairo:
      import Image
      img = Image.new("RGBA",size,"white")
      canvas = Canvas(img)
    else:
      from spingCanvas import Canvas
      canvas = Canvas(size=size,name='MolToImageFile')
      img = canvas._image
  drawer = MolDrawing(canvas)

  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
    
  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
  
  drawer.wedgeDashedBonds=wedgeBonds
  drawer.AddMol(mol,**kwargs)
  canvas.flush()

  return img

def MolToFile(mol,fileName,size=(300,300),kekulize=True, wedgeBonds=True,
              imageType=None,**kwargs):
  """ Generates a drawing of a molecule and writes it to a file
  """
  # original contribution from Uwe Hoffmann
  if not fileName:
    raise ValueError,'no fileName provided'
  if not mol:
    raise ValueError,'Null molecule provided'

  
  if imageType is None:
    imageType=os.path.splitext(fileName)[1][1:]

  useAGG,useCairo,Canvas = _getCanvas()
  if useCairo or useAGG:
    canvas = Canvas(size=size,imageType=imageType,
                              fileName=fileName)
  else:
    canvas = Canvas(size=size,name=fileName,imageType=imageType)
  drawer = MolDrawing(canvas)
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
    
  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
  
  drawer.wedgeDashedBonds=wedgeBonds
  drawer.AddMol(mol,**kwargs)
  if useCAIRO or useAgg:
    canvas.flush()
  else:
    canvas.save()

def MolToImageFile(mol,filename,size=(300,300),kekulize=True, wedgeBonds=True,
                   **kwargs):
  """  DEPRECATED:  please use MolToFile instead

  """
  img = MolToImage(mol,size=size,kekulize=kekulize,wedgeBonds=wedgeBonds,**kwargs)
  img.save(filename)
    
tkRoot=None
tkLabel=None
tkPI=None
def ShowMol(mol,size=(300,300),kekulize=True,wedgeBonds=True,
            title='RDKit Molecule',**kwargs):
  """ Generates a picture of a molecule and displays it in a Tkinter window
  """
  global tkRoot,tkLabel,tkPI
  import Tkinter
  import ImageTk

  img = MolToImage(mol,size,kekulize,wedgeBonds,**kwargs)

  if not tkRoot:
    tkRoot = Tkinter.Tk()
    tkRoot.title(title)
    tkPI = ImageTk.PhotoImage(img)
    tkLabel = Tkinter.Label(tkRoot,image=tkPI)
    tkLabel.place(x=0,y=0,width=img.size[0],height=img.size[1])
  else:
    tkPI.paste(img)
  tkRoot.geometry('%dx%d'%(img.size))

    
