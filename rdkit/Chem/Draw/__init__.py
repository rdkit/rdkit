# $Id$
#
# Copyright (C) 2006-2010 Greg Landrum
#  All Rights Reserved
#
import os.path

def MolToImage(mol,size=(300,300),kekulize=True, wedgeBonds=True,highlightAtoms=[]):
  """ returns a PIL image containing a drawing of the molecule
  """
  if not mol:
    raise ValueError,'Null molecule provided'
  import MolDrawing
  try:
    from aggdraw import Draw
    import Image
    MolDrawing.registerCanvas('agg')
    Canvas = Draw
    useAGG=True
  except:
    useAGG=False
    try:
      import cairo
      import Image
      MolDrawing.registerCanvas('cairo')
      from cairoCanvas import Canvas
      useCAIRO=True
    except:
      useCAIRO=False
      from rdkit.sping.PIL.pidPIL import PILCanvas as Canvas
      canvas = Canvas(size=size,name='MolToImageFile')
      img = canvas._image
      MolDrawing.registerCanvas('sping')
      drawer = MolDrawing.MolDrawing(canvas)
  if useAGG or useCAIRO:
    img = Image.new("RGBA",size,"white")
    canvas = Canvas(img)
    if useAGG:
      canvas.setantialias(True)
    drawer = MolDrawing.MolDrawing(canvas)

  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
    
  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
  
  drawer.wedgeDashedBonds=wedgeBonds
  drawer.AddMol(mol,highlightAtoms=highlightAtoms)
  canvas.flush()

  return img

def MolToFile(mol,fileName,size=(300,300),kekulize=True, wedgeBonds=True,
              highlightAtoms=[],imageType=None):
  """ Generates a drawing of a molecule and writes it to a file
  """
  # original contribution from Uwe Hoffmann
  if not fileName:
    raise ValueError,'no fileName provided'
  if not mol:
    raise ValueError,'Null molecule provided'

  import MolDrawing
  if imageType is None:
    imageType=os.path.splitext(fileName)[1][1:]
  try:
    import cairo
    MolDrawing.registerCanvas('cairo')
    canvas=cairoCanvas.Canvas(size=size,imageType=imageType,
                              fileName=fileName)
    useCAIRO=True
  except ImportError:
    useCAIRO=False
    MolDrawing.registerCanvas('sping')
    if imageType=="pdf":
      from rdkit.sping.PDF.pidPDF import PDFCanvas as Canvas
    elif imageType=="ps":
      from rdkit.sping.PS.pidPS import PSCanvas as Canvas
    elif imageType=="svg":
      from rdkit.sping.SVG.pidSVG import SVGCanvas as Canvas
    elif imageType=="png":
      from rdkit.sping.PIL.pidPIL import PILCanvas as Canvas
    canvas = Canvas(size=size,name=fileName)
  drawer = MolDrawing.MolDrawing(canvas)
  if kekulize:
    from rdkit import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
    
  if not mol.GetNumConformers():
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
  
  drawer.wedgeDashedBonds=wedgeBonds
  drawer.AddMol(mol,highlightAtoms=highlightAtoms)
  if useCAIRO:
    canvas.flush()
  else:
    canvas.save()

def MolToImageFile(mol,filename,size=(300,300),kekulize=True, wedgeBonds=True,
                   highlightAtoms=[]):
  """  DEPRECATED:  please use MolToFile instead

  """
  img = MolToImage(mol,size=size,kekulize=kekulize,wedgeBonds=wedgeBonds,highlightAtoms=highlightAtoms)
  img.save(filename)
    
tkRoot=None
tkLabel=None
tkPI=None
def ShowMol(mol,size=(300,300),kekulize=True,wedgeBonds=True,
            title='RDKit Molecule'):
  """ Generates a picture of a molecule and displays it in a Tkinter window
  """
  global tkRoot,tkLabel,tkPI
  import Tkinter
  import ImageTk

  img = MolToImage(mol,size,kekulize,wedgeBonds)

  if not tkRoot:
    tkRoot = Tkinter.Tk()
    tkRoot.title(title)
    tkPI = ImageTk.PhotoImage(img)
    tkLabel = Tkinter.Label(tkRoot,image=tkPI)
    tkLabel.place(x=0,y=0,width=img.size[0],height=img.size[1])
  else:
    tkPI.paste(img)
  tkRoot.geometry('%dx%d'%(img.size))

    
