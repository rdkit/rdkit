# $Id$
#
# Copyright (C) 2006-2008 Greg Landrum
#  All Rights Reserved
#


def MolToImage(mol,size=(300,300),kekulize=True, wedgeBonds=True):
  if not mol:
    raise ValueError,'Null molecule provided'
  import MolDrawing
  try:
    from aggdraw import Draw
    from PIL import Image
    MolDrawing.registerCanvas('agg')
    useAGG=True
  except:
    import traceback
    traceback.print_exc()
    useAGG=False
    from sping.PIL.pidPIL import PILCanvas as Canvas
    canvas = Canvas(size=size,name='MolToImageFile')
    img = canvas._image
    MolDrawing.registerCanvas('sping')
    drawer = MolDrawing.MolDrawing(canvas)
  if useAGG:
    img = Image.new("RGBA",size,"white")
    canvas = Draw(img)
    canvas.setantialias(True)
    drawer = MolDrawing.MolDrawing(canvas)

  if kekulize:
    import Chem
    mol = Chem.Mol(mol.ToBinary())
    Chem.Kekulize(mol)
    
  if not mol.GetNumConformers():
    from Chem import AllChem
    AllChem.Compute2DCoords(mol)
  
  drawer.wedgeDashedBonds=wedgeBonds
  drawer.AddMol(mol)
  canvas.flush()

  return img

def MolToImageFile(mol,filename,size=(300,300),kekulize=True, wedgeBonds=True):
  img = MolToImage(mol,size=size,kekulize=kekulize,wedgeBonds=wedgeBonds)
  img.save(filename)
    
tkRoot=None
tkLabel=None
tkPI=None
def ShowMol(mol,size=(300,300),kekulize=True,wedgeBonds=True,
            title='RDKit Molecule'):
  global tkRoot,tkLabel,tkPI
  import Tkinter
  from PIL import ImageTk

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

    
