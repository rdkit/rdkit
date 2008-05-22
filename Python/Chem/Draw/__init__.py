# $Id$
#
# Copyright (C) 2006-2008 Greg Landrum
#  All Rights Reserved
#


def MolToImageFile(mol,filename,size=(300,300),kekulize=True, wedgeBonds=True):
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
    canvas = Canvas(size=(300,300),name=filename)
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

  if useAGG:
    canvas.flush()
    img.save(filename)
  else:
    canvas.save()
    
