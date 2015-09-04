# $Id$
#
#  Copyright (C) 2011  greg Landrum 
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for molecule drawing
"""
from rdkit import RDConfig
import unittest,os,tempfile
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.RDLogger import logger
logger = logger()

class TestCase(unittest.TestCase):
  def setUp(self):
    self.mol = Chem.MolFromSmiles('c1c(C[15NH3+])ccnc1[C@](Cl)(Br)[C@](Cl)(Br)F')

  def testCairoFile(self):
    try:
      from rdkit.Chem.Draw.cairoCanvas import Canvas
    except ImportError:
      logger.info("Skipping cairo test")
      return
    os.environ['RDKIT_CANVAS']='cairo'

    foo,fn=tempfile.mkstemp(suffix='.png')
    foo=None
    self.assertEqual(os.path.getsize(fn),0)

    Draw.MolToFile(self.mol,fn)

    self.assertNotEqual(os.path.getsize(fn),0)
    try:
      os.unlink(fn)
    except:
      pass
    
  def testAggFile(self):
    try:
      from rdkit.Chem.Draw.aggCanvas import Canvas
    except ImportError:
      logger.info("Skipping agg test")
      return
    os.environ['RDKIT_CANVAS']='agg'

    foo,fn=tempfile.mkstemp(suffix='.png')
    foo=None
    self.assertEqual(os.path.getsize(fn),0)

    Draw.MolToFile(self.mol,fn)

    self.assertNotEqual(os.path.getsize(fn),0)
    try:
      os.unlink(fn)
    except:
      pass

  def testSpingFile(self):
    try:
      from rdkit.Chem.Draw.spingCanvas import Canvas
    except ImportError:
      logger.info("Skipping sping test")
      return
    os.environ['RDKIT_CANVAS']='sping'

    foo,fn=tempfile.mkstemp(suffix='.png')
    foo=None
    self.assertEqual(os.path.getsize(fn),0)

    Draw.MolToFile(self.mol,fn)

    self.assertNotEqual(os.path.getsize(fn),0)
    try:
      os.unlink(fn)
    except:
      pass

  def testCairoImage(self):
    try:
      from rdkit.Chem.Draw.cairoCanvas import Canvas
    except ImportError:
      return
    os.environ['RDKIT_CANVAS']='cairo'

    img=Draw.MolToImage(self.mol,size=(300,300))
    self.assertTrue(img)
    self.assertEqual(img.size[0],300)
    self.assertEqual(img.size[1],300)
    
  def testAggImage(self):
    try:
      from rdkit.Chem.Draw.aggCanvas import Canvas
    except ImportError:
      return
    os.environ['RDKIT_CANVAS']='agg'
    img=Draw.MolToImage(self.mol,size=(300,300))
    self.assertTrue(img)
    self.assertEqual(img.size[0],300)
    self.assertEqual(img.size[1],300)

  def testSpingImage(self):
    try:
      from rdkit.Chem.Draw.spingCanvas import Canvas
    except ImportError:
      return
    os.environ['RDKIT_CANVAS']='sping'
    img=Draw.MolToImage(self.mol,size=(300,300))
    self.assertTrue(img)
    self.assertEqual(img.size[0],300)
    self.assertEqual(img.size[1],300)

  def testQtImage(self):
    import sys
    try:
      from PySide import QtGui
      from rdkit.Chem.Draw.qtCanvas import Canvas
    except ImportError:
      return
    app = QtGui.QApplication(sys.argv)
    img = Draw.MolToQPixmap(self.mol, size=(300, 300))
    self.assertTrue(img)
    self.assertEqual(img.size().height(), 300)
    self.assertEqual(img.size().width(), 300)

  def testCairoImageDash(self):
    try:
      from rdkit.Chem.Draw.cairoCanvas import Canvas
    except ImportError:
      return
    os.environ['RDKIT_CANVAS']='cairo'

    img=Draw.MolToImage(self.mol,size=(300,300),kekulize=False)
    self.assertTrue(img)
    self.assertEqual(img.size[0],300)
    self.assertEqual(img.size[1],300)
    
  def testAggImageDash(self):
    try:
      from rdkit.Chem.Draw.aggCanvas import Canvas
    except ImportError:
      return
    os.environ['RDKIT_CANVAS']='agg'
    img=Draw.MolToImage(self.mol,size=(300,300),kekulize=False)
    self.assertTrue(img)
    self.assertEqual(img.size[0],300)
    self.assertEqual(img.size[1],300)

  def testSpingImageDash(self):
    try:
      from rdkit.Chem.Draw.spingCanvas import Canvas
    except ImportError:
      return
    os.environ['RDKIT_CANVAS']='sping'
    img=Draw.MolToImage(self.mol,size=(300,300),kekulize=False)
    self.assertTrue(img)
    self.assertEqual(img.size[0],300)
    self.assertEqual(img.size[1],300)

  def testGithubIssue54(self):
    try:
      from rdkit.Chem.Draw.spingCanvas import Canvas
    except ImportError:
      return
    os.environ['RDKIT_CANVAS']='sping'
    mol = Chem.MolFromSmiles('c1([O])ccc(O)cc1')
    img = Draw.MolToImage(mol)
    self.assertTrue(img)
    
  def testGithubIssue86(self):
    mol = Chem.MolFromSmiles('F[C@H](Cl)Br')
    for b in mol.GetBonds():
      self.assertEqual(b.GetBondDir(),Chem.BondDir.NONE)
    img = Draw.MolToImage(mol,kekulize=False)
    self.assertTrue(img)
    for b in mol.GetBonds():
      self.assertEqual(b.GetBondDir(),Chem.BondDir.NONE)

    Chem.WedgeMolBonds(mol,mol.GetConformer())
    obds = [x.GetBondDir() for x in mol.GetBonds()]
    self.assertEqual(obds.count(Chem.BondDir.NONE),2)
    img = Draw.MolToImage(mol,kekulize=False)
    self.assertTrue(img)
    nbds = [x.GetBondDir() for x in mol.GetBonds()]
    self.assertEqual(obds,nbds)

    


    
if __name__ == '__main__':
  unittest.main()

