#
#  Copyright (C) 2011-2017  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os
import sys
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import Draw

try:
  from rdkit.Chem.Draw import IPythonConsole
except ImportError:
  IPythonConsole = None
try:
  from rdkit.Chem.Draw import cairoCanvas
except ImportError:
  cairoCanvas = None
try:
  from rdkit.Chem.Draw import spingCanvas
except ImportError:
  spingCanvas = None
try:
  from rdkit.Chem.Draw import aggCanvas
except ImportError:
  aggCanvas = None
try:
  from rdkit.Chem.Draw import qtCanvas
except ImportError:
  qtCanvas = None


class TestCase(unittest.TestCase):
  showAllImages = False

  def test_interactive(self):
    # We avoid checking in the code with development flag set
    self.assertFalse(self.showAllImages)

  def setUp(self):
    if IPythonConsole is not None and Draw.MolsToGridImage == IPythonConsole.ShowMols:
      IPythonConsole.UninstallIPythonRenderer()
    self.mol = Chem.MolFromSmiles('c1c(C[15NH3+])ccnc1[C@](Cl)(Br)[C@](Cl)(Br)F')

  def _testMolToFile(self):
    try:
      fhdl, fn = tempfile.mkstemp(suffix='.png')
      # mkstemp returns a file handle that we don't need; close it
      os.close(fhdl)
      fhdl = None
      self.assertEqual(os.path.getsize(fn), 0)
      Draw.MolToFile(self.mol, fn)
      self.assertNotEqual(os.path.getsize(fn), 0)
    finally:
      os.remove(fn)

  @unittest.skipIf(cairoCanvas is None, 'Skipping cairo test')
  def testCairoFile(self):
    os.environ['RDKIT_CANVAS'] = 'cairo'
    self._testMolToFile()

  @unittest.skipIf(aggCanvas is None, 'Skipping agg test')
  def testAggFile(self):
    os.environ['RDKIT_CANVAS'] = 'agg'
    self._testMolToFile()

  @unittest.skipIf(spingCanvas is None, 'Skipping sping test')
  def testSpingFile(self):
    os.environ['RDKIT_CANVAS'] = 'sping'
    self._testMolToFile()

  def _testMolToImage(self, mol=None, kekulize=True, options=None, showImage=False, **kwargs):
    mol = mol or self.mol
    img = Draw.MolToImage(mol, size=(300, 300), kekulize=kekulize, options=options, **kwargs)
    self.assertTrue(img)
    self.assertEqual(img.size[0], 300)
    self.assertEqual(img.size[1], 300)
    if self.showAllImages or showImage:
      img.show()

  @unittest.skipIf(cairoCanvas is None, 'Skipping cairo test')
  def testCairoImage(self):
    os.environ['RDKIT_CANVAS'] = 'cairo'
    self._testMolToImage()

  @unittest.skipIf(aggCanvas is None, 'Skipping agg test')
  def testAggImage(self):
    os.environ['RDKIT_CANVAS'] = 'agg'
    self._testMolToImage()

  @unittest.skipIf(spingCanvas is None, 'Skipping sping test')
  def testSpingImage(self):
    os.environ['RDKIT_CANVAS'] = 'sping'
    self._testMolToImage()

  @unittest.skipIf(qtCanvas is None, 'Skipping Qt test')
  def testQtImage(self):
    try:
      from PySide import QtGui
      _ = QtGui.QApplication(sys.argv)
    except ImportError:
      from PyQt5 import QtGui
      _ = QtGui.QGuiApplication(sys.argv)
    img = Draw.MolToQPixmap(self.mol, size=(300, 300))
    self.assertTrue(img)
    self.assertEqual(img.size().height(), 300)
    self.assertEqual(img.size().width(), 300)
    # img.save('/tmp/D_me.png')

  @unittest.skipIf(cairoCanvas is None, 'Skipping cairo test')
  def testCairoImageDash(self):
    os.environ['RDKIT_CANVAS'] = 'cairo'
    self._testMolToImage(kekulize=False)

  @unittest.skipIf(aggCanvas is None, 'Skipping agg test')
  def testAggImageDash(self):
    os.environ['RDKIT_CANVAS'] = 'agg'
    self._testMolToImage(kekulize=False)

  @unittest.skipIf(spingCanvas is None, 'Skipping sping test')
  def testSpingImageDash(self):
    os.environ['RDKIT_CANVAS'] = 'sping'
    self._testMolToImage(kekulize=False, showImage=False)

  @unittest.skipIf(spingCanvas is None, 'Skipping sping test')
  def testGithubIssue54(self):
    # Assert that radicals depict with PIL
    os.environ['RDKIT_CANVAS'] = 'sping'
    mol = Chem.MolFromSmiles('c1([O])ccc(O)cc1')
    img = Draw.MolToImage(mol)
    self.assertTrue(img)
    # img.show()

  def testSpecialCases(self):
    options = Draw.DrawingOptions()
    options.atomLabelDeuteriumTritium = True
    self._testMolToImage(mol=Chem.MolFromSmiles('[2H][C@]([3H])(C)F'), options=options)
    # shared rings
    self._testMolToImage(mol=Chem.MolFromSmiles('c1cccc2cc(cccc3)c3cc21'))
    self._testMolToImage(mol=Chem.MolFromSmiles('C1=CC=CC=CC=C1'))
    self._testMolToImage(mol=Chem.MolFromSmiles('C=C=C'))
    self._testMolToImage(mol=Chem.MolFromSmiles('CC#N'), showImage=False)
    self._testMolToImage(mol=Chem.MolFromSmiles('[CH2-][C-2]C[CH3+][CH5+2]'))
    self._testMolToImage(mol=Chem.MolFromSmiles('[Na+].[OH-]'))
    self._testMolToImage(mol=Chem.MolFromSmiles('c1ccccc1c1ccccc1'),
                         highlightAtoms=(0, 1, 2, 3, 4, 5, 6))
    self._testMolToImage(mol=Chem.MolFromSmiles('c1ccccc1c1ccccc1'),
                         highlightBonds=(0, 2, 4, 6, 8, 10))
    self._testMolToImage(mol=Chem.MolFromSmiles('c1ccccc1c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)'))

  def testGithubIssue86(self):
    # Assert that drawing code doesn't modify wedge bonds
    mol = Chem.MolFromSmiles('F[C@H](Cl)Br')
    for b in mol.GetBonds():
      self.assertEqual(b.GetBondDir(), Chem.BondDir.NONE)
    img = Draw.MolToImage(mol, kekulize=False)
    self.assertTrue(img)
    # img.show()
    for b in mol.GetBonds():
      self.assertEqual(b.GetBondDir(), Chem.BondDir.NONE)

    Chem.WedgeMolBonds(mol, mol.GetConformer())
    obds = [x.GetBondDir() for x in mol.GetBonds()]
    self.assertEqual(obds.count(Chem.BondDir.NONE), 2)
    img = Draw.MolToImage(mol, kekulize=False)
    self.assertTrue(img)
    # img.show()
    nbds = [x.GetBondDir() for x in mol.GetBonds()]
    self.assertEqual(obds, nbds)

  def testGridSVG(self):
    mols = [Chem.MolFromSmiles('NC(C)C(=O)' * x) for x in range(10)]
    legends = ['mol-%d' % x for x in range(len(mols))]
    svg = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=3, subImgSize=(200, 200),
                               useSVG=True)
    self.assertTrue(svg.find("width='600px' height='800px'") > -1)
    svg = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=4, subImgSize=(200, 200),
                               useSVG=True)
    self.assertTrue(svg.find("width='800px' height='600px'") > -1)
    svg = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=3, subImgSize=(300, 300),
                               useSVG=True)
    self.assertTrue(svg.find("width='900px' height='1200px'") > -1)

  def testDrawMorgan(self):
    from rdkit.Chem import rdMolDescriptors
    m = Chem.MolFromSmiles('c1ccccc1CC1CC1')
    bi = {}
    fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(m,radius=2,bitInfo=bi)
    self.assertTrue(872 in bi)

    svg1 = Draw.DrawMorganBit(m,872,bi)
    aid,r = bi[872][0]
    svg2 = Draw.DrawMorganEnv(m,aid,r)
    self.assertEqual(svg1,svg2)
    self.assertTrue("style='fill:#CCCCCC;" in svg1)
    self.assertTrue("style='fill:#E5E533;" in svg1)
    self.assertTrue("style='fill:#9999E5;" in svg1)

    svg1 = Draw.DrawMorganBit(m,872,bi,centerColor=None)
    aid,r = bi[872][0]
    svg2 = Draw.DrawMorganEnv(m,aid,r,centerColor=None)
    self.assertEqual(svg1,svg2)
    self.assertTrue("style='fill:#CCCCCC;" in svg1)
    self.assertTrue("style='fill:#E5E533;" in svg1)
    self.assertFalse("style='fill:#9999E5;" in svg1)
    with self.assertRaises(KeyError):
        Draw.DrawMorganBit(m,32,bi)

  def testDrawRDKit(self):
    m = Chem.MolFromSmiles('c1ccccc1CC1CC1')
    bi = {}
    rdkfp = Chem.RDKFingerprint(m,maxPath=5,bitInfo=bi)
    self.assertTrue(1553 in bi)
    svg1 = Draw.DrawRDKitBit(m,1553,bi)
    path = bi[1553][0]
    svg2 = Draw.DrawRDKitEnv(m,path)
    self.assertEqual(svg1,svg2)
    self.assertTrue("style='fill:#E5E533;" in svg1)
    self.assertFalse("style='fill:#CCCCCC;" in svg1)
    self.assertFalse("style='fill:#9999E5;" in svg1)
    with self.assertRaises(KeyError):
        Draw.DrawRDKitBit(m,32,bi)




if __name__ == '__main__':
  unittest.main()
