#
#  Copyright (C) 2011-2021  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import re
import os
import pathlib
import sys
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdMolDescriptors

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

  def setUp(self):
    if IPythonConsole is not None and Draw.MolsToGridImage == IPythonConsole.ShowMols:
      IPythonConsole.UninstallIPythonRenderer()
    self.mol = Chem.MolFromSmiles('c1c(C[15NH3+])ccnc1[C@](Cl)(Br)[C@](Cl)(Br)F')

  def test_interactive(self):
    # We avoid checking in the code with development flag set
    self.assertFalse(self.showAllImages)

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
    from rdkit.Chem.Draw.rdMolDraw2DQt import rdkitQtVersion
    if rdkitQtVersion.startswith('6'):
      try:
        from PyQt6 import QtGui
      except ImportError:
        # PySide version numbers leapt at Qt6
        from PySide6 import QtGui
    else:
      try:
        from PyQt5 import QtGui
      except ImportError:
        try:
          from PySide import QtGui
        except ImportError:
          # PySide2 supports Qt >= 5.12
          from PySide2 import QtGui

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

    rdDepictor.Compute2DCoords(mol)
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
    self.assertIn("width='600px' height='800px'", svg)
    svg = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=4, subImgSize=(200, 200),
                               useSVG=True)
    self.assertIn("width='800px' height='600px'", svg)
    svg = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=3, subImgSize=(300, 300),
                               useSVG=True)
    self.assertIn("width='900px' height='1200px'", svg)
    self.assertNotIn("class='note'", svg)
    dopts = Draw.rdMolDraw2D.MolDrawOptions()
    dopts.addAtomIndices = True
    svg = Draw.MolsToGridImage(mols, legends=legends, molsPerRow=3, subImgSize=(300, 300),
                               useSVG=True, drawOptions=dopts)
    self.assertIn("class='note'", svg)

  def testDrawMorgan(self):
    m = Chem.MolFromSmiles('c1ccccc1CC1CC1')
    bi = {}
    _ = rdMolDescriptors.GetMorganFingerprintAsBitVect(m, radius=2, bitInfo=bi)
    self.assertTrue(872 in bi)

    svg1 = Draw.DrawMorganBit(m, 872, bi)
    aid, r = bi[872][0]
    svg2 = Draw.DrawMorganEnv(m, aid, r)
    self.assertEqual(svg1, svg2)
    self.assertTrue("style='fill:#CCCCCC;" in svg1)
    self.assertTrue("style='fill:#E5E533;" in svg1)
    self.assertTrue("style='fill:#9999E5;" in svg1)

    svg1 = Draw.DrawMorganBit(m, 872, bi, centerColor=None)
    aid, r = bi[872][0]
    svg2 = Draw.DrawMorganEnv(m, aid, r, centerColor=None)
    self.assertEqual(svg1, svg2)
    self.assertTrue("style='fill:#CCCCCC;" in svg1)
    self.assertTrue("style='fill:#E5E533;" in svg1)
    self.assertFalse("style='fill:#9999E5;" in svg1)
    with self.assertRaises(KeyError):
      Draw.DrawMorganBit(m, 32, bi)

    if hasattr(Draw, 'MolDraw2DCairo'):
      # Github #3796: make sure we aren't trying to generate metadata:
      png = Draw.DrawMorganBit(m, 872, bi, useSVG=False)
      self.assertIn(b'PNG', png)
      self.assertIsNone(Chem.MolFromPNGString(png))

  def testDrawRDKit(self):
    m = Chem.MolFromSmiles('c1ccccc1CC1CC1')
    bi = {}
    _ = Chem.RDKFingerprint(m, maxPath=5, bitInfo=bi)
    self.assertTrue(1553 in bi)
    svg1 = Draw.DrawRDKitBit(m, 1553, bi)
    path = bi[1553][0]
    svg2 = Draw.DrawRDKitEnv(m, path)
    self.assertEqual(svg1, svg2)
    self.assertTrue("style='fill:#E5E533;" in svg1)
    self.assertFalse("style='fill:#CCCCCC;" in svg1)
    self.assertFalse("style='fill:#9999E5;" in svg1)
    with self.assertRaises(KeyError):
      Draw.DrawRDKitBit(m, 32, bi)

    if hasattr(Draw, 'MolDraw2DCairo'):
      # Github #3796: make sure we aren't trying to generate metadata:
      png = Draw.DrawRDKitBit(m, 1553, bi, useSVG=False)
      self.assertIn(b'PNG', png)
      self.assertIsNone(Chem.MolFromPNGString(png))

  def testDrawReaction(self):
    # this shouldn't throw an exception...
    rxn = AllChem.ReactionFromSmarts(
      "[c;H1:3]1:[c:4]:[c:5]:[c;H1:6]:[c:7]2:[nH:8]:[c:9]:[c;H1:1]:[c:2]:1:2.O=[C:10]1[#6;H2:11][#6;H2:12][N:13][#6;H2:14][#6;H2:15]1>>[#6;H2:12]3[#6;H1:11]=[C:10]([c:1]1:[c:9]:[n:8]:[c:7]2:[c:6]:[c:5]:[c:4]:[c:3]:[c:2]:1:2)[#6;H2:15][#6;H2:14][N:13]3"
    )
    _ = Draw.ReactionToImage(rxn)

  def testGithub3762(self):
    m = Chem.MolFromSmiles('CC(=O)O')
    ats = [1, 2, 3]
    svg1 = Draw._moltoSVG(m, (250, 200), ats, "", False)
    svg2 = Draw._moltoSVG(m, (250, 200), ats, "", False, highlightBonds=[])
    with open('testGithub_3762_1.svg', 'w') as f:
      f.write(svg1)
    with open('testGithub_3762_2.svg', 'w') as f:
      f.write(svg2)

    re_str = r"path class='bond-\d atom-\d atom-\d' d='M \d+.\d+,\d+.\d+ L \d+.\d+,\d+.\d+ L \d+.\d+,\d+.\d+ L \d+.\d+,\d+.\d+ Z' style='fill:#FF7F7F;"
    patt = re.compile(re_str)
    self.assertEqual(len(patt.findall(svg1)), 2)
    self.assertEqual(len(patt.findall(svg2)), 0)
    
    pathlib.Path('testGithub_3762_1.svg').unlink()
    pathlib.Path('testGithub_3762_2.svg').unlink()
    
  def testGithub5863(self):
    smiles = "C[C@]12C[C@H](O)[C@H]3[C@@H](CCC4=CC(=O)C=C[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO"
    mol = Chem.MolFromSmiles(smiles)

    info = {}
    rdMolDescriptors.GetMorganFingerprint(mol, radius=2, bitInfo=info, useChirality=True)
    bitId = 1236726849
    # this should run without generating an exception:
    Draw.DrawMorganBit(mol, bitId, info)


if __name__ == '__main__':
  unittest.main()
