#
#  Copyright (C) 2011-2021  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import os
import pathlib
import re
import sys
import tempfile
import unittest

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

try:
  from rdkit.Chem.Draw import IPythonConsole
except ImportError:
  IPythonConsole = None


class TestCase(unittest.TestCase):
  showAllImages = False

  def setUp(self):
    if IPythonConsole is not None and Draw.MolsToGridImage == IPythonConsole.ShowMols:
      IPythonConsole.UninstallIPythonRenderer()
    self.mol = Chem.MolFromSmiles('c1c(C[15NH3+])ccnc1[C@](Cl)(Br)[C@](Cl)(Br)F')

    # For testMolsMatrixToLinear and testMolsMatrixToLinear (which test MolsMatrixToGridImage and its helper _MolsNestedToLinear)
    s = "NC(C)C(=O)"
    mol = Chem.MolFromSmiles(s)
    natoms = mol.GetNumAtoms()
    nbonds = mol.GetNumBonds()

    # Set up matrix with oligomer count for the molecules
    # Should produce this grid:
    # NC(C)C(=O)
    #                                NC(C)C(=O)NC(C)C(=O)
    # NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)                      NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)NC(C)C(=O)
    repeats = [[1], [0, 2], [3, 0, 4]]

    # Create molecule if there are 1 or more oligomers;
    # otherwise, use None for molecule because drawing is distorted if use Chem.MolFromSmiles("")
    self.molsMatrix = [[Chem.MolFromSmiles(s * count) if count else None for count in row]
                       for row in repeats]

    self.legendsMatrix = [[str(count) + " unit(s)" for count in row] for row in repeats]

    def ithItemList(nunits, itemsPerUnit, i=0):
      return [((n * itemsPerUnit) + i) for n in range(nunits)]

    self.highlightAtomListsMatrix = [[ithItemList(count, natoms, 0) for count in row]
                                     for row in repeats]

    # Another bond is created when molecule is oligomerized, so to keep the bond type consistent,
    #   make items per unit one more than the number of bonds
    self.highlightBondListsMatrix = [[ithItemList(count, nbonds + 1, 1) for count in row]
                                     for row in repeats]

    # Parametrize tests: In addition to supplying molsMatrix, supply 0-3 other matrices
    # col labels: legendsMatrix, highlightAtomListsMatrix, highlightBondListsMatrix
    # Zero other matrices: 1 parameter set
    self.paramSets = [
      (None, None, None),
      # One other matrix: 3 parameter sets
      (self.legendsMatrix, None, None),
      (None, self.highlightAtomListsMatrix, None),
      (None, None, self.highlightBondListsMatrix),
      # Two other matrices: 3 parameter sets
      (self.legendsMatrix, self.highlightAtomListsMatrix, None),
      (self.legendsMatrix, None, self.highlightBondListsMatrix),
      (None, self.highlightAtomListsMatrix, self.highlightBondListsMatrix),
      # All three other matrices: 1 parameter set
      (self.legendsMatrix, self.highlightAtomListsMatrix, self.highlightBondListsMatrix),
    ]

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

  @unittest.skipIf(not hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Skipping cairo test')
  def testCairoFile(self):
    self._testMolToFile()

  def _testMolToImage(self, mol=None, kekulize=True, options=None, showImage=False, **kwargs):
    mol = mol or self.mol
    img = Draw.MolToImage(mol, size=(300, 300), kekulize=kekulize, options=options, **kwargs)
    self.assertTrue(img)
    self.assertEqual(img.size[0], 300)
    self.assertEqual(img.size[1], 300)
    if self.showAllImages or showImage:
      img.show()

  @unittest.skipIf(not hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Skipping cairo test')
  def testCairoImage(self):
    self._testMolToImage()

  @unittest.skipIf(not hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Skipping cairo test')
  def testCairoImageDash(self):
    os.environ['RDKIT_CANVAS'] = 'cairo'
    self._testMolToImage(kekulize=False)

  @unittest.skipIf(not hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Skipping test requiring cairo')
  def testSpecialCases(self):
    options = rdMolDraw2D.MolDrawOptions()
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

  @unittest.skipIf(not hasattr(rdMolDraw2D, 'MolDraw2DCairo'), 'Skipping test requiring cairo')
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

  def testMolsMatrixToLinear(self):
    mols, molsPerRow, legends, highlightAtomLists, highlightBondLists = Draw._MolsNestedToLinear(
      self.molsMatrix, self.legendsMatrix, self.highlightAtomListsMatrix,
      self.highlightBondListsMatrix)

    nrows = len(self.molsMatrix)

    def _nestedOrder(self, molsMatrix, legendsMatrix, highlightAtomListsMatrix,
                     highlightBondListsMatrix):
      for r, row in enumerate(molsMatrix):
        for c, item in enumerate(row):
          linearIndex = (r * molsPerRow) + c
          # Test that items in 2D list are in correct position in 1D list
          self.assertTrue(mols[linearIndex] == item)
          # Test that 1D list items are not lists
          self.assertFalse(isinstance(item, list))

      # Other three matrices (legends; atom and bond highlighting) need not be supplied;
      #   only test each if it's supplied
      if legendsMatrix is not None:
        for r, row in enumerate(legendsMatrix):
          for c, item in enumerate(row):
            linearIndex = (r * molsPerRow) + c
            # Test that items in 2D list are in correct position in 1D list
            self.assertTrue(legends[linearIndex] == item)
            # Test that 1D list items are not lists
            self.assertFalse(isinstance(item, list))

      if highlightAtomListsMatrix is not None:
        for r, row in enumerate(highlightAtomListsMatrix):
          for c, item in enumerate(row):
            linearIndex = (r * molsPerRow) + c
            # Test that items in 2D list are in correct position in 1D list
            self.assertTrue(highlightAtomLists[linearIndex] == item)
            # For highlight parameters, entries are lists, so check that sub-items are not lists
            for subitem in item:
              self.assertFalse(isinstance(subitem, list))

      if highlightBondListsMatrix is not None:
        for r, row in enumerate(highlightBondListsMatrix):
          for c, item in enumerate(row):
            linearIndex = (r * molsPerRow) + c
            # Test that items in 2D list are in correct position in 1D list
            self.assertTrue(highlightBondLists[linearIndex] == item)
            # For highlight parameters, entries are lists, so check that sub-items are not lists
            for subitem in item:
              self.assertFalse(isinstance(subitem, list))

      # Test that 1D list has the correct length
      self.assertTrue(len(mols) == nrows * molsPerRow)

    # Parametrize tests: In addition to supplying molsMatrix, supply 0-3 other matrices
    for paramSet in self.paramSets:
      _nestedOrder(self, self.molsMatrix, *paramSet)

    ## Test that exceptions are thrown appropriately

    # Test that supplying a non-nested list raises a ValueError
    # Set up non-nested lists = first sublist of nested lists
    molsNotNested = self.molsMatrix[0]
    legendsNotNested = self.legendsMatrix[0]
    highlightAtomListsNotNested = self.highlightAtomListsMatrix[0]
    highlightBondListsNotNested = self.highlightBondListsMatrix[0]

    with self.assertRaises(TypeError):
      Draw._MolsNestedToLinear(molsNotNested)

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix, legendsMatrix=legendsNotNested)

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix,
                               highlightAtomListsMatrix=highlightAtomListsNotNested)

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix,
                               highlightBondListsMatrix=highlightBondListsNotNested)

    # Test that raises ValueError if other matrices aren't same size (# rows) as molsMatrix
    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix, legendsMatrix=self.legendsMatrix[0:1])

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix,
                               highlightAtomListsMatrix=self.highlightAtomListsMatrix[0:1])

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix,
                               highlightBondListsMatrix=self.highlightBondListsMatrix[0:1])

    # Test that raises ValueError if other matrices' rows aren't same length as molsMatrix's corresponding row
    # Remove last element from first row of each other matrix
    legendsMatrixShortRow0 = [self.legendsMatrix[0][0:-1]] + [self.legendsMatrix[1:]]
    highlightAtomListsMatrixShortRow0 = [self.highlightAtomListsMatrix[0][0:-1]
                                         ] + [self.highlightAtomListsMatrix[1:]]
    highlightBondListsMatrixShortRow0 = [self.highlightBondListsMatrix[0][0:-1]
                                         ] + [self.highlightBondListsMatrix[1:]]

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix, legendsMatrix=legendsMatrixShortRow0)

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix,
                               highlightAtomListsMatrix=highlightAtomListsMatrixShortRow0)

    with self.assertRaises(ValueError):
      Draw._MolsNestedToLinear(molsMatrix=self.molsMatrix,
                               highlightBondListsMatrix=highlightBondListsMatrixShortRow0)

  def testMolsMatrixToGridImage(self):
    subImgSize = (200, 200)

    kwargsValue = Draw.rdMolDraw2D.MolDrawOptions()
    kwargsValue.addAtomIndices = True

    # Tests are that running Draw.MolsMatrixToGridImage doesn't give an error
    # with any combination of parameters supplied (or not)
    # by parametrizing tests: In addition to supplying molsMatrix, supply 0-3 other matrices
    for (legendsMatrix, highlightAtomListsMatrix, highlightBondListsMatrix) in self.paramSets:
      for useSVG in (True, False):
        if not (useSVG or hasattr(rdMolDraw2D, 'MolDraw2DCairo')):
          continue
        for returnPNG in (True, False):
          dwgSubImgSizeNokwargs = Draw.MolsMatrixToGridImage(
            self.molsMatrix, subImgSize, legendsMatrix=legendsMatrix,
            highlightAtomListsMatrix=highlightAtomListsMatrix,
            highlightBondListsMatrix=highlightBondListsMatrix, useSVG=useSVG, returnPNG=returnPNG)
          if useSVG:
            self.assertIn('<svg', dwgSubImgSizeNokwargs)
          dwgSubImgSizeKwargs = Draw.MolsMatrixToGridImage(
            self.molsMatrix, subImgSize, legendsMatrix=legendsMatrix,
            highlightAtomListsMatrix=highlightAtomListsMatrix,
            highlightBondListsMatrix=highlightBondListsMatrix, useSVG=useSVG, returnPNG=returnPNG,
            drawOptions=kwargsValue)
          if useSVG:
            self.assertIn('<svg', dwgSubImgSizeKwargs)
          dwgNosubImgSizeNokwargs = Draw.MolsMatrixToGridImage(
            self.molsMatrix, legendsMatrix=legendsMatrix,
            highlightAtomListsMatrix=highlightAtomListsMatrix,
            highlightBondListsMatrix=highlightBondListsMatrix, useSVG=useSVG, returnPNG=returnPNG)
          if useSVG:
            self.assertIn('<svg', dwgNosubImgSizeNokwargs)
          dwgNosubImgSizeKwargs = Draw.MolsMatrixToGridImage(
            self.molsMatrix, legendsMatrix=legendsMatrix,
            highlightAtomListsMatrix=highlightAtomListsMatrix,
            highlightBondListsMatrix=highlightBondListsMatrix, useSVG=useSVG, returnPNG=returnPNG,
            drawOptions=kwargsValue)
          if useSVG:
            self.assertIn('<svg', dwgNosubImgSizeKwargs)

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

    if hasattr(Draw, 'MolDraw2DCairo') and hasattr(Draw, 'MolFromPNGString'):
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

    if hasattr(Draw, 'MolDraw2DCairo') and hasattr(Draw, 'MolFromPNGString'):
      # Github #3796: make sure we aren't trying to generate metadata:
      png = Draw.DrawRDKitBit(m, 1553, bi, useSVG=False)
      self.assertIn(b'PNG', png)
      self.assertIsNone(Chem.MolFromPNGString(png))

  def testDrawReaction(self):
    # this shouldn't throw an exception...
    rxn = AllChem.ReactionFromSmarts(
      "[c;H1:3]1:[c:4]:[c:5]:[c;H1:6]:[c:7]2:[nH:8]:[c:9]:[c;H1:1]:[c:2]:1:2.O=[C:10]1[#6;H2:11][#6;H2:12][N:13][#6;H2:14][#6;H2:15]1>>[#6;H2:12]3[#6;H1:11]=[C:10]([c:1]1:[c:9]:[n:8]:[c:7]2:[c:6]:[c:5]:[c:4]:[c:3]:[c:2]:1:2)[#6;H2:15][#6;H2:14][N:13]3"
    )
    _ = Draw.ReactionToImage(rxn, useSVG=True)

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
