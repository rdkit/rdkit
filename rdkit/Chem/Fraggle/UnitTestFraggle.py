# $Id$
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import doctest
import unittest

from rdkit import Chem
from rdkit.Chem.Fraggle import FraggleSim
from rdkit.Chem.Fraggle.FraggleSim import select_fragments


def load_tests(loader, tests, ignore):  # pylint: disable=unused-argument
  """ Add the Doctests from the module """
  tests.addTests(
    doctest.DocTestSuite(FraggleSim, optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE))
  return tests


class TestCase(unittest.TestCase):

  def test_generate_fraggle_fragmentation(self):
    mol = Chem.MolFromSmiles('COc1cc(CN2CCC(CC2)NC(=O)c2cncc(C)c2)c(OC)c2ccccc12')

    frags = FraggleSim.generate_fraggle_fragmentation(mol)
    self.assertEqual(len(frags), 16)
    expected = (
      '[*]C(=O)NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '[*]C(=O)c1cncc(C)c1.[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '[*]C(=O)c1cncc(C)c1.[*]c1cc(OC)c2ccccc2c1OC', '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '[*]C(=O)c1cncc(C)c1.[*]Cc1cc(OC)c2ccccc2c1OC',
      '[*]Cc1cc(OC)c2ccccc2c1OC.[*]NC(=O)c1cncc(C)c1', '[*]Cc1cc(OC)c2ccccc2c1OC.[*]c1cncc(C)c1',
      '[*]NC(=O)c1cncc(C)c1.[*]c1cc(OC)c2ccccc2c1OC', '[*]NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '[*]NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1cncc(C)c1',
      '[*]c1c(CN2CCC(NC(=O)c3cncc(C)c3)CC2)cc(OC)c2ccccc12',
      '[*]c1c(OC)cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c1[*]',
      '[*]c1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12',
      '[*]N1CCC(NC(=O)c2cncc(C)c2)CC1.[*]c1cc(OC)c2ccccc2c1OC',
      '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1cncc(C)c1', '[*]c1cc(OC)c2ccccc2c1OC.[*]c1cncc(C)c1')
    for smi in frags:
      self.assertTrue(smi in expected)

    # Test case for fragments that contain a cyclic and acyclic component
    mol = Chem.MolFromSmiles('c12c(CCC)cccc2cccc1')
    frags = FraggleSim.generate_fraggle_fragmentation(mol)
    expected = ['[*]CCC.[*]c1ccccc1[*]', '[*]Cc1cccc2ccccc12', '[*]c1cccc(CCC)c1[*]',
                '[*]c1cccc2ccccc12', '[*]c1ccccc1[*]']
    for smi in frags:
      self.assertTrue(smi in expected)

  def testFragmentation2(self):
    mol = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')

    frags = FraggleSim.generate_fraggle_fragmentation(mol)
    self.assertEqual(len(frags), 13)

    expected = (
      '[*]C(=O)c1ccccc1.[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '[*]C(=O)c1ccccc1.[*]Cc1cc(OC)c2ccccc2c1OC', '[*]C(=O)c1ccccc1.[*]c1cc(OC)c2ccccc2c1OC',
      '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1ccccc1',
      '[*]Cc1cc(OC)c2ccccc2c1OC.[*]NC(=O)c1ccccc1', '[*]Cc1cc(OC)c2ccccc2c1OC.[*]c1ccccc1',
      '[*]N1CCC(NC(=O)c2ccccc2)CC1.[*]c1cc(OC)c2ccccc2c1OC',
      '[*]NC(=O)c1ccccc1.[*]c1cc(OC)c2ccccc2c1OC',
      '[*]NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.[*]c1ccccc1',
      '[*]c1c(CN2CCC(NC(=O)c3ccccc3)CC2)cc(OC)c2ccccc12',
      '[*]c1c(OC)cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c1[*]',
      '[*]c1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12', '[*]c1cc(OC)c2ccccc2c1OC.[*]c1ccccc1')
    for smi in frags:
      self.assertTrue(smi in expected)

  def test_select_fragments(self):
    self.assertRaises(NotImplementedError, select_fragments, 'CCCC', 'InvalidFragmentationType', 10)

  def test_select_fragments_acyclic(self):
    # acyclic fragments: returns fragments with more than two heavy atoms
    smiles = '[*]C.[*]O[*].[*]c1cccc2ccccc12'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 12), '[*]c1cccc2ccccc12')
    # ignore the fragments if less than 60% of original molecule
    smiles = '[*]CC.[*]O[*].[*]c1ccccc1'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 9), '[*]c1ccccc1')
    smiles = '[*]CC.[*]OCC[*].[*]c1ccccc1'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 11), '[*]c1ccccc1')
    smiles = '[*]CC.[*]CCCO[*].[*]c1ccccc1'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 12), None)
    # ignore fragments that are smaller than 2 atoms
    smiles = '[*]C.[*]O[*].[*]c1ccccc1'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 8), '[*]c1ccccc1')
    smiles = '[*]CC.[*]O[*].[*]c1ccccc1'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 9), '[*]c1ccccc1')
    smiles = '[*]CCC.[*]O[*].[*]c1ccccc1'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 10), '[*]CCC.[*]c1ccccc1')
    # Only small fragments
    smiles = '[*]CC.[*]CC.[*]CC'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 6), None)
    smiles = '[*]CCC.[*]CCC.[*]CCC'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_ACYCLIC, 6), smiles)

  def test_select_fragments_cyclic(self):
    smiles = '[*]c1ccccc1[*].[*]cccc[*]'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_CYCLIC, 10), '[*]c1ccccc1[*]')
    smiles = '[*]c1ccccc1[*].[*]cccc[*]'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_CYCLIC, 10), '[*]c1ccccc1[*]')
    smiles = '[*]c1ccccc1[*].[*]C.[*]C'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_CYCLIC, 8), None)
    # Fragment too small
    smiles = '[*]c1ccccc1[*].[*]c(CCCCCCCCCC)ccc[*]'
    self.assertEqual(select_fragments(smiles, FraggleSim.FTYPE_CYCLIC, 20), None)

  def test_select_fragments_cyclic_acyclic(self):
    smiles = '[*]c1ccccc1[*].[*]CCC.[*]C'
    self.assertEqual(
      select_fragments(smiles, FraggleSim.FTYPE_CYCLIC_ACYCLIC, 10), '[*]c1ccccc1[*].[*]CCC')
    smiles = '[*]CCC.[*]C.[*]c1ccccc1[*]'
    self.assertEqual(
      select_fragments(smiles, FraggleSim.FTYPE_CYCLIC_ACYCLIC, 10), '[*]CCC.[*]c1ccccc1[*]')
    smiles = '[*]CCC.[*]CCCCCCCCCCC[*].[*]c1ccccc1[*]'
    self.assertEqual(
      select_fragments(smiles, FraggleSim.FTYPE_CYCLIC_ACYCLIC,
                       Chem.MolFromSmiles(smiles).GetNumAtoms()), None)

  def test_is_ring_cut_valid(self):
    self.assertEqual(FraggleSim.is_ring_cut_valid('[*]CCC[*]'), (False, 0))
    self.assertEqual(FraggleSim.is_ring_cut_valid('[*]C1CC1[*]'), (True, 5))
    self.assertEqual(FraggleSim.is_ring_cut_valid('[*]c1ccccc1[*]'), (True, 8))
    self.assertEqual(FraggleSim.is_ring_cut_valid('[*]cccc[*]'), (False, 0))

  def test_GetFraggleSimilarity(self):
    q = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12')
    m = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.980, places=2)
    self.assertEqual(match, '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1')

    m = Chem.MolFromSmiles('COc1cc(CN2CCC(Nc3nc4ccccc4s3)CC2)c(OC)c2ccccc12')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.794, places=2)
    self.assertEqual(match, '[*]C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1')

    q = Chem.MolFromSmiles('COc1ccccc1')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.347, places=2)
    self.assertEqual(match, '[*]c1ccccc1')

    m = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.266, places=2)
    self.assertEqual(match, '[*]c1ccccc1')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
