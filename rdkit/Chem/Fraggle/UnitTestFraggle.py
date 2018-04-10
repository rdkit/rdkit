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
from rdkit import rdBase
from rdkit.Chem.Fraggle import FraggleSim
from rdkit.Chem.Fraggle.FraggleSim import select_fragments


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(
    doctest.DocTestSuite(FraggleSim, optionflags=doctest.ELLIPSIS + doctest.NORMALIZE_WHITESPACE))
  return tests


def _of(smiles):
  """ Order the fragments alphabetically. If smiles is None, returns None """
  if smiles is None:
    return None
  return '.'.join(sorted(smiles.split('.')))


class TestCase(unittest.TestCase):

  def test_generate_fraggle_fragmentation(self):
    mol = Chem.MolFromSmiles('COc1cc(CN2CCC(CC2)NC(=O)c2cncc(C)c2)c(OC)c2ccccc12')

    frags = FraggleSim.generate_fraggle_fragmentation(mol)
    self.assertEqual(len(frags), 16)
    expected = (
      '*C(=O)NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*C(=O)c1cncc(C)c1.*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*C(=O)c1cncc(C)c1.*c1cc(OC)c2ccccc2c1OC', '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*C(=O)c1cncc(C)c1.*Cc1cc(OC)c2ccccc2c1OC',
      '*Cc1cc(OC)c2ccccc2c1OC.*NC(=O)c1cncc(C)c1', '*Cc1cc(OC)c2ccccc2c1OC.*c1cncc(C)c1',
      '*NC(=O)c1cncc(C)c1.*c1cc(OC)c2ccccc2c1OC', '*NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.*c1cncc(C)c1',
      '*c1c(CN2CCC(NC(=O)c3cncc(C)c3)CC2)cc(OC)c2ccccc12',
      '*c1c(OC)cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c1*',
      '*c1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12',
      '*N1CCC(NC(=O)c2cncc(C)c2)CC1.*c1cc(OC)c2ccccc2c1OC',
      '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.*c1cncc(C)c1', '*c1cc(OC)c2ccccc2c1OC.*c1cncc(C)c1')
    expected = [_of(s) for s in expected]
    for smi in frags:
      self.assertTrue(_of(smi) in expected)

    # Test case for fragments that contain a cyclic and acyclic component
    mol = Chem.MolFromSmiles('c12c(CCC)cccc2cccc1')
    frags = FraggleSim.generate_fraggle_fragmentation(mol)
    expected = ['*CCC.*c1ccccc1*', '*Cc1cccc2ccccc12', '*c1cccc(CCC)c1*',
                '*c1cccc2ccccc12', '*c1ccccc1*']
    expected = [_of(s) for s in expected]
    for smi in frags:
      self.assertTrue(_of(smi) in expected)

  def testFragmentation2(self):
    mol = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')

    frags = FraggleSim.generate_fraggle_fragmentation(mol)
    self.assertEqual(len(frags), 13)

    expected = (
      '*C(=O)c1ccccc1.*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1',
      '*C(=O)c1ccccc1.*Cc1cc(OC)c2ccccc2c1OC', '*C(=O)c1ccccc1.*c1cc(OC)c2ccccc2c1OC',
      '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.*c1ccccc1',
      '*Cc1cc(OC)c2ccccc2c1OC.*NC(=O)c1ccccc1', '*Cc1cc(OC)c2ccccc2c1OC.*c1ccccc1',
      '*N1CCC(NC(=O)c2ccccc2)CC1.*c1cc(OC)c2ccccc2c1OC',
      '*NC(=O)c1ccccc1.*c1cc(OC)c2ccccc2c1OC',
      '*NC1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1.*c1ccccc1',
      '*c1c(CN2CCC(NC(=O)c3ccccc3)CC2)cc(OC)c2ccccc12',
      '*c1c(OC)cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c1*',
      '*c1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12', '*c1cc(OC)c2ccccc2c1OC.*c1ccccc1')
    expected = [_of(s) for s in expected]
    for smi in frags:
      self.assertIn(_of(smi), expected)

  def test_select_fragments(self):
    self.assertRaises(NotImplementedError, select_fragments, 'CCCC', 'InvalidFragmentationType', 10)

  def test_select_fragments_acyclic(self):
    # acyclic fragments: returns fragments with more than two heavy atoms
    self._testSelectFragment('*C.*O*.*c1cccc2ccccc12', FraggleSim.FTYPE_ACYCLIC, 12,
                             '*c1cccc2ccccc12')
    # ignore the fragments if less than 60% of original molecule
    self._testSelectFragment('*CC.*O*.*c1ccccc1', FraggleSim.FTYPE_ACYCLIC, 9,
                             '*c1ccccc1')
    self._testSelectFragment('*CC.*OCC*.*c1ccccc1', FraggleSim.FTYPE_ACYCLIC, 11,
                             '*c1ccccc1')
    self._testSelectFragment('*CC.*OCCC*.*c1ccccc1', FraggleSim.FTYPE_ACYCLIC, 12, None)
    # ignore fragments that are smaller than 2 atoms
    self._testSelectFragment('*C.*O*.*c1ccccc1', FraggleSim.FTYPE_ACYCLIC, 8, '*c1ccccc1')
    self._testSelectFragment('*CC.*O*.*c1ccccc1', FraggleSim.FTYPE_ACYCLIC, 9,
                             '*c1ccccc1')
    self._testSelectFragment('*CCC.*O*.*c1ccccc1', FraggleSim.FTYPE_ACYCLIC, 10,
                             '*CCC.*c1ccccc1')
    # Only small fragments
    self._testSelectFragment('*CC.*CC.*CC', FraggleSim.FTYPE_ACYCLIC, 6, None)
    self._testSelectFragment('*CCC.*CCC.*CCC', FraggleSim.FTYPE_ACYCLIC, 6,
                             '*CCC.*CCC.*CCC')

  def test_select_fragments_cyclic(self):
    self._testSelectFragment('*c1ccccc1*.*cccc*', FraggleSim.FTYPE_CYCLIC, 10,
                             '*c1ccccc1*')
    self._testSelectFragment('*c1ccccc1*.*C.*C', FraggleSim.FTYPE_CYCLIC, 8, None)
    # Fragment too small
    self._testSelectFragment('*c1ccccc1*.*c(CCCCCCCCCC)ccc*', FraggleSim.FTYPE_CYCLIC, 20,
                             None)

  def test_select_fragments_cyclic_acyclic(self):
    self._testSelectFragment('*c1ccccc1*.*CCC.*C', FraggleSim.FTYPE_CYCLIC_ACYCLIC, 10,
                             '*c1ccccc1*.*CCC')
    self._testSelectFragment('*CCC.*C.*c1ccccc1*', FraggleSim.FTYPE_CYCLIC_ACYCLIC, 10,
                             '*CCC.*c1ccccc1*')
    self._testSelectFragment('*CCC.*CCCCCCCCCCC*.*c1ccccc1*',
                             FraggleSim.FTYPE_CYCLIC_ACYCLIC, 25, None)

  def _testSelectFragment(self, smiles, fragmentType, numAtoms, expected):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    self.assertEqual(_of(select_fragments(fragments, fragmentType, numAtoms)), _of(expected))

  def test_isValidRingCut(self):
    rdBase.DisableLog('rdApp.error')
    self.assertEqual(FraggleSim.isValidRingCut(Chem.MolFromSmiles('*CCC*')), False)
    self.assertEqual(FraggleSim.isValidRingCut(Chem.MolFromSmiles('*C1CC1*')), True)
    self.assertEqual(FraggleSim.isValidRingCut(Chem.MolFromSmiles('*c1ccccc1*')), True)
    self.assertEqual(
      FraggleSim.isValidRingCut(Chem.MolFromSmiles('*cccc*', sanitize=False)), False)
    rdBase.EnableLog('rdApp.error')

  def test_GetFraggleSimilarity(self):
    q = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3cncc(C)c3)CC2)c(OC)c2ccccc12')
    m = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.980, places=2)
    self.assertEqual(match, '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1')

    m = Chem.MolFromSmiles('COc1cc(CN2CCC(Nc3nc4ccccc4s3)CC2)c(OC)c2ccccc12')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.794, places=2)
    self.assertEqual(match, '*C1CCN(Cc2cc(OC)c3ccccc3c2OC)CC1')

    q = Chem.MolFromSmiles('COc1ccccc1')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.347, places=2)
    self.assertEqual(match, '*c1ccccc1')

    m = Chem.MolFromSmiles('COc1cc(CN2CCC(NC(=O)c3ccccc3)CC2)c(OC)c2ccccc12')
    sim, match = FraggleSim.GetFraggleSimilarity(q, m)
    self.assertAlmostEqual(sim, 0.266, places=2)
    self.assertEqual(match, '*c1ccccc1')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
