#
#  Copyright (C) 2003  Greg Landrum and Rational Discovery LLC
#
"""unit testing code for the SimilarityScreeners

"""
import unittest

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import SimilarityScreener


def fingerprinter(mol):
  return Chem.RDKFingerprint(mol, minPath=2, maxPath=7, fpSize=2048)


class TestCase(unittest.TestCase):

  def test_SimilarityScreener(self):
    probe = fingerprinter(Chem.MolFromSmiles('C1OCCCC1'))

    screener = SimilarityScreener.SimilarityScreener()
    self.assertEqual(screener.probe, None)
    screener.SetProbe(probe)
    self.assertEqual(screener.probe, probe)

    screener = SimilarityScreener.SimilarityScreener(probe=probe, fingerprinter=fingerprinter)
    self.assertEqual(screener.probe, probe)

    mol = Chem.MolFromSmiles('CCCN')
    self.assertEqual(fingerprinter(mol), screener.GetSingleFingerprint(mol))

    screener.Reset()

  def test1_TopNScreener(self):
    smis = ['C1CCCCC1', 'C1OCCCC1', 'C1NCCCC1', 'c1ccccc1', 'C1C(C)CCCC1', 'C1C(C)C(C)CCC1']
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(smis), delimiter=",", smilesColumn=0,
                                           nameColumn=-1, titleLine=0)
    metric = DataStructs.TanimotoSimilarity
    mol = Chem.MolFromSmiles('C1OCCCC1')
    probe = fingerprinter(mol)

    screener = SimilarityScreener.TopNScreener(3, probe=probe, metric=metric,
                                               fingerprinter=fingerprinter, dataSource=suppl)
    self.assertIsNone(screener.topN)
    matches1 = [x for x in screener]
    self.assertIsNotNone(screener.topN)

    self.assertEqual(len(matches1), 3)
    matches2 = [x for x in screener]
    self.assertEqual(len(matches2), 3)
    self.assertEqual(matches1, matches2)

    self.assertEqual(probe, screener.GetSingleFingerprint(mol))

    # Getting the length also triggers the execution of the screen
    screener = SimilarityScreener.TopNScreener(3, probe=probe, metric=metric,
                                               fingerprinter=fingerprinter, dataSource=suppl)
    self.assertIsNone(screener.topN)
    self.assertEqual(len(screener), 3)
    self.assertIsNotNone(screener.topN)

    # as does accessing elements by index
    screener = SimilarityScreener.TopNScreener(3, probe=probe, metric=metric,
                                               fingerprinter=fingerprinter, dataSource=suppl)
    self.assertIsNone(screener.topN)
    screener[1]
    self.assertIsNotNone(screener.topN)

  def test2_ThresholdScreener(self):
    smis = ['C1CCCCC1', 'C1OCCCC1', 'C1NCCCC1', 'c1ccccc1', 'C1C(C)CCCC1', 'C1C(C)C(C)CCC1']
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(smis), delimiter=",", smilesColumn=0,
                                           nameColumn=-1, titleLine=0)

    metric = DataStructs.TanimotoSimilarity
    probe = fingerprinter(Chem.MolFromSmiles('C1OCCCC1'))

    screener = SimilarityScreener.ThresholdScreener(0.09, probe=probe, metric=metric,
                                                    fingerprinter=fingerprinter, dataSource=suppl)
    matches1 = [x[0] for x in screener]
    self.assertEqual(len(matches1), 5)
    matches2 = [x[0] for x in screener]
    self.assertEqual(len(matches2), 5)
    self.assertEqual(matches1, matches2)

  def test3_ThresholdScreener_folding(self):
    smis = ['C1CCCCC1', 'C1OCCCC1', 'C1NCCCC1', 'c1ccccc1', 'C1C(C)CCCC1', 'C1C(C)C(C)CCC1']
    suppl = Chem.SmilesMolSupplierFromText('\n'.join(smis), delimiter=",", smilesColumn=0,
                                           nameColumn=-1, titleLine=0)

    metric = DataStructs.TanimotoSimilarity
    probe = Chem.RDKFingerprint(Chem.MolFromSmiles('C1OCCCC1'), minPath=2, maxPath=7, fpSize=4096)

    screener = SimilarityScreener.ThresholdScreener(0.09, probe=probe, metric=metric,
                                                    fingerprinter=fingerprinter, dataSource=suppl)
    matches1 = [x[0] for x in screener]
    self.assertEqual(len(matches1), 5)
    matches2 = [x[0] for x in screener]
    self.assertEqual(len(matches2), 5)
    self.assertEqual(matches1, matches2)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
