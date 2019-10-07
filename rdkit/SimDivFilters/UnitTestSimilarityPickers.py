#
#  Copyright (C) 2007  greg Landrum
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
from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.SimDivFilters import SimilarityPickers


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(SimilarityPickers, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def getTestData(self):
    dbName = RDConfig.RDTestDatabase
    conn = DbConnect(dbName, 'simple_mols1')
    mols = []
    for smi, ID in conn.GetData():
      mol = Chem.MolFromSmiles(str(smi))
      mol.SetProp('_Name', str(ID))
      mols.append(mol)

    # Calculate fingerprints
    probefps = []
    for mol in mols:
      fp = Chem.RDKFingerprint(mol)
      fp._id = mol.GetProp('_Name')
      probefps.append(fp)
    return probefps

  def test_GenericPicker(self):
    picker = SimilarityPickers.GenericPicker()
    self.assertRaises(NotImplementedError, picker.MakePicks)

  def test_TopNOverallPicker(self):
    probefps = self.getTestData()

    # Find topN matches
    mol = Chem.MolFromSmiles('COC')
    probeFp = Chem.RDKFingerprint(mol)

    # The picker is initialized lazy; calculation is triggered either by len
    picker = SimilarityPickers.TopNOverallPicker(numToPick=2, probeFps=[probeFp], dataSet=probefps)
    self.assertEqual(picker._picks, None)
    self.assertEqual(len(picker), 2)
    self.assertNotEqual(picker._picks, None)
    # or by addressing the elements
    picker = SimilarityPickers.TopNOverallPicker(numToPick=2, probeFps=[probeFp], dataSet=probefps)
    self.assertEqual(picker._picks, None)
    fp, score = picker[0]
    self.assertEqual(fp._id, 'ether-1')
    self.assertEqual(score, 1.0)
    self.assertNotEqual(picker._picks, None)

  def test_SpreadPicker(self):
    probefps = self.getTestData()

    mol = Chem.MolFromSmiles('COC')
    probeFp = Chem.RDKFingerprint(mol)

    # The picker is initialized lazy; calculation is triggered either by len
    picker = SimilarityPickers.SpreadPicker(numToPick=2, probeFps=[probeFp], dataSet=probefps)
    self.assertEqual(picker._picks, None)
    self.assertEqual(len(picker), 2)

    # or by addressing the elements
    picker = SimilarityPickers.SpreadPicker(numToPick=2, probeFps=[probeFp], dataSet=probefps)
    self.assertEqual(picker._picks, None)
    fp, score = picker[0]
    self.assertEqual(fp._id, 'ether-1')
    self.assertEqual(score, 1.0)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
