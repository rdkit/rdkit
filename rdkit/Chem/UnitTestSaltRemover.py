import doctest
import unittest

import os
from rdkit import Chem

import Chem.SaltRemover
from Chem.SaltRemover import SaltRemover, InputFormat

def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(Chem.SaltRemover, optionflags=doctest.ELLIPSIS))
  return tests

class TestCase(unittest.TestCase):
  def test_withSmiles(self):
    remover = SaltRemover(defnData="[Na+]\nCC(=O)O", defnFormat=InputFormat.SMILES)
    self.assertEqual(len(remover.salts), 2)
    mol = Chem.MolFromSmiles('CC(=O)O.[Na+]')
    res = remover.StripMol(mol)
    self.assertEqual(res.GetNumAtoms(), 0)

  def test_withSdfFile(self):
    testFile = os.sep.join(
      [os.path.dirname(os.path.abspath(__file__)), 'test_data', 'witch-salts.sdf'])
    remover = SaltRemover(defnFilename=testFile, defnFormat=InputFormat.MOL)
    self.assertEqual(len(remover.salts), 240)
    m = Chem.MolFromSmiles("Cc1onc(-c2ccccc2)c1C([O-])=NC1C(=O)N2C1SC(C)(C)C2C(=O)O.O.[Na+]")
    tuple = remover.StripMolWithDeleted(m)
    self.assertEqual(Chem.MolToSmiles(tuple.mol), 'Cc1onc(-c2ccccc2)c1C([O-])=NC1C(=O)N2C1SC(C)(C)C2C(=O)O.O')
    self.assertEqual(len(tuple.deleted), 1)
    self.assertEqual(Chem.MolToSmiles(tuple.deleted[0]), '[Na+]')

  def test_withSmiFile(self):
    testFile = os.sep.join(
      [os.path.dirname(os.path.abspath(__file__)), 'test_data', 'c6h6-cdk.smi'])
    remover = SaltRemover(defnFilename=testFile, defnFormat=InputFormat.SMILES)
    self.assertEqual(len(remover.salts), 216)

if __name__ == '__main__':  # pragma: nocover
  unittest.main()
