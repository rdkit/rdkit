import doctest
import os
import unittest

from rdkit import Chem
from rdkit.Chem.SaltRemover import InputFormat, SaltRemover


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
    self.assertEqual(Chem.MolToSmiles(tuple.mol),
                     'Cc1onc(-c2ccccc2)c1C([O-])=NC1C(=O)N2C1SC(C)(C)C2C(=O)O.O')
    self.assertEqual(len(tuple.deleted), 1)
    self.assertEqual(Chem.MolToSmiles(tuple.deleted[0]), '[Na+]')

  def test_withSmiFile(self):
    testFile = os.sep.join(
      [os.path.dirname(os.path.abspath(__file__)), 'test_data', 'c6h6-cdk.smi'])
    remover = SaltRemover(defnFilename=testFile, defnFormat=InputFormat.SMILES)
    self.assertEqual(len(remover.salts), 216)

  def test_withDontRemoveEverything(self):
    testFile = os.sep.join(
      [os.path.dirname(os.path.abspath(__file__)), 'test_data', 'witch-salts.sdf'])
    remover = SaltRemover(defnFilename=testFile, defnFormat=InputFormat.MOL)
    m = Chem.MolFromSmiles('Cc1ccccc1')
    mol, deleted = remover.StripMolWithDeleted(m, dontRemoveEverything=True)
    # List should be empty
    self.assertFalse(deleted)
    self.assertEqual(m, mol)

  def test_SmilesVsSmarts(self):
    # SMARTS
    remover = SaltRemover(defnData="[Cl,Br]")
    mol = Chem.MolFromSmiles('CN(Br)Cl.Cl')
    res = remover.StripMol(mol)
    self.assertEqual(res.GetNumAtoms(), 4)
    self.assertEqual(Chem.MolToSmiles(res), 'CN(Cl)Br')
    mol = Chem.MolFromSmiles('CN(C)C.Cl.Br')
    res, deleted = remover.StripMolWithDeleted(mol)
    self.assertEqual(Chem.MolToSmiles(res), 'CN(C)C')
    # Because we read in SMARTS, we should output as well. Otherwise, we will have
    # mismatches
    self.assertListEqual([Chem.MolToSmarts(m) for m in deleted], ['[Cl,Br]'])
    # SMILES
    remover = SaltRemover(defnData="Cl", defnFormat=InputFormat.SMILES)
    mol = Chem.MolFromSmiles('CN(Br)Cl.Cl')
    res = remover.StripMol(mol)
    self.assertEqual(res.GetNumAtoms(), 4)
    self.assertEqual(Chem.MolToSmiles(res), 'CN(Cl)Br')

  def test_github_4550(self):
    m = Chem.MolFromSmiles('Cl.C[N]1=CC=CC=C1', sanitize=False)
    self.assertEqual(m.GetNumAtoms(), 8)
    saltstrip = SaltRemover()
    res = saltstrip.StripMol(m, sanitize=False)
    self.assertEqual(Chem.MolToSmiles(res), 'CN1=CC=CC=C1')

  def test_github_7327(self):
    m = Chem.MolFromSmiles('C=CC=O')
    assert m

    saltstrip = SaltRemover()
    m = saltstrip.StripMol(m)

    # No atoms removed
    assert m.GetNumAtoms() == 4

    # Rotatable bond definition from mmpdb
    q = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]')
    assert q

    assert m.GetSubstructMatches(q) == ((1, 2), )


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
