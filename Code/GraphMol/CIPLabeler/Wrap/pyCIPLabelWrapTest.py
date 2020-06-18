import unittest

from rdkit import Chem


class TestCase(unittest.TestCase):

  def testLabelAtomsList(self):
    mol = Chem.MolFromSmiles("C[C@H](Cl)CC[C@H](Cl)C")
    self.assertIsNotNone(mol)

    atom1 = mol.GetAtomWithIdx(1)
    atom5 = mol.GetAtomWithIdx(5)

    self.assertTrue(atom1.HasProp("_CIPCode"))
    self.assertTrue(atom5.HasProp("_CIPCode"))

    atom1.ClearProp("_CIPCode")
    atom5.ClearProp("_CIPCode")

    Chem.rdCIPLabeler.AssignCIPLabels(mol, [1], None)

    self.assertEqual(atom1.GetProp("_CIPCode"), "S")
    self.assertFalse(atom5.HasProp("_CIPCode"))

  def testLabelAtomsSet(self):
    mol = Chem.MolFromSmiles("C[C@H](Cl)CC[C@H](Cl)C")
    self.assertIsNotNone(mol)

    atom1 = mol.GetAtomWithIdx(1)
    atom5 = mol.GetAtomWithIdx(5)

    self.assertTrue(atom1.HasProp("_CIPCode"))
    self.assertTrue(atom5.HasProp("_CIPCode"))

    atom1.ClearProp("_CIPCode")
    atom5.ClearProp("_CIPCode")

    Chem.rdCIPLabeler.AssignCIPLabels(mol, {1}, None)

    self.assertEqual(atom1.GetProp("_CIPCode"), "S")
    self.assertFalse(atom5.HasProp("_CIPCode"))

  def testLabelBondsList(self):
    mol = Chem.MolFromSmiles(r"C\C=C\C=C/C")
    self.assertIsNotNone(mol)

    bond1 = mol.GetBondWithIdx(1)
    bond3 = mol.GetBondWithIdx(3)

    self.assertEqual(bond1.GetBondType(), Chem.BondType.DOUBLE)
    self.assertEqual(bond3.GetBondType(), Chem.BondType.DOUBLE)

    self.assertFalse(bond1.HasProp("_CIPCode"))
    self.assertFalse(bond3.HasProp("_CIPCode"))

    Chem.rdCIPLabeler.AssignCIPLabels(mol, None, [3])

    self.assertFalse(bond1.HasProp("_CIPCode"))
    self.assertEqual(bond3.GetProp("_CIPCode"), "Z")


if __name__ == '__main__':
  unittest.main()
