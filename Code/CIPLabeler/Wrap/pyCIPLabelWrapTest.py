import unittest

from rdkit import Chem


class TestCase(unittest.TestCase):

  def testLabelAtomsList(self):
    mol = Chem.MolFromSmiles("C[C@H](Cl)CC[C@H](Cl)C")
    assert mol

    atom1 = mol.GetAtomWithIdx(1)
    atom5 = mol.GetAtomWithIdx(5)

    assert atom1.HasProp("_CIPCode") == True
    assert atom5.HasProp("_CIPCode") == True

    atom1.ClearProp("_CIPCode")
    atom5.ClearProp("_CIPCode")

    Chem.rdCIPLabeler.AssignCIPLabels(mol, [1], None)

    assert atom1.GetProp("_CIPCode") == "S"
    assert atom5.HasProp("_CIPCode") == False

  def testLabelAtomsSet(self):
    mol = Chem.MolFromSmiles("C[C@H](Cl)CC[C@H](Cl)C")
    assert mol

    atom1 = mol.GetAtomWithIdx(1)
    atom5 = mol.GetAtomWithIdx(5)

    assert atom1.HasProp("_CIPCode") == True
    assert atom5.HasProp("_CIPCode") == True

    atom1.ClearProp("_CIPCode")
    atom5.ClearProp("_CIPCode")

    Chem.rdCIPLabeler.AssignCIPLabels(mol, {1}, None)

    assert atom1.GetProp("_CIPCode") == "S"
    assert atom5.HasProp("_CIPCode") == False

  def testLabelBondsList(self):
    mol = Chem.MolFromSmiles(r"C\C=C\C=C/C")
    assert mol

    bond1 = mol.GetBondWithIdx(1)
    bond3 = mol.GetBondWithIdx(3)

    assert bond1.GetBondType() == Chem.BondType.DOUBLE
    assert bond3.GetBondType() == Chem.BondType.DOUBLE

    assert bond1.HasProp("_CIPCode") == False
    assert bond3.HasProp("_CIPCode") == False

    Chem.rdCIPLabeler.AssignCIPLabels(mol, None, [3])

    assert bond1.HasProp("_CIPCode") == False
    assert bond3.GetProp("_CIPCode") == "Z"


if __name__ == '__main__':
  unittest.main()
