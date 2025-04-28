import unittest

from rdkit import Chem

class TestCase(unittest.TestCase):

  def testLabelAtomsList(self):
    mol = Chem.MolFromSmiles("C[C@H](Cl)CC[C@H](Cl)C")
    self.assertIsNotNone(mol)

    atom1 = mol.GetAtomWithIdx(1)
    atom5 = mol.GetAtomWithIdx(5)

    self.assertTrue(atom1.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertTrue(atom5.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)

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

    self.assertTrue(atom1.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertTrue(atom5.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED)

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

  def doOneAtropIomerMandP(self,  inputSmiles , expected):

    ps = Chem.SmilesParserParams()
    ps.allowCXSMILES = True
    ps.parseName = False
    ps.sanitize = True
    ps.removeHs = False

    mol = Chem.MolFromSmiles(inputSmiles, ps)

    self.assertIsNotNone(mol)
    Chem.rdCIPLabeler.AssignCIPLabels(mol)

    cipsCodes =""
    for bondIndex in range(mol.GetNumBonds()):

      bond = mol.GetBondWithIdx(bondIndex)

      if (bond.HasProp("_CIPCode")):
        cipsCodes += str(bondIndex) + bond.GetProp("_CIPCode") + ":" 

    self.assertEqual(cipsCodes, expected)

  def testAtropIsomer(self):
    mol = "FC1=C(C2=C(C)C(N3C(=O)C4=C(N(C)C3=O)C(F)=CC=C4)=CC=C2)C2=C(NC3=C2CC[C@H](C(O)(C)C)C3)C(C(=O)N)=C1 |(2.1158,0.5489,;1.4029,0.9642,;0.6554,0.5402,;0.6459,-0.2846,;1.3556,-0.7053,;2.0747,-0.3011,;1.346,-1.5302,;2.3682,-1.8618,;2.9748,-1.3027,;2.794,-0.4978,;3.7623,-1.5486,;3.9431,-2.3535,;3.3364,-2.9126,;3.5173,-3.7175,;2.549,-2.6667,;1.9423,-3.2258,;4.7594,-2.6222,;4.9309,-3.4292,;5.3958,-2.0448,;5.2074,-1.2063,;4.3851,-0.9565,;0.6268,-1.9345,;-0.0827,-1.5137,;-0.0732,-0.6889,;-0.0819,0.9813,;-0.0819,1.8063,;-0.8665,2.0612,;-1.3515,1.3938,;-0.8665,0.7264,;-1.2039,-0.0639,;-2.0578,-0.1602,;-2.5629,0.5349,;-3.3837,0.4518,;-4.2045,0.3687,;-3.4668,1.2725,;-3.3006,-0.369,;-2.2074,1.3172,;0.6555,2.2474,;0.6459,3.0723,;-0.0731,3.4765,;1.3556,3.4931,;1.403,1.8235,),wU:7.14wD:31.35|"
    self.doOneAtropIomerMandP(mol, "6P:")

    mol = "C1(N2C(C)=CC=C2Br)=C(C)C(C)=C(N2C(C)=CC=C2Br)C(C)=C1C |(-0.0002,1.5403,;-0.0002,3.0805,;-1.334,3.8508,;-2.6678,3.0807,;-1.334,5.391,;1.3338,5.391,;1.3338,3.8508,;2.6676,3.0807,;-1.3338,0.7702,;-2.6678,1.5403,;-1.3338,-0.7702,;-2.6678,-1.5401,;-0.0002,-1.5403,;-0.0002,-3.0805,;1.3338,-3.8508,;2.6676,-3.0805,;1.3338,-5.391,;-1.334,-5.391,;-1.334,-3.8508,;-2.6678,-3.0805,;1.3338,-0.7702,;2.6678,-1.5403,;1.3338,0.7702,;2.6678,1.5404,),wU:1.6,13.14|"
    self.doOneAtropIomerMandP(mol, "0m:12m:")

    mol = "N1(n2c(C)ccc2Br)C(=O)[C@H](C)[C@H](C)C1=O |(-11.1517,1.8306,;-11.1517,3.3708,;-12.4855,4.1411,;-13.8193,3.371,;-12.4855,5.6813,;-9.8177,5.6813,;-9.8177,4.1411,;-8.4839,3.371,;-12.3975,0.9252,;-13.8622,1.4011,;-11.9217,-0.5394,;-12.8269,-1.7852,;-10.3817,-0.5394,;-9.4765,-1.7852,;-9.9059,0.9252,;-8.4413,1.4011,),wU:0.8,10.11,12.13|"
    self.doOneAtropIomerMandP(mol, "0p:")

if __name__ == '__main__':
  unittest.main()
