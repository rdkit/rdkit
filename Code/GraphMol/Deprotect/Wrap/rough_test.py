import unittest
from rdkit import Chem
from rdkit.Chem import rdDeprotect as rd


class TestCase(unittest.TestCase):

  def testDeprotect(self):
    smiles = "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C"
    m = Chem.MolFromSmiles(smiles)
    m2 = rd.Deprotect(m)
    self.assertEqual(Chem.MolToSmiles(m2), Chem.CanonSmiles("NCc1ccccc1N"))
    self.assertEqual(list(m2.GetPropsAsDict()["DEPROTECTIONS"]), ["Boc", "Boc"])
    self.assertEqual(m2.GetPropsAsDict()["DEPROTECTION_COUNT"], 2)

  def testDeprotectInPlace(self):
    smiles = "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C"
    m = Chem.MolFromSmiles(smiles)
    self.assertTrue(rd.DeprotectInPlace(m))
    self.assertEqual(Chem.MolToSmiles(m), Chem.CanonSmiles("NCc1ccccc1N"))
    self.assertEqual(list(m.GetPropsAsDict()["DEPROTECTIONS"]), ["Boc", "Boc"])
    self.assertEqual(m.GetPropsAsDict()["DEPROTECTION_COUNT"], 2)

  def testDeprotectVector(self):
    reaction_class = "amine"
    reaction_smarts = "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]"
    abbreviation = "Boc"
    full_name = "tert-butyloxycarbonyl"
    data = rd.DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)

    smiles = "N(C(=O)OC(C)(C)C)Cc1ccccc1NC(=O)OC(C)(C)C"
    m = Chem.MolFromSmiles(smiles)
    m2 = rd.Deprotect(m, [])
    self.assertEqual(Chem.MolToSmiles(m), Chem.MolToSmiles(m2))
    self.assertEqual(m2.GetPropsAsDict()["DEPROTECTION_COUNT"], 0)

    m2 = rd.Deprotect(m, [data])
    self.assertEqual(Chem.MolToSmiles(m2), Chem.CanonSmiles("NCc1ccccc1N"))
    self.assertEqual(list(m2.GetPropsAsDict()["DEPROTECTIONS"]), ["Boc", "Boc"])
    self.assertEqual(m2.GetPropsAsDict()["DEPROTECTION_COUNT"], 2)

  def test_examples(self):
    count = 0
    for data in rd.GetDeprotections():
      if data.example:
        start, end = data.example.split(">>")
        m = Chem.MolFromSmiles(start)
        m2 = rd.Deprotect(m, [data])
        print("Testing", data.full_name)
        self.assertEqual(Chem.MolToSmiles(m2), Chem.CanonSmiles(end))
        count += 1
    assert count

  def test_examples_inplace(self):
    count = 0
    for data in rd.GetDeprotections():
      if data.example:
        start, end = data.example.split(">>")
        m = Chem.MolFromSmiles(start)
        rd.DeprotectInPlace(m, [data])
        print("Testing", data.full_name)
        self.assertEqual(Chem.MolToSmiles(m), Chem.CanonSmiles(end))
        count += 1
    assert count


if __name__ == "__main__":
  unittest.main()
