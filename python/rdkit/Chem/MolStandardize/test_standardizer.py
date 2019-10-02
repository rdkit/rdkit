import unittest
from rdkit import Chem
from rdkit.Chem.MolStandardize.standardize import Standardizer


class FakeStandardizer(Standardizer):
    def normalize(self):
        def fake_normalize(y):
            props = y.GetPropsAsDict()
            for k, v in props:
                y.ClearProp(k)
            return y
        return fake_normalize

class TestCase(unittest.TestCase):

    def testPreserveProps(self):
        PROP_NAME = "MyProp"
        PROP_VALUE = "foo"
        standardizer = FakeStandardizer()
        m = Chem.MolFromSmiles("C")
        m.SetProp(PROP_NAME, PROP_VALUE)

        standardized_mol = standardizer.standardize(m)
        self.assertTrue(standardized_mol.HasProp(PROP_NAME))
        self.assertEqual(PROP_VALUE, standardized_mol.GetProp(PROP_NAME))
