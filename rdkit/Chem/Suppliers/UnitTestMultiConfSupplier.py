import unittest
import os
from rdkit import RDConfig
from rdkit.Chem.Suppliers.MultiConfSupplier import MultiConfSupplier


class TestCase(unittest.TestCase):
    def testMultiConfParse(self):
        path = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data')
        # check ValueError raised when no 3D coordinates supplied
        with self.assertRaises(ValueError):
            suppl = MultiConfSupplier(os.path.join(path, 'read_multi_conf_error.sdf'))
            for mol in suppl:
                print(mol)

        # test function returns expected outputs
        truth = [62, 7, 200]

        # check default
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'))]], truth)

        # check propertyName = SMILES
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'), "SMILES")]], truth)

        # check propertyName = countLine
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'), "countLine")]],
                         truth)

        # check propertyName using GetProp
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'), "ID")]], truth)


if __name__ == '__main__':
    unittest.main()
