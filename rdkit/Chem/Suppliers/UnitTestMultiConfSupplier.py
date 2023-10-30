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
        # TODO maybe change these to just use the file containing energy values
        # check default
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'))]], truth)

        # check propertyName = SMILES
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'), "SMILES")]], truth)

        # check propertyName = countLine
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'), "countLine")]],
                         truth)

        # check propertyName using GetProp
        self.assertEqual([m.GetNumConformers() for m in [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'), "ID")]], truth)

        # check correct CID:Property Value pairing
        mols = [x for x in MultiConfSupplier(os.path.join(path, 'read_multi_conf.sdf'), suppliedConfId='CID')]
        self.assertEqual(mols[0].GetConformer(4).GetProp('Energy'), '75.56787014654348')

if __name__ == '__main__':
    unittest.main()
