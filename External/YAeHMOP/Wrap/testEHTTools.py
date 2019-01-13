#  Copyright (c) 2019 Greg Landrum
#  All rights reserved.
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
from rdkit import RDConfig
import os
import sys
import unittest
from rdkit import DataStructs, Chem
from rdkit.Chem import rdEHTTools


class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        mol = Chem.MolFromMolFile(os.path.join(
            RDConfig.RDBaseDir, "External", "YAeHMOP", "test_data", "benzene.mol"), removeHs=False)
        self.assertEqual(mol.GetNumAtoms(), 12)
        self.assertTrue(rdEHTTools.RunMol(mol))
        self.assertAlmostEqual(mol.GetAtomWithIdx(1).GetDoubleProp("_EHTCharge"), -0.026, places=3)
        self.assertAlmostEqual(mol.GetAtomWithIdx(7).GetDoubleProp("_EHTCharge"), 0.026, places=3)


if __name__ == '__main__':
    unittest.main()
