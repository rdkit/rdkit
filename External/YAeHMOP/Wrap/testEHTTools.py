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

    def test2(self):
        mol = Chem.MolFromMolFile(os.path.join(
            RDConfig.RDBaseDir, "External", "YAeHMOP", "test_data", "benzene.mol"), removeHs=False)
        self.assertEqual(mol.GetNumAtoms(), 12)
        self.assertTrue(rdEHTTools.RunMol(mol))
        cm = rdEHTTools.GetChargeMatrix(mol)
        self.assertEqual(cm.shape, (12, 12))
        self.assertAlmostEqual(cm[0][0], 0.161, places=3)
        self.assertAlmostEqual(cm[1][0], 0.118, places=3)
        self.assertAlmostEqual(cm[11][10], 0.004, places=3)

    def test3(self):
        mol = Chem.MolFromMolFile(os.path.join(
            RDConfig.RDBaseDir, "External", "YAeHMOP", "test_data", "benzene.mol"), removeHs=False)
        self.assertEqual(mol.GetNumAtoms(), 12)
        self.assertTrue(rdEHTTools.RunMol(mol))
        opm = rdEHTTools.GetOverlapPopulationMatrix(mol)
        self.assertEqual(opm.shape, (int(12 * 13 / 2),))
        self.assertAlmostEqual(opm[0], 2.7035, 3)
        self.assertAlmostEqual(opm[2], 2.7035, 3)
        self.assertAlmostEqual(opm[3], -0.0785, 3)


if __name__ == '__main__':
    unittest.main()
