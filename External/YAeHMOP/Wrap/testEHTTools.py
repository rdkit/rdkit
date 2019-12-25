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
        ok, res = rdEHTTools.RunMol(mol)
        self.assertTrue(ok)
        chgs = res.GetAtomicCharges()
        self.assertAlmostEqual(chgs[1], -0.026, places=3)
        self.assertAlmostEqual(chgs[7], 0.026, places=3)

    def test2(self):
        mol = Chem.MolFromMolFile(os.path.join(
            RDConfig.RDBaseDir, "External", "YAeHMOP", "test_data", "benzene.mol"), removeHs=False)
        self.assertEqual(mol.GetNumAtoms(), 12)
        ok, res = rdEHTTools.RunMol(mol)
        self.assertTrue(ok)
        cm = res.GetReducedChargeMatrix()
        self.assertEqual(cm.shape, (12, res.numOrbitals))
        for i in range(6):
            self.assertAlmostEqual(cm[i][0], 0.161, places=3)
            self.assertAlmostEqual(cm[i+6][0], 0.005, places=3)
            self.assertAlmostEqual(cm[i][6], 0.1066, places=3)
            self.assertAlmostEqual(cm[i+6][6], 0.060, places=3)
            self.assertAlmostEqual(cm[i][9], 0.167, places=3)
            self.assertAlmostEqual(cm[i+6][9], 0.000, places=3)
        

    def test3(self):
        mol = Chem.MolFromMolFile(os.path.join(
            RDConfig.RDBaseDir, "External", "YAeHMOP", "test_data", "benzene.mol"), removeHs=False)
        self.assertEqual(mol.GetNumAtoms(), 12)
        ok, res = rdEHTTools.RunMol(mol)
        self.assertTrue(ok)
        opm = res.GetReducedOverlapPopulationMatrix()
        self.assertEqual(opm.shape, (int(12 * 13 / 2),))
        self.assertAlmostEqual(opm[0], 2.7035, 3)
        self.assertAlmostEqual(opm[2], 2.7035, 3)
        self.assertAlmostEqual(opm[3], -0.0785, 3)

    def test1(self):
        mol = Chem.MolFromMolFile(os.path.join(
            RDConfig.RDBaseDir, "External", "YAeHMOP", "test_data", "benzene.mol"), removeHs=False)
        self.assertEqual(mol.GetNumAtoms(), 12)
        ok, res = rdEHTTools.RunMol(mol)
        self.assertTrue(ok)
        self.assertAlmostEqual(res.fermiEnergy, -12.804, places=3)
        self.assertAlmostEqual(res.totalEnergy, -535.026, places=3)


if __name__ == '__main__':
    unittest.main()
