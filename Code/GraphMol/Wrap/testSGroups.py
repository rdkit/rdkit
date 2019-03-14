#
#  Copyright (C) 2019  Greg Landrum
#         All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

from rdkit import Chem
from rdkit.Chem import RDConfig
import os
import sys
import unittest


class TestCase(unittest.TestCase):
    def setUp(self):
        mb = """
  ACCLDraw04231812452D

 13 13  0  0  0  0  0  0  0  0999 V2000
    5.0960   -4.3327    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4016   -4.3321    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0960   -5.6718    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2539   -6.3314    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2510   -2.3329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    7.4056   -1.6662    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0963   -1.6662    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    6.2510   -3.6663    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5505   -6.3291    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    8.5505   -7.6624    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    9.7052   -5.6625    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    7.4016   -5.6658    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.9416   -2.3329    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  3  4  2  0  0  0  0
  1  3  1  0  0  0  0
  5  8  1  0  0  0  0
  5  7  2  0  0  0  0
  5  6  1  0  0  0  0
  1  8  2  0  0  0  0
  2  8  1  0  0  0  0
  9 12  1  0  0  0  0
  9 11  2  0  0  0  0
  9 10  1  0  0  0  0
  4 12  1  0  0  0  0
  2 12  2  0  0  0  0
  7 13  1  0  0  0  0
M  STY  2   1 DAT   2 DAT
M  SLB  2   1   1   2   2
M  SAL   1  1  10
M  SDT   1 pH                                                    
M  SDD   1     0.0000    0.0000    DR    ALL  1       6
M  SED   1 4.6
M  SAL   2  2   5   7
M  SBL   2  1   4
M  SDT   2 Stereo                                                
M  SDD   2     0.0000    0.0000    DR    ALL  1       6
M  SED   2 E/Z unknown
M  END"""
        self.m1 = Chem.MolFromMolBlock(mb)

    def testBasics(self):
        self.assertTrue(self.m1 is not None)
        sgs = Chem.GetMolSGroups(self.m1)
        self.assertEqual(len(sgs), 2)
        self.assertTrue(sgs[0].HasProp("TYPE"))
        self.assertTrue(sgs[1].HasProp("TYPE"))
        self.assertEqual(sgs[0].GetProp("TYPE"), "DAT")
        self.assertEqual(sgs[1].GetProp("TYPE"), "DAT")

        self.assertTrue(sgs[0].HasProp("FIELDNAME"))
        self.assertEqual(sgs[0].GetProp("FIELDNAME"), "pH")

        self.assertEqual(sorted(sgs[0].GetPropNames()), [
                         'DATAFIELDS', 'FIELDDISP', 'FIELDINFO', 'FIELDNAME', 'FIELDTYPE', 'ID', 'QUERYOP', 'QUERYTYPE', 'TYPE'])
        dd = sgs[0].GetPropsAsDict()
        self.assertTrue("TYPE" in dd)
        self.assertEqual(dd["TYPE"], "DAT")
        self.assertTrue("FIELDNAME" in dd)
        self.assertEqual(dd["FIELDNAME"], "pH")

        Chem.ClearMolSGroups(self.m1)
        self.assertEqual(len(Chem.GetMolSGroups(self.m1)), 0)

    def testLifetime(self):
        self.assertTrue(self.m1 is not None)
        mcpy = Chem.Mol(self.m1)
        smi = Chem.MolToSmiles(mcpy)
        sgs = Chem.GetMolSGroups(mcpy)
        self.assertEqual(len(sgs), 2)
        mcpy = None
        parent = sgs[0].GetOwningMol()
        self.assertEqual(smi, Chem.MolToSmiles(parent))
        sgl = list(sgs)
        sgs = None
        parent = sgl[0].GetOwningMol()
        self.assertEqual(smi, Chem.MolToSmiles(parent))


if __name__ == '__main__':
    print("Testing SGroups wrapper")
    unittest.main()
