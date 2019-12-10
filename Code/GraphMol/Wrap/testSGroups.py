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
from rdkit import Geometry
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
    sgs = Chem.GetMolSubstanceGroups(self.m1)
    self.assertEqual(len(sgs), 2)
    self.assertTrue(sgs[0].HasProp("TYPE"))
    self.assertTrue(sgs[1].HasProp("TYPE"))
    self.assertEqual(sgs[0].GetProp("TYPE"), "DAT")
    self.assertEqual(sgs[1].GetProp("TYPE"), "DAT")

    self.assertTrue(sgs[0].HasProp("FIELDNAME"))
    self.assertEqual(sgs[0].GetProp("FIELDNAME"), "pH")

    self.assertEqual(sorted(sgs[0].GetPropNames()), [
      'DATAFIELDS', 'FIELDDISP', 'FIELDINFO', 'FIELDNAME', 'FIELDTYPE', 'ID', 'QUERYOP',
      'QUERYTYPE', 'TYPE'
    ])
    dd = sgs[0].GetPropsAsDict()
    self.assertTrue("TYPE" in dd)
    self.assertEqual(dd["TYPE"], "DAT")
    self.assertTrue("FIELDNAME" in dd)
    self.assertEqual(dd["FIELDNAME"], "pH")

    Chem.ClearMolSubstanceGroups(self.m1)
    self.assertEqual(len(Chem.GetMolSubstanceGroups(self.m1)), 0)

  def testLifetime(self):
    self.assertTrue(self.m1 is not None)
    mcpy = Chem.Mol(self.m1)
    smi = Chem.MolToSmiles(mcpy)
    sgs = Chem.GetMolSubstanceGroups(mcpy)
    self.assertEqual(len(sgs), 2)
    mcpy = None
    parent = sgs[0].GetOwningMol()
    self.assertEqual(smi, Chem.MolToSmiles(parent))
    sgl = list(sgs)
    sgs = None
    parent = sgl[0].GetOwningMol()
    self.assertEqual(smi, Chem.MolToSmiles(parent))

  def testSetProp(self):
    self.assertTrue(self.m1 is not None)
    mcpy = Chem.Mol(self.m1)
    sg = Chem.GetMolSubstanceGroupWithIdx(mcpy, 0)
    sg.SetProp("foo", "bar")
    sgs2 = Chem.GetMolSubstanceGroups(mcpy)
    pd = sgs2[0].GetPropsAsDict()
    self.assertTrue('foo' in pd)
    self.assertEqual(pd['foo'], 'bar')

  def testCreateSGroup(self):
    mol = Chem.MolFromMolBlock('''
  Mrv1810 10061910532D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 O -14.8632 3.7053 0 0
M  V30 2 C -13.3232 3.7053 0 0
M  V30 3 O -11.7832 3.7053 0 0
M  V30 4 * -10.2453 3.6247 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 END BOND
M  V30 END CTAB
M  END''')
    self.assertTrue(mol is not None)
    mcpy = Chem.Mol(mol)
    sg = Chem.CreateMolSubstanceGroup(mcpy, "SRU")
    self.assertEqual(sg.GetProp("TYPE"), "SRU")
    sg = Chem.GetMolSubstanceGroupWithIdx(mcpy, 0)
    sg.AddBracket((Geometry.Point3D(-11.1, 4.6, 0), Geometry.Point3D(-11.1, 2.7, 0)))
    sg.AddBracket((Geometry.Point3D(-13.9, 2.7, 0), Geometry.Point3D(-13.9, 4.6, 0)))
    sg.AddAtomWithIdx(1)
    sg.AddAtomWithIdx(2)
    sg.AddBondWithIdx(0)
    sg.AddBondWithIdx(2)
    sg.SetProp("CONNECT", "HT")
    sg.SetProp("LABEL", "n")
    mb = Chem.MolToMolBlock(mcpy, forceV3000=True)
    self.assertNotEqual(mb.find('V30 1 SRU'), -1)
    self.assertNotEqual(mb.find('BRKXYZ'), -1)
    self.assertNotEqual(mb.find('CONNECT=HT'), -1)


if __name__ == '__main__':
  print("Testing SubstanceGroups wrapper")
  unittest.main()
