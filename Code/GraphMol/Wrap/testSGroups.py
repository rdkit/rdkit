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
    mb3000 = '''
  Mrv2102 09042105562D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 13 13 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 5.8858 -5.0042 0 0
M  V30 2 C 8.5487 -5.0035 0 0
M  V30 3 C 5.8858 -6.5508 0 0
M  V30 4 C 7.2232 -7.3127 0 0
M  V30 5 C 7.2198 -2.6945 0 0
M  V30 6 C 8.5534 -1.9244 0 0
M  V30 7 C 5.8862 -1.9244 0 0
M  V30 8 C 7.2198 -4.2345 0 0
M  V30 9 C 9.8757 -7.31 0 0
M  V30 10 O 9.8757 -8.85 0 0
M  V30 11 O 11.2094 -6.5401 0 0
M  V30 12 C 8.5487 -6.5439 0 0
M  V30 13 C 4.5525 -2.6945 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 3 4
M  V30 2 1 1 3
M  V30 3 1 5 8
M  V30 4 2 5 7
M  V30 5 1 5 6
M  V30 6 2 1 8
M  V30 7 1 2 8
M  V30 8 1 9 12
M  V30 9 2 9 11
M  V30 10 1 9 10
M  V30 11 1 4 12
M  V30 12 2 2 12
M  V30 13 1 7 13
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 10) FIELDNAME=pH -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  1       6" FIELDDATA=4.6
M  V30 2 DAT 0 ATOMS=(2 5 7) FIELDNAME=Stereo -
M  V30 FIELDDISP="    0.0000    0.0000    DR    ALL  1       6" -
M  V30 FIELDDATA="E/Z unknown"
M  V30 END SGROUP
M  V30 END CTAB
M  END
'''
    self.m2 = Chem.MolFromMolBlock(mb3000)

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

    self.assertEqual(sorted(sgs[0].GetPropNames()),
                     ['DATAFIELDS', 'FIELDDISP', 'FIELDNAME', 'ID', 'TYPE', 'index'])
    dd = sgs[0].GetPropsAsDict()
    self.assertTrue("TYPE" in dd)
    self.assertEqual(dd["TYPE"], "DAT")
    self.assertTrue("FIELDNAME" in dd)
    self.assertEqual(dd["FIELDNAME"], "pH")

    sgs[0].ClearProp("FIELDNAME")
    self.assertFalse(sgs[0].HasProp("FIELDNAME"))

    # The property doesn't exist anymore, but this should not fail
    sgs[0].ClearProp("FIELDNAME")
    self.assertFalse(sgs[0].HasProp("FIELDNAME"))

    Chem.ClearMolSubstanceGroups(self.m1)
    self.assertEqual(len(Chem.GetMolSubstanceGroups(self.m1)), 0)

  def testBasics3000(self):
    self.assertTrue(self.m2 is not None)
    sgs = Chem.GetMolSubstanceGroups(self.m2)
    self.assertEqual(len(sgs), 2)
    self.assertTrue(sgs[0].HasProp("TYPE"))
    self.assertTrue(sgs[1].HasProp("TYPE"))
    self.assertEqual(sgs[0].GetProp("TYPE"), "DAT")
    self.assertEqual(sgs[1].GetProp("TYPE"), "DAT")

    self.assertTrue(sgs[0].HasProp("FIELDNAME"))
    self.assertEqual(sgs[0].GetProp("FIELDNAME"), "pH")

    self.assertEqual(sorted(sgs[0].GetPropNames()),
                     ['DATAFIELDS', 'FIELDDISP', 'FIELDNAME', 'TYPE', 'index'])
    dd = sgs[0].GetPropsAsDict()
    self.assertTrue("TYPE" in dd)
    self.assertEqual(dd["TYPE"], "DAT")
    self.assertTrue("FIELDNAME" in dd)
    self.assertEqual(dd["FIELDNAME"], "pH")

    Chem.ClearMolSubstanceGroups(self.m2)
    self.assertEqual(len(Chem.GetMolSubstanceGroups(self.m2)), 0)

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

  def testXBONDS(self):
    mol = Chem.MolFromMolBlock('''
  Mrv1824 06192020192D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -1.25 0.7484 0 0
M  V30 2 C -2.5837 -0.0216 0 0
M  V30 3 C -2.5837 -1.5617 0 0
M  V30 4 C -1.25 -2.3317 0 0
M  V30 5 C 0.0837 -1.5617 0 0
M  V30 6 C 0.0837 -0.0216 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 1 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(1 1) XBONDS=(2 1 6) XBHEAD=(2 6 1) XBCORR=(4 6 6 1 1) -
M  V30 BRKXYZ=(9 -2.02 -0.0216 0 -2.02 1.5184 0 0 0 0) BRKXYZ=(9 -0.48 1.5184 -
M  V30 0 -0.48 -0.0216 0 0 0 0) CONNECT=HT LABEL="1-3"
M  V30 END SGROUP
M  V30 END CTAB
M  END
''')
    self.assertIsNotNone(mol)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 1)
    sg = sgs[0]
    self.assertEqual(sg.GetProp('TYPE'), 'SRU')
    v = sg.GetUnsignedVectProp('XBHEAD')
    self.assertEqual(len(v), 2)
    self.assertEqual(list(v), [5, 0])
    v = sg.GetUnsignedVectProp('XBCORR')
    self.assertEqual(len(v), 4)
    self.assertEqual(list(v), [5, 5, 0, 0])

  def testDataFields(self):
    mol = Chem.MolFromMolBlock('''
  Mrv2007 06242013252D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.1782 0.0827 0 0
M  V30 2 C -6.5119 -0.6873 0 0
M  V30 3 C -3.8446 -0.6873 0 0
M  V30 4 O -2.5109 0.0826 0 0
M  V30 5 C -1.0732 -0.4692 0 0
M  V30 6 C -5.1782 1.6227 0 0
M  V30 7 C -6.5119 2.3927 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 1 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 1 6
M  V30 6 1 7 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(2 1 6) FIELDNAME=property -
M  V30 FIELDDISP="   -5.1782    0.0827    DA    ALL  0       0" -
M  V30 FIELDDATA=val2 FIELDDATA=val1
M  V30 END SGROUP
M  V30 END CTAB
M  END
''')
    self.assertIsNotNone(mol)
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 1)
    sg = sgs[0]
    self.assertEqual(sg.GetProp('TYPE'), 'DAT')
    pv = sg.GetStringVectProp('DATAFIELDS')
    self.assertEqual(list(pv), ['val2', 'val1'])

  def testCopying(self):
    mol = Chem.MolFromMolBlock('''
  Mrv2014 07312005252D          
 
  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 3 0 0
M  V30 BEGIN ATOM
M  V30 1 * -12.75 11.5 0 0
M  V30 2 O -11.4163 12.27 0 0
M  V30 3 C -10.0826 11.5 0 0
M  V30 4 C -8.749 12.27 0 0
M  V30 5 O -10.0826 9.96 0 0
M  V30 6 N -7.4153 11.5 0 0
M  V30 7 C -6.0816 12.27 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 2 3 5
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SRU 0 ATOMS=(3 2 3 5) XBONDS=(2 1 3) BRKXYZ=(9 -9.9955 12.6173 0 -
M  V30 -9.0715 11.0169 0 0 0 0) BRKXYZ=(9 -11.5035 11.1527 0 -12.4275 12.7531 -
M  V30 0 0 0 0) CONNECT=HT LABEL=n
M  V30 2 DAT 0 ATOMS=(1 6) FIELDNAME=foo_data -
M  V30 FIELDDISP="   -7.4153   11.5000    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=bar
M  V30 3 DAT 0 ATOMS=(1 7) FIELDNAME=bar_data -
M  V30 FIELDDISP="   -6.0816   12.2700    DAU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=baz
M  V30 END SGROUP
M  V30 END CTAB
M  END''')
    self.assertEqual(len(Chem.GetMolSubstanceGroups(mol)), 3)
    mol2 = Chem.Mol(mol)
    Chem.ClearMolSubstanceGroups(mol2)
    self.assertEqual(len(Chem.GetMolSubstanceGroups(mol2)), 0)
    sgs = Chem.GetMolSubstanceGroups(mol)
    Chem.AddMolSubstanceGroup(mol2, sgs[0])
    Chem.AddMolSubstanceGroup(mol2, sgs[2])
    self.assertEqual(len(Chem.GetMolSubstanceGroups(mol2)), 2)
    molb = Chem.MolToV3KMolBlock(mol2)
    self.assertEqual(molb.find("foo_data"), -1)
    self.assertGreater(molb.find("M  V30 2 DAT 0 ATOMS=(1 7) FIELDNAME=bar_data"), 0)

    # we can also use this to copy SGroups:
    sgs2 = Chem.GetMolSubstanceGroups(mol2)
    newsg = Chem.AddMolSubstanceGroup(mol2, sgs[1])
    newsg.SetProp("FIELDNAME", "blah_data")
    molb = Chem.MolToV3KMolBlock(mol2)
    self.assertGreater(molb.find("M  V30 3 DAT 0 ATOMS=(1 6) FIELDNAME=blah_data"), 0)

  def testCStates(self):
    mol = Chem.MolFromMolBlock('''
  ACCLDraw08282007542D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 9 9 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.759 -6.2839 0 0 
M  V30 2 C 10.8013 -6.2833 0 0 
M  V30 3 C 9.7821 -5.6936 0 0 
M  V30 4 C 10.8013 -7.4648 0 0 
M  V30 5 C 8.759 -7.4701 0 0 
M  V30 6 C 9.7847 -8.0543 0 0 
M  V30 7 C 11.8245 -5.6926 0 0 
M  V30 8 O 12.8482 -6.2834 0 0 
M  V30 9 O 11.8245 -4.5104 0 0 
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 6 4 
M  V30 2 2 5 6 
M  V30 3 1 2 3 
M  V30 4 1 1 5 
M  V30 5 2 4 2 
M  V30 6 2 3 1 
M  V30 7 2 7 9 
M  V30 8 1 7 8 
M  V30 9 1 2 7 
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 1 ATOMS=(3 7 8 9) XBONDS=(1 9) CSTATE=(4 9 -1.02 -0.59 0) LABEL=-
M  V30 CO2H 
M  V30 END SGROUP
M  V30 END CTAB
M  END
''')
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 1)
    sg0 = sgs[0]
    pd = sg0.GetPropsAsDict()
    self.assertTrue('TYPE' in pd)
    self.assertEqual(pd['TYPE'], 'SUP')
    cstates = sg0.GetCStates()
    self.assertEqual(len(cstates), 1)
    cs0 = cstates[0]
    self.assertEqual(cs0.bondIdx, 8)
    self.assertAlmostEqual(cs0.vector.x, -1.02, 2)
    self.assertAlmostEqual(cs0.vector.y, -0.59, 2)

  def testBrackets(self):
    mol = Chem.MolFromMolBlock('''
  Mrv2014 08282011142D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 5 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -15.7083 2 0 0
M  V30 2 C -14.3747 2.77 0 0
M  V30 3 O -13.041 2 0 0
M  V30 4 C -11.7073 2.77 0 0
M  V30 5 C -10.3736 2 0 0
M  V30 6 C -9.0399 2.77 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 MON 0 ATOMS=(2 3 4) BRKXYZ=(9 -13.811 1.23 0 -13.811 3.54 0 0 0 0) -
M  V30 BRKXYZ=(9 -10.9373 3.54 0 -10.9373 1.23 0 0 0 0)
M  V30 END SGROUP
M  V30 END CTAB
M  END''')
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 1)
    sg0 = sgs[0]
    pd = sg0.GetPropsAsDict()
    self.assertTrue('TYPE' in pd)
    self.assertEqual(pd['TYPE'], 'MON')
    brackets = sg0.GetBrackets()
    self.assertEqual(len(brackets), 2)
    b = brackets[0]
    self.assertEqual(len(b), 3)
    self.assertAlmostEqual(b[0].x, -13.811, 3)
    self.assertAlmostEqual(b[0].y, 1.230, 3)
    self.assertAlmostEqual(b[1].x, -13.811, 3)
    self.assertAlmostEqual(b[1].y, 3.540, 3)
    self.assertAlmostEqual(b[2].x, 0, 3)
    self.assertAlmostEqual(b[2].y, 0, 3)

    b = brackets[1]
    self.assertEqual(len(b), 3)
    self.assertAlmostEqual(b[0].x, -10.937, 3)
    self.assertAlmostEqual(b[0].y, 3.54, 3)
    self.assertAlmostEqual(b[1].x, -10.937, 3)
    self.assertAlmostEqual(b[1].y, 1.23, 3)
    self.assertAlmostEqual(b[2].x, 0, 3)
    self.assertAlmostEqual(b[2].y, 0, 3)

  def testAttachPoints(self):
    mol = Chem.MolFromMolBlock('''
  Mrv2014 09012006262D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.0833 0.0833 0 0
M  V30 2 C -3.7497 0.8533 0 0
M  V30 3 O -2.416 0.0833 0 0
M  V30 4 O -3.7497 2.3933 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 2 4
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(3 2 3 4) SAP=(3 2 1 1) XBONDS=(1 1) LABEL=CO2H ESTATE=E
M  V30 END SGROUP
M  V30 END CTAB
M  END''')
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 1)
    sg0 = sgs[0]
    pd = sg0.GetPropsAsDict()
    self.assertTrue('TYPE' in pd)
    self.assertEqual(pd['TYPE'], 'SUP')
    aps = sg0.GetAttachPoints()
    self.assertEqual(len(aps), 1)
    self.assertEqual(aps[0].aIdx, 1)
    self.assertEqual(aps[0].lvIdx, 0)
    self.assertEqual(aps[0].id, '1')

    mol = Chem.MolFromMolBlock('''
  Mrv2014 09012006262D          

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 4 3 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.0833 0.0833 0 0
M  V30 2 C -3.7497 0.8533 0 0
M  V30 3 O -2.416 0.0833 0 0
M  V30 4 O -3.7497 2.3933 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 2 2 4
M  V30 END BOND
M  V30 END CTAB
M  END''')
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 0)

    sg = Chem.CreateMolSubstanceGroup(mol, "SUP")
    sg.AddAtomWithIdx(1)
    sg.AddAtomWithIdx(2)
    sg.AddAtomWithIdx(3)
    sg.AddBondWithIdx(0)
    sg.SetProp('LABEL', 'CO2H')
    sg.AddAttachPoint(1, 0, '1')
    molb = Chem.MolToV3KMolBlock(mol)
    self.assertGreater(
      molb.find('M  V30 1 SUP 0 ATOMS=(3 2 3 4) XBONDS=(1 1) LABEL=CO2H SAP=(3 2 1 1)'), 0)

  def testClearValues(self):
    mol = Chem.MolFromMolBlock('''example
 -ISIS-  10171405052D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 14 15 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 6.4292 -1.1916 0 0 CFG=3
M  V30 2 C 7.0125 -0.6042 0 0
M  V30 3 N 6.4292 -0.0250999 0 0
M  V30 4 C 5.8416 -0.6042 0 0
M  V30 5 C 5.8416 -1.7708 0 0
M  V30 6 N 6.4292 -2.3584 0 0 CFG=3
M  V30 7 C 7.0125 -1.7708 0 0
M  V30 8 O 5.7166 -3.5875 0 0
M  V30 9 C 5.7166 -4.4125 0 0 CFG=3
M  V30 10 C 4.8875 -4.4125 0 0
M  V30 11 C 6.5376 -4.4166 0 0
M  V30 12 C 5.7166 -5.2376 0 0
M  V30 13 C 6.4292 -3.175 0 0
M  V30 14 O 7.1375 -3.5875 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 1
M  V30 5 1 1 5
M  V30 6 1 5 6
M  V30 7 1 6 7
M  V30 8 1 7 1
M  V30 9 1 6 13
M  V30 10 1 8 9
M  V30 11 1 9 10
M  V30 12 1 9 11
M  V30 13 1 9 12
M  V30 14 2 13 14
M  V30 15 1 8 13
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 SUP 0 ATOMS=(7 8 9 10 11 12 13 14) XBONDS=(1 9) BRKXYZ=(9 6.24 -2.9 0 -
M  V30 6.24 -2.9 0 0 0 0) CSTATE=(4 9 0 0.82 0) LABEL=Boc SAP=(3 13 6 1)
M  V30 END SGROUP
M  V30 END CTAB
M  END''')
    sgs = Chem.GetMolSubstanceGroups(mol)
    self.assertEqual(len(sgs), 1)
    self.assertEqual(len(sgs[0].GetBrackets()), 1)
    sgs[0].ClearBrackets()
    self.assertEqual(len(sgs[0].GetBrackets()), 0)

    self.assertEqual(len(sgs[0].GetCStates()), 1)
    sgs[0].ClearCStates()
    self.assertEqual(len(sgs[0].GetCStates()), 0)

    self.assertEqual(len(sgs[0].GetAttachPoints()), 1)
    sgs[0].ClearAttachPoints()
    self.assertEqual(len(sgs[0].GetAttachPoints()), 0)

  def testCreateDataSGroup(self):
    mol = Chem.MolFromSmiles('CC(=O)O')
    sg = Chem.CreateMolDataSubstanceGroup(mol, "pKa", "4.5")
    sg.SetAtoms([1, 2, 3])
    sg = Chem.GetMolSubstanceGroups(mol)[0]
    self.assertEqual(sg.GetProp("TYPE"), "DAT")
    self.assertEqual(sg.GetProp("FIELDNAME"), "pKa")
    self.assertEqual(sg.GetStringVectProp("DATAFIELDS")[0], "4.5")


if __name__ == '__main__':
  print("Testing SubstanceGroups wrapper")
  unittest.main()
