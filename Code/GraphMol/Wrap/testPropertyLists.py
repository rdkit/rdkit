#
#  Copyright (C) 2019  Greg Landrum
#         All Rights Reserved
#
from rdkit import RDConfig, rdBase
from rdkit import Chem
from io import BytesIO
import unittest


class TestCase(unittest.TestCase):
    def setUp(self):
        self.sdf = b"""
     RDKit  2D

  3  3  0  0  0  0  0  0  0  0999 V2000
    0.8660    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4330    0.7500    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4330   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  1  1  0
M  END
>  <atom.dprop.PartialCharge>  (1) 
0.008 -0.314 0.008

>  <atom.iprop.NumHeavyNeighbors>  (1) 
2 2 2

>  <atom.prop.AtomLabel>  (1) 
C1 N2 C3

>  <atom.bprop.IsCarbon>  (1) 
1 0 1

>  <atom.prop.PartiallyMissing>  (1) 
one n/a three

>  <atom.iprop.PartiallyMissingInt>  (1) 
[?] 2 2 ?
$$$$"""

    def testForwardSupplier(self):
        sio = BytesIO(self.sdf)
        suppl = Chem.ForwardSDMolSupplier(sio)
        suppl.SetProcessPropertyLists(False)
        m = next(suppl)
        self.assertTrue(m.HasProp("atom.prop.AtomLabel"))
        self.assertFalse(m.GetAtomWithIdx(0).HasProp("AtomLabel"))

        sio = BytesIO(self.sdf)
        suppl = Chem.ForwardSDMolSupplier(sio)
        self.assertTrue(suppl.GetProcessPropertyLists())
        m = next(suppl)
        self.assertTrue(m.HasProp("atom.prop.AtomLabel"))
        self.assertTrue(m.GetAtomWithIdx(0).HasProp("AtomLabel"))

    def testSupplier(self):
        suppl = Chem.SDMolSupplier()
        suppl.SetData(self.sdf)
        suppl.SetProcessPropertyLists(False)
        m = suppl[0]
        self.assertFalse(suppl.GetProcessPropertyLists())
        self.assertTrue(m.HasProp("atom.prop.AtomLabel"))
        self.assertFalse(m.GetAtomWithIdx(0).HasProp("AtomLabel"))

        suppl.SetProcessPropertyLists(True)
        m = suppl[0]
        self.assertTrue(m.HasProp("atom.prop.AtomLabel"))
        self.assertTrue(m.GetAtomWithIdx(0).HasProp("AtomLabel"))

    def testCreateLists(self):
        suppl = Chem.SDMolSupplier()
        suppl.SetData(self.sdf)
        m = suppl[0]
        self.assertTrue(m.GetAtomWithIdx(0).HasProp("NumHeavyNeighbors"))
        m.ClearProp("atom.iprop.NumHeavyNeighbors")
        self.assertFalse(m.HasProp("atom.iprop.NumHeavyNeighbors"))
        Chem.CreateAtomIntPropertyList(m,"NumHeavyNeighbors")
        self.assertTrue(m.HasProp("atom.iprop.NumHeavyNeighbors"))

        self.assertTrue(m.GetAtomWithIdx(0).HasProp("PartialCharge"))
        m.ClearProp("atom.dprop.PartialCharge")
        self.assertFalse(m.HasProp("atom.dprop.PartialCharge"))
        Chem.CreateAtomDoublePropertyList(m,"PartialCharge")
        self.assertTrue(m.HasProp("atom.dprop.PartialCharge"))

        self.assertTrue(m.GetAtomWithIdx(0).HasProp("IsCarbon"))
        m.ClearProp("atom.bprop.IsCarbon")
        self.assertFalse(m.HasProp("atom.bprop.IsCarbon"))
        Chem.CreateAtomBoolPropertyList(m,"IsCarbon")
        self.assertTrue(m.HasProp("atom.bprop.IsCarbon"))

        self.assertTrue(m.GetAtomWithIdx(0).HasProp("PartiallyMissing"))
        m.ClearProp("atom.prop.PartiallyMissing")
        self.assertFalse(m.HasProp("atom.prop.PartiallyMissing"))
        Chem.CreateAtomStringPropertyList(m,"PartiallyMissing")
        self.assertTrue(m.HasProp("atom.prop.PartiallyMissing"))
        self.assertEqual(m.GetProp("atom.prop.PartiallyMissing"),"one n/a three")
        Chem.CreateAtomStringPropertyList(m,"PartiallyMissing",missingValueMarker="?")
        self.assertTrue(m.HasProp("atom.prop.PartiallyMissing"))
        self.assertEqual(m.GetProp("atom.prop.PartiallyMissing"),"[?] one ? three")

    def testGithubPR4160(self):
        # this shouldn't fail with a bad any cast anymore
        from rdkit import Chem
        m = Chem.MolFromSmiles("CC")
        for a in m.GetAtoms():
             a.SetIntProp("foo", 1)
        Chem.CreateAtomIntPropertyList(m, "foo")
 


if __name__ == '__main__':
  unittest.main()

