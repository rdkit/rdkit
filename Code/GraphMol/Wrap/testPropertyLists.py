#
#  Copyright (C) 2019  Greg Landrum
#         All Rights Reserved
#
import unittest
from io import BytesIO

from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem.rdmolops import _TestSetProps


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

>  <bond.iprop.Number>  (1) 
3 2 1

$$$$"""

  def testForwardSupplier(self):
    sio = BytesIO(self.sdf)
    suppl = Chem.ForwardSDMolSupplier(sio)
    suppl.SetProcessPropertyLists(False)
    m = next(suppl)
    self.assertTrue(m.HasProp("atom.prop.AtomLabel"))
    self.assertFalse(m.GetAtomWithIdx(0).HasProp("AtomLabel"))
    self.assertTrue(m.HasProp("bond.iprop.Number"))
    self.assertFalse(m.GetBondWithIdx(0).HasProp("Number"))

    sio = BytesIO(self.sdf)
    suppl = Chem.ForwardSDMolSupplier(sio)
    self.assertTrue(suppl.GetProcessPropertyLists())
    m = next(suppl)
    self.assertTrue(m.HasProp("atom.prop.AtomLabel"))
    self.assertTrue(m.GetAtomWithIdx(0).HasProp("AtomLabel"))

    self.assertTrue(m.HasProp("bond.iprop.Number"))
    self.assertTrue(m.GetBondWithIdx(0).HasProp("Number"))
    self.assertTrue('Number' in m.GetBondWithIdx(0).GetPropsAsDict())
    self.assertEqual(m.GetBondWithIdx(0).GetIntProp("Number"), 3)

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
    Chem.CreateAtomIntPropertyList(m, "NumHeavyNeighbors")
    self.assertTrue(m.HasProp("atom.iprop.NumHeavyNeighbors"))

    self.assertTrue(m.GetAtomWithIdx(0).HasProp("PartialCharge"))
    m.ClearProp("atom.dprop.PartialCharge")
    self.assertFalse(m.HasProp("atom.dprop.PartialCharge"))
    Chem.CreateAtomDoublePropertyList(m, "PartialCharge")
    self.assertTrue(m.HasProp("atom.dprop.PartialCharge"))

    self.assertTrue(m.GetAtomWithIdx(0).HasProp("IsCarbon"))
    m.ClearProp("atom.bprop.IsCarbon")
    self.assertFalse(m.HasProp("atom.bprop.IsCarbon"))
    Chem.CreateAtomBoolPropertyList(m, "IsCarbon")
    self.assertTrue(m.HasProp("atom.bprop.IsCarbon"))

    self.assertTrue(m.GetAtomWithIdx(0).HasProp("PartiallyMissing"))
    m.ClearProp("atom.prop.PartiallyMissing")
    self.assertFalse(m.HasProp("atom.prop.PartiallyMissing"))
    Chem.CreateAtomStringPropertyList(m, "PartiallyMissing")
    self.assertTrue(m.HasProp("atom.prop.PartiallyMissing"))
    self.assertEqual(m.GetProp("atom.prop.PartiallyMissing"), "one n/a three")
    Chem.CreateAtomStringPropertyList(m, "PartiallyMissing", missingValueMarker="?")
    self.assertTrue(m.HasProp("atom.prop.PartiallyMissing"))
    self.assertEqual(m.GetProp("atom.prop.PartiallyMissing"), "[?] one ? three")

  def testGithubPR4160(self):
    # this shouldn't fail with a bad any cast anymore
    from rdkit import Chem
    m = Chem.MolFromSmiles("CC")
    for a in m.GetAtoms():
      a.SetIntProp("foo", 1)
    Chem.CreateAtomIntPropertyList(m, "foo")

  def testSetProps(self):
    from rdkit import Chem
    m = Chem.MolFromSmiles("CC")

    conf = Chem.Conformer(m.GetNumAtoms())
    for i in range(m.GetNumAtoms()):
      conf.SetAtomPosition(i, (0., 0., 0.))
    m.AddConformer(conf)

    _TestSetProps(m)
    default_expected = {
      'bool': True,
      'uint': 4294967295,
      'double': 3.14159,
      'svint': [0, 1, 2, -2],
      'svuint': [0, 1, 2, 4294967294],
      'svdouble': [0.0, 1.0, 2.0],
      'svstring': ['The', 'RDKit']
    }

    def check(ob, prefix):
      expected = {prefix + k: v for k, v in default_expected.items()}
      d = ob.GetPropsAsDict(False, False)
      for k, v in d.items():
        if 'sv' in k:
          d[k] = list(v)
      assert d == expected, repr((d, expected))
      for k, v in expected.items():
        v2 = ob.GetProp(k, True)
        if 'sv' in k:
          v2 = list(v2)
        assert v2 == v, repr(k, v2, v)
        assert type(ob.GetProp(k)) == str
      return len(d) > 0

    assert check(m, "mol_")
    for atom in m.GetAtoms():
      check(atom, f"atom_{atom.GetIdx()}")

    for bond in m.GetBonds():
      check(bond, f"bond_{bond.GetIdx()}")

    for idx, conf in enumerate(m.GetConformers()):
      check(conf, f"conf_{idx}")


if __name__ == '__main__':
  unittest.main()
