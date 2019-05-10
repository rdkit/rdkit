
from rdkit import RDConfig
import unittest
from rdkit import Chem
from rdkit.Chem import rdMMPA


def natoms(tpl):
  return tuple(x.GetNumAtoms() if x is not None else 0 for x in tpl)


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    m = Chem.MolFromSmiles('c1ccccc1OC')
    frags = rdMMPA.FragmentMol(m)
    self.assertEqual(len(frags), 3)
    for frag in frags:
      self.assertEqual(len(frag), 2)
    frags = sorted(frags, key=natoms)
    self.assertEqual(frags[0][0], None)
    self.assertEqual(frags[1][0], None)
    self.assertNotEqual(frags[2][0], None)
    self.assertNotEqual(frags[0][1], None)
    self.assertNotEqual(frags[1][1], None)
    self.assertNotEqual(frags[2][1], None)

    self.assertEqual(frags[0][1].GetNumAtoms(), m.GetNumAtoms() + 2)
    self.assertEqual(frags[1][1].GetNumAtoms(), m.GetNumAtoms() + 2)

    fs = Chem.GetMolFrags(frags[0][1], asMols=True)
    self.assertEqual(len(fs), 2)
    self.assertEqual(Chem.MolToSmiles(fs[0], True), 'c1ccc([*:1])cc1')
    self.assertEqual(Chem.MolToSmiles(fs[1], True), 'CO[*:1]')

    fs = Chem.GetMolFrags(frags[1][1], asMols=True)
    self.assertEqual(len(fs), 2)
    self.assertEqual(Chem.MolToSmiles(fs[0], True), 'c1ccc(O[*:1])cc1')
    self.assertEqual(Chem.MolToSmiles(fs[1], True), 'C[*:1]')

    fs = Chem.GetMolFrags(frags[2][0], asMols=True)
    self.assertEqual(len(fs), 1)
    self.assertEqual(Chem.MolToSmiles(fs[0], True), 'O([*:1])[*:2]')
    fs = Chem.GetMolFrags(frags[2][1], asMols=True)
    self.assertEqual(len(fs), 2)
    self.assertEqual(Chem.MolToSmiles(fs[0], True), 'c1ccc([*:2])cc1')
    self.assertEqual(Chem.MolToSmiles(fs[1], True), 'C[*:1]')

  def test2(self):
    m = Chem.MolFromSmiles('c1ccccc1OC')
    frags = rdMMPA.FragmentMol(m, resultsAsMols=False)
    self.assertEqual(len(frags), 3)
    for frag in frags:
      self.assertEqual(len(frag), 2)
    frags = sorted(frags)
    self.assertEqual(frags[0][0], '')
    self.assertEqual(frags[1][0], '')
    self.assertNotEqual(frags[2][0], '')
    self.assertNotEqual(frags[0][1], '')
    self.assertNotEqual(frags[1][1], '')
    self.assertNotEqual(frags[2][1], '')

    self.assertEqual(frags[0][1], 'CO[*:1].c1ccc([*:1])cc1')
    self.assertEqual(frags[1][1], 'C[*:1].c1ccc(O[*:1])cc1')
    self.assertEqual(frags[2][0], 'O([*:1])[*:2]')
    self.assertEqual(frags[2][1], 'C[*:1].c1ccc([*:2])cc1')

  def test3(self):
    m = Chem.MolFromSmiles('c1ccccc1OC')
    frags = rdMMPA.FragmentMol(m, resultsAsMols=False, pattern='cO')
    self.assertEqual(len(frags), 1)
    for frag in frags:
      self.assertEqual(len(frag), 2)
    frags = sorted(frags)
    self.assertEqual(frags[0][0], '')
    self.assertNotEqual(frags[0][1], '')

    self.assertEqual(frags[0][1], 'CO[*:1].c1ccc([*:1])cc1')

  def test4(self):
    m = Chem.MolFromSmiles('Cc1ccccc1NC(=O)C(C)[NH+]1CCCC1')  # ZINC00000051
    frags = rdMMPA.FragmentMol(m, resultsAsMols=False)
    #for frag in sorted(frags):
    #    print(frag)
    cores = set(x[0] for x in frags)
    self.assertTrue('C([*:1])([*:2])[*:3]' in cores)
    # FIX: this needs to be investigated, it's not currently passing
    #self.assertTrue('O=C(N[*:3])C([*:1])[*:2]' in cores)
    self.assertEqual(len(frags), 18)
    for frag in frags:
      self.assertEqual(len(frag), 2)

  def test5(self):
    m = Chem.MolFromSmiles(
      "CC[C@H](C)[C@@H](C(=O)N[C@H]1CSSC[C@H]2C(=O)NCC(=O)N3CCC[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H](C(=O)N[C@H](C(=O)N4CCC[C@H]4C(=O)N[C@H](C(=O)N2)C)CC(=O)N)NC1=O)C(=O)N)CO)Cc5ccc(cc5)O)CCCC[NH3+])N")  # ALPHA-CONOTOXIN SI
    frags = rdMMPA.FragmentMol(m, resultsAsMols=False)
    self.assertFalse(len(frags))
    frags = rdMMPA.FragmentMol(m, maxCuts=2, maxCutBonds=21, resultsAsMols=False)
    self.assertEqual(len(frags), 231)

  def test6(self):
    m = Chem.MolFromSmiles(
      "CC[C@H](C)[C@@H](C(=O)N[C@H]1CSSC[C@H]2C(=O)NCC(=O)N3CCC[C@H]3C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](CSSC[C@@H](C(=O)N[C@H](C(=O)N4CCC[C@H]4C(=O)N[C@H](C(=O)N2)C)CC(=O)N)NC1=O)C(=O)N)CO)Cc5ccc(cc5)O)CCCC[NH3+])N")  # ALPHA-CONOTOXIN SI
    frags = rdMMPA.FragmentMol(m, resultsAsMols=False)
    self.assertFalse(len(frags))
    frags1 = rdMMPA.FragmentMol(m, minCuts=1, maxCuts=1, maxCutBonds=21,
                                resultsAsMols=False)
    frags2 = rdMMPA.FragmentMol(m, minCuts=2, maxCuts=2, maxCutBonds=21,
                                resultsAsMols=False)

    frags = rdMMPA.FragmentMol(m, maxCuts=2, maxCutBonds=21, resultsAsMols=False)

    self.assertEqual(set(frags1+frags2), set(frags))
    self.assertEqual(set(frags1).intersection(set(frags2)), set())


  def test7(self):
    m = Chem.MolFromSmiles("Oc1ccccc1N")

    frags1 = rdMMPA.FragmentMol(m, minCuts=1, maxCuts=1, maxCutBonds=21,
                                resultsAsMols=False)

    frags2 = rdMMPA.FragmentMol(m, minCuts=2, maxCuts=2, maxCutBonds=21,
                                resultsAsMols=False)

    frags = rdMMPA.FragmentMol(m, maxCuts=2, maxCutBonds=21, resultsAsMols=False)

    self.assertEqual(set(frags1+frags2), set(frags))


  def test8(self):
    m = Chem.MolFromSmiles('Cc1ccccc1NC(=O)C(C)[NH+]1CCCC1')  # ZINC00000051
    sm = Chem.MolFromSmarts("[#6+0;!$(*=,#[!#6])]!@!=!#[*]")
    matching_atoms = m.GetSubstructMatches(sm)
    bonds = []
    for a,b in matching_atoms:
      bond = m.GetBondBetweenAtoms(a,b)
      bonds.append(bond.GetIdx())

    frags = rdMMPA.FragmentMol(m, resultsAsMols=False)
    frags2 =  rdMMPA.FragmentMol(m, bonds, resultsAsMols=False)
    frags3 =  rdMMPA.FragmentMol(m, tuple(bonds), resultsAsMols=False)
    self.assertEqual(frags, frags2)
    self.assertEqual(frags2, frags3)

  def test9(self):
    m = Chem.MolFromSmiles("Oc1ccccc1N")

    try:
      
      frags1 = rdMMPA.FragmentMol(m, minCuts=1, maxCuts=0, maxCutBonds=21,
                                  resultsAsMols=False)
      self.assertTrue(False) # should not get here
    except ValueError as e:
      self.assertEqual(str(e), "supplied maxCuts is less than minCuts")

    try:
      
      frags1 = rdMMPA.FragmentMol(m, minCuts=0, maxCuts=0, maxCutBonds=21,
                                  resultsAsMols=False)
      self.assertTrue(False) # should not get here
    except ValueError as e:
      self.assertEqual(str(e), "minCuts must be greater than 0")
    


if __name__ == "__main__":
  unittest.main()
