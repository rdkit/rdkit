import unittest
import time

from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.TestRunner import isDebugBuild

def load_smiles(text):
  mols = []
  for line in text.strip().splitlines():
    smiles = line.split()[0]
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    mols.append(mol)
  return mols


_ignore = object()


class MCSTestCase(unittest.TestCase):

  def assert_search(self, smiles, numAtoms, numBonds, smarts=_ignore, verbose=False, **kwargs):
    result = rdFMCS.FindMCS(smiles, **kwargs)
    if verbose:
      print("RESULT: ",result.smartsString)
    self.assert_result(result, canceled=False, numAtoms=numAtoms, numBonds=numBonds, smarts=smarts)

  def assert_result(self, result, canceled=_ignore, numAtoms=_ignore, numBonds=_ignore,
                    smarts=_ignore):
    if canceled is not _ignore:
      self.assertEqual(result.canceled, canceled)
    if numAtoms is not _ignore:
      self.assertEqual(result.numAtoms, numAtoms)
    if numBonds is not _ignore:
      self.assertEqual(result.numBonds, numBonds)
    if smarts is not _ignore:
      self.assertEqual(result.smartsString, smarts)


simple_mols = load_smiles("""
c1ccccc1O phenol
CO methanol""")


# class TestMinAtoms(MCSTestCase):
#
#   def test_min_atoms_2(self):
#     self.assert_search(simple_mols, 2, 1, minNumAtoms=2)
#
#   def test_min_atoms_3(self):
#     self.assert_search(simple_mols, -1, -1, smarts=None, minNumAtoms=3)
#
#   def test_min_atoms_1(self):
#     try:
#       result = rdFMCS.FindMCS(simple_mols, minNumAtoms=1)
#     except ValueError:
#       pass
#     else:
#       raise AssertionError("should have raised an exception")


maximize_mols = load_smiles("""
C12CCC1CC2OCCCCCCC 2-rings-and-chain-with-O
C12CCC1CC2SCCCCCCC 2-rings-and-chain-with-S
""")


class TextMaximize(MCSTestCase):
  # C12CCC1CC2OCCCCCCC 2-rings-and-chain-with-O
  # C12CCC1CC2SCCCCCCC 2-rings-and-chain-with-S
  def test_maximize_default(self):
    # default maximizes the number of bonds
    self.assert_search(maximize_mols, 6, 7)

  def _test_maximize_atoms(self):
    # FIX: This fails to maximize the number of atoms
    self.assert_search(maximize_mols, 7, 6, maximizeBonds=False)

  def test_maximize_bonds(self):
    self.assert_search(maximize_mols, 6, 7, maximizeBonds=True)


atomtype_mols = load_smiles("""
c1ccccc1O phenol
CCCCCCOn1cccc1 different-answers-depending-on-type
""")


class TestAtomTypes(MCSTestCase):
  # The tests compare:
  #   c1ccccc1O
  #   CCCCCCOn1cccc1
  def test_atom_compare_default(self):
    self.assert_search(atomtype_mols, 4, 3, smarts='[#6](:[#6]:[#6]):[#6]',
                       bondCompare=rdFMCS.BondCompare.CompareOrderExact)

  def test_atom_compare_elements(self):
    self.assert_search(atomtype_mols, 4, 3, smarts='[#6](:[#6]:[#6]):[#6]',
                       atomCompare=rdFMCS.AtomCompare.CompareElements,
                       bondCompare=rdFMCS.BondCompare.CompareOrderExact)

  def test_atom_compare_any(self):
    # Note: bond aromaticities must still match!
    # 'cccccO' matches 'ccccnO'
    self.assert_search(atomtype_mols, 6, 5,
                       atomCompare=rdFMCS.AtomCompare.CompareAny,
                       bondCompare=rdFMCS.BondCompare.CompareOrderExact)

  def test_atom_compare_any_bond_compare_any(self):
    # Linear chain of 7 atoms
    self.assert_search(atomtype_mols, 7, 6,
                       atomCompare=rdFMCS.AtomCompare.CompareAny,
                        bondCompare=rdFMCS.BondCompare.CompareAny)

  def test_bond_compare_any(self):
    # Linear chain of 7 atoms
    self.assert_search(atomtype_mols, 7, 6,
                       bondCompare=rdFMCS.BondCompare.CompareAny)


isotope_mols = load_smiles("""
C1C[0N]CC[5C]1[1C][2C][2C][3C] C1223
C1CPCC[4C]1[2C][2C][1C][3C] C2213
""")


class TestIsotopes(MCSTestCase):
  # C1C[0N]CC[5C]1[1C][2C][2C][3C] C1223
  # C1CPCC[4C]1[2C][2C][1C][3C] C2213
  def test_without_isotope(self):
    # The entire system, except the N/P in the ring
    self.assert_search(isotope_mols, numAtoms=9, numBonds=8)

  def test_isotopes(self):
    # 5 atoms of class '0' in the ring
    self.assert_search(isotope_mols, 5, 4,
                       atomCompare=rdFMCS.AtomCompare.CompareIsotopes)

  def test_isotope_complete_ring_only(self):
    # the 122 in the chain
    self.assert_search(isotope_mols, 3, 2,
                       atomCompare=rdFMCS.AtomCompare.CompareIsotopes,
                       completeRingsOnly=True)


bondtype_mols = load_smiles("""
C1CCCCC1OC#CC#CC#CC#CC#CC first
c1ccccc1ONCCCCCCCCCC second
""")


class TestBondTypes(MCSTestCase):
  # C1CCCCC1OC#CC#CC#CC#CC#CC
  # c1ccccc1ONCCCCCCCCCC second
  def test_bond_compare_default(self):
    self.assert_search(bondtype_mols, 7, 7)

  def test_bond_compare_bondtypes(self):
    # Match the 'CCCCCC' part of the first ring, with the second's tail
    self.assert_search(bondtype_mols, 6, 5,
                       bondCompare=rdFMCS.BondCompare.CompareOrderExact)

  def test_bond_compare_any(self):
    # the CC#CC chain matches the CCCC tail
    self.assert_search(bondtype_mols, 10, 9,
                       bondCompare=rdFMCS.BondCompare.CompareAny)

  def test_atom_compare_elements_bond_compare_any(self):
    self.assert_search(bondtype_mols, 10, 9,
                       atomCompare=rdFMCS.AtomCompare.CompareElements,
                       bondCompare=rdFMCS.BondCompare.CompareAny)

  def test_atom_compare_any_bond_compare_any(self):
    # complete match!
    self.assert_search(bondtype_mols, 18, 18,
                       atomCompare=rdFMCS.AtomCompare.CompareAny,
                       bondCompare=rdFMCS.BondCompare.CompareAny)


valence_mols = load_smiles("""
CCCCCCCCN
CCCC(C)CCCC
""")


class TestValences(MCSTestCase):

  def test_valence_compare_default(self):
    # match 'CCCCCCCC'
    self.assert_search(valence_mols, 8, 7)

  def _test_valence_compare_valence(self):
    # FIX: matchValence isn't used?
    # match 'CCCC'
    self.assert_search(valence_mols, 4, 3, matchValences=True)

  def _test_valence_compare_valence(self):
    # FIX: matchValence isn't used?
    # match 'CCCCN' to '[CH-]CCCC' (but in reverse)
    self.assert_search(valence_mols, 5, 4, matchValences=True,
                       atomCompare=rdFMCS.AtomCompare.CompareAny)


ring_mols = load_smiles("""
C12CCCC(N2)CCCC1 6-and-7-bridge-rings-with-N
C1CCCCN1 6-ring
C1CCCCCN1 7-ring
C1CCCCCCCC1 9-ring
NC1CCCCCC1 N+7-ring
C1CC1CCCCCC 3-ring-with-tail
C12CCCC(O2)CCCC1 6-and-7-bridge-rings-with-O
""")


def SELECT(mols, *offsets):
  return [mols[offset - 1] for offset in offsets]


class TestRingMatchesRingOnly(MCSTestCase):
  # C12CCCC(N2)CCCC1 6-and-7-bridge-rings-with-N
  # C1CCCCN1 6-ring
  # C1CCCCCN1 7-ring
  # C1CCCCCCCC1 9-ring
  # NC1CCCCCC1 N+7-ring
  # C1CC1CCCCCC 3-ring-with-tail
  # C12CCCC(O2)CCCC1 6-and-7-bridge-rings-with-O
  def test_default(self):
    # Should match 'CCCCC'
    self.assert_search(ring_mols, 5, 4)

  def test_ring_only(self):
    # Should match "CCC"
    self.assert_search(ring_mols, 3, 2, ringMatchesRingOnly=True)

  def test_ring_only_select_1_2(self):
    # Should match "C1CCCCCN1"
    self.assert_search(SELECT(ring_mols, 1, 2), 6, 6, ringMatchesRingOnly=True)

  def test_ring_only_select_1_3(self):
    # Should match "C1CCCCCCN1"
    self.assert_search(SELECT(ring_mols, 1, 3), 7, 7, ringMatchesRingOnly=True)

  def test_ring_only_select_1_4(self):
    # Should match "C1CCCCCCCC1"
    self.assert_search(SELECT(ring_mols, 1, 4), 9, 9, ringMatchesRingOnly=True)

  def test_select_1_5(self):
    # Should match "NCCCCCC"
    self.assert_search(SELECT(ring_mols, 1, 5), 8, 7, ringMatchesRingOnly=False)

  def test_ring_only_select_1_5(self):
    # Should match "CCCCCC"
    self.assert_search(SELECT(ring_mols, 1, 5), 7, 6, ringMatchesRingOnly=True)

  def test_select_1_6(self):
    # Should match "CCCCCCCCC" by breaking one of the 3-carbon ring bonds
    self.assert_search(SELECT(ring_mols, 1, 6), 9, 8)

  def test_ring_only_select_1_6(self):
    # Should match "CCC" from the three atom ring
    self.assert_search(SELECT(ring_mols, 1, 6), 3, 2, ringMatchesRingOnly=True)

  def test_ring_only_select_1_7(self):
    # Should match the outer ring "C1CCCCCCCC1"
    self.assert_search(SELECT(ring_mols, 1, 7), 9, 9)

  def test_ring_only_select_1_7_any_atoms(self):
    # Should match everything
    self.assert_search(SELECT(ring_mols, 1, 7), 10, 11, ringMatchesRingOnly=True,
                       atomCompare=rdFMCS.AtomCompare.CompareAny)


class TestCompleteRingsOnly(MCSTestCase):
  # C12CCCC(N2)CCCC1 6-and-7-bridge-rings-with-N
  # C1CCCCN1 6-ring
  # C1CCCCCN1 7-ring
  # C1CCCCCCCC1 9-ring
  # NC1CCCCCC1 N+7-ring
  # C1CC1CCCCCC 3-ring-with-tail
  # C12CCCC(O2)CCCC1 6-and-7-bridge-rings-with-O
  def test_ring_only(self):
    # No match: "CCC" is not in a ring
    self.assert_search(ring_mols, 0, 0, completeRingsOnly=True)

  def test_ring_only_select_1_2(self):
    # Should match "C1CCCCCN1"
    self.assert_search(SELECT(ring_mols, 1, 2), 6, 6, completeRingsOnly=True)

  def test_ring_only_select_1_3(self):
    # Should match "C1CCCCCCN1"
    self.assert_search(SELECT(ring_mols, 1, 3), 7, 7, completeRingsOnly=True)

  def _test_ring_only_select_1_4(self):
    # FIX: This does not match and should
    # Should match "C1CCCCCCCC1"
    self.assert_search(SELECT(ring_mols, 1, 4), 9, 9, completeRingsOnly=True)

  def test_ring_only_select_1_5(self):
    # No match: "CCCCCC" is not in a ring
    self.assert_search(SELECT(ring_mols, 1, 5), 0, 0, completeRingsOnly=True)

  def _test_ring_only_select_1_7(self):
    # FIX: This does not match and should
    # Should match the outer ring "C1CCCCCCCC1"
    self.assert_search(SELECT(ring_mols, 1, 7), 9, 9, completeRingsOnly=True)

  def test_ring_only_select_1_7_any_atoms(self):
    # Should match everything
    self.assert_search(SELECT(ring_mols, 1, 7), 10, 11, completeRingsOnly=True,
                       atomCompare=rdFMCS.AtomCompare.CompareAny)


  def test_ring_to_nonring_bond(self):
    # Should allow the cO in phenol to match the CO in the other structure
    self.assert_search(atomtype_mols, 2, 1, completeRingsOnly=True)


lengthy_mols = [Chem.MolFromSmiles("Nc1ccccc1" * 20), Chem.MolFromSmiles("Nc1ccccccccc1" * 20)]


class TestTimeout(MCSTestCase):
  # This should take over two minutes to process. Give it 1 seconds.
  @unittest.skipIf(isDebugBuild(), "Timeout test will fail on debug builds.")
  def test_timeout(self):
    t1 = time.time()
    result = rdFMCS.FindMCS(lengthy_mols, timeout=1)
    self.assert_result(result, canceled=True)
    self.assertTrue(result.numAtoms > 1)
    self.assertTrue(result.numBonds >= result.numAtoms - 1, (result.numAtoms, result.numBonds))
    t2 = time.time()
    self.assertTrue(t2 - t1 < 2, t2 - t1)

  # Check for non-negative values
  def test_timeout_negative(self):
    try:
      rdFMCS.FindMCS(lengthy_mols, timeout=-1)
    except OverflowError:
      pass
    else:
      raise AssertionError("bad range check for timeout")


if __name__ == "__main__":
  unittest.main()
