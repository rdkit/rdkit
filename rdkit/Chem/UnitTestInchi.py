#
#  Copyright (c) 2011-2025, Novartis Institutes for BioMedical Research Inc.
#   and other RDKit contributors
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import gzip
import io
import os
import pickle
import re
import unittest

from rdkit import RDConfig, RDLogger
from rdkit.Chem import (INCHI_AVAILABLE, ForwardSDMolSupplier, MolFromMolBlock, MolFromSmiles,
                        MolToMolBlock, MolToSmiles, SanitizeMol, rdDepictor)

if INCHI_AVAILABLE:
  from rdkit.Chem import (InchiReadWriteError, InchiToInchiKey, MolBlockToInchi, MolFromInchi,
                          MolFromInchiAndAuxInfo, MolToInchi, MolToInchiAndAuxInfo, MolToInchiKey,
                          GetInchiVersion)

COLOR_RED = '\033[31m'
COLOR_GREEN = '\033[32m'
COLOR_RESET = '\033[0m'


def inchiDiffPrefix(inchi1, inchi2):
  inchi1 = inchi1.split('/')
  inchi2 = inchi2.split('/')
  for i in range(len(inchi1) + 1):
    if i == len(inchi1):
      break
    if i == len(inchi2) or inchi1[i] != inchi2[i]:
      break
  if len(inchi1) >= i:
    return inchi1[i][0]
  else:
    return inchi2[i][0]


def inchiDiff(inchi1, inchi2):
  inchi1 = inchi1.split('/')
  inchi2 = inchi2.split('/')
  for i in range(len(inchi1) + 1):
    if i == len(inchi1):
      break
    if i == len(inchi2) or inchi1[i] != inchi2[i]:
      break
  return ('/'.join(inchi1[:i]) + COLOR_RED + '/' + '/'.join(inchi1[i:]) + '\n' + COLOR_RESET +
          '/'.join(inchi2[:i]) + COLOR_RED + '/' + '/'.join(inchi2[i:]) + COLOR_RESET)


@unittest.skipUnless(INCHI_AVAILABLE, 'Inchi support not available')
class RegressionTest(unittest.TestCase):

  def testPrechloricAcid(self):
    examples = (
      ('OCl(=O)(=O)=O', 'InChI=1S/ClHO4/c2-1(3,4)5/h(H,2,3,4,5)'),
      ('CC1=CC2=NCC(CN2C=C1)C(=O)c3ccc4cc(C)ccc4c3.OCl(=O)(=O)=O',
       'InChI=1S/C21H20N2O.ClHO4/c1-14-3-4-17-11-18(6-5-16(17)9-14)21(24)19-12-22-20-10-15(2)7-8-23(20)13-19;2-1(3,4)5/h3-11,19H,12-13H2,1-2H3;(H,2,3,4,5)'
       ),
      ('CNc1ccc2nc3ccccc3[n+](C)c2c1.[O-]Cl(=O)(=O)=O',
       'InChI=1S/C14H13N3.ClHO4/c1-15-10-7-8-12-14(9-10)17(2)13-6-4-3-5-11(13)16-12;2-1(3,4)5/h3-9H,1-2H3;(H,2,3,4,5)'
       ),
    )
    for smiles, expected in examples:
      m = MolFromSmiles(smiles)
      inchi = MolToInchi(m)
      self.assertEqual(inchi, expected)


@unittest.skipUnless(INCHI_AVAILABLE, 'Inchi support not available')
class TestCase(unittest.TestCase):

  def setUp(self):
    self.dataset = dict()
    self.dataset_inchi = dict()
    inf = gzip.open(os.path.join(RDConfig.RDCodeDir, 'Chem/test_data', 'pubchem-hard-set.sdf.gz'),
                    'r')
    self.dataset['problematic'] = ForwardSDMolSupplier(inf, sanitize=False, removeHs=False)
    with open(os.path.join(RDConfig.RDCodeDir, 'Chem/test_data', 'pubchem-hard-set.inchi'),
              'r') as intF:
      buf = intF.read().replace('\r\n', '\n').encode('latin1')
      intF.close()
    with io.BytesIO(buf) as inF:
      pkl = inF.read()
    self.dataset_inchi['problematic'] = pickle.loads(pkl, encoding='latin1')
    # disable logging
    RDLogger.DisableLog('rdApp.warning')

  def tearDown(self):
    RDLogger.EnableLog('rdApp.warning')
    RDLogger.EnableLog('rdApp.error')

  def test0InchiWritePubChem(self):
    for fp, f in self.dataset.items():
      inchi_db = self.dataset_inchi[fp]
      same, diff, reasonable = 0, 0, 0
      for m in f:
        if m is None:  # pragma: nocover
          continue
        ref_inchi = inchi_db[m.GetProp('PUBCHEM_COMPOUND_CID')]
        x, y = MolToInchi(m), ref_inchi
        if x != y:
          # print("---------------")
          # print(m.GetProp('PUBCHEM_COMPOUND_CID'))
          # print(MolToSmiles(m))
          # print(y)
          # print(x)
          if re.search(r'.[1-9]?ClO4', x) is not None:
            reasonable += 1
            continue
          SanitizeMol(m)
          if filter(lambda i: i >= 8, [len(r) for r in m.GetRingInfo().AtomRings()]):
            reasonable += 1
            continue
          # THERE ARE NO EXAMPLES FOR THE FOLLOWING (no coverage)
          # if it is because RDKit does not think the bond is stereo
          z = MolToInchi(MolFromMolBlock(MolToMolBlock(m)))
          if y != z and inchiDiffPrefix(y, z) == 'b':
            reasonable += 1
            continue
          # some warning
          try:
            MolToInchi(m, treatWarningAsError=True)
          except InchiReadWriteError as inst:
            _, error = inst.args
            if 'Metal' in error:
              reasonable += 1
              continue

          diff += 1
          print('InChI mismatch for PubChem Compound ' + m.GetProp('PUBCHEM_COMPOUND_CID'))
          print(MolToSmiles(m, True))
          print(inchiDiff(x, y))
          print()

        else:
          same += 1

      fmt = "\n{0}InChI write Summary: {1} identical, {2} suffix variance, {3} reasonable{4}"
      print(fmt.format(COLOR_GREEN, same, diff, reasonable, COLOR_RESET))
      self.assertEqual(same, 1162)
      self.assertEqual(diff, 0)
      self.assertEqual(reasonable, 19)

  def test1InchiReadPubChem(self):
    for f in self.dataset.values():
      same, diff, reasonable = 0, 0, 0
      for m in f:
        if m is None:  # pragma: nocover
          continue
        x = MolToInchi(m)
        y = None
        RDLogger.DisableLog('rdApp.error')
        mol = MolFromInchi(x)
        RDLogger.EnableLog('rdApp.error')
        if mol is not None:
          y = MolToInchi(MolFromSmiles(MolToSmiles(mol, isomericSmiles=True)))
        if y is None:
          # metal involved?
          try:
            MolToInchi(m, treatWarningAsError=True)
          except InchiReadWriteError as inst:
            _, error = inst.args
            if 'Metal' in error or \
                    'Charges were rearranged' in error:
              reasonable += 1
              continue
          # THERE ARE NO EXAMPLES FOR THE FOLLOWING (no coverage)
          # RDKit does not like the SMILES? use MolBlock instead
          inchiMol = MolFromInchi(x)
          if inchiMol:
            rdDepictor.Compute2DCoords(inchiMol)
            z = MolToInchi(MolFromMolBlock(MolToMolBlock(inchiMol)))
            if x == z:
              reasonable += 1
              continue
          # InChI messed up the radical?
          unsanitizedInchiMol = MolFromInchi(x, sanitize=False)
          if sum([
              a.GetNumRadicalElectrons() * a.GetAtomicNum() for a in m.GetAtoms()
              if a.GetNumRadicalElectrons() != 0
          ]) != sum([
              a.GetNumRadicalElectrons() * a.GetAtomicNum() for a in unsanitizedInchiMol.GetAtoms()
              if a.GetNumRadicalElectrons() != 0
          ]):
            reasonable += 1
            continue

          diff += 1
          cid = m.GetProp('PUBCHEM_COMPOUND_CID')
          print(COLOR_GREEN + 'Empty mol for PubChem Compound ' + cid + '\n' + COLOR_RESET)
          continue
        if x != y:
          # if there was warning in the first place, then this is
          # tolerable
          try:
            MolToInchi(m, treatWarningAsError=True)
            MolFromInchi(x, treatWarningAsError=True)
          except InchiReadWriteError as inst:
            reasonable += 1
            continue
          # or if there are big rings
          SanitizeMol(m)
          if filter(lambda i: i >= 8, [len(r) for r in m.GetRingInfo().AtomRings()]):
            reasonable += 1
            continue
          # THERE ARE NO EXAMPLES FOR THE FOLLOWING (no coverage)
          # or if RDKit loses bond stereo
          s = MolToSmiles(m, True)
          if MolToSmiles(MolFromSmiles(s), True) != s:
            reasonable += 1
            continue
          # or if it is RDKit SMILES writer unhappy about the mol
          inchiMol = MolFromInchi(x)
          rdDepictor.Compute2DCoords(inchiMol)
          z = MolToInchi(MolFromMolBlock(MolToMolBlock(inchiMol)))
          if x == z:
            reasonable += 1
            continue

          diff += 1
          print(COLOR_GREEN + 'Molecule mismatch for PubChem Compound ' + cid + COLOR_RESET)
          print(inchiDiff(x, y))
          print()
        else:
          same += 1
      fmt = "\n{0}InChI read Summary: {1} identical, {2} variance, {3} reasonable variance{4}"
      print(fmt.format(COLOR_GREEN, same, diff, reasonable, COLOR_RESET))
      self.assertEqual(same, 688)
      self.assertEqual(diff, 0)
      self.assertEqual(reasonable, 493)

  def test2InchiOptions(self):
    m = MolFromSmiles("CC=C(N)C")
    inchi1 = MolToInchi(m).split('/', 1)[1]
    inchi2 = MolToInchi(m, "/SUU").split('/', 1)[1]
    self.assertEqual(inchi1 + '/b4-3?', inchi2)

  def test3InchiKey(self):
    inchi = 'InChI=1S/C9H12/c1-2-6-9-7-4-3-5-8-9/h3-5,7-8H,2,6H2,1H3'
    self.assertEqual(InchiToInchiKey(inchi), 'ODLMAHJVESYWTB-UHFFFAOYSA-N')

  def test4MolToInchiKey(self):
    m = MolFromSmiles("CC=C(N)C")
    inchi = MolToInchi(m)
    k1 = InchiToInchiKey(inchi)
    k2 = MolToInchiKey(m)
    self.assertEqual(k1, k2)

  def test5MolBlockToInchi(self):
    mb = """
  Mrv1824 02111920092D          

  6  6  0  0  0  0            999 V2000
   -5.5134    3.5259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    3.1134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5134    1.8759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    3.1134    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  2  0  0  0  0
M  END"""
    inchi = MolBlockToInchi(mb)
    self.assertEqual(inchi, "InChI=1S/C5H8O/c1-2-4-6-5-3-1/h1-2H,3-5H2")
    # make sure that options work
    mb2 = """
  Mrv1824 02121905282D          

 10 11  0  0  0  0            999 V2000
   -4.6875   -1.1393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4020   -1.5518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.4020   -2.3768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.6875   -2.7893    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9730   -2.3768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.9730   -1.5518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.2586   -2.7893    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5441   -1.5518    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5441   -2.3768    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -1.9608   -0.9684    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0  0  0  0
  2  3  1  0  0  0  0
  3  4  2  0  0  0  0
  4  5  1  0  0  0  0
  5  6  2  0  0  0  0
  1  6  1  0  0  0  0
  7  9  1  0  0  0  0
  6  8  1  0  0  0  0
  7  5  1  0  0  0  0
  8 10  1  0  0  0  0
  8  9  2  0  0  0  0
M  END"""
    inchi2 = MolBlockToInchi(mb2, options="/FixedH")
    self.assertEqual(inchi2, "InChI=1/C8H8N2/c1-6-7-4-2-3-5-8(7)10-9-6/h2-5H,1H3,(H,9,10)/f/h10H")

  def test6GetInchiVersion(self):
    version = GetInchiVersion()
    self.assertIsInstance(version, str)
    self.assertGreaterEqual(version, "1.07.2")


@unittest.skipUnless(INCHI_AVAILABLE, 'Inchi support not available')
class TestMolFromInchiAndAuxInfo(unittest.TestCase):

  def test0RoundTripAtomOrder(self):
    """Verify that round-tripping through InChI+AuxInfo preserves atom ordering."""
    from rdkit.Chem import MolFromSmiles, MolToSmiles
    smiles_list = ['c1ccccc1O', 'CC(=O)O', 'C(=O)(N)C', 'c1cc(O)ccc1N']
    for smi in smiles_list:
      mol = MolFromSmiles(smi)
      inchi, aux = MolToInchiAndAuxInfo(mol)
      mol2 = MolFromInchiAndAuxInfo(inchi, aux)
      self.assertIsNotNone(mol2)
      # The atom ordering should be restored, so atom symbols in order should match
      orig_atoms = [a.GetAtomicNum() for a in mol.GetAtoms()]
      new_atoms = [a.GetAtomicNum() for a in mol2.GetAtoms()]
      self.assertEqual(orig_atoms, new_atoms,
                       f"Atom order mismatch for {smi}: {orig_atoms} vs {new_atoms}")

  def test1StereoPreservation(self):
    """Verify stereochemistry is preserved through round-trip."""
    from rdkit.Chem import MolFromSmiles
    smi = '[C@@H](O)(F)Cl'
    mol = MolFromSmiles(smi)
    inchi, aux = MolToInchiAndAuxInfo(mol)
    mol2 = MolFromInchiAndAuxInfo(inchi, aux)
    self.assertIsNotNone(mol2)
    inchi2 = MolToInchi(mol2)
    self.assertEqual(inchi, inchi2)

  def test2NoneAuxInfo(self):
    """MolFromInchiAndAuxInfo with None auxinfo should still return a molecule."""
    inchi = 'InChI=1S/CH4/h1H4'
    mol = MolFromInchiAndAuxInfo(inchi, None)
    self.assertIsNotNone(mol)
    self.assertEqual(mol.GetNumAtoms(), 1)

  def test3EmptyAuxInfo(self):
    """MolFromInchiAndAuxInfo with empty auxinfo should still return a molecule."""
    inchi = 'InChI=1S/CH4/h1H4'
    mol = MolFromInchiAndAuxInfo(inchi, '')
    self.assertIsNotNone(mol)

  def test4InvalidInchi(self):
    """MolFromInchiAndAuxInfo with invalid InChI should return None."""
    mol = MolFromInchiAndAuxInfo('not_an_inchi', 'AuxInfo=1/0/N:1/')
    self.assertIsNone(mol)

  def test5MultiFragmentAuxInfo(self):
    """Test with a molecule that produces multi-fragment AuxInfo (semicolons)."""
    from rdkit.Chem import MolFromSmiles
    smi = 'CC.OO'
    mol = MolFromSmiles(smi)
    inchi, aux = MolToInchiAndAuxInfo(mol)
    mol2 = MolFromInchiAndAuxInfo(inchi, aux)
    self.assertIsNotNone(mol2)
    self.assertEqual(mol2.GetNumAtoms(), mol.GetNumAtoms())

  def test6CoordinateRestoration(self):
    """Round-trip a mol block with 2D coords and verify conformer is restored and
    trailing semicolons in /rC: are handled correctly."""
    from rdkit.Chem import MolFromMolBlock
    mb = """
  Mrv1824 02111920092D

  6  6  0  0  0  0            999 V2000
   -5.5134    3.5259    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    3.1134    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.2279    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -5.5134    1.8759    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    2.2884    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7989    3.1134    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  4  1  0  0  0  0
  4  5  1  0  0  0  0
  5  6  1  0  0  0  0
  1  6  1  0  0  0  0
  2  3  2  0  0  0  0
M  END"""
    mol = MolFromMolBlock(mb)
    self.assertIsNotNone(mol)
    inchi, aux = MolToInchiAndAuxInfo(mol)
    mol2 = MolFromInchiAndAuxInfo(inchi, aux)
    self.assertIsNotNone(mol2)
    self.assertEqual(mol2.GetNumConformers(), 1)
    # Verify round-trip produces populated /rC: layer
    inchi2, aux2 = MolToInchiAndAuxInfo(mol2)
    self.assertIn('/rC:', aux2)
    self.assertNotIn('/rC:;', aux2)
    # Real AuxInfo has trailing semicolons in /rC: — verify our parser handles them
    rc_match = re.search(r'/rC:([^/]+)', aux)
    self.assertIsNotNone(rc_match, "AuxInfo should contain /rC: layer")
    self.assertTrue(rc_match.group(1).endswith(';'),
                    "Real AuxInfo /rC: should end with trailing semicolon")

  def test7EmptyCoordinates(self):
    """AuxInfo with empty /rC: or no /rC: should produce no conformer, no crash."""
    from rdkit.Chem import MolFromSmiles
    mol = MolFromSmiles('C')
    inchi, aux = MolToInchiAndAuxInfo(mol)
    # Manually strip or replace /rC: to simulate empty coords
    aux_empty = re.sub(r'/rC:[^/]*', '/rC:;', aux)
    mol2 = MolFromInchiAndAuxInfo(inchi, aux_empty)
    self.assertIsNotNone(mol2)
    self.assertEqual(mol2.GetNumConformers(), 0)

  def test7bAllZeroCoordinates(self):
    """AuxInfo with all-zero /rC: entries (,,;) should produce no conformer."""
    from rdkit.Chem import MolFromSmiles
    mol = MolFromSmiles('CCO')
    inchi, aux = MolToInchiAndAuxInfo(mol)
    # Replace /rC: content with all-zero entries (,,; format from InChI spec)
    num_atoms = mol.GetNumAtoms()
    zero_coords = ';'.join([',,'] * num_atoms) + ';'
    aux_zeros = re.sub(r'/rC:[^/]*', '/rC:' + zero_coords, aux)
    mol2 = MolFromInchiAndAuxInfo(inchi, aux_zeros)
    self.assertIsNotNone(mol2)
    self.assertEqual(mol2.GetNumConformers(), 0)

  def test8CoordinateAtomOrderMatch(self):
    """Verify coordinates are on the correct atoms after reordering."""
    from rdkit.Chem import MolFromMolBlock
    mb = """
  Mrv1824 02111920092D

  3  2  0  0  0  0            999 V2000
    1.0000    2.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    3.0000    4.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    5.0000    6.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END"""
    mol = MolFromMolBlock(mb)
    self.assertIsNotNone(mol)
    orig_data = list(zip([a.GetSymbol() for a in mol.GetAtoms()], mol.GetConformer().GetPositions()))
    inchi, aux = MolToInchiAndAuxInfo(mol)
    mol2 = MolFromInchiAndAuxInfo(inchi, aux)
    self.assertIsNotNone(mol2)
    self.assertEqual(mol2.GetNumConformers(), 1)
    # Each atom should have the same symbol and (approximately) same coordinates
    for i, a in enumerate(mol2.GetAtoms()):
      orig_sym, orig_pos = orig_data[i]
      self.assertEqual(a.GetSymbol(), orig_sym)
      pos = mol2.GetConformer().GetAtomPosition(i)
      for j in range(3):
        self.assertAlmostEqual([pos.x, pos.y, pos.z][j], orig_pos[j], places=4,
                               msg=f"Coord mismatch for atom {i} ({a.GetSymbol()})")


if __name__ == '__main__':  # pragma: nocover
  # only run the test if InChI is available
  if INCHI_AVAILABLE:
    unittest.main()
