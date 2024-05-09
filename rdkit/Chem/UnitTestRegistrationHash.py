"""
Tests for rdkit.Chem.RegistrationHash

Focus is on molecules that have differing SMILES, but
are actually the same.
"""

import os
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import RegistrationHash
from rdkit.Chem.RegistrationHash import HashLayer


def hash_sdf(molblock: str, escape=None, data_field_names=None, enable_tautomer_hash_v2=False):
  """
    Gets the layers of the SDF, and generates the hash based on only the layers passed in

    :param molblock: The molblock to hash
    :param escape: optional field which can contain arbitrary information
    :param data_field_names: optional sequence of names of SGroup DAT fields which
       will be included in the hash.

    :return: A dict with the hash & the layers
    """
  mol = Chem.MolFromMolBlock(molblock)
  return RegistrationHash.GetMolLayers(mol, escape=escape, data_field_names=data_field_names,
                                       enable_tautomer_hash_v2=enable_tautomer_hash_v2)


class CanonicalizerTest(unittest.TestCase):
  maxDiff = 2000  # hash diffs can be long!

  def test_example_structure(self):

    structure = """
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 21 23 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.00838 1.55562 0 0
M  V30 2 C -2.00918 0.555619 0 0
M  V30 3 C -2.87558 0.056219 0 0
M  V30 4 C -2.87638 -0.943781 0 0
M  V30 5 C -3.74278 -1.44298 0 0
M  V30 6 N -2.01078 -1.44438 0 0
M  V30 7 C -1.14438 -0.945181 0 0
M  V30 8 C -0.278781 -1.44578 0 0
M  V30 9 C 0.587619 -0.946381 0 0
M  V30 10 N 1.45322 -1.44718 0 0
M  V30 11 C 2.31962 -0.947781 0 0
M  V30 12 O 3.18522 -1.44858 0 0
M  V30 13 C 2.32042 0.0522191 0 0
M  V30 14 C 1.45482 0.552819 0 0
M  V30 15 C 1.45562 1.55282 0 0
M  V30 16 F 2.45562 1.55202 0 0
M  V30 17 F 1.95642 2.41842 0 0
M  V30 18 F 0.590019 2.05362 0 0
M  V30 19 C 0.588419 0.0536191 0 0
M  V30 20 C -0.277181 0.554219 0 0
M  V30 21 C -1.14358 0.054819 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 2 1 CFG=3
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5 CFG=3
M  V30 5 1 4 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 2 9 10
M  V30 10 1 10 11
M  V30 11 1 11 12
M  V30 12 2 11 13
M  V30 13 1 13 14
M  V30 14 1 14 15
M  V30 15 1 15 16
M  V30 16 1 15 17
M  V30 17 1 15 18
M  V30 18 2 14 19
M  V30 19 1 19 20
M  V30 20 2 20 21
M  V30 21 1 21 2
M  V30 22 1 21 7
M  V30 23 1 19 9
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEREL1 ATOMS=(2 2 4)
M  V30 END COLLECTION
M  V30 END CTAB
M  END

$$$$
"""

    layers = hash_sdf(structure)
    expected_layers = {
      HashLayer.CANONICAL_SMILES:
      "C[C@H]1C[C@@H](C)c2cc3c(C(F)(F)F)cc(O)nc3cc2N1 |o1:1,3|",
      HashLayer.ESCAPE:
      "",
      HashLayer.FORMULA:
      "C15H15F3N2O",
      HashLayer.NO_STEREO_SMILES:
      "CC1CC(C)c2cc3c(C(F)(F)F)cc(O)nc3cc2N1",
      HashLayer.NO_STEREO_TAUTOMER_HASH:
      "CC1CC(C)[C]2[CH][C]3[C]([CH][C]2[N]1)[N][C]([O])[CH][C]3C(F)(F)F_2_0",
      HashLayer.SGROUP_DATA:
      "[]",
      HashLayer.TAUTOMER_HASH:
      "C[C@H]1C[C@@H](C)[C]2[CH][C]3[C]([CH][C]2[N]1)[N][C]([O])[CH][C]3C(F)(F)F_2_0 |o1:1,3|",
    }
    self.assertEqual(layers, expected_layers)

    # should have the same NO_STEREO_SMILES the structure with stereo removed
    mol = Chem.MolFromMolBlock(structure)
    Chem.rdmolops.RemoveStereochemistry(mol)
    stereo_insensitive_layers = hash_sdf(Chem.MolToMolBlock(mol))
    self.assertEqual(stereo_insensitive_layers[HashLayer.NO_STEREO_SMILES],
                     expected_layers[HashLayer.NO_STEREO_SMILES])

  def test_XBHEAD_XBCORR(self):
    # rdkit/Code/GraphMol/FileParsers/sgroup_test_data/repeat_groups_query1.mol
    structure = """
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
"""

    expected_layers = {
      HashLayer.CANONICAL_SMILES: "C1CCCCC1",
      HashLayer.ESCAPE: "",
      HashLayer.FORMULA: "C6H12",
      HashLayer.NO_STEREO_SMILES: "C1CCCCC1",
      HashLayer.NO_STEREO_TAUTOMER_HASH: "C1CCCCC1_0_0",
      HashLayer.SGROUP_DATA:
      '[{"type": "SRU", "atoms": [2], "bonds": [[1, 2], [2, 3]], "index": 1, "connect": "HT", "label": "1-3", "XBHEAD": [[0, 5], [1, 0]], "XBCORR": [[[0, 5], [0, 5]], [[1, 0], [1, 0]]]}]',
      HashLayer.TAUTOMER_HASH: "C1CCCCC1_0_0",
    }

    layers = hash_sdf(structure)
    self.assertEqual(layers, expected_layers)

    # test for the same hash if the order in the XBHEAD block is changed
    structure.replace("XBHEAD=(2 6 1)", "XBHEAD=(2 1 6)")

    layers = hash_sdf(structure)
    self.assertEqual(layers, expected_layers)

  def test_data_sgroups(self):
    # see SHARED-7879
    structure = """
   Mrv1808 05312108392D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 7 6 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C -5.6458 2.7017 0 0
M  V30 2 C -4.3121 1.9317 0 0
M  V30 3 C -2.9785 2.7017 0 0
M  V30 4 C -1.6448 1.9317 0 0
M  V30 5 C -0.3111 2.7017 0 0
M  V30 6 C 1.0226 1.9317 0 0
M  V30 7 C 2.3563 2.7017 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 1 5 6
M  V30 6 1 6 7
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 7) FIELDNAME="lambda max" FIELDINFO=nm -
M  V30 FIELDDISP="    1.4300   -0.7700    DRU   ALL  0       0" -
M  V30 MRV_FIELDDISP=0 FIELDDATA=250
M  V30 END SGROUP
M  V30 END CTAB
M  END
"""

    expected_layers = {
      HashLayer.CANONICAL_SMILES: 'CCCCCCC',
      HashLayer.ESCAPE: '',
      HashLayer.FORMULA: 'C7H16',
      HashLayer.NO_STEREO_SMILES: 'CCCCCCC',
      HashLayer.NO_STEREO_TAUTOMER_HASH: 'CCCCCCC_0_0',
      HashLayer.SGROUP_DATA: '[]',
      HashLayer.TAUTOMER_HASH: 'CCCCCCC_0_0'
    }

    expected_layers_v2 = {
      HashLayer.CANONICAL_SMILES: 'CCCCCCC',
      HashLayer.ESCAPE: '',
      HashLayer.FORMULA: 'C7H16',
      HashLayer.NO_STEREO_SMILES: 'CCCCCCC',
      HashLayer.NO_STEREO_TAUTOMER_HASH: '[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]_0_0',
      HashLayer.SGROUP_DATA: '[]',
      HashLayer.TAUTOMER_HASH: '[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH3]_0_0'
    }

    layers = hash_sdf(structure)
    self.assertEqual(layers, expected_layers)

    layers_v2 = hash_sdf(structure, enable_tautomer_hash_v2=True)
    self.assertEqual(layers_v2, expected_layers_v2)

  def test_stereo_imine_hash(self):
    stereo_mol = Chem.MolFromSmiles(r'[H]\N=C\C', sanitize=False)
    stereo = RegistrationHash.GetMolLayers(stereo_mol)

    no_stereo_mol = Chem.MolFromSmiles('N=CC', sanitize=False)
    no_stereo = RegistrationHash.GetMolLayers(no_stereo_mol)

    hydrogen_as_is_mol = Chem.MolFromSmiles('N=CC[H]', sanitize=False)
    hydrogen_as_is = RegistrationHash.GetMolLayers(hydrogen_as_is_mol)

    self.assertNotEqual(stereo[HashLayer.CANONICAL_SMILES], no_stereo[HashLayer.CANONICAL_SMILES])
    self.assertEqual(no_stereo[HashLayer.CANONICAL_SMILES],
                     hydrogen_as_is[HashLayer.CANONICAL_SMILES])

    self.assertEqual(stereo[HashLayer.NO_STEREO_SMILES], no_stereo[HashLayer.NO_STEREO_SMILES])
    self.assertEqual(no_stereo[HashLayer.NO_STEREO_SMILES],
                     hydrogen_as_is[HashLayer.NO_STEREO_SMILES])

    self.assertEqual(stereo[HashLayer.NO_STEREO_TAUTOMER_HASH],
                     no_stereo[HashLayer.NO_STEREO_TAUTOMER_HASH])
    self.assertEqual(no_stereo[HashLayer.NO_STEREO_TAUTOMER_HASH],
                     hydrogen_as_is[HashLayer.NO_STEREO_TAUTOMER_HASH])

    self.assertEqual(stereo[HashLayer.TAUTOMER_HASH], no_stereo[HashLayer.TAUTOMER_HASH])
    self.assertEqual(no_stereo[HashLayer.TAUTOMER_HASH], hydrogen_as_is[HashLayer.TAUTOMER_HASH])

  def test_empty(self):
    """Can the hasher operate on an empty molecule without traceback?"""
    mol = Chem.Mol()
    RegistrationHash.GetMolLayers(mol)

  def test_enhanced_stereo_canonicalizer(self):
    """Can we correctly canonicalize a molecule with enhanced stereo?"""
    groups = (
      # basics
      ('C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)F |a:3,&1:1,7,&2:5,r|',
       'C[C@H](O)[C@H](C)[C@@H](C)[C@H](C)F |a:3,&1:1,7,&2:5,r|'),

      # renumbering the groups
      ('C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)F |a:3,&1:1,7,&2:5,r|',
       'C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)F |a:3,&2:1,7,&1:5,r|'),

      # AND/OR canonicalization
      ('C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,o1:7,&2:1,&3:5,r|',
       'C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,&2:7,o1:1,&3:5,r|'),

      # Dan's example in SHARED-7811 (corrected)
      (
        r'C[C@@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@H](C)C1 |o1:1,7,o2:5,11|',  # s1
        r'C[C@H](O)[C@H](C)CC1C[C@H](C)C[C@@H](C)C1 |o1:1,8;o2:3,11|',  # s2
        r'C[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,5,o2:7,11|'  # s3 fixed!
        r'C[C@@H](O)[C@H](C)CC1C[C@H](C)C[C@H](C)C1 |o1:1,11;o2:3,8|',  # s4
      ),
    )
    for mols_smis in groups:
      csmis = set()
      for smi in mols_smis:
        mol = Chem.MolFromSmiles(smi)
        csmi = Chem.MolToCXSmiles(mol)
        self.assertIsNotNone(csmi)
        csmis.add(csmi)
      self.assertEqual(len(csmis), 1)

  def test_enhanced_stereo_canonicalizer_non_matching(self):
    # Uncorrected example in SHARED-7811. These are NOT equivalent,
    # since the groups are defined on non-interconvertible pairs of atomss
    s1 = r'C[C@@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@H](C)C1 |o1:1,7,o2:5,11|'
    s3 = r'C[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,11,o2:5,7|'

    csmis = set()
    for smi in (s1, s3):
      mol = Chem.MolFromSmiles(smi)
      csmi = Chem.MolToCXSmiles(mol)
      csmis.add(csmi)
    self.assertEqual(len(csmis), 2)

  def test_enhanced_stereo_regex(self):
    cxsmileses = ((
      'C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)F |a:3,&1:1,7,&2:5,r|',
      'C[C@@H](O)[C@H](C)[C@@H](C)[C@@H](C)O |&3:3,o1:7,&2:1,&3:5,r|',
      'C[C@@H](O)[C@H](C)CC1C[C@H](C)C[C@H](C)C1 |o1:1,11;o2:3,8|',
    ))
    for cxsmiles in cxsmileses:
      self.assertEqual(cxsmiles.count(':'),
                       len(RegistrationHash.ENHANCED_STEREO_GROUP_REGEX.findall(cxsmiles)))

  def test_shared_7984(self):
    smi = '[2H]C([2H])=C([2H])C(=O)O'
    no_stereo_mol = Chem.MolFromSmiles(smi)
    self.assertEqual(no_stereo_mol.GetNumAtoms(), 8)

    stereo_mol = Chem.Mol(no_stereo_mol)

    # We consider a double bond with these features equal to a STEREONONE one.
    stereo_mol.GetBondWithIdx(2).SetStereo(Chem.BondStereo.STEREOANY)
    stereo_mol.GetBondWithIdx(2).SetBondDir(Chem.BondDir.EITHERDOUBLE)

    no_stereo_hash = RegistrationHash.GetMolLayers(no_stereo_mol)
    stereo_hash = RegistrationHash.GetMolLayers(stereo_mol)

    self.assertEqual(no_stereo_hash, stereo_hash)

  def test_non_matching_pseudoatom_label(self):

    pol_sdf = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Pol 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"""

    mod_sdf = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Mod 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"""

    carbon_sdf = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"""

    pol_hash = RegistrationHash.GetMolHash(hash_sdf(pol_sdf))
    mod_hash = RegistrationHash.GetMolHash(hash_sdf(mod_sdf))
    carbon_hash = RegistrationHash.GetMolHash(hash_sdf(carbon_sdf))

    hashes = {pol_hash, mod_hash, carbon_hash}
    self.assertEqual(len(hashes), 3)

  def test_matching_pseudoatom_label(self):

    pol_sdf = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 Pol 0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"""

    pol_sdf_2 = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    2.5981   -0.0000    0.0000 Pol 0  0  0  0  0  0  0  0  0  0  0  0
    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END"""

    pol_hash = RegistrationHash.GetMolHash(hash_sdf(pol_sdf))
    pol_hash_2 = RegistrationHash.GetMolHash(hash_sdf(pol_sdf_2))

    hashes = {pol_hash, pol_hash_2}
    self.assertEqual(len(hashes), 1)

  def test_hash_schemes(self):

    def get_hash(m, scheme):
      hash = RegistrationHash.GetMolHash(RegistrationHash.GetMolLayers(m), scheme)
      self.assertEqual(len(hash), 40)  # SHA-1 hashes are 40 hex digits long
      return hash

    def same_hash(m1, m2, scheme):
      return get_hash(m1, scheme) == get_hash(m2, scheme)

    mol1 = Chem.MolFromSmiles("C[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,11,o2:5,7|")
    self.assertEqual(mol1.GetNumAtoms(), 14)

    all_layers_hash = get_hash(mol1, RegistrationHash.HashScheme.ALL_LAYERS)
    stereo_insensitive_layers_hash = get_hash(mol1,
                                              RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS)
    tautomer_insensitive_layers_hash = get_hash(
      mol1, RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS)

    self.assertNotEqual(all_layers_hash, tautomer_insensitive_layers_hash)
    self.assertNotEqual(all_layers_hash, stereo_insensitive_layers_hash)
    self.assertNotEqual(stereo_insensitive_layers_hash, tautomer_insensitive_layers_hash)

    # Compare stereoisomers
    mol2 = Chem.MolFromSmiles("CC(O)C(C)CC1CC(C)CC(C)C1")
    self.assertFalse(same_hash(mol1, mol2, RegistrationHash.HashScheme.ALL_LAYERS))
    self.assertTrue(same_hash(mol1, mol2, RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS))
    self.assertFalse(same_hash(mol1, mol2, RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS))

    # Compare tautomers
    mol1 = Chem.MolFromSmiles("N1C=NC2=CN=CN=C12")
    mol2 = Chem.MolFromSmiles("N1C=NC2=NC=NC=C12")
    self.assertFalse(same_hash(mol1, mol2, RegistrationHash.HashScheme.ALL_LAYERS))
    self.assertFalse(same_hash(mol1, mol2, RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS))
    self.assertTrue(same_hash(mol1, mol2, RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS))

    # Compare with and without sgroup data we canonicalize on (DAT/SRU/COP)
    smiles = "CNCC(=O)OC"
    mol1 = Chem.MolFromSmiles(smiles)
    mol2 = Chem.MolFromSmiles(f"{smiles} |Sg:n:5,3,4::ht|")
    # ...where all hash schemes account for sgroup data
    self.assertFalse(same_hash(mol1, mol2, RegistrationHash.HashScheme.ALL_LAYERS))
    self.assertFalse(same_hash(mol1, mol2, RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS))
    self.assertFalse(same_hash(mol1, mol2, RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS))

  def test_no_stereo_hash(self):
    """
        Check that none of the layers involved in the STEREO_INSENSITIVE_LAYERS HashScheme
        includes any SMILES stereo markers
        """
    smi = r'C\C=C/[C@@H](C)Cl'
    mol = Chem.MolFromSmiles(smi)
    self.assertEqual(mol.GetNumAtoms(), 6)

    layers = RegistrationHash.GetMolLayers(mol)

    for layer in RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS.value:
      self.assertNotIn('@', layers[layer])
      self.assertNotIn('/', layers[layer])
      self.assertNotIn('\\', layers[layer])

  def test_atom_map_numbers(self):
    """
        Check that atom map numbers are stripped in the molhash (SHARED-8007)
        """
    mol_w_atom_map = Chem.MolFromSmiles('Oc1cc([*:1])ccn1')
    mol = Chem.MolFromSmiles('Oc1cc([*])ccn1')
    atom_map_hash_layers = RegistrationHash.GetMolLayers(mol_w_atom_map)
    hash_layers = RegistrationHash.GetMolLayers(mol)
    assert atom_map_hash_layers[HashLayer.CANONICAL_SMILES] == hash_layers[
      HashLayer.CANONICAL_SMILES]

  def test_compare_stereoisomers(self):
    stereo_mol = Chem.MolFromSmiles("C[C@H]1CC(C[C@@H](C)[C@@H](C)O)C[C@@H](C)C1 |o1:1,11,o2:5,7|")
    no_stereo_mol = Chem.MolFromSmiles("CC(O)C(C)CC1CC(C)CC(C)C1")

    stereo_layers = RegistrationHash.GetMolLayers(stereo_mol)
    stereo_insensitive_layers = RegistrationHash.GetMolLayers(no_stereo_mol)
    self.assertEqual(
      RegistrationHash.GetMolHash(stereo_layers,
                                  RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS),
      RegistrationHash.GetMolHash(stereo_insensitive_layers,
                                  RegistrationHash.HashScheme.STEREO_INSENSITIVE_LAYERS))
    self.assertEqual(stereo_layers[HashLayer.NO_STEREO_SMILES], "CC1CC(C)CC(CC(C)C(C)O)C1")
    self.assertEqual(stereo_insensitive_layers[HashLayer.NO_STEREO_SMILES],
                     "CC1CC(C)CC(CC(C)C(C)O)C1")

  def _compare_layers_with_expected_diffs(self, new_layers, old_layers, expected_diffs):
    for k in old_layers:
      if k not in expected_diffs:
        self.assertEqual(new_layers[k], old_layers[k])
      else:
        self.assertEqual(new_layers[k], expected_diffs[k])

  def test_sgroup_data_hashes(self):
    structure1 = '''
  Mrv2108 11252113552D

  0  0  0     0  0            999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 30 34 2 0 0
M  V30 BEGIN ATOM
M  V30 1 C 8.78 -9.2456 0 0
M  V30 2 C 7.2925 -8.8469 0 0
M  V30 3 C 6.2033 -9.9358 0 0
M  V30 4 C 6.6021 -11.4235 0 0
M  V30 5 C 8.0896 -11.8222 0 0
M  V30 6 C 9.1786 -10.733 0 0
M  V30 7 C 10.6662 -11.1317 0 0
M  V30 8 C 11.7552 -10.0426 0 0
M  V30 9 C 11.3565 -8.5552 0 0
M  V30 10 C 9.8688 -8.1565 0 0
M  V30 11 C 13.2428 -10.441 0 0
M  V30 12 C 14.3317 -9.3519 0 0 CFG=2
M  V30 13 C 15.8524 -9.5927 0 0
M  V30 14 C 16.5517 -8.2207 0 0
M  V30 15 C 15.463 -7.1316 0 0
M  V30 16 C 14.0905 -7.8307 0 0
M  V30 17 C 16.5515 -10.9649 0 0
M  V30 18 C 8.4882 -13.3096 0 0
M  V30 19 C 5.5131 -12.5124 0 0
M  V30 20 C 5.9117 -14.0001 0 0
M  V30 21 C 4.8227 -15.0888 0 0
M  V30 22 C 3.3352 -14.6902 0 0
M  V30 23 C 2.9365 -13.2027 0 0
M  V30 24 C 4.0255 -12.1139 0 0
M  V30 25 C 3.6269 -10.6262 0 0
M  V30 26 Cl 4.7157 -9.5372 0 0
M  V30 27 C 5.2472 -16.5691 0 0
M  V30 28 C 6.7414 -16.9417 0 0
M  V30 29 C 7.8112 -15.8339 0 0
M  V30 30 C 7.3867 -14.3535 0 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 2 1 2
M  V30 2 1 2 3
M  V30 3 2 3 4
M  V30 4 1 4 5
M  V30 5 2 5 6
M  V30 6 1 6 1
M  V30 7 1 6 7
M  V30 8 2 7 8
M  V30 9 1 8 9
M  V30 10 2 9 10
M  V30 11 1 10 1
M  V30 12 1 12 11 CFG=1
M  V30 13 1 12 13
M  V30 14 1 13 14
M  V30 15 1 14 15
M  V30 16 1 15 16
M  V30 17 1 16 12
M  V30 18 1 13 17
M  V30 19 1 5 18
M  V30 20 1 4 19
M  V30 21 1 19 20
M  V30 22 2 20 21
M  V30 23 1 21 22
M  V30 24 2 22 23
M  V30 25 1 23 24
M  V30 26 2 24 19
M  V30 27 1 24 25
M  V30 28 1 3 26
M  V30 29 1 21 27
M  V30 30 1 20 30
M  V30 31 2 29 30
M  V30 32 1 28 29
M  V30 33 2 27 28
M  V30 34 1 8 11
M  V30 END BOND
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 12) FIELDNAME=stereolabel -
M  V30 FIELDDISP="   14.8391  -10.2748    DA    ALL  1       5" FIELDDATA=R
M  V30 2 DAT 0 ATOMS=(2 4 19) FIELDNAME=STEREOLABEL -
M  V30 FIELDDISP="    6.4455  -12.6102    DA    ALL  1       5" FIELDDATA=R
M  V30 END SGROUP
M  V30 END CTAB
M  END
'''
    structure2 = structure1.replace('FIELDDATA=R', 'FIELDDATA=S')
    structure3 = structure1.replace('FIELDDATA=R', 'FIELDDATA=S', 1)

    # default behavior is that these are all the same:
    layers1 = hash_sdf(structure1)
    layers2 = hash_sdf(structure2)
    expected_layers1 = {
      HashLayer.CANONICAL_SMILES:
      "Cc1ccc2ccccc2c1-c1c(Cl)cc2ccc(C[C@H]3CCCC3C)cc2c1C",
      HashLayer.ESCAPE:
      "",
      HashLayer.FORMULA:
      "C29H29Cl",
      HashLayer.NO_STEREO_SMILES:
      "Cc1ccc2ccccc2c1-c1c(Cl)cc2ccc(CC3CCCC3C)cc2c1C",
      HashLayer.NO_STEREO_TAUTOMER_HASH:
      "C[C]1[CH][CH][C]2[CH][CH][CH][CH][C]2[C]1[C]1[C](Cl)[CH][C]2[CH][CH][C](CC3CCCC3C)[CH][C]2[C]1C_0_0",
      HashLayer.SGROUP_DATA:
      "[]",
      HashLayer.TAUTOMER_HASH:
      "C[C]1[CH][CH][C]2[CH][CH][CH][CH][C]2[C]1[C]1[C](Cl)[CH][C]2[CH][CH][C](C[C@H]3CCCC3C)[CH][C]2[C]1C_0_0",
    }
    self.assertEqual(layers1, expected_layers1)
    self.assertEqual(layers2, layers1)

    # using both stereo labels they are all different from each other:
    layers1a = hash_sdf(structure1, data_field_names=['STEREOLABEL', 'stereolabel'])
    expected_diffs = {
      HashLayer.SGROUP_DATA:
      '[{"fieldName": "STEREOLABEL", "atom": [10, 11], "bonds": [], "value": "R"}, {"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "R"}]',
    }
    self._compare_layers_with_expected_diffs(layers1a, layers1, expected_diffs)

    layers2a = hash_sdf(structure2, data_field_names=['STEREOLABEL', 'stereolabel'])
    expected_diffs = {
      HashLayer.SGROUP_DATA:
      '[{"fieldName": "STEREOLABEL", "atom": [10, 11], "bonds": [], "value": "S"}, {"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "S"}]',
    }
    self._compare_layers_with_expected_diffs(layers2a, layers1, expected_diffs)

    layers3a = hash_sdf(structure3, data_field_names=['STEREOLABEL', 'stereolabel'])
    expected_diffs = {
      HashLayer.SGROUP_DATA:
      '[{"fieldName": "STEREOLABEL", "atom": [10, 11], "bonds": [], "value": "R"}, {"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "S"}]',
    }
    self._compare_layers_with_expected_diffs(layers3a, layers1, expected_diffs)

    # ensure that we can pick the data field to use:
    layers1b = hash_sdf(structure1, data_field_names=['stereolabel'])
    expected_diffs = {
      HashLayer.SGROUP_DATA:
      '[{"fieldName": "stereolabel", "atom": [20], "bonds": [], "value": "R"}]',
    }
    self._compare_layers_with_expected_diffs(layers1b, layers1, expected_diffs)

  def test_escape_layer(self):
    mol = Chem.MolFromSmiles("c1ccccc1")

    layers = RegistrationHash.GetMolLayers(mol)
    layers_with_escape = RegistrationHash.GetMolLayers(mol, escape="make it unique")
    assert layers_with_escape[HashLayer.ESCAPE] == "make it unique"
    for k in layers_with_escape:
      if k != HashLayer.ESCAPE:
        assert layers_with_escape[k] == layers[k]

  def test_shared_8516(self):
    no_explicit_h_mol = Chem.MolFromSmiles('C[C@H](N)C(=O)O', sanitize=False)
    with_explicit_h_mol = Chem.MolFromSmiles('[H][C@@](C)(N)C(=O)O', sanitize=False)
    no_explicit_h_layers = RegistrationHash.GetMolLayers(no_explicit_h_mol)
    with_explicit_h_layers = RegistrationHash.GetMolLayers(with_explicit_h_mol)
    self.assertEqual(no_explicit_h_layers, with_explicit_h_layers)

  def test_mol_layers_distinguish_carbon_isotopes(self):
    benzene = Chem.MolFromSmiles("c1ccccc1")
    benzene_with_c13 = Chem.MolFromSmiles("[13c]1ccccc1")

    layers1 = RegistrationHash.GetMolLayers(benzene)
    layers2 = RegistrationHash.GetMolLayers(benzene_with_c13)

    self.assertNotEqual(layers1, layers2)

  def test_mol_layers_distinguish_non_tautomeric_hydrogen_isotopes(self):
    benzene_with_deuterium = Chem.MolFromSmiles("[2H]c1ccccc1")
    layers = RegistrationHash.GetMolLayers(benzene_with_deuterium)

    expected_layers = {
      HashLayer.CANONICAL_SMILES: "[2H]c1ccccc1",
      HashLayer.ESCAPE: "",
      HashLayer.FORMULA: "C6H6",
      HashLayer.NO_STEREO_SMILES: "[2H]c1ccccc1",
      HashLayer.NO_STEREO_TAUTOMER_HASH: "[2H][C]1[CH][CH][CH][CH][CH]1_0_0",
      HashLayer.SGROUP_DATA: "[]",
      HashLayer.TAUTOMER_HASH: "[2H][C]1[CH][CH][CH][CH][CH]1_0_0"
    }

    assert layers == expected_layers

  def test_mol_layers_distinguish_tautomeric_hydrogen_isotopes(self):
    mol1 = Chem.MolFromSmiles("[2H]N1C=NC2=CN=CN=C12")
    mol2 = Chem.MolFromSmiles("[2H]N1C=NC2=NC=NC=C12")

    layers1 = RegistrationHash.GetMolLayers(mol1)
    layers2 = RegistrationHash.GetMolLayers(mol2)

    expected_layers1 = {
      HashLayer.CANONICAL_SMILES: "[2H]n1cnc2cncnc21",
      HashLayer.ESCAPE: "",
      HashLayer.FORMULA: "C5H4N4",
      HashLayer.NO_STEREO_SMILES: "[2H]n1cnc2cncnc21",
      HashLayer.NO_STEREO_TAUTOMER_HASH: "[2H]N1[CH][N][C]2[CH][N][CH][N][C]21_0_0",
      HashLayer.SGROUP_DATA: "[]",
      HashLayer.TAUTOMER_HASH: "[2H]N1[CH][N][C]2[CH][N][CH][N][C]21_0_0"
    }

    expected_layers2 = {
      HashLayer.CANONICAL_SMILES: "[2H]n1cnc2ncncc21",
      HashLayer.ESCAPE: "",
      HashLayer.FORMULA: "C5H4N4",
      HashLayer.NO_STEREO_SMILES: "[2H]n1cnc2ncncc21",
      HashLayer.NO_STEREO_TAUTOMER_HASH: "[2H]N1[CH][N][C]2[N][CH][N][CH][C]21_0_0",
      HashLayer.SGROUP_DATA: "[]",
      HashLayer.TAUTOMER_HASH: "[2H]N1[CH][N][C]2[N][CH][N][CH][C]21_0_0"
    }

    self.assertEqual(layers1, expected_layers1)
    self.assertEqual(layers2, expected_layers2)

  def testIodine(self):
    """Does stereo group canonicalization mess up isotopes?"""
    mol = Chem.MolFromSmiles('CC[C@H](C)[999C@H](C)O  |o1:2,4|')
    layers = RegistrationHash.GetMolLayers(mol)
    self.assertEqual(layers[HashLayer.CANONICAL_SMILES], 'CC[C@H](C)[999C@H](C)O |o1:2,4|')

  def testBadBondDir(self):
    """"""
    molBlock = """
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 12 12 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 5.192672 -7.044349 0.000000 0
M  V30 2 C 5.272538 -5.675621 0.000000 0
M  V30 3 C 6.415595 -5.133214 0.000000 0
M  V30 4 C 7.487292 -5.854749 0.000000 0
M  V30 5 C 6.572156 -3.950482 0.000000 0
M  V30 6 C 5.935220 -2.698719 0.000000 0
M  V30 7 C 4.531603 -2.325867 0.000000 0
M  V30 8 C 4.159340 -0.962544 0.000000 0
M  V30 9 C 2.686034 -0.657857 0.000000 0
M  V30 10 C 5.166201 0.053779 0.000000 0
M  V30 11 C 6.621450 -0.316643 0.000000 0
M  V30 12 C 6.973858 -1.673599 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 3 1 3 4
M  V30 4 2 3 5
M  V30 5 1 6 5 CFG=3
M  V30 6 1 6 7
M  V30 7 1 7 8
M  V30 8 2 8 9
M  V30 9 1 8 10
M  V30 10 1 11 10
M  V30 11 1 6 12
M  V30 12 1 11 12
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STERAC1 ATOMS=(1 6)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
$$$$
"""
    mol = Chem.MolFromMolBlock(molBlock)
    # Simulate reapplying the wedging from the original molBlock
    mol.GetBondBetweenAtoms(4, 5).SetBondDir(Chem.BondDir.BEGINWEDGE)
    # this shouldn't throw an exception
    RegistrationHash.GetMolLayers(mol)

  def test_tautomer_v2_hash(self):
    """
        Check that using v2 of the tautomer hash works
    """
    enol = Chem.MolFromSmiles('CC=CO')
    keto = Chem.MolFromSmiles('CCC=O')

    # Default, v1 of the tautomer hash:
    enol_layers = RegistrationHash.GetMolLayers(enol)
    keto_layers = RegistrationHash.GetMolLayers(keto)

    for layer in (RegistrationHash.HashLayer.TAUTOMER_HASH,
                  RegistrationHash.HashLayer.NO_STEREO_TAUTOMER_HASH):
      self.assertNotEqual(enol_layers[layer], keto_layers[layer])

    self.assertNotEqual(
      RegistrationHash.GetMolHash(
        enol_layers, hash_scheme=RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS),
      RegistrationHash.GetMolHash(
        keto_layers, hash_scheme=RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS))

    # v2 of the tautomer hash:
    enol_layers = RegistrationHash.GetMolLayers(enol, enable_tautomer_hash_v2=True)
    self.assertIn(RegistrationHash.HashLayer.TAUTOMER_HASH, enol_layers)

    keto_layers = RegistrationHash.GetMolLayers(keto, enable_tautomer_hash_v2=True)

    for layer in (RegistrationHash.HashLayer.TAUTOMER_HASH,
                  RegistrationHash.HashLayer.NO_STEREO_TAUTOMER_HASH):
      self.assertEqual(enol_layers[layer], keto_layers[layer])

    self.assertEqual(
      RegistrationHash.GetMolHash(
        enol_layers, hash_scheme=RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS),
      RegistrationHash.GetMolHash(
        keto_layers, hash_scheme=RegistrationHash.HashScheme.TAUTOMER_INSENSITIVE_LAYERS))

  def testAtropisomer(self):
    atropTestFile = os.path.join(
      RDConfig.RDBaseDir, 'Code/GraphMol/FileParsers/test_data/atropisomers/RP-6306_atrop1.sdf')
    mol = Chem.MolFromMolFile(atropTestFile)

    bond = mol.GetBondWithIdx(3)
    self.assertEqual(bond.GetStereo(), Chem.BondStereo.STEREOATROPCW)

    layersCw = RegistrationHash.GetMolLayers(mol)

    smiles1, smiExt1 = layersCw[HashLayer.CANONICAL_SMILES].split()
    taut1, tautExt1 = layersCw[HashLayer.TAUTOMER_HASH].split()
    self.assertEqual(smiles1, 'Cc1cc2c(C(N)=O)c(N)n(-c3c(C)ccc(O)c3C)c2nc1C')
    self.assertEqual(
      taut1,
      'C[C]1[CH][C]2[C]([C]([N])[O])[C]([N])N([C]3[C](C)[CH][CH][C]([O])[C]3C)[C]2[N][C]1C_5_0')
    self.assertEqual(smiExt1, '|wD:10.9|')
    self.assertEqual(smiExt1, tautExt1)

    # Now look at the other atropisomer
    bond = mol.GetBondWithIdx(3)
    bond.SetStereo(Chem.BondStereo.STEREOATROPCCW)

    # Hashes should not match
    layersCcw = RegistrationHash.GetMolLayers(mol)
    self.assertNotEqual(layersCw, layersCcw)

    smiles2, smiExt2 = layersCcw[HashLayer.CANONICAL_SMILES].split()
    taut2, tautExt2 = layersCcw[HashLayer.TAUTOMER_HASH].split()

    # SMILES and tautomer hash should be the same, but the extensions must be different
    self.assertEqual(smiles2, smiles1)
    self.assertEqual(taut2, taut1)
    self.assertEqual(smiExt2, '|wD:10.20|')  # same atom, but "down" wedge on a different bond
    self.assertEqual(smiExt2, tautExt2)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
