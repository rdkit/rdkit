#  Copyright (c) 2025 David Cosgrove and other RDKit contributors
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
#       products derived from this software without specific prior written
#       permission.
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

# These tests are just to check that the Python wrappers are working
# ok.  The bulk of the tests are in the C++ code.

import unittest

from rdkit import Chem
from rdkit.Chem import rdEnumerateStereoisomers

class TestCase(unittest.TestCase):

  def testEnumerateStereoisomersBasic(self):
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C')
    smiles = set()
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    enum = rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts)
    while True:
        iso = enum.next()
        if iso is None:
            break
        smiles.add(Chem.MolToSmiles(iso))
    self.assertEqual(len(smiles), 4)

    # Now do the same again, but iterating over the StereoisomerEnumerator
    enum2 = rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts)
    smiles2 = {Chem.MolToSmiles(iso) for iso in enum2}

    self.assertEqual(len(smiles2), 4)
    self.assertEqual(smiles, smiles2)

  def testEnumerateStereoisomersOnlyStereoGroups(self):
    # The last chiral center in the SMILES (the [CH]) is intentionally
    # left unspecified. This atom won't ever be enumerated, since
    # it's never included in any of the stereo groups.
    base = Chem.MolFromSmiles('C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)C[CH](O)C')
    centers = [idx for idx, _ in Chem.FindMolChiralCenters(base, includeCIP=False)]
    self.assertEqual(len(centers), 4)  # The unspecified center is not detected!

    def mol_with_groups(group_specs):
      mol = Chem.RWMol(base)
      groups = [
        Chem.rdchem.CreateStereoGroup(group_type, mol, atom_indices)
        for group_type, atom_indices in group_specs
      ]
      mol.SetStereoGroups(groups)
      return mol

    def enumerate_cxsmiles(mol):
      options = rdEnumerateStereoisomers.StereoEnumerationOptions()
      options.onlyStereoGroups = True
      isomers = [isomer for isomer in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, options)]

      # Do not canonicalize the SMILES to preserve the atom ordering
      # from the input and make it easier to compare enumerated parities.
      params = Chem.SmilesWriteParams()
      params.canonical = False
      params.rootedAtAtom = 0
      return isomers, {Chem.MolToCXSmiles(isomer, params) for isomer in isomers}

    # Absolute groups are fixed. A racemic or relative group can be inverted as
    # a whole; centers outside the group remain fixed.
    # Once enumerated, the racemic and relative groups should be stripped,
    # since the atoms in them have already been enumerated.
    for group_type, atom_indices, expected_cxsmiles in (
        (Chem.StereoGroupType.STEREO_ABSOLUTE, centers[0:2], {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C |a:1,4|'
        }),
        (Chem.StereoGroupType.STEREO_AND, centers[0:2], {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C',
          'C[C@@H](F)C[C@@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C'
        }),
        (Chem.StereoGroupType.STEREO_OR, centers[0:2], {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C',
          'C[C@@H](F)C[C@@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C'
        })):
      mol = mol_with_groups([(group_type, atom_indices)])
      isomers, cxsmiles = enumerate_cxsmiles(mol)
      self.assertEqual(len(isomers), len(expected_cxsmiles))
      self.assertEqual(cxsmiles, expected_cxsmiles)

    # Multiple group types are handled independently. The ABS center and the
    # ungrouped center stay fixed while the AND and OR groups can each flip.
    mol = mol_with_groups([
      (Chem.StereoGroupType.STEREO_ABSOLUTE, [centers[0]]),
      (Chem.StereoGroupType.STEREO_AND, centers[1:3]),
      (Chem.StereoGroupType.STEREO_OR, [centers[3]]),
    ])
    isomers, cxsmiles = enumerate_cxsmiles(mol)
    expected_cxsmiles = {
      'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C',
      'C[C@H](F)C[C@@H](Cl)C[C@@H](Br)C[C@H](I)CC(O)C',
      'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@@H](I)CC(O)C',
      'C[C@H](F)C[C@@H](Cl)C[C@@H](Br)C[C@@H](I)CC(O)C',
    }
    self.assertEqual(len(isomers), 4)
    self.assertEqual(cxsmiles, expected_cxsmiles)

  def testEnumerateStereoisomersEmbeddingNotRejectingInput(self):
    mol = Chem.MolFromSmiles('C1C/C=C/CCC[C@H]1C(=O)Nc2ccc(cc2)C[C@@H](C(=O)[O-])[NH3+]')
    canon_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    smiles = set()
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.tryEmbedding=True
    enum = rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts)
    while True:
        iso = enum.next()
        if iso is None:
            break
        smiles.add(Chem.MolToSmiles(iso, isomericSmiles=True))
    self.assertEqual(len(smiles), 1)
    self.assertIn(canon_smiles, smiles)

  def testEnumerateStereoisomersWithCrazyNumberOfCenters(self):
    # insanely large numbers of isomers aren't a problem
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C' * 101)
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.maxIsomers=13
    smiles = set()
    enum = rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts)
    while True:
      iso = enum.next()
      if iso is None:
          break
      self.assertEqual(iso.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(iso, isomericSmiles=True))
    self.assertEqual(len(smiles), 13)

  def testEnumerateStereoisomersWithRingPatterns(self):

    def _getIsomers(mol, opts):
      smiles = list()
      enum = rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts)

      # This is the maximum number of stereoisomers the mol can be enumerated
      # into. Depending on the options we use, some of these might get discarded
      # because they are impossible in 3D space, but this number should
      # always be the same, independently of the options, since it only
      # depends on the number of stereo features on the molecule.
      self.assertEqual(enum.GetStereoisomerCount(), 16)

      while True:
        iso = enum.next()
        if iso is None:
            break
        smiles.append(iso)
      return smiles

    mol = Chem.MolFromSmiles("BrC=CC1OC(C2)(F)C2(Cl)C1")

    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()

    # Check the "old" behavior first
    opts.useRingSystemFilter = False  # the default is to use them
    allIsomers = _getIsomers(mol, opts)

    self.assertEqual(len(allIsomers), 16)

    # Now check the behavior when patterns are enabled
    opts.useRingSystemFilter = True  # restore the default value
    validIsomers = _getIsomers(mol, opts)

    self.assertEqual(len(validIsomers), 8)


if __name__ == '__main__':
  unittest.main()

