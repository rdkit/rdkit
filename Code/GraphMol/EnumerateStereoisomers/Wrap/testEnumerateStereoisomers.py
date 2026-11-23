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

import os
from rdkit import Chem
from rdkit.Chem import rdEnumerateStereoisomers

class TestCase(unittest.TestCase):

  def testEnumerateStereoisomersBasic(self):
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C')
    smiles = set()
    enum = rdEnumerateStereoisomers.StereoisomerEnumerator(mol)
    while True:
        iso = enum.next()
        if iso is None:
            break
        smiles.add(Chem.MolToSmiles(iso))
    self.assertEqual(len(smiles), 4)

    # Now do the same again, but iterating over the StereoisomerEnumerator
    enum2 = rdEnumerateStereoisomers.StereoisomerEnumerator(mol)
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

    undefinedCenterIdx = 13
    undefinedCenter = base.GetAtomWithIdx(undefinedCenterIdx)
    self.assertEqual(undefinedCenter.GetTotalDegree(), 4)
    self.assertEqual(undefinedCenter.GetChiralTag(), Chem.ChiralType.CHI_UNSPECIFIED)

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
    firstTwoTetrahedrals = centers[0:2]
    lastAndUndefinedTetrahedrals = [centers[-1], undefinedCenterIdx]
    for group_type, atom_indices, expected_cxsmiles in (
        (Chem.StereoGroupType.STEREO_ABSOLUTE, firstTwoTetrahedrals, {
          # In this case, nothing is enumerated, because the group is absolute.
          # No flippers are detected, so we just return the original molecule
          # (including the ABS stereo group).
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C |a:1,4|'
        }),
        (Chem.StereoGroupType.STEREO_AND, firstTwoTetrahedrals, {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C',
          'C[C@@H](F)C[C@@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C'
        }),
        (Chem.StereoGroupType.STEREO_OR, firstTwoTetrahedrals, {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C',
          'C[C@@H](F)C[C@@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C'
        }),

        # Undefined centers are not enumerated when in a stereo groups.
        (Chem.StereoGroupType.STEREO_ABSOLUTE, lastAndUndefinedTetrahedrals, {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C |a:10|'
        }),
        (Chem.StereoGroupType.STEREO_AND, lastAndUndefinedTetrahedrals, {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C',
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@@H](I)CC(O)C'
        }),
        (Chem.StereoGroupType.STEREO_OR, lastAndUndefinedTetrahedrals, {
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@H](I)CC(O)C',
          'C[C@H](F)C[C@H](Cl)C[C@H](Br)C[C@@H](I)CC(O)C'
        }),
        ):
      mol = mol_with_groups([(group_type, atom_indices)])
      isomers, cxsmiles = enumerate_cxsmiles(mol)
      self.assertEqual(len(isomers), len(expected_cxsmiles))
      self.assertEqual(cxsmiles, expected_cxsmiles)

    # Multiple group types are handled independently. The ABS center and the
    # ungrouped center stay fixed while the AND and OR groups can each flip.
    mol = mol_with_groups([
      (Chem.StereoGroupType.STEREO_ABSOLUTE, [centers[0]]),
      (Chem.StereoGroupType.STEREO_AND, centers[1:3]),
      (Chem.StereoGroupType.STEREO_OR, lastAndUndefinedTetrahedrals),
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

  def testMaxIsomersOptionDefault(self):
    # Check that there's a maxIsomers value even when
    # not explicitly specified. We don't check the actual
    # enumeration of this many isomers because it takes a while.
    # Other tests do check that this limit is correctly enforced.
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    self.assertEqual(opts.maxIsomers, 16384)

  def testEnumerateStereoisomersWithCrazyNumberOfCenters(self):
    # insanely large numbers of isomers aren't a problem
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C' * 31)
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

    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
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

  def testGithub4144(self):
    """ github4144: EnumerateStereiosomers clearing ring info """
    m = Chem.MolFromSmiles('CSCc1cnc(C=Nn2c(C)nc3sc4c(c3c2=O)CCCCC4)s1')
    sssr = [list(x) for x in Chem.GetSymmSSSR(m)]
    sms = rdEnumerateStereoisomers.StereoisomerEnumerator(m)
    for sm in sms:
      sssr2 = [list(x) for x in Chem.GetSymmSSSR(sm)]
      self.assertEqual(sssr, sssr2)

  def testWigglyBondsAndEnumerateStereoisomers(self):

    m = Chem.MolFromSmiles('CC=CC |w:1.0|')
    sms = [Chem.MolToSmiles(x) for x in rdEnumerateStereoisomers.StereoisomerEnumerator(m)]
    self.assertEqual(len(sms), 2)
    self.assertIn('C/C=C/C', sms)
    self.assertIn('C/C=C\\C', sms)

    m = Chem.MolFromSmiles('CC(F)(Cl)(Br) |w:1.0|')
    sms = [Chem.MolToSmiles(x) for x in rdEnumerateStereoisomers.StereoisomerEnumerator(m)]
    self.assertEqual(len(sms), 2)
    self.assertIn('C[C@](F)(Cl)Br', sms)
    self.assertIn('C[C@@](F)(Cl)Br', sms)

  def testEnumerateStereoisomersRandomSamplingShouldBeDeterministicAndPortable(self):
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C=C(Br)C(I)N')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.maxIsomers=2
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    expected = set([
      'C/C(F)=C\\[C@@H](Cl)/C=C(\\Br)[C@H](N)I',
      'C/C(F)=C\\[C@H](Cl)/C=C(\\Br)[C@H](N)I',
    ])
    self.assertEqual(smiles, expected)

  def testEnumerateStereoisomersMaxIsomersShouldBeReturnedEvenWithTryEmbedding(self):
    mol = Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.tryEmbedding=True
    opts.maxIsomers=8
    isomers = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      isomers.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(len(isomers), 8)

  def testEnumerateStereoisomersTryEmbeddingShouldNotInfiniteLoopWhenMaxIsomersIsLargerThanActual(
      self):
    mol = Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.tryEmbedding=True
    opts.maxIsomers=1024
    isomers = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      isomers.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(len(isomers), 8)

  def testEnumerateStereoisomersRandomSeeding(self):
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.randomSeed=-1
    opts.maxIsomers=3
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C')
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))

    expected = {
      'C/C(F)=C\\[C@@H](C)Cl',
      'C/C(F)=C/[C@H](C)Cl',
      'C/C(F)=C\\[C@H](C)Cl',
      'C/C(F)=C/[C@@H](C)Cl'
    }

    self.assertEqual(len(smiles), 3)
    self.assertTrue(expected.issuperset(smiles))

    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.randomSeed=0xBEEF
    mol = Chem.MolFromSmiles('c1ccc2c(c1)C(=O)N(C2=O)C3CCC(=O)NC3=O')
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(
      smiles,
      set([
        'O=C1CC[C@@H](N2C(=O)c3ccccc3C2=O)C(=O)N1',
        'O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1',
      ]))

    mol = Chem.MolFromSmiles('CCCC(=C(CCl)C(C)CBr)[C@H](F)C(C)C')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    smiles = list(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in rdEnumerateStereoisomers.StereoisomerEnumerator(mol))

    self.assertEqual(smiles, [
      'CCC/C(=C(\\CCl)[C@H](C)CBr)[C@H](F)C(C)C',
      'CCC/C(=C(/CCl)[C@H](C)CBr)[C@H](F)C(C)C',
      'CCC/C(=C(\\CCl)[C@@H](C)CBr)[C@H](F)C(C)C',
      'CCC/C(=C(/CCl)[C@@H](C)CBr)[C@H](F)C(C)C',
    ])

  def testEnumerateStereoisomersOnlyUnassigned(self):
    # shouldn't enumerate anything
    fully_assigned = Chem.MolFromSmiles('C/C(F)=C/[C@@H](C)Cl')
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(fully_assigned):
      self.assertEqual(m.HasProp('_MolFileChiralFlag'), False)
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(smiles, set(['C/C(F)=C/[C@@H](C)Cl']))

    # should only enumerate the bond stereo
    partially_assigned = Chem.MolFromSmiles('CC(F)=C[C@@H](C)Cl')
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(partially_assigned):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(smiles, set([
      'C/C(F)=C/[C@@H](C)Cl',
      'C/C(F)=C\\[C@@H](C)Cl',
    ]))

    # should enumerate everything
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.onlyUnassigned=False
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(fully_assigned, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(
      smiles,
      set([
        'C/C(F)=C\\[C@@H](C)Cl',
        'C/C(F)=C\\[C@H](C)Cl',
        'C/C(F)=C/[C@H](C)Cl',
        'C/C(F)=C/[C@@H](C)Cl',
      ]))

  def testEnumerateStereoisomersOnlyUnique(self):
    mol = Chem.MolFromSmiles('FC(Cl)C(Cl)F')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.unique=False
    smiles = []
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.append(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(len(smiles), 4)
    self.assertEqual(len(set(smiles)), 3)

    mol = Chem.MolFromSmiles('FC(Cl)C(Cl)F')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.unique=True
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(
      smiles, set(['F[C@@H](Cl)[C@@H](F)Cl', 'F[C@H](Cl)[C@@H](F)Cl', 'F[C@H](Cl)[C@H](F)Cl']))

    mol = Chem.MolFromSmiles('CC=CC=CC')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.unique=False
    smiles = []
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.append(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(len(smiles), 4)
    self.assertEqual(len(set(smiles)), 3)

    mol = Chem.MolFromSmiles('CC=CC=CC')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.unique=True
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(smiles, set(['C/C=C/C=C/C', 'C/C=C\\C=C\\C', 'C/C=C\\C=C/C']))

    mol = Chem.MolFromSmiles('FC(Cl)C=CC=CC(F)Cl')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.unique=False
    smiles = []
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.append(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(len(smiles), 16)
    self.assertEqual(len(set(smiles)), 10)

    mol = Chem.MolFromSmiles('FC(Cl)C=CC=CC(F)Cl')
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.unique=True
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(
      smiles,
      set(
        sorted([
          'F[C@H](Cl)/C=C/C=C/[C@@H](F)Cl', 'F[C@H](Cl)/C=C/C=C/[C@H](F)Cl',
          'F[C@@H](Cl)/C=C/C=C/[C@@H](F)Cl', 'F[C@@H](Cl)/C=C\\C=C\\[C@@H](F)Cl',
          'F[C@H](Cl)/C=C\\C=C/[C@@H](F)Cl', 'F[C@H](Cl)/C=C/C=C\\[C@@H](F)Cl',
          'F[C@H](Cl)/C=C\\C=C\\[C@H](F)Cl', 'F[C@@H](Cl)/C=C\\C=C/[C@@H](F)Cl',
          'F[C@H](Cl)/C=C\\C=C/[C@H](F)Cl', 'F[C@H](Cl)/C=C\\C=C\\[C@@H](F)Cl'
        ])))

  def testEnumerateStereoisomersOnlyEnhancedStereo(self):
    rdbase = os.environ["RDBASE"]
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/two_centers_or.mol')
    mol = Chem.MolFromMolFile(filename)
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m))
    # switches the centers linked by an "OR", but not the absolute group
    self.assertEqual(smiles, {r'C[C@H]([C@@H](C)F)[C@@H](C)Br', r'C[C@@H]([C@H](C)F)[C@@H](C)Br'})

    # we need the SMILES without the influence of the stereo groups:
    mol = Chem.RWMol(mol)
    mol.SetStereoGroups([])
    original_smiles = Chem.MolToSmiles(mol)
    self.assertIn(original_smiles, smiles)

  def testNoExtrasIfEnumeratingAllWithEnhancedStereo(self):
    """
    If the onlyUnassigned option is False, make sure that enhanced stereo
    groups aren't double-counted.
    """
    rdbase = os.environ["RDBASE"]
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/two_centers_or.mol')
    mol = Chem.MolFromMolFile(filename)

    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.onlyUnassigned=False
    opts.unique=False
    smiles = []
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.append(Chem.MolToSmiles(m))
    self.assertEqual(len(smiles), len(set(smiles)))
    self.assertEqual(len(smiles), 2**3)

  def testIssue2890(self):
    mol = Chem.MolFromSmiles('CC=CC')
    mol.GetBondWithIdx(1).SetStereo(Chem.rdchem.BondStereo.STEREOANY)
    isomers = []
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      isomers.append(m)
    self.assertEqual(len(isomers), 2)

  def testIssue3231(self):
    mol = Chem.MolFromSmiles(
      'C[C@H](OC1=C(N)N=CC(C2=CN(C3C[C@H](C)NCC3)N=C2)=C1)C4=C(Cl)C=CC(F)=C4Cl')
    Chem.AssignStereochemistry(mol, force=True, flagPossibleStereoCenters=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    self.assertEqual(centers, [(1, 'S'), (12, '?'), (14, 'S')])
    enumsi_opt = rdEnumerateStereoisomers.StereoEnumerationOptions()
    enumsi_opt.maxIsomers=20
    enumsi_opt.onlyUnassigned=False
    isomers = list(rdEnumerateStereoisomers.StereoisomerEnumerator(mol, enumsi_opt))
    chi_cents = []
    for iso in isomers:
      self.assertEqual(iso.GetProp('_MolFileChiralFlag'), '1')
      Chem.AssignStereochemistry(iso)
      chi_cents.append(Chem.FindMolChiralCenters(iso, includeUnassigned=True))
    self.assertEqual(sorted(chi_cents), [
      [(1, 'R'), (12, 'R'), (14, 'R')],
      [(1, 'R'), (12, 'R'), (14, 'S')],
      [(1, 'R'), (12, 'S'), (14, 'R')],
      [(1, 'R'), (12, 'S'), (14, 'S')],
      [(1, 'S'), (12, 'R'), (14, 'R')],
      [(1, 'S'), (12, 'R'), (14, 'S')],
      [(1, 'S'), (12, 'S'), (14, 'R')],
      [(1, 'S'), (12, 'S'), (14, 'S')],
    ])

  def testIssue3505(self):
    m = Chem.MolFromSmiles('CCC(C)Br')
    mols = list(rdEnumerateStereoisomers.StereoisomerEnumerator(m))
    self.assertEqual(len(mols), 2)
    for mol in mols:
      self.assertEqual(mol.GetProp('_MolFileChiralFlag'), '1')
      at = mol.GetAtomWithIdx(2)
      self.assertIn(at.GetChiralTag(),
                    [Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW])
      self.assertTrue(at.HasProp("_ChiralityPossible"))

  def testEnumerateEitherDoubleStereo(self):
    """ EnumerateStereoisomers from MOL with explicit either cis/trans bond """
    rdbase = os.environ["RDBASE"]
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/simple_either.mol')
    mol = Chem.MolFromMolFile(filename)
    smiles = set()
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(m, isomericSmiles=True))
    self.assertEqual(smiles, {"C/C=C/C", "C/C=C\\C"})

  def testTryEmbeddingManyChirals(self):
    smiles = "C1" + "C(Cl)(Br)" * 40 + "C1"
    mol = Chem.MolFromSmiles(smiles)
    opts = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opts.tryEmbedding=True
    opts.maxIsomers=2
    isomers = []
    for m in rdEnumerateStereoisomers.StereoisomerEnumerator(mol, opts):
      self.assertEqual(m.GetProp('_MolFileChiralFlag'), '1')
      isomers.append(m)
    self.assertEqual(len(isomers), 2)

  def testGithub6045(self):
    mol = Chem.MolFromSmiles('O[C@H](Br)[C@H](F)C |&1:1,3|')
    prods = list(rdEnumerateStereoisomers.StereoisomerEnumerator(mol))
    self.assertEqual(len(prods), 2)
    for prod in prods:
      self.assertEqual(prod.GetProp('_MolFileChiralFlag'), '1')
      self.assertEqual(len(prod.GetStereoGroups()), 0)

  def testGithub7516(self):
    m = Chem.MolFromSmiles('CC1CC(C)C1')
    sis = list(rdEnumerateStereoisomers.StereoisomerEnumerator(m))
    self.assertEqual(len(sis), 2)
    self.assertEqual(Chem.MolToSmiles(sis[0]), 'C[C@H]1C[C@H](C)C1')
    self.assertEqual(Chem.MolToSmiles(sis[1]), 'C[C@H]1C[C@@H](C)C1')

    m = Chem.MolFromSmiles('COC(=O)C1CC(NC(N)=O)C1')
    sis = list(rdEnumerateStereoisomers.StereoisomerEnumerator(m))
    self.assertEqual(len(sis), 2)

    m = Chem.MolFromSmiles('O=C(NC1CC2[NH+](C(C1)CC2)Cc3ccccc3)N')
    opt = rdEnumerateStereoisomers.StereoEnumerationOptions()
    opt.useRingSystemFilter=False
    sis = list(rdEnumerateStereoisomers.StereoisomerEnumerator(m, opt))
    self.assertEqual(len(sis), 8)
    opt.useRingSystemFilter=True
    sis = list(rdEnumerateStereoisomers.StereoisomerEnumerator(m, opt))
    self.assertEqual(len(sis), 4)

  def testGithub7608(self):
    m = Chem.MolFromSmiles('N=C(N1CCC1)NCCOc2ccc(C(F)(F)F)cc2')
    sis = list(rdEnumerateStereoisomers.StereoisomerEnumerator(m))
    self.assertEqual(len(sis), 1)



if __name__ == '__main__':
  unittest.main()

