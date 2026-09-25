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
    enum = rdEnumerateStereoisomers.StereoisomerEnumerator(mol)
    while True:
        iso = enum.next()
        if iso is None:
            break
        smiles.add(Chem.MolToSmiles(iso))
    self.assertEqual(len(smiles), 4)

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
          break;
      self.assertEqual(iso.GetProp('_MolFileChiralFlag'), '1')
      smiles.add(Chem.MolToSmiles(iso, isomericSmiles=True))
    self.assertEqual(len(smiles), 13)


if __name__ == '__main__':
  unittest.main()

