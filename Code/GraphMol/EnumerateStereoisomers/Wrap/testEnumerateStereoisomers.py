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

  def atestEnumerateStereoisomersBasic(self):
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


  def atestEnumerateStereoisomersWithCrazyNumberOfCenters(self):
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

  def testIsomerSets(self):
    mol = Chem.MolFromSmiles("C[C@@H]1N[C@H](C)[C@@H]([C@H](C)[C@@H]1C)C1[C@@H](C)O[C@@H](C)[C@@H](C)[C@H]1C |a:5,o1:1,8,o2:14,16,&1:11,18,&2:3,6,r|")
    enum = rdEnumerateStereoisomers.StereoisomerEnumerator(mol)
    sets = enum.GetStereoisomerSets()
    self.assertEqual(enum.GetStereoisomerCount() , 32)
    self.assertEqual(sets.numIsomerSets , 4)
    self.assertEqual(sets.numIsomersInSet , 8)
    self.assertEqual(sets.numIsomers , 32)

    isomer_sets = []
    for i in range(sets.numIsomerSets):
      isomer_sets.append(set())
      
    for isomer in enum:
      s = isomer.GetIntProp("isomerSet")
      isomer_sets[s].add(Chem.MolToSmiles(isomer))

    assert isomer_sets == [
      {'C[C@@H]1[C@@H](C)[C@@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@@H](C)[C@@H](C)O[C@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)N[C@H](C)[C@@H]1[C@H]1[C@@H](C)[C@@H](C)[C@@H](C)O[C@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)O[C@H](C)[C@H]1[C@@H]1[C@@H](C)[C@@H](C)[C@@H](C)N[C@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)O[C@H](C)[C@H]1[C@@H]1[C@H](C)[C@@H](C)[C@@H](C)N[C@@H]1C',
       'C[C@@H]1[C@H](C)[C@@H]([C@H]2[C@@H](C)[C@@H](C)[C@@H](C)O[C@H]2C)[C@H](C)N[C@@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@H](C)[C@@H](C)[C@@H](C)O[C@@H]1C',
       'C[C@@H]1[C@H](C)[C@@H]([C@@H]2[C@@H](C)[C@@H](C)[C@@H](C)O[C@H]2C)[C@H](C)N[C@@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)O[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@@H](C)[C@@H](C)N[C@H]1C'},
      {'C[C@@H]1[C@H](C)[C@H](C)N[C@H](C)[C@@H]1[C@H]1[C@@H](C)[C@@H](C)[C@@H](C)O[C@H]1C',
       'C[C@@H]1[C@H](C)[C@@H]([C@@H]2[C@@H](C)[C@H](C)[C@H](C)N[C@H]2C)[C@H](C)O[C@@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)O[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)N[C@H]1C',
       'C[C@@H]1[C@H](C)[C@H]([C@@H]2[C@@H](C)[C@H](C)[C@H](C)N[C@H]2C)[C@H](C)O[C@@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@H](C)[C@@H](C)[C@@H](C)O[C@@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@@H](C)[C@@H](C)O[C@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)O[C@H](C)[C@H]1[C@@H]1[C@H](C)[C@H](C)[C@H](C)N[C@@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)O[C@H](C)[C@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)N[C@H]1C'},
      {'C[C@@H]1[C@H](C)[C@@H]([C@@H]2[C@@H](C)[C@H](C)[C@H](C)O[C@H]2C)[C@H](C)N[C@@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)O[C@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)N[C@H](C)[C@@H]1[C@H]1[C@@H](C)[C@H](C)[C@H](C)O[C@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)O[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@@H](C)[C@@H](C)N[C@H]1C',
       'C[C@@H]1[C@H](C)[C@@H]([C@H]2[C@@H](C)[C@H](C)[C@H](C)O[C@H]2C)[C@H](C)N[C@@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)O[C@H](C)[C@H]1[C@@H]1[C@H](C)[C@@H](C)[C@@H](C)N[C@@H]1C',
       'C[C@@H]1[C@@H](C)[C@@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@H](C)[C@H](C)[C@H](C)O[C@@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)O[C@H](C)[C@H]1[C@@H]1[C@@H](C)[C@@H](C)[C@@H](C)N[C@H]1C'},
      {'C[C@@H]1[C@H](C)[C@H](C)O[C@H](C)[C@H]1[C@@H]1[C@H](C)[C@H](C)[C@H](C)N[C@@H]1C',
       'C[C@H]1[C@H](C)[C@H](C)O[C@@H](C)[C@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)N[C@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)N[C@H](C)[C@@H]1[C@H]1[C@@H](C)[C@H](C)[C@H](C)O[C@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)O[C@H](C)[C@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)N[C@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@H](C)[C@H](C)[C@H](C)O[C@@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)O[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)N[C@H]1C',
       'C[C@H]1[C@H](C)[C@H](C)N[C@@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)O[C@H]1C',
       'C[C@@H]1[C@H](C)[C@H](C)N[C@H](C)[C@@H]1[C@@H]1[C@@H](C)[C@H](C)[C@H](C)O[C@H]1C'}]
    


    


if __name__ == '__main__':
  unittest.main()
