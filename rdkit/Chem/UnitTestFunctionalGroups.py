#
#  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
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
#  Created by Greg Landrum October 2006
#


import os.path
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import FunctionalGroups


class TestCase(unittest.TestCase):

  def test1Basics(self):
    txt = """
AcidChloride\tC(=O)Cl\tAcid Chloride
  AcidChloride.Benzoyl\tC(=O)(Cl)c1ccccc1\tBenzoyl
Amine\tN\tAmine
  Amine.Primary\t[N;H2]\tPrimary
    Amine.Primary.Aromatic\t[N;H2][a]\tPrimary Aromatic
  Amine.Aromatic\tN[a]\tAromatic
"""
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy(data=txt)
    self.assertTrue(hierarchy)
    self.assertEqual(len(hierarchy), 2)
    self.assertEqual(len(hierarchy[0]), 2)
    self.assertEqual(len(hierarchy[1]), 4)
    self.assertEqual(hierarchy[0].name, 'Acid Chloride')
    self.assertEqual(hierarchy[0].children[0].name, 'Benzoyl')
    self.assertEqual(hierarchy[0].label, 'AcidChloride')
    self.assertEqual(hierarchy[0].rxnSmarts, '')
    m = Chem.MolFromSmiles('ClC(=O)CCCNc1ccccc1')
    fp = FunctionalGroups.CreateMolFingerprint(m, hierarchy)
    self.assertEqual(fp, [1, 0, 1, 0, 0, 1])

    m = Chem.MolFromSmiles('OC(=O)CCC')
    fp = FunctionalGroups.CreateMolFingerprint(m, hierarchy)
    self.assertEqual(fp, [0, 0, 0, 0, 0, 0])

    # make sure we get the same hierarchy on the second call:
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy(data=txt)
    self.assertTrue(hierarchy)
    self.assertEqual(len(hierarchy), 2)
    self.assertEqual(len(hierarchy[0]), 2)
    self.assertEqual(len(hierarchy[1]), 4)
    self.assertEqual(hierarchy[0].name, 'Acid Chloride')
    self.assertEqual(hierarchy[0].children[0].name, 'Benzoyl')
    self.assertEqual(hierarchy[0].label, 'AcidChloride')
    self.assertEqual(hierarchy[0].rxnSmarts, '')

    # if we edit this hierarchy it doesn't affect the global one:
    hierarchy.pop(0)
    self.assertEqual(len(hierarchy[0]), 4)
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy(data=txt)
    self.assertTrue(hierarchy)
    self.assertEqual(len(hierarchy), 2)
    self.assertEqual(len(hierarchy[0]), 2)
    self.assertEqual(len(hierarchy[1]), 4)
    self.assertEqual(hierarchy[0].name, 'Acid Chloride')
    self.assertEqual(hierarchy[0].children[0].name, 'Benzoyl')
    self.assertEqual(hierarchy[0].label, 'AcidChloride')
    self.assertEqual(hierarchy[0].rxnSmarts, '')

    # and if we edit the global one and don't force, we get the edited one:
    FunctionalGroups.hierarchy.pop(0)
    self.assertEqual(len(FunctionalGroups.hierarchy[0]), 4)
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy(data=txt)
    self.assertTrue(hierarchy)
    self.assertEqual(len(hierarchy), 1)
    self.assertEqual(len(hierarchy[0]), 4)

    # but a force gets us back:
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy(data=txt, force=True)
    self.assertEqual(len(hierarchy), 2)
    self.assertEqual(len(hierarchy[0]), 2)
    self.assertEqual(len(hierarchy[1]), 4)

  def test2Comments(self):
    txt = """
AcidChloride\tC(=O)Cl\tAcid Chloride
  AcidChloride.Benzoyl\tC(=O)(Cl)c1ccccc1\tBenzoyl
Amine\tN\tAmine
  Amine.Primary\t[N;H2]\tPrimary
    //Amine.Primary.Aromatic\t[N;H2][a]\tPrimary Aromatic
  Amine.Aromatic\tN[a]\tAromatic
"""
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy(data=txt)
    self.assertTrue(hierarchy)
    self.assertEqual(len(hierarchy), 2)
    self.assertEqual(len(hierarchy[0]), 2)
    self.assertEqual(len(hierarchy[1]), 3)

  def test3Reactions(self):
    txt = """BoronicAcid\t[$(B-!@[#6])](O)(O)\tBoronic Acid\t[#6:1][B:2]([O:3])[O:4]>>[#6:1].[B:2]([O:3])[O:4]
  BoronicAcid.Aromatic\t[$(B-!@c)](O)(O)\tAromatic\t[c:1][B:2]([O:3])[O:4]>>[c:1].[B:2]([O:3])[O:4]
  BoronicAcid.Aliphatic\t[$(B-!@C)](O)(O)\tAliphatic\t[C:1][B:2]([O:3])[O:4]>>[C:1].[B:2]([O:3])[O:4]
  """
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy(data=txt)
    self.assertTrue(hierarchy)
    self.assertEqual(len(hierarchy), 1)
    self.assertEqual(len(hierarchy[0].children), 2)
    self.assertNotEqual(hierarchy[0].rxnSmarts, '')
    self.assertNotEqual(hierarchy[0].children[0].rxnSmarts, '')

  def test4Hs(self):
    hierarchy = FunctionalGroups.BuildFuncGroupHierarchy()

    inName = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'NCI_5K_TPSA.csv')
    with open(inName, 'r') as inF:
      ms = [Chem.MolFromSmiles(x.split(',')[0]) for x in inF if x[0] != '#']
    for m in ms:
      mh = Chem.AddHs(m)
      fp = FunctionalGroups.CreateMolFingerprint(m, hierarchy)
      fph = FunctionalGroups.CreateMolFingerprint(mh, hierarchy)
      if fp != fph:
        print(Chem.MolToSmiles(m))
        print(fp.ToBitString())
        print(fph.ToBitString())
      self.assertEqual(fp, fph)


if __name__ == '__main__':
  unittest.main()
