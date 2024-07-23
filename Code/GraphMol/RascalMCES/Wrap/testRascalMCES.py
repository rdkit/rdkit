#  Copyright (c) 2023 David Cosgrove and other RDKit contributors
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
import os
import unittest

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdRascalMCES


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    mol1 = Chem.MolFromSmiles("c1ccccc1Cl")
    mol2 = Chem.MolFromSmiles("c1ccccc1F")
    opts = rdRascalMCES.RascalOptions()

    results = rdRascalMCES.FindMCES(mol1, mol2, opts)
    self.assertEqual(len(results), 1)
    self.assertEqual(results[0].smartsString, 'c1:c:c:c:c:c:1')
    self.assertEqual(len(results[0].bondMatches()), 6)
    self.assertEqual(len(results[0].atomMatches()), 6)

  def test2(self):
    # Test single largest fragment extraction from results
    ad1 = Chem.MolFromSmiles("CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934")
    ad2 = Chem.MolFromSmiles("N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336")

    opts = rdRascalMCES.RascalOptions()
    results = rdRascalMCES.FindMCES(ad1, ad2, opts)
    self.assertEqual(len(results), 1)
    self.assertEqual(results[0].smartsString,
                     'N(-C)-c1:c:c:c(-CC(=O)-NCCCCCCCCCC):c:c:1.NC12CC3CC(-C1)-CC(-C2)-C3')
    results[0].largestFragmentOnly()
    self.assertEqual(results[0].smartsString, 'N(-C)-c1:c:c:c(-CC(=O)-NCCCCCCCCCC):c:c:1')

  def test3(self):
    # Test not specifying options
    mol1 = Chem.MolFromSmiles("c1ccccc1Cl")
    mol2 = Chem.MolFromSmiles("c1ccccc1F")

    results = rdRascalMCES.FindMCES(mol1, mol2)
    self.assertEqual(len(results), 1)
    self.assertEqual(results[0].smartsString, 'c1:c:c:c:c:c:1')
    self.assertEqual(len(results[0].bondMatches()), 6)
    self.assertEqual(len(results[0].atomMatches()), 6)

  def test4(self):
    # Test setting non-default option
    mol1 = Chem.MolFromSmiles('Oc1cccc2C(=O)C=CC(=O)c12')
    mol2 = Chem.MolFromSmiles('O1C(=O)C=Cc2cc(OC)c(O)cc12')
    results = rdRascalMCES.FindMCES(mol1, mol2)
    self.assertEqual(len(results), 0)

    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.5
    results = rdRascalMCES.FindMCES(mol1, mol2, opts)
    self.assertEqual(len(results), 1)

  def test5(self):
    # Test setting non-default option singleLargestFrag
    ad1 = Chem.MolFromSmiles("CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934")
    ad2 = Chem.MolFromSmiles("N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336")

    opts = rdRascalMCES.RascalOptions()
    opts.singleLargestFrag = True
    results = rdRascalMCES.FindMCES(ad1, ad2, opts)
    self.assertEqual(len(results), 1)
    self.assertEqual(results[0].smartsString, 'CCCCCCCCCCNC12CC3CC(-C1)-CC(-C2)-C3')

  def test6(self):
    # Test the threshold and examine the tier1 and tier2 similarities.
    ad1 = Chem.MolFromSmiles("CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934")
    ad2 = Chem.MolFromSmiles("N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336")

    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.95
    results = rdRascalMCES.FindMCES(ad1, ad2, opts)
    self.assertEqual(len(results), 0)
    opts.returnEmptyMCES = True
    results = rdRascalMCES.FindMCES(ad1, ad2, opts)
    self.assertEqual(len(results), 1)

  def testRascalCluster(self):
    cdk2_file = Path(os.environ['RDBASE']) / 'Contrib' / 'Fastcluster' / 'cdk2.smi'
    suppl = Chem.SmilesMolSupplier(str(cdk2_file), '\t', 1, 0, False)
    mols = [mol for mol in suppl]
    clusters = rdRascalMCES.RascalCluster(mols)
    self.assertEqual(len(clusters), 8)
    expClusters = [7, 7, 6, 2, 2, 2, 2, 20]
    for clus, expClusSize in zip(clusters, expClusters):
      self.assertEqual(expClusSize, len(clus))

    clusOpts = rdRascalMCES.RascalClusterOptions()
    clusOpts.similarityCutoff = 0.6
    clusters = rdRascalMCES.RascalCluster(mols, clusOpts)
    expClusters = [9, 8, 6, 2, 2, 2, 2, 2, 2, 2, 11]
    for clus, expClusSize in zip(clusters, expClusters):
      self.assertEqual(expClusSize, len(clus))

  def testRascalButinaCluster(self):
    cdk2_file = Path(os.environ['RDBASE']) / 'Contrib' / 'Fastcluster' / 'cdk2.smi'
    suppl = Chem.SmilesMolSupplier(str(cdk2_file), '\t', 1, 0, False)
    mols = [mol for mol in suppl]
    clusters = rdRascalMCES.RascalButinaCluster(mols)
    self.assertEqual(len(clusters), 29)
    expClusters = [
      6, 6, 6, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
    ]
    for clus, expClusSize in zip(clusters, expClusters):
      self.assertEqual(expClusSize, len(clus))

  def testMaxBondMatchPairs(self):
    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.0
    opts.returnEmptyMCES = True
    opts.singleLargestFrag = True
    opts.allBestMCESs = True
    opts.completeAromaticRings = False
    opts.timeout = -1

    too_long_1 =  Chem.MolFromSmiles('CCCC=CCCCC=CCCCCCCCCCCCCCCCCCCCCC1CNCCC1')
    too_long_2 =  Chem.MolFromSmiles('CCCC=CCCCCC=CCCCCCCCCCCCCCCCCCCCC1CNCCC1')
    results = rdRascalMCES.FindMCES(too_long_1, too_long_2, opts)
    self.assertEqual(len(results[0].bondMatches()), 0)

    opts.maxBondMatchPairs = 1200
    results = rdRascalMCES.FindMCES(too_long_1, too_long_2, opts)
    self.assertEqual(len(results[0].bondMatches()), 26)

  def testExactConnectionsMatch(self):
    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.5
    opts.allBestMCESs = True
    mol1 = Chem.MolFromSmiles('c1ccccc1C1CCC(C(C)C)C1')
    mol2 = Chem.MolFromSmiles('c1ccccc1C(C)C')
    results = rdRascalMCES.FindMCES(mol1, mol2, opts)
    self.assertEqual(results[0].numFragments, 1)
    self.assertEqual(results[0].smartsString, 'c1:c:c:c:c:c:1-C(-C)-C')

    opts.exactConnectionsMatch = True
    results = rdRascalMCES.FindMCES(mol1, mol2, opts)
    self.assertEqual(results[0].numFragments, 2)
    self.assertEqual(results[0].smartsString,
                     '[#6&a&D2]1:[#6&a&D2]:[#6&a&D2]:[#6&a&D2]:[#6&a&D2]:[#6&a&D3]:1.[#6&A&D3](-[#6&A&D1])-[#6&A&D1]')

  def testEquivalentAtoms(self):
    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.5
    opts.equivalentAtoms = "[F,Cl,Br,I]"
    mol1 = Chem.MolFromSmiles('c1ccccc1F')
    mol2 = Chem.MolFromSmiles('c1ccccc1Br')
    results = rdRascalMCES.FindMCES(mol1, mol2, opts)
    self.assertEqual(results[0].numFragments, 1)
    self.assertEqual(results[0].smartsString, 'c1:c:c:c:c:c:1-[F,Cl,Br,I]')

  def testEquivalentBonds(self):
    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.5
    opts.ignoreBondOrders = True
    mol1 = Chem.MolFromSmiles('CC=CC')
    mol2 = Chem.MolFromSmiles('CCCC')
    results = rdRascalMCES.FindMCES(mol1, mol2, opts)
    self.assertEqual(results[0].numFragments, 1)
    self.assertEqual(results[0].smartsString, 'C~C~C~C')
    
    
if __name__ == "__main__":
  unittest.main()
