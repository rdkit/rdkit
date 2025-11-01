#  Copyright (c) 2024 David Cosgrove and other RDKit contributors
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
import tempfile
import unittest

from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import (rdSynthonSpaceSearch, rdFingerprintGenerator,
                        rdRascalMCES, rdGeneralizedSubstruct, rdMolDescriptors)


class TestCase(unittest.TestCase):

  def setUp(self):
    self.sssDir = Path(os.environ["RDBASE"]) / "Code" / "GraphMol" / "SynthonSpaceSearch" / "data"

  def testSubstructSearch(self):
    fName = self.sssDir / "idorsia_toy_space_a.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = 10
    query = Chem.MolFromSmarts("c1ccccc1C(=O)N1CCCC1")
    results = synthonspace.SubstructureSearch(query, params=params)
    self.assertEqual(10, len(results.GetHitMolecules()))
    smParams = Chem.SubstructMatchParameters()
    results = synthonspace.SubstructureSearch(query,
                                              substructMatchParams=smParams,
                                              params=params)
    self.assertEqual(10, len(results.GetHitMolecules()))

    # callback returns None, stil get all results
    mols = []
    synthonspace.SubstructureSearchIncremental(
            query,
            lambda results: mols.extend(results),
            substructMatchParams=smParams,
            params=params)
    self.assertEqual(10, len(mols))

    # callback returns True, get one chunk
    mols = []
    params.toTryChunkSize = 2
    def callback_returns_true(chunk):
        mols.extend(chunk)
        return True
    synthonspace.SubstructureSearchIncremental(
            query,
            callback_returns_true,
            substructMatchParams=smParams,
            params=params)
    self.assertEqual(2, len(mols))

    # Exceptions thrown in the callback propagate back here
    mols = []
    params.toTryChunkSize = 2
    def callback_raises(chunk):
        mols.extend(chunk)
        raise StopIteration

    try:
        synthonspace.SubstructureSearchIncremental(
                query,
                callback_raises,
                substructMatchParams=smParams,
                params=params)
    except StopIteration:
        pass
    else:
        assert False, "Expected exception"
    self.assertEqual(2, len(mols))



  def testFingerprintSearch(self):
    fName = self.sssDir / "idorsia_toy_space_a.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)
    self.assertEqual(6, synthonspace.GetNumReactions())
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = -1
    params.similarityCutoff = 0.45
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048, useBondOrder=True)
    results = synthonspace.FingerprintSearch(
      Chem.MolFromSmiles("O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"), fpgen, params)
    self.assertEqual(278, len(results.GetHitMolecules()))
    mols = []
    synthonspace.FingerprintSearchIncremental(
      Chem.MolFromSmiles("O=C(Nc1c(CNC=O)cc[s]1)c1nccnc1"),
      fpgen,
      lambda results: mols.extend(results),
      params)
    self.assertEqual(278, len(mols))
    

  def testEnumerate(self):
    fName = self.sssDir / "amide_space.txt"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadTextFile(fName)
    with tempfile.NamedTemporaryFile() as tmp:
      synthonspace.WriteEnumeratedFile(tmp.name)

  def testTimeOut(self):
    fName = self.sssDir / "idorsia_toy_space_a.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.timeOut = 1
    params.maxHits = -1
    params.similarityCutoff = 0.3
    params.fragSimilarityAdjuster = 0.3
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048, useBondOrder=True)
    results = synthonspace.FingerprintSearch(
      Chem.MolFromSmiles("c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"), fpgen, params)
    self.assertTrue(results.GetTimedOut())

  def testSynthonError(self):
    fName = self.sssDir / "amide_space_error.txt"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    self.assertRaises(RuntimeError, synthonspace.ReadTextFile, fName)

    fName = self.sssDir / "synthon_error.txt"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    self.assertRaises(RuntimeError, synthonspace.ReadTextFile, fName)

  def testRascalSearch(self):
    fName = self.sssDir / "Syntons_5567.csv"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadTextFile(fName)
    self.assertEqual(10, synthonspace.GetNumReactions())
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = 10
    rascalOpts = rdRascalMCES.RascalOptions()
    results = synthonspace.RascalSearch(Chem.MolFromSmiles("c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"),
                                             rascalOpts, params)
    self.assertEqual(10, len(results.GetHitMolecules()))

  def testExtendedSubsructureSearch(self):
    fName = self.sssDir / "extended_query.csv"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadTextFile(fName)
    self.assertEqual(7, synthonspace.GetNumReactions())
    m = Chem.MolFromSmarts('[#6]-*.c1nc2cccnc2n1 |m:1:3.10|')
    xqry = rdGeneralizedSubstruct.CreateExtendedQueryMol(m)
    results = synthonspace.SubstructureSearch(xqry)
    self.assertEqual(12, len(results.GetHitMolecules()))

  def testHeavyAtomCutoffs(self):
    fName = self.sssDir / "idorsia_toy_space_a.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = -1
    params.maxHitHeavyAtoms = 30
    params.minHitHeavyAtoms = 25
    q = Chem.MolFromSmarts("c1ccccc1C(=O)N1CCCC1")
    results = synthonspace.SubstructureSearch(q)
    self.assertEqual(220, len(results.GetHitMolecules()))
    results = synthonspace.SubstructureSearch(q, params=params)
    self.assertEqual(141, len(results.GetHitMolecules()))

  def testMolWeightCutoffs(self):
    fName = self.sssDir / "idorsia_toy_space_a.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = -1
    params.maxHitMolWt = 450
    params.minHitMolWt = 350
    q = Chem.MolFromSmarts("c1ccccc1C(=O)N1CCCC1")
    results = synthonspace.SubstructureSearch(q)
    self.assertEqual(220, len(results.GetHitMolecules()))
    results = synthonspace.SubstructureSearch(q, params=params)
    self.assertEqual(185, len(results.GetHitMolecules()))
    
  def testChiralAtomtCutoffs(self):
    fName = self.sssDir / "idorsia_toy_space_a.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = -1
    params.maxHitChiralAtoms = 2
    params.minHitChiralAtoms = 1
    q = Chem.MolFromSmarts("Cc1nc(-c2ccnn2CC(C))n([C@@H]2CCOC2)n1")
    results = synthonspace.SubstructureSearch(q)
    self.assertEqual(20, len(results.GetHitMolecules()))
    results = synthonspace.SubstructureSearch(q, params=params)
    self.assertEqual(10, len(results.GetHitMolecules()))
    
    
if __name__ == "__main__":
  unittest.main()
