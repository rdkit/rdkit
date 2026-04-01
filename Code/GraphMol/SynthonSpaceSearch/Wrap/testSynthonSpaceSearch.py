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
                        rdRascalMCES, rdGeneralizedSubstruct, rdMolDescriptors,
                        rdDistGeom)


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

  def testBestHit(self):
    fName = self.sssDir / "small_freedom_shapes.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)
    
    query = Chem.MolFromSmiles("CCC(C(=O)NCc1ccco1)N(Cc1sccc1C)C(C)C |(1.19967,-2.26511,-1.8853;-0.0674677,-1.53728,-1.54329;-0.0395195,-0.735921,-0.2957;0.978688,0.316536,-0.356493;0.610079,1.54481,-0.391485;2.37979,0.114038,-0.377756;3.31574,1.24811,-0.443092;4.723,0.774447,-0.458657;5.58509,0.534718,0.592748;6.76773,0.108395,-0.00740365;6.56938,0.109237,-1.38929;5.34763,0.509833,-1.60401;-1.33892,-0.260032,0.0922388;-2.08764,0.476264,-0.832263;-3.39845,0.92709,-0.383151;-4.99177,0.223473,-0.828611;-5.94415,1.34626,0.133161;-5.00918,2.1798,0.729651;-3.7235,1.96002,0.462253;-2.60368,2.81164,1.05297;-1.99138,-1.01276,1.09008;-1.10607,-0.965532,2.34165;-2.26411,-2.45544,0.780694)|")
    ssparams = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    ssparams.maxHits = -1
    ssparams.numThreads = 1
    ssparams.numConformers = 200
    ssparams.confRMSThreshold = 0.25
    ssparams.randomSeed = 0xdac
    
    ssparams.fragSimilarityAdjuster = 0.1
    ssparams.approxSimilarityAdjuster = 0.1
    ssparams.similarityCutoff = 1.0
    ssparams.timeOut = 0
    
    results = synthonspace.ShapeSearch(query, ssparams)
    print(f"Number of hits : {len(results.GetHitMolecules())}")
    self.assertEqual(len(results.GetHitMolecules()), 0)

    bestHit = results.GetBestHit()
    self.assertIsNotNone(bestHit)
    self.assertLess(float(bestHit.GetProp('Similarity')), 1.0)

  def testUserConfGen(self):
    def makeConformers(smiles, numConfs):
      mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
      rdDistGeom.EmbedMultipleConfs(mol, 10, randomSeed=0xdac)
      mol = Chem.RemoveHs(mol)
      return mol

    fName = self.sssDir / "amide_space.txt"

    buildParams = rdSynthonSpaceSearch.ShapeBuildParams()
    buildParams.numThreads = 1
    buildParams.setUserConformerGenerator(makeConformers)
    
    synthonspace2 = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace2.ReadTextFile(fName)
    synthonspace2.BuildSynthonShapes(buildParams)
    # There's not an easy test for this, but we can at least
    # check there's there correct number of reactions and products.
    self.assertEqual(synthonspace2.GetNumReactions(), 1)
    self.assertEqual(synthonspace2.GetNumProducts(), 12)

    
if __name__ == "__main__":
  unittest.main()
