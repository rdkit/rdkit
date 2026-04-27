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
                        rdDistGeom, rdGaussianShape)


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

  def testExcludedVolume(self):
    excVolMol = Chem.MolFromSmiles("C.C=O.CC.CC.CC(N)=O.CCC.CCC.CCC.CCC(=O)N[C@@H](C)C(=O)NCC=O.CCC(C)C.CCC(C)C.CCCC(=O)N[C@H](CN)CO.CCCCNC=O.CCCN[C@H](C=O)[C@@H](C)O.CCNC(=N)N.CNCC(N)=O.CNCCN.N.N.N=C(N)N.NC[C@H](CC(=O)O)NC=O.O.O.c1c[nH]cn1.CCC(C)C |(23.735,-5.551,4.831;10.558,-20.222,6.049;10.782,-19.593,7.082;25.756,-7.961,1.508;24.4,-8.595,1.992;7.488,-19.618,10.529;8.946,-19.661,10.098;16.988,-18.147,6.134;16.576,-18.196,7.564;15.432,-17.604,7.854;17.24,-18.815,8.408;15.229,-10.345,10.204;14.148,-9.909,9.204;13.16,-11.048,8.967;22.155,-3.456,-2.463;20.831,-2.784,-2.822;19.665,-3.765,-2.706;13.808,-17.366,17.137;12.337,-17.174,16.743;12.24,-16.466,15.411;18.52,-8.706,8.462;19.494,-9.564,7.605;19.097,-9.56,6.108;17.957,-9.849,5.775;20.037,-9.291,5.223;19.811,-9.168,3.789;20.952,-8.374,3.183;19.66,-10.496,3.038;20.161,-11.523,3.493;19.033,-10.416,1.856;18.853,-11.525,0.927;17.429,-11.999,0.744;16.5,-11.418,1.302;20.813,-11.741,-2.325;22.348,-11.427,-2.65;23.303,-11.616,-1.434;24.76,-11.817,-1.924;23.177,-10.46,-0.471;15.598,-14.535,14.266;14.867,-13.332,13.467;13.366,-13.102,13.804;13.147,-12.641,15.286;12.71,-12.136,12.79;18.012,-12.649,11.567;19.478,-12.385,11.172;20.187,-13.687,10.696;19.563,-14.252,9.404;18.599,-15.004,9.477;20.098,-13.86,8.229;19.648,-14.29,6.885;19.589,-15.803,6.809;18.634,-16.353,6.048;20.625,-13.781,5.83;20.759,-12.365,5.889;10.588,-11.836,3.288;11.37,-12.826,2.444;12.677,-13.19,3.171;13.48,-14.365,2.585;14.129,-13.976,1.342;15.412,-14.212,1.087;16.186,-14.727,1.907;7.421,-11.667,6.162;6.675,-12.99,6.581;6.654,-13.078,8.098;7.823,-13.191,8.731;7.885,-13.096,10.19;8.84,-11.926,10.495;9.77,-11.686,9.72;8.317,-14.401,10.877;7.872,-15.676,10.12;9.729,-14.357,11.106;13.973,-20.975,6.768;14.488,-21.058,8.204;13.543,-20.563,9.217;12.621,-21.314,9.805;11.859,-20.81,10.766;12.435,-22.565,9.423;14.659,-3.563,5.5;14.605,-4.881,5.293;13.35,-5.606,5.104;12.661,-5.673,6.457;11.536,-6.4,6.552;13.138,-5.069,7.429;21.154,-5.393,6.387;20.372,-4.647,5.612;19.099,-5.165,5.119;17.966,-4.53,5.886;17.138,-3.732,5.2;12.496,-8.081,8.665;12.367,-16.371,3.412;9.009,-17.463,14.046;9.374,-18.5,14.785;8.763,-18.749,15.936;10.373,-19.272,14.373;20.368,-1.798,-0.556;20.957,-1.177,0.463;20.595,-1.63,1.881;19.074,-1.555,2.174;18.319,-2.781,1.667;18.177,-2.926,0.433;17.852,-3.574,2.5;21.396,-0.863,2.83;21.647,-1.266,4.08;21.129,-2.256,4.588;26.022,-4.63,2.497;26.497,-10.537,-0.853;16.554,-20.899,12.518;16.264,-20.015,11.532;15.125,-19.344,11.93;14.762,-19.853,13.111;15.582,-20.814,13.498;27.715,-4.512,-2.12;27.605,-5.904,-2.016;28.301,-6.74,-2.887;29.098,-6.158,-3.87;28.228,-8.243,-2.741),wD:25.17,35.26,65.51,68.55,96.75,wU:40.30,49.38|")
    excVol = rdGaussianShape.ShapeInput(excVolMol, -1, rdGaussianShape.ShapeInputOptions(),
                                        rdGaussianShape.ShapeOverlayOptions())
      
    ssparams = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    ssparams.excludedVolume = excVol
    ssparams.fragSimilarityAdjuster = 0.2
    ssparams.approxSimilarityAdjuster = 0.2
    ssparams.numConformers = 100
    ssparams.confRMSThreshold = 0.5
    ssparams.randomSeed = 0xdac
    ssparams.bestHit = True
    ssparams.similarityCutoff = 0.5
    ssparams.maxMeanExcludedVolume = 5.0

    ovlyOptions = rdGaussianShape.ShapeOverlayOptions()
    ovlyOptions.simAlpha = 0.95
    ovlyOptions.simBeta = 0.05
    ssparams.shapeOverlayOptions = ovlyOptions
    
    fName = self.sssDir / "4ala_shapes.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)

    comb_4aji_4aj1 = Chem.MolFromSmiles("CNc1nc2ccc(NC(C)=O)cc2s1.COc1ccc(CC(C(=O)O)C(=O)O)cc1OC |(23.956,-7.951,-4.055;24.081,-7.139,-2.841;23.068,-6.904,-1.967;23.237,-6.265,-0.838;22.03,-6.164,-0.135;21.811,-5.548,1.089;20.55,-5.531,1.66;19.465,-6.131,1.009;18.155,-6.121,1.583;17.318,-7.197,1.736;16.262,-7.032,2.809;17.43,-8.221,1.058;19.663,-6.743,-0.207;20.934,-6.758,-0.767;21.422,-7.467,-2.285;16.331,-12.235,6.625;15.378,-13.287,6.54;14.89,-13.803,7.711;15.669,-14.364,8.712;15.068,-14.88,9.87;13.688,-14.822,10.03;13.008,-15.402,11.226;12.32,-16.745,10.932;13.341,-17.807,10.61;14.399,-17.852,11.169;12.961,-18.68,9.696;11.481,-17.202,12.095;11.812,-18.168,12.768;10.397,-16.444,12.327;12.917,-14.275,9.014;13.494,-13.757,7.872;12.793,-13.19,6.852;11.466,-12.683,7.147)|")

    hits = synthonspace.ShapeSearch(comb_4aji_4aj1, ssparams)
    self.assertEqual(len(hits.GetHitMolecules()), 1)
    

  def testPossibleHitsWrite(self):
    fName = self.sssDir / "amide_space_shapes.spc"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadDBFile(fName)


    ssparams = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    ssparams.fragSimilarityAdjuster = 0.2
    ssparams.approxSimilarityAdjuster = 0.2
    ssparams.numConformers = 100
    ssparams.confRMSThreshold = 0.5
    ssparams.randomSeed = 0xdac
    ssparams.bestHit = True
    ssparams.similarityCutoff = 0.8
    ssparams.maxMeanExcludedVolume = 5.0
    ssparams.possibleHitsFile = "amide_space_shapes_poss_hits.txt"
    ssparams.writePossibleHitsAndStop = True

    query = Chem.MolFromSmiles("O=C(c1ccccc1)N1CCCC1 |(0.0443291,-1.81486,-1.76886;0.0506321,-0.858174,-0.921491;1.37975,-0.430412,-0.483603;2.18964,-1.35506,0.144714;3.47088,-1.00454,0.585539;3.93803,0.297573,0.388032;3.1267,1.22739,-0.242406;1.85597,0.849751,-0.670531;-1.14837,-0.261434,-0.446583;-1.26073,0.836916,0.520219;-2.73583,1.04666,0.696614;-3.34033,-0.283345,0.290893;-2.46516,-0.679843,-0.874401)|")
    hits = synthonspace.ShapeSearch(query, ssparams)
    self.assertEqual(len(hits.GetHitMolecules()), 0)
    phf = Path(ssparams.possibleHitsFile)
    self.assertTrue(phf.exists())

    hits = synthonspace.ShapeSearch(query, ssparams, 0, -1)
    self.assertEqual(len(hits.GetHitMolecules()), 3)
    phf.unlink()
    

if __name__ == "__main__":
  unittest.main()
