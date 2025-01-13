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

from rdkit import Chem
from rdkit.Chem import rdSynthonSpaceSearch, rdFingerprintGenerator


class TestCase(unittest.TestCase):

  def setUp(self):
    self.sssDir = Path(os.environ["RDBASE"]) / "Code" / "GraphMol" / "SynthonSpaceSearch" / "data"

  def testSubstructSearch(self):
    fName = self.sssDir / "Syntons_5567.csv"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadTextFile(fName)
    self.assertEqual(10, synthonspace.GetNumReactions())
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = 10
    results = synthonspace.SubstructureSearch(Chem.MolFromSmarts("c1ccccc1C(=O)N1CCCC1"), params)
    self.assertEqual(10, len(results.GetHitMolecules()))

  def testFingerprintSearch(self):
    fName = self.sssDir / "Syntons_5567.csv"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadTextFile(fName)
    self.assertEqual(10, synthonspace.GetNumReactions())
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.maxHits = 10
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048, useBondOrder=True)
    results = synthonspace.FingerprintSearch(
      Chem.MolFromSmiles("c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"), fpgen, params)
    self.assertEqual(10, len(results.GetHitMolecules()))

  def testEnumerate(self):
    fName = self.sssDir / "amide_space.txt"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadTextFile(fName)
    with tempfile.NamedTemporaryFile() as tmp:
      synthonspace.WriteEnumeratedFile(tmp.name)

  def testTimeOut(self):
    fName = self.sssDir / "Syntons_5567.csv"
    synthonspace = rdSynthonSpaceSearch.SynthonSpace()
    synthonspace.ReadTextFile(fName)
    self.assertEqual(10, synthonspace.GetNumReactions())
    params = rdSynthonSpaceSearch.SynthonSpaceSearchParams()
    params.timeOut = 50
    params.similarityCutoff = 0.3
    params.fragSimilarityAdjuster = 0.3
    fpgen = rdFingerprintGenerator.GetRDKitFPGenerator(fpSize=2048, useBondOrder=True)
    results = synthonspace.FingerprintSearch(
      Chem.MolFromSmiles("c12ccc(C)cc1[nH]nc2C(=O)NCc1cncs1"), fpgen, params)
    self.assertFalse(results.GetTimedOut())

    params.timeOut = 1
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


if __name__ == "__main__":
  unittest.main()
