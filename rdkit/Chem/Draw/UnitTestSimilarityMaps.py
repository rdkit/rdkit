# $Id$
#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
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
# Created by Sereina Riniker, Aug 2013


""" unit testing code for molecule drawing
"""
from rdkit import RDConfig
import unittest,os,tempfile
from rdkit import Chem
from rdkit.Chem import Draw
try:
  from rdkit.Chem.Draw import SimilarityMaps as sm
except ImportError:
  sm = None
from rdkit.RDLogger import logger
logger = logger()

class TestCase(unittest.TestCase):
  def setUp(self):
    self.mol1 = Chem.MolFromSmiles('c1ccccc1')
    self.mol2 = Chem.MolFromSmiles('c1ccncc1')

  def testSimilarityMap(self):
    # Morgan2 BV
    refWeights = [0.5, 0.5, 0.5, -0.5, 0.5, 0.5]
    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'))
    for w,r in zip(weights, refWeights): self.assertEqual(w, r)

    fig, maxWeight = sm.GetSimilarityMapForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'))
    self.assertEqual(maxWeight, 0.5)
    
    weights, maxWeight = sm.GetStandardizedWeights(weights)
    self.assertEqual(maxWeight, 0.5)
    refWeights = [1.0, 1.0, 1.0, -1.0, 1.0, 1.0]
    for w,r in zip(weights, refWeights): self.assertEqual(w, r)

    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetMorganFingerprint(m, i, fpType='count'))
    self.assertTrue(weights[3] < 0)
    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetMorganFingerprint(m, i, fpType='bv', useFeatures=True))
    self.assertTrue(weights[3] < 0)

    # hashed AP BV
    refWeights = [0.09523, 0.17366, 0.17366, -0.23809, 0.17366, 0.17366]
    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetAPFingerprint(m, i, fpType='bv', nBits=1024))
    for w,r in zip(weights, refWeights): self.assertAlmostEqual(w, r, 4)

    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetAPFingerprint(m, i, fpType='normal'))
    self.assertTrue(weights[3] < 0)
    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetAPFingerprint(m, i, fpType='hashed'))
    self.assertTrue(weights[3] < 0)
    
    # hashed TT BV
    refWeights = [0.5, 0.5, -0.16666, -0.5, -0.16666, 0.5]
    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetTTFingerprint(m, i, fpType='bv', nBits=1024, nBitsPerEntry=1))
    for w,r in zip(weights, refWeights): self.assertAlmostEqual(w, r, 4)

    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetTTFingerprint(m, i, fpType='normal'))
    self.assertTrue(weights[3] < 0)
    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetTTFingerprint(m, i, fpType='hashed'))
    self.assertTrue(weights[3] < 0)

    # RDK fingerprint BV
    refWeights = [0.42105, 0.42105, 0.42105, -0.32895, 0.42105, 0.42105]
    weights = sm.GetAtomicWeightsForFingerprint(self.mol1, self.mol2, lambda m, i: sm.GetRDKFingerprint(m, i, nBits=1024, nBitsPerHash=1))
    for w,r in zip(weights, refWeights): self.assertAlmostEqual(w, r, 4)

    
if __name__ == '__main__':
  try:
    import matplotlib
    from rdkit.Chem.Draw.mplCanvas import Canvas
  except ImportError:
    pass
  except RuntimeError:  # happens with GTK can't initialize
    pass
  else:
    unittest.main()
