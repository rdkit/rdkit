#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  Copyright (c) 2021, Greg Landrum
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

import sys
import unittest

from rdkit import Chem
from rdkit.RDLogger import logger

try:
  import matplotlib
except ImportError:
  matplotlib = None

from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import SimilarityMaps as sm

logger = logger()


class TestCase(unittest.TestCase):

  def setUp(self):
    self.mol1 = Chem.MolFromSmiles('c1ccccc1')
    self.mol2 = Chem.MolFromSmiles('c1ccncc1')

  @unittest.skipUnless(matplotlib, 'Matplotlib required')
  def testSimilarityMap(self):
    # Morgan2 BV
    refWeights = [0.5, 0.5, 0.5, -0.5, 0.5, 0.5]
    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'))
    for w, r in zip(weights, refWeights):
      self.assertEqual(w, r)

    d2d = rdMolDraw2D.MolDraw2DSVG(400, 400)
    _, maxWeight = sm.GetSimilarityMapForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'), d2d)
    self.assertEqual(maxWeight, 0.5)

    weights, maxWeight = sm.GetStandardizedWeights(weights)
    self.assertEqual(maxWeight, 0.5)
    refWeights = [1.0, 1.0, 1.0, -1.0, 1.0, 1.0]
    for w, r in zip(weights, refWeights):
      self.assertEqual(w, r)

    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetMorganFingerprint(m, i, fpType='count'))
    self.assertTrue(weights[3] < 0)
    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2,
      lambda m, i: sm.GetMorganFingerprint(m, i, fpType='bv', useFeatures=True))
    self.assertTrue(weights[3] < 0)

    # hashed AP BV
    refWeights = [0.09523, 0.17366, 0.17366, -0.23809, 0.17366, 0.17366]
    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetAPFingerprint(m, i, fpType='bv', nBits=1024))
    for w, r in zip(weights, refWeights):
      self.assertAlmostEqual(w, r, 4)

    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetAPFingerprint(m, i, fpType='normal'))
    self.assertTrue(weights[3] < 0)
    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetAPFingerprint(m, i, fpType='hashed'))
    self.assertTrue(weights[3] < 0)

    # hashed TT BV
    refWeights = [0.5, 0.5, -0.16666, -0.5, -0.16666, 0.5]
    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2,
      lambda m, i: sm.GetTTFingerprint(m, i, fpType='bv', nBits=1024, nBitsPerEntry=1))
    for w, r in zip(weights, refWeights):
      self.assertAlmostEqual(w, r, 4)

    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetTTFingerprint(m, i, fpType='normal'))
    self.assertTrue(weights[3] < 0)
    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetTTFingerprint(m, i, fpType='hashed'))
    self.assertTrue(weights[3] < 0)

    # RDK fingerprint BV
    refWeights = [0.42105, 0.42105, 0.42105, -0.32895, 0.42105, 0.42105]
    weights = sm.GetAtomicWeightsForFingerprint(
      self.mol1, self.mol2, lambda m, i: sm.GetRDKFingerprint(m, i, nBits=1024, nBitsPerHash=1))
    for w, r in zip(weights, refWeights):
      self.assertAlmostEqual(w, r, 4)

  @unittest.skipUnless(matplotlib, 'Matplotlib required')
  def testSimilarityMapKWArgs(self):
    # Morgan2 BV
    m1 = Chem.MolFromSmiles('CC[C@](F)(Cl)c1ccccc1')
    m2 = Chem.MolFromSmiles('CC[C@@](F)(Cl)c1ccccc1')
    weights = sm.GetAtomicWeightsForFingerprint(
      m1, m2, lambda m, i: sm.GetAPFingerprint(m, atomId=i, includeChirality=False))
    for w in weights:
      self.assertAlmostEqual(w, 0.100, 4)
    weights = sm.GetAtomicWeightsForFingerprint(
      m1, m2, lambda m, i: sm.GetAPFingerprint(m, atomId=i, includeChirality=True))
    for i, w in enumerate(weights):
      if i != 2:
        self.assertAlmostEqual(w, 0.098, 3)
      else:
        self.assertAlmostEqual(w, -0.082, 3)

    weights = sm.GetAtomicWeightsForFingerprint(
      m1, m2, lambda m, i: sm.GetTTFingerprint(m, atomId=i, includeChirality=False))
    for w in weights:
      self.assertTrue(w > 0.0)
    weights = sm.GetAtomicWeightsForFingerprint(
      m1, m2, lambda m, i: sm.GetTTFingerprint(m, atomId=i, includeChirality=True))
    for i, w in enumerate(weights):
      if i > 4:
        self.assertTrue(w > 0.0)
      else:
        self.assertTrue(w < 0.0)

    weights = sm.GetAtomicWeightsForFingerprint(
      m1, m2, lambda m, i: sm.GetMorganFingerprint(m, radius=1, atomId=i, useChirality=False))
    weights2 = sm.GetAtomicWeightsForFingerprint(
      m1, m2, lambda m, i: sm.GetMorganFingerprint(m, radius=1, atomId=i, useChirality=True))
    # testing explicit values here seems silly, just check that the contribution of the
    # chiral center drops:
    self.assertTrue(weights[2] > weights2[2])

  def testSimilarityMapsMolDraw2D(self):
    # nothing really sensible to test here, just make sure things run
    mol = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')
    refmol = Chem.MolFromSmiles('CCCN(CCCCN1CCN(c2ccccc2OC)CC1)Cc1ccc2ccccc2c1')
    d = Draw.MolDraw2DSVG(400, 400)
    d.ClearDrawing()
    _, maxWeight = sm.GetSimilarityMapForFingerprint(
      refmol, mol, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'), draw2d=d)
    d.FinishDrawing()
    with open('similarityMap1_out.svg', 'w+') as outf:
      outf.write(d.GetDrawingText())

    # Github #2904: make sure we can provide our own colormap as a list:
    colors = [(0, 1, 0, 0.5), (1, 1, 1), (0, 0, 1, 0.5)]
    d = Draw.MolDraw2DSVG(400, 400)
    d.ClearDrawing()
    _, maxWeight = sm.GetSimilarityMapForFingerprint(
      refmol, mol, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'), draw2d=d,
      colorMap=colors)
    d.FinishDrawing()
    with open('similarityMap1_out2.svg', 'w+') as outf:
      outf.write(d.GetDrawingText())

    # Github #2904: make sure we can provide our own colormap as a matplotlib colormap:
    try:
      from matplotlib import cm
      d = Draw.MolDraw2DSVG(400, 400)
      d.ClearDrawing()
      _, maxWeight = sm.GetSimilarityMapForFingerprint(
        refmol, mol, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'), draw2d=d,
        colorMap=cm.PiYG)
      d.FinishDrawing()
      with open('similarityMap1_out3.svg', 'w+') as outf:
        outf.write(d.GetDrawingText())
    except ImportError:
      pass

  @unittest.skipUnless(matplotlib, 'Matplotlib required')
  def testGithub4763(self):
    mol = Chem.MolFromSmiles('COc1cccc2cc(C(=O)NCCCCN3CCN(c4cccc5nccnc54)CC3)oc21')
    refmol = Chem.MolFromSmiles('CCCN(CCCCN1CCN(c2ccccc2OC)CC1)Cc1ccc2ccccc2c1')
    d = Draw.MolDraw2DSVG(400, 400)
    d.ClearDrawing()
    _, maxWeight = sm.GetSimilarityMapForFingerprint(
      refmol, mol, lambda m, i: sm.GetMorganFingerprint(m, i, radius=2, fpType='bv'), draw2d=d,
      colorMap="coolwarm")
    d.FinishDrawing()
    svg = d.GetDrawingText()
    with open('github4763.svg', 'w+') as outf:
      outf.write(svg)
    self.assertFalse('fill:#FBFCFB7F' in svg)
    self.assertTrue('fill:#DDDCDB' in svg)


if __name__ == '__main__':
  unittest.main()
