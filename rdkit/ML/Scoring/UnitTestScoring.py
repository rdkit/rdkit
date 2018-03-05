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
""" unit testing code for the Scoring functionality

"""

import math
import unittest

from rdkit.ML.Scoring import Scoring


class TestCase(unittest.TestCase):

  def setUp(self):
    # generate a scored list of 400 molecules with 20 actives
    self.numActives = 20
    self.numDecoys = 380
    self.numMol = self.numActives + self.numDecoys
    act = [[1] for _ in range(0, self.numActives)]
    dcy = [[0] for _ in range(0, self.numDecoys)]
    self.scoreBestCase = act + dcy
    self.scoreWorstCase = dcy + act
    self.scoreEmptyList = []
    self.scoreAllActives = [[1] for _ in range(0, self.numMol)]
    self.scoreAllDecoys = [[0] for _ in range(0, self.numMol)]
    self.index = 0  # where the active/inactive information lies
    # test the 5% fraction
    self.fractions = [float(self.numActives) / self.numMol]
    self.fracSmall = [self.fractions[0] / 100]
    # exponential weight
    self.alpha = 100.0
    # accuracy for float number comparison
    self.acc = 4

  def test1_enrichment_factor(self):
    # """ test enrichment factor """
    # best case
    enrich = Scoring.CalcEnrichment(self.scoreBestCase, self.index, self.fractions)
    self.assertAlmostEqual(enrich[0], float(self.numActives), self.acc)
    # worst case
    enrich = Scoring.CalcEnrichment(self.scoreWorstCase, self.index, self.fractions)
    self.assertAlmostEqual(enrich[0], 0.0, self.acc)
    # empty list
    self.assertRaises(ValueError, Scoring.CalcEnrichment, self.scoreEmptyList, self.index,
                      self.fractions)
    # all actives
    enrich = Scoring.CalcEnrichment(self.scoreAllActives, self.index, self.fractions)
    self.assertAlmostEqual(enrich[0], 1.0, self.acc)
    # all decoys
    enrich = Scoring.CalcEnrichment(self.scoreAllDecoys, self.index, self.fractions)
    self.assertEqual(enrich[0], 0.0)
    # fraction * numMol is smaller than 1
    enrich = Scoring.CalcEnrichment(self.scoreBestCase, self.index, self.fracSmall)
    self.assertAlmostEqual(enrich[0], float(self.numActives), self.acc)
    # fraction list is empty
    self.assertRaises(ValueError, Scoring.CalcEnrichment, self.scoreBestCase, self.index, [])
    # fraction == 0.0
    enrich = Scoring.CalcEnrichment(self.scoreBestCase, self.index, [0.0])
    self.assertAlmostEqual(enrich[0], float(self.numActives), self.acc)
    # fraction < 0
    self.assertRaises(ValueError, Scoring.CalcEnrichment, self.scoreBestCase, self.index, [-0.05])
    # fraction > 1
    self.assertRaises(ValueError, Scoring.CalcEnrichment, self.scoreBestCase, self.index, [1.5])

  def test2_RIE(self):
    # """ test RIE """
    ratio = float(self.numActives) / self.numMol
    # best case
    RIEmax = ((1 - math.exp(-self.alpha * ratio)) / (1 - math.exp(-self.alpha))) / ratio
    rie = Scoring.CalcRIE(self.scoreBestCase, self.index, self.alpha)
    self.assertAlmostEqual(rie, RIEmax, self.acc)
    # worst case
    RIEmin = ((1 - math.exp(self.alpha * ratio)) / (1 - math.exp(self.alpha))) / ratio
    rie = Scoring.CalcRIE(self.scoreWorstCase, self.index, self.alpha)
    self.assertAlmostEqual(rie, RIEmin, self.acc)
    # empty list
    self.assertRaises(ValueError, Scoring.CalcRIE, self.scoreEmptyList, self.index, self.alpha)
    # alpha == 0
    self.assertRaises(ValueError, Scoring.CalcRIE, self.scoreBestCase, self.index, 0.0)
    # all decoys
    rie = Scoring.CalcRIE(self.scoreAllDecoys, self.index, self.alpha)
    self.assertEqual(rie, 0.0)

  def test3_AUC(self):
    # """ test area under the curve (AUC) of ROC """
    # best case
    auc = Scoring.CalcAUC(self.scoreBestCase, self.index)
    self.assertAlmostEqual(auc, 1.0, self.acc)
    # worst case
    auc = Scoring.CalcAUC(self.scoreWorstCase, self.index)
    self.assertAlmostEqual(auc, 0.0, self.acc)
    # empty list
    self.assertRaises(ValueError, Scoring.CalcAUC, self.scoreEmptyList, self.index)
    # all actives
    auc = Scoring.CalcAUC(self.scoreAllActives, self.index)
    self.assertAlmostEqual(auc, 0.0, self.acc)
    # all decoys
    auc = Scoring.CalcAUC(self.scoreAllDecoys, self.index)
    self.assertAlmostEqual(auc, 0.0, self.acc)

  def test4_BEDROC(self):
    # """ test BEDROC """
    # best case
    bedroc = Scoring.CalcBEDROC(self.scoreBestCase, self.index, self.alpha)
    self.assertAlmostEqual(bedroc, 1.0, self.acc)
    # worst case
    bedroc = Scoring.CalcBEDROC(self.scoreWorstCase, self.index, self.alpha)
    self.assertAlmostEqual(bedroc, 0.0, self.acc)
    # empty list
    self.assertRaises(ValueError, Scoring.CalcBEDROC, self.scoreEmptyList, self.index, self.alpha)
    # alpha == 0.0
    self.assertRaises(ValueError, Scoring.CalcBEDROC, self.scoreBestCase, self.index, 0.0)
    # all actives
    bedroc = Scoring.CalcBEDROC(self.scoreAllActives, self.index, self.alpha)
    self.assertEqual(bedroc, 1.0)
    # all decoys
    bedroc = Scoring.CalcBEDROC(self.scoreAllDecoys, self.index, self.alpha)
    self.assertEqual(bedroc, 0.0)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
