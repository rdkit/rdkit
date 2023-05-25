r"""
$Id$

Scoring - Calculate rank statistics

Created by Sereina Riniker, October 2012
after a file from Peter Gedeck, Greg Landrum

\param scores: ordered list with descending similarity containing
               active/inactive information
\param col: column index in scores where active/inactive information is stored
\param fractions: list of fractions at which the value shall be calculated
\param alpha: exponential weight
"""

import math
from collections import namedtuple


def CalcROC(scores, col):
  """ Determines a ROC curve """
  numMol = len(scores)
  if numMol == 0:
    raise ValueError('score list is empty')
  TPR = [0] * numMol  # True positive rate: TP/(TP+FN)
  FPR = [0] * numMol  # False positive rate: FP/(TN+FP)
  numActives = 0
  numInactives = 0

  # loop over score list
  for i in range(numMol):
    if scores[i][col]:
      numActives += 1
    else:
      numInactives += 1
    TPR[i] = numActives  # TP
    FPR[i] = numInactives  # FP

  # normalize, check that there are actives and inactives
  if numActives > 0:
    TPR = [1.0 * i / numActives for i in TPR]
  if numInactives > 0:
    FPR = [1.0 * i / numInactives for i in FPR]

  RocCurve = namedtuple('RocCurve', ['FPR', 'TPR'])
  return RocCurve(FPR=FPR, TPR=TPR)


def CalcAUC(scores, col):
  """ Determines the area under the ROC curve """
  # determine the ROC curve
  roc = CalcROC(scores, col)
  FPR = roc.FPR
  TPR = roc.TPR

  numMol = len(scores)
  AUC = 0

  # loop over score list
  for i in range(0, numMol - 1):
    AUC += (FPR[i + 1] - FPR[i]) * (TPR[i + 1] + TPR[i])

  return 0.5 * AUC


def _RIEHelper(scores, col, alpha):
  numMol = len(scores)
  alpha = float(alpha)
  if numMol == 0:
    raise ValueError('score list is empty')
  if alpha <= 0.0:
    raise ValueError('alpha must be greater than zero')

  denom = 1.0 / numMol * ((1 - math.exp(-alpha)) / (math.exp(alpha / numMol) - 1))
  numActives = 0
  sum_exp = 0

  # loop over score list
  for i in range(numMol):
    active = scores[i][col]
    if active:
      numActives += 1
      sum_exp += math.exp(-(alpha * (i + 1)) / numMol)

  if numActives > 0:  # check that there are actives
    RIE = sum_exp / (numActives * denom)
  else:
    RIE = 0.0

  return RIE, numActives


def CalcRIE(scores, col, alpha):
  """ RIE original definded here:
    Sheridan, R.P., Singh, S.B., Fluder, E.M. & Kearsley, S.K.
    Protocols for Bridging the Peptide to Nonpeptide Gap in Topological Similarity Searches.
    J. Chem. Inf. Comp. Sci. 41, 1395-1406 (2001).
    """
  RIE, _ = _RIEHelper(scores, col, alpha)
  return RIE


def CalcBEDROC(scores, col, alpha):
  """ BEDROC original defined here:
    Truchon, J. & Bayly, C.I.
    Evaluating Virtual Screening Methods: Good and Bad Metric for the "Early Recognition"
    Problem. J. Chem. Inf. Model. 47, 488-508 (2007).
    ** Arguments**

      - scores: 2d list or numpy array
             0th index representing sample
             scores must be in sorted order with low indexes "better"
             scores[sample_id] = vector of sample data
      -  col: int
             Index of sample data which reflects true label of a sample
             scores[sample_id][col] = True iff that sample is active
      -  alpha: float
             hyper parameter from the initial paper for how much to enrich the top
     **Returns**
       float BedROC score
    """
  # calculate RIE
  RIE, numActives = _RIEHelper(scores, col, alpha)

  if numActives > 0:
    numMol = len(scores)
    ratio = 1.0 * numActives / numMol
    RIEmax = (1 - math.exp(-alpha * ratio)) / (ratio * (1 - math.exp(-alpha)))
    RIEmin = (1 - math.exp(alpha * ratio)) / (ratio * (1 - math.exp(alpha)))

    if RIEmax != RIEmin:
      BEDROC = (RIE - RIEmin) / (RIEmax - RIEmin)
    else:  # numActives = numMol
      BEDROC = 1.0
  else:
    BEDROC = 0.0

  return BEDROC


def CalcEnrichment(scores, col, fractions):
  """ Determines the enrichment factor for a set of fractions """
  numMol = len(scores)
  if numMol == 0:
    raise ValueError('score list is empty')
  if len(fractions) == 0:
    raise ValueError('fraction list is empty')
  for i in fractions:
    if i > 1 or i < 0:
      raise ValueError('fractions must be between [0,1]')

  numPerFrac = [math.ceil(numMol * f) for f in fractions]
  numPerFrac.append(numMol)
  numActives = 0
  enrich = []

  # loop over score list
  for i in range(numMol):
    if i > (numPerFrac[0] - 1) and i > 0:
      enrich.append(1.0 * numActives * numMol / i)
      numPerFrac.pop(0)
    active = scores[i][col]
    if active:
      numActives += 1

  if numActives > 0:  # check that there are actives
    enrich = [e / numActives for e in enrich]
  else:
    enrich = [0.0] * len(fractions)
  return enrich


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
