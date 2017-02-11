"""
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
from __future__ import print_function
import math
import numpy as np
from operator import itemgetter


def SortByColumns(scores, classCol, reverse=False):
  """
  Safeguard against wrong ROCs (& AUCs & RIEs and whatnot) due to ordering artifacts.

  Like `sorted` but preshuffles and explicitly ignores class attributes.

  This avoids having artificially good/bad scores because:
    - `scores` are ordered by class (col) and python stable (tim)sort
      keeps the order on ties; this is specially bad if our model
      returns (nan) scores, or any other constant.
    - related, we leak the class in sorting; it is easy to introduce
      this kind of bugs as one would tend to put together scores and
      classes when ranking.

  Parameters
  ----------
  scores : list, tuple, np.array or the likes
    A 2-dimensional matrix of scores.
    num_rows=len(scores), num_columns=len(scores[0])

  classCol : int or int collection
    The columns to ignore.

  reverse : bool, default False
    If True, scores are sorted in ascending order
    (e.g. the lower a distance is, the higher is the score).
    If False, use descending.
    (e.g. the higher a probability estimate is, the higher the score).

  Examples
  --------
    This is an extreme example with:
      - every example tied
      - the data pre-sorted by class
    >>> scores = [
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 0],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ...     [0.2, 0.3, 1],
    ... ]

    The correct value of the AUC for these scores is 0.5.
    (in fact, the scores are completely uninformative).
    Using `SortByColumns` make it explicit how correct sorting
    must be done to get the correct result.
    >>> shuffled_sorted_scores = SortByColumns(scores, -1)
    >>> print('Scores are uninformative, so AUC must be',
    ...       CalcAUC(shuffled_sorted_scores, -1))
    Scores are uninformative, so AUC must be 0.5

    If we simply sort using scores, data preordering will make the
    AUC estimation wrong.
    >>> sorted_scores = sorted(scores, key=lambda x: (x[0], x[1]))
    >>> print('If we let initial ordering to leak, AUC will wrongly be',
    ...       CalcAUC(sorted_scores, -1))
    If we let initial ordering to leak, AUC will wrongly be 0.0

    If we forget to remove the class when sorting, class membership leaking
    will render the AUC estimation wrong.
    >>> class_leaked_scores = sorted(scores, reverse=True)
    >>> print('If we let the class to leak when sorting, AUC will wrongly be',
    ...       CalcAUC(class_leaked_scores, -1))
    If we let the class to leak when sorting, AUC will wrongly be 1.0

  See also
  --------
    Github PR: https://github.com/rdkit/rdkit/pull/1306
    Discusses some performance penalties incurred when using this function
    as currently implemented
  """
  # Shuffle
  scores = list(scores)                     # No side effects
  if 0 == len(scores):
    raise ValueError('scores cannot be empty')
  np.random.RandomState(0).shuffle(scores)  # Shuffle (deterministically)

  # Allow to pass collections of columns to ignore
  if isinstance(classCol, int):
    classCol = (classCol,)

  # Sort without taking the class into account
  # This works for any of {list, tuple, numpy array, ...}
  numCols = len(scores[0])
  classCol = set(c if c >= 0 else c + numCols for c in classCol)
  allButClass = itemgetter(*[i for i in range(numCols) if i not in classCol])
  return sorted(scores, key=allButClass, reverse=reverse)

  # Could these safeguards be introduced in ROC/BEDROC etc.
  # without changing the API but for adding keyword arguments?


def CalcROC(scores, col):
  """ Determines a ROC curve """
  numMol = len(scores)
  if numMol == 0:
    raise ValueError('score list is empty')
  TPR = [0] * numMol  # True positive rate: TP/(TP+FP)
  TNR = [0] * numMol  # True negative rate: TN/(TN+FN)
  numActives = 0
  numInactives = 0

  # loop over score list
  for i in range(numMol):
    if scores[i][col]:
      numActives += 1
    else:
      numInactives += 1
    TPR[i] = numActives  # TP
    TNR[i] = numInactives  # TN

  # normalize, check that there are actives and inactives
  if numActives > 0:
    TPR = [1.0 * i / numActives for i in TPR]
  if numInactives > 0:
    TNR = [1.0 * i / numInactives for i in TNR]

  return [TNR, TPR]


def CalcAUC(scores, col):
  """ Determines the area under the ROC curve """
  # determine the ROC curve
  roc = CalcROC(scores, col)
  TNR = roc[0]
  TPR = roc[1]

  numMol = len(scores)
  AUC = 0

  # loop over score list
  for i in range(0, numMol - 1):
    AUC += (TNR[i + 1] - TNR[i]) * (TPR[i + 1] + TPR[i])

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
