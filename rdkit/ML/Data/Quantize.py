# $Id$
#
#  Copyright (C) 2001-2008  Greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
""" Automatic search for quantization bounds

This uses the expected informational gain to determine where quantization bounds should
lie.

**Notes**:

  - bounds are less than, so if the bounds are [1.,2.],
    [0.9,1.,1.1,2.,2.2] -> [0,1,1,2,2]

"""

import numpy

from rdkit.ML.InfoTheory import entropy

try:
  from rdkit.ML.Data import cQuantize
except ImportError:
  hascQuantize = 0
else:
  hascQuantize = 1

_float_tol = 1e-8


def feq(v1, v2, tol=_float_tol):
  """ floating point equality with a tolerance factor

      **Arguments**

        - v1: a float

        - v2: a float

        - tol: the tolerance for comparison

      **Returns**

        0 or 1

    """
  return abs(v1 - v2) < tol


def FindVarQuantBound(vals, results, nPossibleRes):
  """ Uses FindVarMultQuantBounds, only here for historic reasons
    """
  bounds, gain = FindVarMultQuantBounds(vals, 1, results, nPossibleRes)
  return (bounds[0], gain)


def _GenVarTable(vals, cuts, starts, results, nPossibleRes):
  """ Primarily intended for internal use

     constructs a variable table for the data passed in
     The table for a given variable records the number of times each possible value
      of that variable appears for each possible result of the function.

     **Arguments**

       - vals: a 1D Numeric array with the values of the variables

       - cuts: a list with the indices of the quantization bounds
         (indices are into _starts_ )

       - starts: a list of potential starting points for quantization bounds

       - results: a 1D Numeric array of integer result codes

       - nPossibleRes: an integer with the number of possible result codes

     **Returns**

       the varTable, a 2D Numeric array which is nVarValues x nPossibleRes

     **Notes**

       - _vals_ should be sorted!

    """
  nVals = len(cuts) + 1
  varTable = numpy.zeros((nVals, nPossibleRes), 'i')
  idx = 0
  for i in range(nVals - 1):
    cut = cuts[i]
    while idx < starts[cut]:
      varTable[i, results[idx]] += 1
      idx += 1
  while idx < len(vals):
    varTable[-1, results[idx]] += 1
    idx += 1
  return varTable


def _PyRecurseOnBounds(vals, cuts, which, starts, results, nPossibleRes, varTable=None):
  """ Primarily intended for internal use

     Recursively finds the best quantization boundaries

     **Arguments**

       - vals: a 1D Numeric array with the values of the variables,
         this should be sorted

       - cuts: a list with the indices of the quantization bounds
         (indices are into _starts_ )

       - which: an integer indicating which bound is being adjusted here
         (and index into _cuts_ )

       - starts: a list of potential starting points for quantization bounds

       - results: a 1D Numeric array of integer result codes

       - nPossibleRes: an integer with the number of possible result codes

     **Returns**

       - a 2-tuple containing:

         1) the best information gain found so far

         2) a list of the quantization bound indices ( _cuts_ for the best case)

     **Notes**

      - this is not even remotely efficient, which is why a C replacement
        was written

    """
  nBounds = len(cuts)
  maxGain = -1e6
  bestCuts = None
  highestCutHere = len(starts) - nBounds + which
  if varTable is None:
    varTable = _GenVarTable(vals, cuts, starts, results, nPossibleRes)
  while cuts[which] <= highestCutHere:
    varTable = _GenVarTable(vals, cuts, starts, results, nPossibleRes)
    gainHere = entropy.InfoGain(varTable)
    if gainHere > maxGain:
      maxGain = gainHere
      bestCuts = cuts[:]
    # recurse on the next vars if needed
    if which < nBounds - 1:
      gainHere, cutsHere = _RecurseOnBounds(vals, cuts[:], which + 1, starts, results, nPossibleRes,
                                            varTable=varTable)
      if gainHere > maxGain:
        maxGain = gainHere
        bestCuts = cutsHere
    # update this cut
    cuts[which] += 1
    for i in range(which + 1, nBounds):
      if cuts[i] == cuts[i - 1]:
        cuts[i] += 1

  return maxGain, bestCuts


def _NewPyRecurseOnBounds(vals, cuts, which, starts, results, nPossibleRes, varTable=None):
  """ Primarily intended for internal use

     Recursively finds the best quantization boundaries

     **Arguments**

       - vals: a 1D Numeric array with the values of the variables,
         this should be sorted

       - cuts: a list with the indices of the quantization bounds
         (indices are into _starts_ )

       - which: an integer indicating which bound is being adjusted here
         (and index into _cuts_ )

       - starts: a list of potential starting points for quantization bounds

       - results: a 1D Numeric array of integer result codes

       - nPossibleRes: an integer with the number of possible result codes

     **Returns**

       - a 2-tuple containing:

         1) the best information gain found so far

         2) a list of the quantization bound indices ( _cuts_ for the best case)

     **Notes**

      - this is not even remotely efficient, which is why a C replacement
        was written

    """
  nBounds = len(cuts)
  maxGain = -1e6
  bestCuts = None
  highestCutHere = len(starts) - nBounds + which
  if varTable is None:
    varTable = _GenVarTable(vals, cuts, starts, results, nPossibleRes)
  while cuts[which] <= highestCutHere:
    gainHere = entropy.InfoGain(varTable)
    if gainHere > maxGain:
      maxGain = gainHere
      bestCuts = cuts[:]
    # recurse on the next vars if needed
    if which < nBounds - 1:
      gainHere, cutsHere = _RecurseOnBounds(vals, cuts[:], which + 1, starts, results, nPossibleRes,
                                            varTable=None)
      if gainHere > maxGain:
        maxGain = gainHere
        bestCuts = cutsHere
    # update this cut
    oldCut = cuts[which]
    cuts[which] += 1
    bot = starts[oldCut]
    if oldCut + 1 < len(starts):
      top = starts[oldCut + 1]
    else:
      top = starts[-1]
    for i in range(bot, top):
      v = results[i]
      varTable[which, v] += 1
      varTable[which + 1, v] -= 1
    for i in range(which + 1, nBounds):
      if cuts[i] == cuts[i - 1]:
        cuts[i] += 1

  return maxGain, bestCuts


def _NewPyFindStartPoints(sortVals, sortResults, nData):
  # --------------------------------
  #
  # find all possible dividing points
  #
  #  There are a couple requirements for a dividing point:
  #    1) the dependent variable (descriptor) must change across it,
  #    2) the result score must change across it
  #
  #  So, in the list [(0,0),(1,0),(1,1),(2,1)]:
  #    we should divide before (1,0) and (2,1)
  #
  # --------------------------------
  startNext = []
  tol = 1e-8
  blockAct = sortResults[0]
  lastBlockAct = None
  lastDiv = None
  i = 1
  while i < nData:
    # move to the end of this block:
    while i < nData and sortVals[i] - sortVals[i - 1] <= tol:
      if sortResults[i] != blockAct:
        # this block is heterogeneous
        blockAct = -1
      i += 1
    if lastBlockAct is None:
      # first time through:
      lastBlockAct = blockAct
      lastDiv = i
    else:
      if blockAct == -1 or lastBlockAct == -1 or blockAct != lastBlockAct:
        startNext.append(lastDiv)
        lastDiv = i
        lastBlockAct = blockAct
      else:
        lastDiv = i
    if i < nData:
      blockAct = sortResults[i]
    i += 1
  # catch the case that the last point also sets a bin:
  if blockAct != lastBlockAct:
    startNext.append(lastDiv)
  return startNext


def FindVarMultQuantBounds(vals, nBounds, results, nPossibleRes):
  """ finds multiple quantization bounds for a single variable

     **Arguments**

       - vals: sequence of variable values (assumed to be floats)

       - nBounds: the number of quantization bounds to find

       - results: a list of result codes (should be integers)

       - nPossibleRes: an integer with the number of possible values of the
         result variable

     **Returns**

       - a 2-tuple containing:

         1) a list of the quantization bounds (floats)

         2) the information gain associated with this quantization


    """
  assert len(vals) == len(results), 'vals/results length mismatch'

  nData = len(vals)
  if nData == 0:
    return [], -1e8

  # sort the variable values:
  svs = list(zip(vals, results))
  svs.sort()
  sortVals, sortResults = zip(*svs)
  startNext = _FindStartPoints(sortVals, sortResults, nData)
  if not len(startNext):
    return [0], 0.0
  if len(startNext) < nBounds:
    nBounds = len(startNext) - 1
  if nBounds == 0:
    nBounds = 1
  initCuts = list(range(nBounds))
  maxGain, bestCuts = _RecurseOnBounds(sortVals, initCuts, 0, startNext, sortResults, nPossibleRes)
  quantBounds = []
  nVs = len(sortVals)
  for cut in bestCuts:
    idx = startNext[cut]
    if idx == nVs:
      quantBounds.append(sortVals[-1])
    elif idx == 0:
      quantBounds.append(sortVals[idx])
    else:
      quantBounds.append((sortVals[idx] + sortVals[idx - 1]) / 2.)

  return quantBounds, maxGain


# hascQuantize=0
if hascQuantize:
  _RecurseOnBounds = cQuantize._RecurseOnBounds
  _FindStartPoints = cQuantize._FindStartPoints
else:
  _RecurseOnBounds = _NewPyRecurseOnBounds
  _FindStartPoints = _NewPyFindStartPoints

if __name__ == '__main__':
  if 1:
    d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 1), (1.4, 0), (1.6, 1), (2., 1), (2.1, 0), (2.1, 0),
         (2.1, 0), (2.2, 1), (2.3, 0)]
    varValues = list(map(lambda x: x[0], d))
    resCodes = list(map(lambda x: x[1], d))
    nPossibleRes = 2
    res = FindVarMultQuantBounds(varValues, 2, resCodes, nPossibleRes)
    print('RES:', res)
    target = ([1.3, 2.05], .34707)
  else:
    d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 1), (1.4, 0), (1.6, 1), (2., 1), (2.1, 0), (2.1, 0),
         (2.1, 0), (2.2, 1), (2.3, 0)]
    varValues = list(map(lambda x: x[0], d))
    resCodes = list(map(lambda x: x[1], d))
    nPossibleRes = 2
    res = FindVarMultQuantBounds(varValues, 1, resCodes, nPossibleRes)
    print(res)
    # sys.exit(1)
    d = [(1.4, 1), (1.4, 0)]

    varValues = list(map(lambda x: x[0], d))
    resCodes = list(map(lambda x: x[1], d))
    nPossibleRes = 2
    res = FindVarMultQuantBounds(varValues, 1, resCodes, nPossibleRes)
    print(res)

    d = [(1.4, 0), (1.4, 0), (1.6, 1)]
    varValues = list(map(lambda x: x[0], d))
    resCodes = list(map(lambda x: x[1], d))
    nPossibleRes = 2
    res = FindVarMultQuantBounds(varValues, 2, resCodes, nPossibleRes)
    print(res)
