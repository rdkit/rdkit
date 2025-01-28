#
#  Copyright (C) 2001,2002,2003  greg Landrum and Rational Discovery LLC
#
""" Functionality for ranking bits using info gains

 **Definitions used in this module**

    - *sequence*: an object capable of containing other objects which supports
      __getitem__() and __len__().  Examples of these include lists, tuples, and
      Numeric arrays.

    - *IntVector*: an object containing integers which supports __getitem__() and
       __len__(). Examples include lists, tuples, Numeric Arrays, and BitVects.


 **NOTE**: Neither *sequences* nor *IntVectors* need to support item assignment.
   It is perfectly acceptable for them to be read-only, so long as they are
   random-access.

"""
import numpy

from rdkit.ML.InfoTheory import entropy


def FormCounts(bitVects, actVals, whichBit, nPossibleActs, nPossibleBitVals=2):
  """ generates the counts matrix for a particular bit

  **Arguments**

    - bitVects: a *sequence* containing *IntVectors*

    - actVals: a *sequence*

    - whichBit: an integer, the bit number to use.

    - nPossibleActs: the (integer) number of possible activity values.

    - nPossibleBitVals: (optional) if specified, this integer provides the maximum
      value attainable by the (increasingly inaccurately named) bits in _bitVects_.

  **Returns**

    a Numeric array with the counts

  **Notes**

    This is really intended for internal use.

  """
  if len(bitVects) != len(actVals):
    raise ValueError('var and activity lists should be the same length')
  res = numpy.zeros((nPossibleBitVals, nPossibleActs), numpy.integer)
  for i in range(len(bitVects)):
    res[bitVects[i][whichBit], actVals[i]] += 1
  return res


def CalcInfoGains(bitVects, actVals, nPossibleActs, nPossibleBitVals=2):
  """  Calculates the information gain for a set of points and activity values

  **Arguments**

    - bitVects: a *sequence* containing *IntVectors*

    - actVals: a *sequence*

    - nPossibleActs: the (integer) number of possible activity values.

    - nPossibleBitVals: (optional) if specified, this integer provides the maximum
      value attainable by the (increasingly inaccurately named) bits in _bitVects_.

   **Returns**

     a list of floats

  """
  if len(bitVects) != len(actVals):
    raise ValueError('var and activity lists should be the same length')
  nBits = len(bitVects[0])
  res = numpy.zeros(nBits, float)

  for bit in range(nBits):
    counts = FormCounts(bitVects, actVals, bit, nPossibleActs, nPossibleBitVals=nPossibleBitVals)
    res[bit] = entropy.InfoGain(counts)
  return res


def RankBits(bitVects, actVals, nPossibleBitVals=2, metricFunc=CalcInfoGains):
  """ Rank a set of bits according to a metric function

  **Arguments**

    - bitVects: a *sequence* containing *IntVectors*

    - actVals: a *sequence*

    - nPossibleBitVals: (optional) if specified, this integer provides the maximum
      value attainable by the (increasingly inaccurately named) bits in _bitVects_.

    - metricFunc: (optional) the metric function to be used.  See _CalcInfoGains()_
      for a description of the signature of this function.

   **Returns**

     A 2-tuple containing:

       - the relative order of the bits (a list of ints)

       - the metric calculated for each bit (a list of floats)

  """
  nPossibleActs = max(actVals) + 1
  metrics = metricFunc(bitVects, actVals, nPossibleActs, nPossibleBitVals=nPossibleBitVals)
  bitOrder = list(numpy.argsort(metrics))
  bitOrder.reverse()
  return bitOrder, metrics


def AnalyzeSparseVects(bitVects, actVals):
  """ #DOC

  **Arguments**

    - bitVects: a *sequence* containing SBVs

    - actVals: a *sequence*

   **Returns**

     a list of floats

   **Notes**

      - these need to be bit vects and binary activities

  """
  nPts = len(bitVects)
  if nPts != len(actVals):
    raise ValueError('var and activity lists should be the same length')
  nBits = bitVects[0].GetSize()

  actives = numpy.zeros(nBits, numpy.integer)
  inactives = numpy.zeros(nBits, numpy.integer)
  nActives, nInactives = 0, 0
  for i in range(nPts):
    sig, act = bitVects[i], actVals[i]
    onBitList = sig.GetOnBits()
    if act:
      for bit in onBitList:
        actives[bit] += 1
      nActives += 1
    else:
      for bit in onBitList:
        inactives[bit] += 1
      nInactives += 1
  resTbl = numpy.zeros((2, 2), numpy.integer)
  res = []
  gains = []
  for bit in range(nBits):
    nAct, nInact = actives[bit], inactives[bit]
    if nAct or nInact:
      resTbl[0, 0] = nAct
      resTbl[1, 0] = nPts - nAct
      resTbl[0, 1] = nInact
      resTbl[1, 1] = nPts - nInact
      gain = entropy.InfoGain(resTbl)
      gains.append(gain)
      res.append((bit, gain, nAct, nInact))
  return res, gains


def SparseRankBits(bitVects, actVals, metricFunc=AnalyzeSparseVects):
  """ Rank a set of bits according to a metric function

  **Arguments**

    - bitVects: a *sequence* containing SBVs

    - actVals: a *sequence*

    - metricFunc: (optional) the metric function to be used.  See _SparseCalcInfoGains()_
      for a description of the signature of this function.

   **Returns**

     A 2-tuple containing:

       - the relative order of the bits (a list of ints)

       - the metric calculated for each bit (a list of floats)

    **Notes**

      - these need to be bit vects and binary activities

  """
  info, metrics = metricFunc(bitVects, actVals)
  bitOrder = list(numpy.argsort(metrics))
  bitOrder.reverse()
  return bitOrder, info
