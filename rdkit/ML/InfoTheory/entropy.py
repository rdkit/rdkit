#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#
""" Informational Entropy functions

  The definitions used are the same as those in Tom Mitchell's
  book "Machine Learning"

"""
import math

import numpy

# try to get the C versions of these routines
try:
  import rdkit.ML.InfoTheory.rdInfoTheory as cEntropy
except ImportError:
  hascEntropy = 0
else:
  hascEntropy = 1

# it's pretty obvious what this is for ;-)
_log2 = math.log(2)


def PyInfoEntropy(results):
  """ Calculates the informational entropy of a set of results.

  **Arguments**

    results is a 1D Numeric array containing the number of times a
    given set hits each possible result.
    For example, if a function has 3 possible results, and the
      variable in question hits them 5, 6 and 1 times each,
      results would be [5,6,1]

  **Returns**

    the informational entropy

  """
  nInstances = float(sum(results))
  if nInstances == 0:
    # to return zero or one... that is the question
    return 0
  probs = results / nInstances

  # -------
  #  NOTE: this is a little hack to allow the use of Numeric
  #   functionality to calculate the informational entropy.
  #    The problem is that the system log function pitches a fit
  #    when you call log(0.0).  We are perfectly happy with that
  #    returning *anything* because we're gonna multiply by 0 anyway.

  # Here's the risky (but marginally faster way to do it:
  #    add a small number to probs and hope it doesn't screw
  #    things up too much.
  # t = probs+1e-10

  # Here's a perfectly safe approach that's a little bit more obfuscated
  #  and a tiny bit slower
  t = numpy.choose(numpy.greater(probs, 0.0), (1, probs))
  return sum(-probs * numpy.log(t) / _log2)


def PyInfoGain(varMat):
  """ calculates the information gain for a variable

    **Arguments**

      varMat is a Numeric array with the number of possible occurrences
        of each result for reach possible value of the given variable.

      So, for a variable which adopts 4 possible values and a result which
        has 3 possible values, varMat would be 4x3

    **Returns**

      The expected information gain
  """
  variableRes = numpy.sum(varMat, 1)  # indexed by variable, Sv in Mitchell's notation
  overallRes = numpy.sum(varMat, 0)  # indexed by result, S in Mitchell's notation

  term2 = 0
  for i in range(len(variableRes)):
    term2 = term2 + variableRes[i] * InfoEntropy(varMat[i])
  tSum = sum(overallRes)
  if tSum != 0.0:
    term2 = 1. / tSum * term2
    gain = InfoEntropy(overallRes) - term2
  else:
    gain = 0
  return gain


# if we have the C versions, use them, otherwise use the python stuff
if hascEntropy:
  InfoEntropy = cEntropy.InfoEntropy
  InfoGain = cEntropy.InfoGain
else:
  InfoEntropy = PyInfoEntropy
  InfoGain = PyInfoGain
