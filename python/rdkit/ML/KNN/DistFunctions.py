# $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#      All Rights Reserved
#

import math


def EuclideanDist(ex1, ex2, attrs):
  """
    >>> v1 = [0,1,0,1]
    >>> v2 = [1,0,1,0]
    >>> EuclideanDist(v1,v2,range(4))
    2.0
    >>> EuclideanDist(v1,v1,range(4))
    0.0
    >>> v2 = [0,0,0,1]
    >>> EuclideanDist(v1,v2,range(4))
    1.0
    >>> v2 = [0,.5,0,.5]
    >>> abs(EuclideanDist(v1,v2,range(4))-1./math.sqrt(2))<1e-4
    1

    """
  dist = 0.0
  for i in attrs:
    dist += (ex1[i] - ex2[i])**2
  dist = math.sqrt(dist)
  return dist


def TanimotoDist(ex1, ex2, attrs):
  """
    >>> v1 = [0,1,0,1]
    >>> v2 = [1,0,1,0]
    >>> TanimotoDist(v1,v2,range(4))
    1.0
    >>> v2 = [1,0,1,1]
    >>> TanimotoDist(v1,v2,range(4))
    0.75
    >>> TanimotoDist(v2,v2,range(4))
    0.0

    # this tests Issue 122
    >>> v3 = [0,0,0,0]
    >>> TanimotoDist(v3,v3,range(4))
    1.0

    """
  inter = 0.0
  unin = 0.0
  for i in attrs:
    if (ex1[i] or ex2[i]):
      unin += 1
      if (ex1[i] and ex2[i]):
        inter += 1
  if (unin != 0.0):
    return (1 - inter / unin)
  else:
    return 1.0


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import sys
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
