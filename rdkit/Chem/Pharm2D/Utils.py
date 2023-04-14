#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" utility functionality for the 2D pharmacophores code

  See Docs/Chem/Pharm2D.triangles.jpg for an illustration of the way
  pharmacophores are broken into triangles and labelled.

  See Docs/Chem/Pharm2D.signatures.jpg for an illustration of bit
  numbering

"""
import itertools

try:
  # his is available in Python 3.8: https://docs.python.org/3/library/math.html
  from math import comb
except ImportError:
  from functools import lru_cache

  def _CheckCombArgument(n: int, k: int) -> None:
    if not isinstance(n, int) or not isinstance(k, int):
      raise ValueError(f"n ({n}) and k ({k}) must be positive integers")
    if n < 0:
      raise ValueError(f"n ({n}) must be a positive integer")
    if k < 0:
      raise ValueError(f"k ({k}) must be a positive integer")

  @lru_cache(maxsize=128)
  def comb(n: int, k: int) -> int:
    """
        Return the number of ways to choose k items from n items without repetition and without order.
        Evaluates to n! / (k! * (n - k)!) when k <= n and evaluates to zero when k > n.
        
        Arguments:
        ----------
            - n (integer): the number of items to choose from 
            - k (integer): the number of items to choose
        
        Returns:
        ----------
            - integer: the number of ways to choose k items from n items without repetition and without order
        
        
        Optimization reference:
        1) https://github.com/python/cpython/blob/main/Modules/mathmodule.c (Line 3381 <-> 3568)
        2) https://github.com/python/cpython/commit/60c320c38e4e95877cde0b1d8562ebd6bc02ac61
        3) https://bugs.python.org/issue37295 
        
        """
    _CheckCombArgument(n, k)
    if k > n:
      return 0
    if k == 0 or n == 0 or k == n:
      return 1
    if n - k > k:
      k = n - k

    # Optimization for small arguments
    if k < 8 or (16 * k < n < 512) or (k < n < 32):
      res = 1
      for i in range(k):
        res *= n - i
        res = res // (i + 1)
      return res

    j = k // 2
    return comb(n, j) * comb(n - j, k - j) // comb(k, j)


#
#  number of points in a scaffold -> sequence of distances (p1, p2) in
#   the scaffold
#
nPointDistDict = {
  2: ((0, 1), ),
  3: ((0, 1), (0, 2), (1, 2)),
  4: ((0, 1), (0, 2), (0, 3), (1, 2), (2, 3)),
  5: ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)),
  6: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5)),
  7: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)),
  8: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (1, 2), (2, 3), (3, 4), (4, 5),
      (5, 6), (6, 7)),
  9: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (1, 2), (2, 3), (3, 4),
      (4, 5), (5, 6), (6, 7), (7, 8)),
  10: ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6), (0, 7), (0, 8), (0, 9), (1, 2), (2, 3),
       (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)),
}

#
#  number of distances in a scaffold -> number of points in the scaffold
#
nDistPointDict = {
  1: 2,
  3: 3,
  5: 4,
  7: 5,
  9: 6,
  11: 7,
  13: 8,
  15: 9,
  17: 10,
}

_trianglesInPharmacophore = {}


def GetTriangles(nPts):
  """ returns a tuple with the distance indices for
     triangles composing an nPts-pharmacophore

    """
  global _trianglesInPharmacophore
  if nPts < 3:
    return []
  res = _trianglesInPharmacophore.get(nPts, [])
  if not res:
    idx1, idx2, idx3 = (0, 1, nPts - 1)
    while idx1 < nPts - 2:
      res.append((idx1, idx2, idx3))
      idx1 += 1
      idx2 += 1
      idx3 += 1
    res = tuple(res)
    _trianglesInPharmacophore[nPts] = res
  return res


def BinsTriangleInequality(d1, d2, d3):
  """ checks the triangle inequality for combinations
      of distance bins.

      the general triangle inequality is:
         d1 + d2 >= d3
      the conservative binned form of this is:
         d1(upper) + d2(upper) >= d3(lower)

    """
  if d1[1] + d2[1] < d3[0]:
    return False
  if d2[1] + d3[1] < d1[0]:
    return False
  if d3[1] + d1[1] < d2[0]:
    return False
  return True


def ScaffoldPasses(combo, bins=None):
  """ checks the scaffold passed in to see if all
    contributing triangles can satisfy the triangle inequality

    the scaffold itself (encoded in combo) is a list of binned distances

    """
  # this is the number of points in the pharmacophore
  tris = GetTriangles(nDistPointDict[len(combo)])
  for tri in tris:
    ds = [bins[combo[x]] for x in tri]
    if not BinsTriangleInequality(ds[0], ds[1], ds[2]):
      return False
  return True


_numCombDict = {}


def NumCombinations(nItems, nSlots):
  """  returns the number of ways to fit nItems into nSlots

      We assume that (x, y) and (y, x) are equivalent, and
      (x, x) is allowed.

      General formula is, for N items and S slots:
        res = (N+S-1)! / ( (N-1)! * S! )

    """
  global _numCombDict
  res = _numCombDict.get((nItems, nSlots), -1)
  if res == -1:
    res = comb(nItems + nSlots - 1, nSlots)
    _numCombDict[(nItems, nSlots)] = res
  return res


_verbose = 0

_countCache = {}


def CountUpTo(nItems, nSlots, vs, idx=0, startAt=0):
  """ Figures out where a given combination of indices would
     occur in the combinatorial explosion generated by _GetIndexCombinations_

     **Arguments**

       - nItems: the number of items to distribute

       - nSlots: the number of slots in which to distribute them

       - vs: a sequence containing the values to find

       - idx: used in the recursion

       - startAt: used in the recursion

    **Returns**

       an integer

    """
  global _countCache
  if _verbose:
    print('  ' * idx, f'CountUpTo({idx})', vs[idx], startAt)
  if idx == 0 and (nItems, nSlots, tuple(vs)) in _countCache:
    return _countCache[(nItems, nSlots, tuple(vs))]
  elif idx >= nSlots:
    accum = 0
  elif idx == nSlots - 1:
    accum = vs[idx] - startAt
  else:
    accum = 0
    # get the digit at idx correct
    for i in range(startAt, vs[idx]):
      nLevsUnder = nSlots - idx - 1
      nValsOver = nItems - i
      numCombs = NumCombinations(nValsOver, nLevsUnder)
      if _verbose:
        print('  ' * idx, ' ', i, nValsOver, nLevsUnder, numCombs)
      accum += numCombs
    accum += CountUpTo(nItems, nSlots, vs, idx + 1, vs[idx])
  if _verbose:
    print('  ' * idx, '>', accum)
  if idx == 0:
    _countCache[(nItems, nSlots, tuple(vs))] = accum
  return accum


_indexCombinations = {}


def GetIndexCombinations(nItems, nSlots, slot=0, lastItemVal=0):
  """ Generates all combinations of nItems in nSlots without including
      duplicates

    **Arguments**

      - nItems: the number of items to distribute

      - nSlots: the number of slots in which to distribute them

      - slot: used in recursion

      - lastItemVal: used in recursion

    **Returns**

      a list of lists

    """
  global _indexCombinations
  if not slot and (nItems, nSlots) in _indexCombinations:
    res = _indexCombinations[(nItems, nSlots)]
  elif slot >= nSlots:
    res = []
  elif slot == nSlots - 1:
    res = [[x] for x in range(lastItemVal, nItems)]
  else:
    res = []
    for x in range(lastItemVal, nItems):
      tmp = GetIndexCombinations(nItems, nSlots, slot + 1, x)
      for entry in tmp:
        res.append([x] + entry)
    if not slot:
      _indexCombinations[(nItems, nSlots)] = res
  return res


def GetAllCombinations(choices, noDups=1, which=0):
  """  Does the combinatorial explosion of the possible combinations
    of the elements of _choices_.

    **Arguments**

      - choices: sequence of sequences with the elements to be enumerated

      - noDups: (optional) if this is nonzero, results with duplicates,
        e.g. (1,1,0), will not be generated

      - which: used in recursion

    **Returns**

      a list of lists

    >>> GetAllCombinations([(0, ), (1, ), (2, )])
    [[0, 1, 2]]
    >>> GetAllCombinations([(0, ), (1, 3), (2, )])
    [[0, 1, 2], [0, 3, 2]]

    >>> GetAllCombinations([(0, 1), (1, 3), (2, )])
    [[0, 1, 2], [0, 3, 2], [1, 3, 2]]

    """
  if which >= len(choices):
    return []
  elif which == len(choices) - 1:
    return [[x] for x in choices[which]]

  res = []
  tmp = GetAllCombinations(choices, noDups=noDups, which=which + 1)
  for thing in choices[which]:
    for other in tmp:
      if not noDups or thing not in other:
        res.append([thing] + other)
  return res


def GetUniqueCombinations(choices, classes, which=0):
  """  Does the combinatorial explosion of the possible combinations
    of the elements of _choices_.

    """
  #   print(choices, classes)
  assert len(choices) == len(classes)
  if which >= len(choices):
    return []
  if which == len(choices) - 1:
    return [[(classes[which], x)] for x in choices[which]]

  res = []
  tmp = GetUniqueCombinations(choices, classes, which=which + 1)
  for thing in choices[which]:
    for other in tmp:
      if not any(x[1] == thing for x in other):
        newL = [(classes[which], thing)] + other
        newL.sort()
        if newL not in res:
          res.append(newL)
  return res


def GetUniqueCombinations_new(choices, classes, which=0):
  """  Does the combinatorial explosion of the possible combinations
    of the elements of _choices_.

    """
  #   print(choices, classes)
  assert len(choices) == len(classes)
  combos = set()
  for choice in itertools.product(*choices):
    # If a choice occurs in more than one of the fields, we ignore this case
    if len(set(choice)) != len(choice):
      continue
    combos.add(tuple(sorted((cls, ch) for cls, ch in zip(classes, choice))))
  return [list(combo) for combo in sorted(combos)]


def UniquifyCombinations(combos):
  """ uniquifies the combinations in the argument

      **Arguments**:

        - combos: a sequence of sequences

      **Returns**

        - a list of tuples containing the unique combos

    """
  resD = {}
  for combo in combos:
    k = combo[:]
    k.sort()
    resD[tuple(k)] = tuple(combo)
  return list(resD.values())


def GetPossibleScaffolds(nPts, bins, useTriangleInequality=True):
  """ gets all realizable scaffolds (passing the triangle inequality) with the
       given number of points and returns them as a list of tuples

    """
  if nPts < 2:
    return 0
  if nPts == 2:
    return [(x, ) for x in range(len(bins))]

  nDists = len(nPointDistDict[nPts])
  combos = GetAllCombinations([range(len(bins))] * nDists, noDups=0)
  return [
    tuple(combo) for combo in combos if not useTriangleInequality or ScaffoldPasses(combo, bins)
  ]


def OrderTriangle(featIndices, dists):
  """
      put the distances for a triangle into canonical order

      It's easy if the features are all different:

      >>> OrderTriangle([0, 2, 4], [1, 2, 3])
      ([0, 2, 4], [1, 2, 3])

      It's trickiest if they are all the same:
      
      >>> OrderTriangle([0, 0, 0], [1, 2, 3])
      ([0, 0, 0], [3, 2, 1])
      >>> OrderTriangle([0, 0, 0], [2, 1, 3])
      ([0, 0, 0], [3, 2, 1])
      >>> OrderTriangle([0, 0, 0], [1, 3, 2])
      ([0, 0, 0], [3, 2, 1])
      >>> OrderTriangle([0, 0, 0], [3, 1, 2])
      ([0, 0, 0], [3, 2, 1])
      >>> OrderTriangle([0, 0, 0], [3, 2, 1])
      ([0, 0, 0], [3, 2, 1])

      >>> OrderTriangle([0, 0, 1], [3, 2, 1])
      ([0, 0, 1], [3, 2, 1])
      >>> OrderTriangle([0, 0, 1], [1, 3, 2])
      ([0, 0, 1], [1, 3, 2])
      >>> OrderTriangle([0, 0, 1], [1, 2, 3])
      ([0, 0, 1], [1, 3, 2])

    """
  if len(featIndices) != 3:
    raise ValueError('bad indices')
  if len(dists) != 3:
    raise ValueError('bad dists')

  fs = set(featIndices)
  if len(fs) == 3:
    return featIndices, dists

  dSums = [0] * 3
  dSums[0] = dists[0] + dists[1]
  dSums[1] = dists[0] + dists[2]
  dSums[2] = dists[1] + dists[2]
  mD = max(dSums)
  if len(fs) == 1:
    if dSums[0] == mD:
      if dists[0] > dists[1]:
        ireorder = (0, 1, 2)
        dreorder = (0, 1, 2)
      else:
        ireorder = (0, 2, 1)
        dreorder = (1, 0, 2)
    elif dSums[1] == mD:
      if dists[0] > dists[2]:
        ireorder = (1, 0, 2)
        dreorder = (0, 2, 1)
      else:
        ireorder = (1, 2, 0)
        dreorder = (2, 0, 1)
    else:
      if dists[1] > dists[2]:
        ireorder = (2, 0, 1)
        dreorder = (1, 2, 0)
      else:
        ireorder = (2, 1, 0)
        dreorder = (2, 1, 0)
  else:
    # two classes
    if featIndices[0] == featIndices[1]:
      if dists[1] > dists[2]:
        ireorder = (0, 1, 2)
        dreorder = (0, 1, 2)
      else:
        ireorder = (1, 0, 2)
        dreorder = (0, 2, 1)
    elif featIndices[0] == featIndices[2]:
      if dists[0] > dists[2]:
        ireorder = (0, 1, 2)
        dreorder = (0, 1, 2)
      else:
        ireorder = (2, 1, 0)
        dreorder = (2, 1, 0)
    else:  # featIndices[1]==featIndices[2]:
      if dists[0] > dists[1]:
        ireorder = (0, 1, 2)
        dreorder = (0, 1, 2)
      else:
        ireorder = (0, 2, 1)
        dreorder = (1, 0, 2)

  dists = [dists[x] for x in dreorder]
  featIndices = [featIndices[x] for x in ireorder]
  return featIndices, dists


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
