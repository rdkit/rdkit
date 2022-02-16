#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the signature utils

"""

import doctest
import itertools
import unittest

from rdkit.Chem.Pharm2D import Utils


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(Utils, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def testCounts(self):
    vals = [
      ((0, 1, 2), 4),
      ((0, 0, 0), 0),
      ((2, 2, 2), 9),
      ((1, 1, 2), 7),
      ((1, 2, 2), 8),
    ]
    for combo, tgt in vals:
      res = Utils.CountUpTo(3, 3, combo)
      assert res == tgt, f'Bad res ({res}) for combo {str((combo, tgt))}'

  def testGetTriangles(self):
    vals = [
      (2, 0, []),
      (3, 1, ((0, 1, 2), )),
      (4, 2, ((0, 1, 3), (1, 2, 4))),
      (5, 3, ((0, 1, 4), (1, 2, 5), (2, 3, 6))),
    ]
    for tpl in vals:
      nPts, cnt, tris = tpl
      r = Utils.GetTriangles(nPts)
      assert len(r) == cnt, f'bad triangle length {len(r)} for probe {str(tpl)}'
      assert r == tris, f'bad triangle list {str(r)} for probe {str(tpl)}'

  def testDistTriangleInequality(self):
    bins = [(1, 2), (2, 3), (5, 6)]
    vals = [
      ((0, 0, 0), 1),
      ((0, 0, 1), 1),
      ((0, 0, 2), 0),
      ((0, 1, 2), 1),
      ((1, 1, 2), 1),
    ]
    for tpl in vals:
      ds, tgt = tpl
      distBins = [bins[x] for x in ds]
      r = Utils.BinsTriangleInequality(distBins[0], distBins[1], distBins[2])
      assert r == tgt, f'bad result {r} for probe {str(tpl)}'

  def testLimitPharmacophores(self):
    bins = [(1, 2), (2, 3), (5, 6)]
    vals = [
      ((0, 0, 0), 1),
      ((0, 0, 1), 1),
      ((0, 0, 2), 0),
      ((0, 1, 2), 1),
      ((1, 1, 2), 1),
    ]
    for tpl in vals:
      ds, tgt = tpl
      r = Utils.ScaffoldPasses(ds, bins)
      assert r == tgt, f'bad result {r} for probe {str(tpl)}'

  def testGetPossiblePharmacophores(self):
    bins = [(1, 2), (2, 3), (5, 6)]
    vals = [
      (2, 3),
      (3, 24),
    ]
    for tpl in vals:
      num, tgt = tpl
      pphores = Utils.GetPossibleScaffolds(num, bins)
      cnt = len(pphores)
      assert cnt == tgt, f'bad pharmacophore count {cnt} for probe {str(tpl)}'
    self.assertEqual(Utils.GetPossibleScaffolds(1, bins), 0)

  def testOrderTriangle(self):
    # Additional tests to complement the doctests
    self.assertRaises(ValueError, Utils.OrderTriangle, [0, 2], [1, 2, 3])
    self.assertRaises(ValueError, Utils.OrderTriangle, [0, 2, 4], [1, 2])
    self.assertEqual(Utils.OrderTriangle([1, 3, 1], [2, 3, 4]), ([1, 3, 1], [4, 3, 2]))
    self.assertEqual(Utils.OrderTriangle([1, 3, 1], [4, 3, 2]), ([1, 3, 1], [4, 3, 2]))

    # If all the features are the same, we want the distances in reverse order
    for dist in itertools.permutations([1, 2, 3], 3):
      self.assertEqual(Utils.OrderTriangle([1, 1, 1], dist), ([1, 1, 1], [3, 2, 1]))

  def testUniquifyCombinations(self):
    combos = [[1, 2, 3], [3, 2, 1]]
    # Last equivalent combination is returned
    self.assertEqual(Utils.UniquifyCombinations(combos), [(3, 2, 1), ])

    combos = [[1], [1], [2]]
    # Last equivalent combination is returned
    self.assertEqual(sorted(Utils.UniquifyCombinations(combos)), sorted([(1, ), (2, )]))

  def _compareCombinations(self, actual, expected):
    if sorted(actual) != sorted(expected):
      print(actual)
    self.assertEqual(sorted(actual), sorted(expected), 'Combinations are different')
    if actual != expected:
      print(actual)
    self.assertEqual(actual, expected, 'Combinations are not in the same order')

  def testGetUniqueCombinations(self):
    x = [[(1, (11, )), (1, (12, )), (2, (31, ))], [(1, (11, )), (1, (12, )), (2, (32, ))],
         [(1, (11, )), (1, (13, )), (2, (31, ))], [(1, (11, )), (1, (13, )), (2, (32, ))],
         [(1, (11, )), (1, (14, )), (2, (31, ))], [(1, (11, )), (1, (14, )), (2, (32, ))],
         [(1, (12, )), (1, (13, )), (2, (31, ))], [(1, (12, )), (1, (13, )), (2, (32, ))],
         [(1, (12, )), (1, (14, )), (2, (31, ))], [(1, (12, )), (1, (14, )), (2, (32, ))],
         [(1, (13, )), (1, (14, )), (2, (31, ))], [(1, (13, )), (1, (14, )), (2, (32, ))]]

    # yapf: disable
    GetUniqueCombinations = Utils.GetUniqueCombinations

    choices = [[(11,), (12,)], [(21,), (22,)], [(31,), (32,)]]
    classes = [1, 2, 3]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(1, (11,)), (2, (21,)), (3, (31,))], [(1, (11,)), (2, (21,)), (3, (32,))],
       [(1, (11,)), (2, (22,)), (3, (31,))], [(1, (11,)), (2, (22,)), (3, (32,))],
       [(1, (12,)), (2, (21,)), (3, (31,))], [(1, (12,)), (2, (21,)), (3, (32,))],
       [(1, (12,)), (2, (22,)), (3, (31,))], [(1, (12,)), (2, (22,)), (3, (32,))]])

    classes = [1, 1, 2]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(1, (11,)), (1, (21,)), (2, (31,))], [(1, (11,)), (1, (21,)), (2, (32,))],
       [(1, (11,)), (1, (22,)), (2, (31,))], [(1, (11,)), (1, (22,)), (2, (32,))],
       [(1, (12,)), (1, (21,)), (2, (31,))], [(1, (12,)), (1, (21,)), (2, (32,))],
       [(1, (12,)), (1, (22,)), (2, (31,))], [(1, (12,)), (1, (22,)), (2, (32,))]])

    choices = [[(11,), (12,)], [(11,), (12,)], [(31,), (32,)]]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(1, (11,)), (1, (12,)), (2, (31,))], [(1, (11,)), (1, (12,)), (2, (32,))]])

    choices = [[(11,), (12,)], [(11,), (12,), (13,)], [(31,), (32,)]]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(1, (11,)), (1, (12,)), (2, (31,))], [(1, (11,)), (1, (12,)), (2, (32,))],
       [(1, (11,)), (1, (13,)), (2, (31,))], [(1, (11,)), (1, (13,)), (2, (32,))],
       [(1, (12,)), (1, (13,)), (2, (31,))], [(1, (12,)), (1, (13,)), (2, (32,))]])

    choices = [[(11,), (12,), (14,)], [(11,), (12,), (13,)], [(31,), (32,)]]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(1, (11,)), (1, (12,)), (2, (31,))], [(1, (11,)), (1, (12,)), (2, (32,))],
       [(1, (11,)), (1, (13,)), (2, (31,))], [(1, (11,)), (1, (13,)), (2, (32,))],
       [(1, (12,)), (1, (13,)), (2, (31,))], [(1, (12,)), (1, (13,)), (2, (32,))],
       [(1, (11,)), (1, (14,)), (2, (31,))], [(1, (11,)), (1, (14,)), (2, (32,))],
       [(1, (12,)), (1, (14,)), (2, (31,))], [(1, (12,)), (1, (14,)), (2, (32,))],
       [(1, (13,)), (1, (14,)), (2, (31,))], [(1, (13,)), (1, (14,)), (2, (32,))]])

    choices = [[(11, 111,), (12,), (14,)], [(11, 111,), (12,), (13,)], [(31,), (32,)]]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(1, (11, 111)), (1, (12,)), (2, (31,))], [(1, (11, 111)), (1, (12,)), (2, (32,))],
       [(1, (11, 111)), (1, (13,)), (2, (31,))], [(1, (11, 111)), (1, (13,)), (2, (32,))],
       [(1, (12,)), (1, (13,)), (2, (31,))], [(1, (12,)), (1, (13,)), (2, (32,))],
       [(1, (11, 111)), (1, (14,)), (2, (31,))], [(1, (11, 111)), (1, (14,)), (2, (32,))],
       [(1, (12,)), (1, (14,)), (2, (31,))], [(1, (12,)), (1, (14,)), (2, (32,))],
       [(1, (13,)), (1, (14,)), (2, (31,))], [(1, (13,)), (1, (14,)), (2, (32,))]])

    choices = [[(11,), (12,), ], [(11,), (13,), ], [(11,), (14,), ]]
    classes = [1, 1, 1]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(1, (11,)), (1, (13,)), (1, (14,))], [(1, (11,)), (1, (12,)), (1, (14,))],
       [(1, (11,)), (1, (12,)), (1, (13,))], [(1, (12,)), (1, (13,)), (1, (14,))]])

    choices = [[(0,), (4,)], [(0,)]]
    classes = [0, 1]
    self._compareCombinations(GetUniqueCombinations(choices, classes), [[(0, (4,)), (1, (0,))]])

    choices = [[], [], []]
    classes = [0, 1, 1]
    self.assertEqual(GetUniqueCombinations(choices, classes), [])

    choices = [[(0,), (4,)], [(0,), (4,)], [(0,), (4,)]]
    classes = [0, 0, 0]
    self.assertEqual(GetUniqueCombinations(choices, classes), [])

    choices = [[(0,), (4,)], [(0,), (4,)], []]
    classes = [0, 0, 1]
    self.assertEqual(GetUniqueCombinations(choices, classes), [])

    choices = [[(4,), (6,), (7,), (10,)], [(2,), (3,), (5,), (13,), (14,)],
               [(2,), (3,), (5,), (13,), (14,)]]
    classes = [0, 1, 1]
    self._compareCombinations(
      GetUniqueCombinations(choices, classes),
      [[(0, (4,)), (1, (2,)), (1, (3,))], [(0, (4,)), (1, (2,)), (1, (5,))],
       [(0, (4,)), (1, (2,)), (1, (13,))], [(0, (4,)), (1, (2,)), (1, (14,))],
       [(0, (4,)), (1, (3,)), (1, (5,))], [(0, (4,)), (1, (3,)), (1, (13,))],
       [(0, (4,)), (1, (3,)), (1, (14,))], [(0, (4,)), (1, (5,)), (1, (13,))],
       [(0, (4,)), (1, (5,)), (1, (14,))], [(0, (4,)), (1, (13,)), (1, (14,))],
       [(0, (6,)), (1, (2,)), (1, (3,))], [(0, (6,)), (1, (2,)), (1, (5,))],
       [(0, (6,)), (1, (2,)), (1, (13,))], [(0, (6,)), (1, (2,)), (1, (14,))],
       [(0, (6,)), (1, (3,)), (1, (5,))], [(0, (6,)), (1, (3,)), (1, (13,))],
       [(0, (6,)), (1, (3,)), (1, (14,))], [(0, (6,)), (1, (5,)), (1, (13,))],
       [(0, (6,)), (1, (5,)), (1, (14,))], [(0, (6,)), (1, (13,)), (1, (14,))],
       [(0, (7,)), (1, (2,)), (1, (3,))], [(0, (7,)), (1, (2,)), (1, (5,))],
       [(0, (7,)), (1, (2,)), (1, (13,))], [(0, (7,)), (1, (2,)), (1, (14,))],
       [(0, (7,)), (1, (3,)), (1, (5,))], [(0, (7,)), (1, (3,)), (1, (13,))],
       [(0, (7,)), (1, (3,)), (1, (14,))], [(0, (7,)), (1, (5,)), (1, (13,))],
       [(0, (7,)), (1, (5,)), (1, (14,))], [(0, (7,)), (1, (13,)), (1, (14,))],
       [(0, (10,)), (1, (2,)), (1, (3,))], [(0, (10,)), (1, (2,)), (1, (5,))],
       [(0, (10,)), (1, (2,)), (1, (13,))], [(0, (10,)), (1, (2,)), (1, (14,))],
       [(0, (10,)), (1, (3,)), (1, (5,))], [(0, (10,)), (1, (3,)), (1, (13,))],
       [(0, (10,)), (1, (3,)), (1, (14,))], [(0, (10,)), (1, (5,)), (1, (13,))],
       [(0, (10,)), (1, (5,)), (1, (14,))], [(0, (10,)), (1, (13,)), (1, (14,))]])
    # yapf: enable


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
