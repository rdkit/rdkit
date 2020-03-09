#
#  Copyright (C) 2001  greg Landrum
#
""" unit testing code for variable quantization

"""
import unittest

from rdkit.ML.Data import Quantize


class TestCase(unittest.TestCase):

    def testOneSplit1(self):
        # """ simple case (clear division) """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 0), (1.4, 0), (1.6, 0), (2., 1), (2.1, 1), (2.2, 1),
             (2.3, 1)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 2
        res = Quantize.FindVarQuantBound(varValues, resCodes, nPossibleRes)
        target = (1.8, 0.97095)
        self.assertEqual(
          [Quantize.feq(x, y, 1e-4) for x, y in zip(res, target)], [1, 1],
          'result comparison failed: %s != %s' % (res, target))

    def testOneSplit2_noise(self):
        # """ some noise """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 0), (1.4, 1), (1.6, 0), (2., 1), (2.1, 1), (2.2, 1),
             (2.3, 1)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 2
        res = Quantize.FindVarQuantBound(varValues, resCodes, nPossibleRes)
        target = (1.8, 0.60999)
        self.assertEqual(
          [Quantize.feq(x, y, 1e-4) for x, y in zip(res, target)], [1, 1],
          'result comparison failed: %s != %s' % (res, target))

    def testOneSplit3(self):
        # """ optimal division not possible """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 2), (1.4, 2), (1.6, 2), (2., 2), (2.1, 1), (2.2, 1),
             (2.3, 1)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 3
        res = Quantize.FindVarQuantBound(varValues, resCodes, nPossibleRes)
        target = (1.3, 0.88129)
        self.assertEqual(
          [Quantize.feq(x, y, 1e-4) for x, y in zip(res, target)], [1, 1],
          'result comparison failed: %s != %s' % (res, target))

    def testOneSplit4_duplicates(self):
        # """ lots of duplicates """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.2, 1), (1.4, 0), (1.4, 0), (1.6, 0), (2., 1), (2.1, 1),
             (2.1, 1), (2.1, 1), (2.1, 1), (2.2, 1), (2.3, 1)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 2
        res = Quantize.FindVarQuantBound(varValues, resCodes, nPossibleRes)
        target = (1.8, 0.68939)
        self.assertEqual(
          [Quantize.feq(x, y, 1e-4) for x, y in zip(res, target)], [1, 1],
          'result comparison failed: %s != %s' % (res, target))

    def testOneSplit5_outOfOrder(self):
        # """ same as testOneSplit1 data, but out of order """
        d = [(1., 0), (1.1, 0), (2.2, 1), (1.2, 0), (1.6, 0), (1.4, 0), (2., 1), (2.1, 1), (1.4, 0),
             (2.3, 1)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 2
        res = Quantize.FindVarQuantBound(varValues, resCodes, nPossibleRes)
        target = (1.8, 0.97095)
        self.assertEqual(
          [Quantize.feq(x, y, 1e-4) for x, y in zip(res, target)], [1, 1],
          'result comparison failed: %s != %s' % (res, target))

    def testMultSplit1_simple_dual(self):
        # """ simple dual split """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 2), (1.4, 2), (1.6, 2), (2., 2), (2.1, 1), (2.1, 1),
             (2.1, 1), (2.2, 1), (2.3, 1)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 3
        res = Quantize.FindVarMultQuantBounds(varValues, 2, resCodes, nPossibleRes)
        target = ([1.3, 2.05], 1.55458)
        self.assertEqual(
          min([Quantize.feq(x, y, 1e-4) for x, y in zip(res[0], target[0])]), 1,
          'split bound comparison failed: %s != %s' % (res[0], target[0]))
        self.assertTrue(
          Quantize.feq(res[1], target[1], 1e-4),
          'InfoGain comparison failed: %s != %s' % (res[1], target[1]))

    def testMultSplit2_outOfOrder(self):
        # """ same test as testMultSplit1, but out of order """
        d = [(1., 0), (2.1, 1), (1.1, 0), (1.2, 0), (1.4, 2), (1.6, 2), (2., 2), (1.4, 2), (2.1, 1),
             (2.2, 1), (2.1, 1), (2.3, 1)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 3
        res = Quantize.FindVarMultQuantBounds(varValues, 2, resCodes, nPossibleRes)
        target = ([1.3, 2.05], 1.55458)
        self.assertTrue(
          Quantize.feq(res[1], target[1], 1e-4),
          'InfoGain comparison failed: %s != %s' % (res[1], target[1]))
        self.assertEqual(
          min([Quantize.feq(x, y, 1e-4) for x, y in zip(res[0], target[0])]), 1,
          'split bound comparison failed: %s != %s' % (res[0], target[0]))

    def testMultSplit3_4results(self):
        # """  4 possible results """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 2), (1.4, 2), (1.6, 2), (2., 2), (2.1, 1), (2.1, 1),
             (2.1, 1), (2.2, 1), (2.3, 1), (3.0, 3), (3.1, 3), (3.2, 3), (3.3, 3)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 4
        res = Quantize.FindVarMultQuantBounds(varValues, 3, resCodes, nPossibleRes)
        target = ([1.30, 2.05, 2.65], 1.97722)
        self.assertTrue(
          Quantize.feq(res[1], target[1], 1e-4),
          'InfoGain comparison failed: %s != %s' % (res[1], target[1]))
        self.assertEqual(
          min([Quantize.feq(x, y, 1e-4) for x, y in zip(res[0], target[0])]), 1,
          'split bound comparison failed: %s != %s' % (res[0], target[0]))

    def testMultSplit4_dualValued_island(self):
        # """ dual valued, with an island """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 1), (1.4, 1), (1.6, 1), (2., 1), (2.1, 0), (2.1, 0),
             (2.1, 0), (2.2, 0), (2.3, 0)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 2
        res = Quantize.FindVarMultQuantBounds(varValues, 2, resCodes, nPossibleRes)
        target = ([1.3, 2.05], .91830)
        self.assertTrue(
          Quantize.feq(res[1], target[1], 1e-4),
          'InfoGain comparison failed: %s != %s' % (res[1], target[1]))
        self.assertEqual(
          min([Quantize.feq(x, y, 1e-4) for x, y in zip(res[0], target[0])]), 1,
          'split bound comparison failed: %s != %s' % (res[0], target[0]))

    def testMultSplit5_dualValued_island_noisy(self):
        # """ dual valued, with an island, a bit noisy """
        d = [(1., 0), (1.1, 0), (1.2, 0), (1.4, 1), (1.4, 0), (1.6, 1), (2., 1), (2.1, 0), (2.1, 0),
             (2.1, 0), (2.2, 1), (2.3, 0)]
        varValues, resCodes = zip(*d)
        nPossibleRes = 2
        res = Quantize.FindVarMultQuantBounds(varValues, 2, resCodes, nPossibleRes)
        target = ([1.3, 2.05], .34707)
        self.assertTrue(
          Quantize.feq(res[1], target[1], 1e-4),
          'InfoGain comparison failed: %s != %s' % (res[1], target[1]))
        self.assertEqual(
          min([Quantize.feq(x, y, 1e-4) for x, y in zip(res[0], target[0])]), 1,
          'split bound comparison failed: %s != %s' % (res[0], target[0]))

    def test9NewSplits(self):
        d = [(0, 0),
             (1, 1),
             (2, 0), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 2], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 2], str(res))

        d = [(0, 1),
             (1, 0),
             (2, 1), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 2], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 2], str(res))

        d = [(0, 0),
             (0, 0),
             (1, 1),
             (1, 1),
             (2, 0),
             (2, 1), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))

        d = [(0, 0),
             (0, 1),
             (1, 1),
             (1, 1),
             (2, 0),
             (2, 1), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))

        d = [(0, 0),
             (0, 0),
             (1, 0),
             (1, 1),
             (2, 0),
             (2, 1), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))

        d = [(0, 0),
             (0, 0),
             (1, 0),
             (1, 0),
             (2, 1),
             (2, 1), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [4], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [4], str(res))

        d = [(0, 0),
             (0, 0),
             (1, 1),
             (1, 1),
             (2, 1),
             (2, 1), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2], str(res))

        d = [(0, 0),
             (0, 0),
             (1, 0),
             (1, 0),
             (2, 0),
             (2, 0), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [], str(res))

        d = [(0, 0),
             (0, 1),
             (1, 0),
             (1, 1),
             (2, 0),
             (2, 0), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [2, 4], str(res))

        d = [(1, 0),
             (2, 1),
             (2, 1),
             (3, 1),
             (3, 1),
             (3, 1),
             (4, 0),
             (4, 1),
             (4, 1), ]
        varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 6], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 6], str(res))

        d = [(1, 1.65175902843, 0), (2, 1.89935600758, 0), (3, 1.89935600758, 1), (4, 1.89935600758, 1),
             (5, 2.7561609745, 1), (6, 2.7561609745, 1), (7, 2.7561609745, 1), (8, 2.7561609745, 1),
             (9, 3.53454303741, 1), (10, 3.53454303741, 1), (11, 3.53454303741, 1),
             (12, 3.53454303741, 1), (13, 3.53454303741, 1)]

        _, varValues, resCodes = zip(*d)
        res = Quantize._NewPyFindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 4], str(res))
        res = Quantize._FindStartPoints(varValues, resCodes, len(d))
        self.assertTrue(res == [1, 4], str(res))

    def testGithubIssue18(self):
        d = [0, 1, 2, 3, 4]
        a = [0, 0, 1, 1, 1]
        _ = Quantize.FindVarMultQuantBounds(d, 1, a, 2)

        d2 = [(x, ) for x in d]
        self.assertRaises(ValueError, lambda: Quantize.FindVarMultQuantBounds(d2, 1, a, 2))
        self.assertRaises(ValueError, lambda: Quantize._FindStartPoints(d2, a, len(d2)))


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
