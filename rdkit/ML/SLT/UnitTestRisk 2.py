#
#  Copyright (C) 2001  greg Landrum
#
""" unit testing code for SLT Risk functions

"""

import unittest
from rdkit.ML.SLT import Risk
import math
import numpy


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dList = [(1, 40), (2, 40), (3, 21), (7, 16), (8, 12), (9, 11), (10, 11)]
    self.nPts = 95
    self.eps = min(4. / math.sqrt(self.nPts), 1.)
    self.eps2 = .10
    self.tol = 1e-4

  def testBurges(self):
    " testing Burges empirical risk bound "
    res = numpy.array([Risk.BurgesRiskBound(x[0], self.nPts, x[1], self.eps2) for x in self.dList])
    target = numpy.array([.7445, .8157, .6698, .7649, .7506, .7658, .7896])
    maxDev = max(abs(res - target))
    assert maxDev < self.tol, 'maxDev too high'

  def testCherkassky(self):
    " testing Cherkassky empirical risk bound "
    res = numpy.array(
      [Risk.CherkasskyRiskBound(x[0], self.nPts, x[1], self.eps) for x in self.dList])
    target = numpy.array([.6654, .7450, .5378, .6329, .6010, .6175, .6501])
    maxDev = max(abs(res - target))
    assert maxDev < self.tol, 'maxDev too high'

  def testChristiani(self):
    " testing Christiani empirical risk bound "
    res = numpy.array(
      [Risk.CristianiRiskBound(x[0], self.nPts, x[1], self.eps) for x in self.dList])
    target = numpy.array([1.5617, 1.7438, 1.4797, 1.7394, 1.7235, 1.7653, 1.8235])
    maxDev = max(abs(res - target))
    assert maxDev < self.tol, 'maxDev too high'

  def test_log2(self):
    self.assertEqual(Risk.log2(1), 0)
    self.assertEqual(Risk.log2(2), 1)
    self.assertEqual(Risk.log2(4), 2)
    self.assertEqual(Risk.log2(0.5), -1)
    self.assertEqual(Risk.log2(0.25), -2)
    self.assertEqual(Risk.log2(0.125), -3)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
