#
#  Copyright (C) 2001  greg Landrum
#

""" unit testing code for SLT Risk functions

"""
from __future__ import print_function
import unittest
from rdkit.ML.SLT import Risk
import math
import numpy    

class TestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(), end='')
    self.dList=[(1, 40),
                (2, 40),
                (3, 21),
                (7, 16),
                (8, 12),
                (9, 11),
                (10, 11)]
    self.nPts = 95
    self.eps = min(4./math.sqrt(self.nPts),1.)
    self.eps2 = .10
    self.tol = 1e-4
    
  def testBurges(self):
    " testing Burges empirical risk bound "
    res = numpy.array(map(lambda x,y=self.nPts,z=self.eps2:Risk.BurgesRiskBound(x[0],y,x[1],z),
                          self.dList))
    res = numpy.array([Risk.BurgesRiskBound(x[0],self.nPts,x[1],self.eps2) 
                       for x in self.dList])
    target = numpy.array([.7445,.8157,.6698,.7649,.7506,.7658,.7896])
    maxDev = max(abs(res-target))
    assert maxDev < self.tol,'maxDev too high' 

  def testCherkassky(self):
    " testing Cherkassky empirical risk bound "
    res = numpy.array([Risk.CherkasskyRiskBound(x[0],self.nPts,x[1],self.eps) 
                       for x in self.dList])
    target = numpy.array([.6654,.7450,.5378,.6329,.6010,.6175,.6501])
    maxDev = max(abs(res-target))
    assert maxDev < self.tol,'maxDev too high'


if __name__ == '__main__':
  unittest.main()

