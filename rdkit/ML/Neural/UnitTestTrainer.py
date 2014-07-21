## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

#
#  Copyright (C) 2000  greg Landrum
#
""" unit tests for the Neural network trainer implementation

   this basically works out **all** of the network code

"""
from __future__ import print_function
import unittest
from rdkit.ML.Neural import Network,Trainers
import numpy
import random
random.seed(23)

class TrainerTestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
    self.trainTol = 0.3
    self.orExamples = [
        [0,0,1,0.1],
        [0,1,1,.9],
        [1,0,1,.9],
        [1,1,1,.9]
        ]
    self.andExamples = [
        [0,0,1,0.1],
        [0,1,1,.1],
        [1,0,1,.1],
        [1,1,1,.9]
        ]
    self.xorExamples = [
        [0,0,1,.1],
        [0,1,1,.9],
        [1,0,1,.9],
        [1,1,1,.1]
        ]
    self.linExamples = [
          [.1,.1],
          [.2,.2],
          [.3,.3],
          [.4,.4],
          [.8,.8]
      ]

  def _trainExamples(self,ex,arch=[3,1]):
    net = Network.Network(arch)
    t = Trainers.BackProp()
    t.TrainOnLine(ex,net,errTol=self.trainTol,useAvgErr=0,
                  silent=1)
    errs = map(lambda x,y=net: abs(x[-1]-y.ClassifyExample(x)),
               ex)
    return errs

  def testBackpropOr(self):
    " testing backprop training on or "
    errs = self._trainExamples(self.orExamples)
    assert max(errs)<self.trainTol,'net did not converge properly on or'
  def testBackpropAnd(self):
    " testing backprop training on and "
    errs = self._trainExamples(self.andExamples)
    assert max(errs)<self.trainTol,'net did not converge properly on and'
  def testBackpropLin(self):
    " testing backprop training on a linear function "
    errs = self._trainExamples(self.linExamples,arch=[1,2,1])
    assert max(errs)<self.trainTol,'net did not converge properly on linear fit'
      
if __name__ == '__main__':
  unittest.main()
