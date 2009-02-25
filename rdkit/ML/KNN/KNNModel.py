# $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#      All Rights Reserved
#

""" Define the class _KNNModel_, used to represent a k-nearest neighbhors model

"""
from rdkit.ML.KNN import DistFunctions
import bisect
class KNNModel(object):
  """ This is a base class used by KNNClassificationModel
  and KNNRegressionModel to represent a k-nearest neighbor predictor. In general
  one of this child classes needs to be instantiated.

  _KNNModel_s can save the following pieces of internal state, accessible via
    standard setter/getter functions - the child object store additional stuff:

    1) _Examples_: a list of examples which have been predicted (either classified
                    or values predicted)

    2) _TrainingExamples_: List of training examples (since this is a KNN model these examples
                           along with the value _k_ below define the model)

    3) _TestExamples_: the list of examples used to test the model

    4) _k_: the number of closest neighbors used for prediction

  """
  def __init__(self, k, attrs, dfunc) :
    self._setup(k, attrs, dfunc)
      
  def _setup(self, k, attrs, dfunc) :
    self._examples = []
    self._trainingExamples = []
    self._testExamples = []
    self._k = k
    self._attrs = attrs
    self._dfunc = dfunc
    self._name = ""

  def GetName(self) :
    return self_name

  def SetName(self, name) :
    self._name = name

  def GetExamples(self) :
    return self._examples

  def SetExamples(self, examples):
    self._examples = examples

  def GetTrainingExamples(self):
    return self._trainingExamples

  def SetTrainingExamples(self,examples):
    self._trainingExamples = examples

  def GetTestExamples(self) :
    return self._testExamples

  def SetTestExamples(self, examples) :
    self._testExamples = examples

  def GetNeighbors(self,example):
    """ Returns the k nearest neighbors of the example

    """
    # we'll maintain a sorted list of nearest neighbors in res:
    res = [(1e8,0)]*self._k
    for trex in self._trainingExamples:
      dist = self._dfunc(trex, example, self._attrs)
      knn = (dist, trex)
      if dist < res[-1][0]:
        loc = bisect.bisect(res,knn)
        res.insert(loc,knn)
        res.pop()
    return res
    
