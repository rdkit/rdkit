# $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#      All Rights Reserved
#
""" Define the class _KNNModel_, used to represent a k-nearest neighbhors model

"""
from rdkit.DataStructs.TopNContainer import TopNContainer


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

  def __init__(self, k, attrs, dfunc, radius=None):
    self._setup(k, attrs, dfunc, radius)

  def _setup(self, k, attrs, dfunc, radius):
    self._examples = []
    self._trainingExamples = []
    self._testExamples = []
    self._k = k
    self._attrs = attrs
    self._dfunc = dfunc
    self._name = ""
    self._radius = radius

  def GetName(self):
    return self._name

  def SetName(self, name):
    self._name = name

  def GetExamples(self):
    return self._examples

  def SetExamples(self, examples):
    self._examples = examples

  def GetTrainingExamples(self):
    return self._trainingExamples

  def SetTrainingExamples(self, examples):
    self._trainingExamples = examples

  def GetTestExamples(self):
    return self._testExamples

  def SetTestExamples(self, examples):
    self._testExamples = examples

  def GetNeighbors(self, example):
    """ Returns the k nearest neighbors of the example

    """
    nbrs = TopNContainer(self._k)
    for trex in self._trainingExamples:
      dist = self._dfunc(trex, example, self._attrs)
      if self._radius is None or dist < self._radius:
        nbrs.Insert(-dist, trex)
    nbrs.reverse()
    return [x for x in nbrs]
