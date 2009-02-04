# $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#      All Rights Reserved
#

""" Define the class _KNNRegressionModel_, used to represent a k-nearest neighbhors
regression model

    Inherits from _KNNModel_
"""

from ML.KNN import KNNModel
from ML.KNN import DistFunctions
import math


class KNNRegressionModel(KNNModel.KNNModel) :
  """ This is used to represent a k-nearest neighbor classifier

  """

  def __init__(self, k, attrs, dfunc) :
    self._setup(k, attrs, dfunc)

    self._badExamples = [] # list of examples incorrectly classified

  def type(self):
    return "Regression Model"

  def SetBadExamples(self, examples) :
    self._badExamples = examples

  def GetBadExamples(self) :
    return self._badExamples

  def NameModel(self, varNames) :
    self.SetName(self.type())

  def PredictExample(self, example, appendExamples=0, weightedAverage=0) :
    """ Generates a prediction for an example by looking at its closest neighbors

    **Arguments**

      - examples: the example to be classified

      - appendExamples: if this is nonzero then the example will be stored on this model

    **Returns**

      - the classification of _example_

    """
    if appendExamples:
      self._examples.append(example)

    # first find the k-closest examples in the traning set
    knnLst = self.GetNeighbors(example)

    accum = 0.0
    for knn in knnLst:
      if not weightedAverage:
        accum += knn[1][-1]
    if not weightedAverage:
      accum /= len(knnLst)

    return accum

