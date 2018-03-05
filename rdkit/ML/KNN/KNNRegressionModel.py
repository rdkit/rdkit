# $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#      All Rights Reserved
#
""" Define the class _KNNRegressionModel_, used to represent a k-nearest neighbhors
regression model

    Inherits from _KNNModel_
"""

from rdkit.ML.KNN import KNNModel


class KNNRegressionModel(KNNModel.KNNModel):
  """ This is used to represent a k-nearest neighbor classifier

  """

  def __init__(self, k, attrs, dfunc, radius=None):
    self._setup(k, attrs, dfunc, radius)

    self._badExamples = []  # list of examples incorrectly classified

  def type(self):
    return "Regression Model"

  def SetBadExamples(self, examples):
    self._badExamples = examples

  def GetBadExamples(self):
    return self._badExamples

  def NameModel(self, varNames):
    self.SetName(self.type())

  def PredictExample(self, example, appendExamples=0, weightedAverage=0, neighborList=None):
    """ Generates a prediction for an example by looking at its closest neighbors

    **Arguments**

      - examples: the example to be classified

      - appendExamples: if this is nonzero then the example will be stored on this model

      - weightedAverage: if provided, the neighbors' contributions to the value will be
                         weighed by their reciprocal square distance

      - neighborList: if provided, will be used to return the list of neighbors

    **Returns**

      - the classification of _example_

    """
    if appendExamples:
      self._examples.append(example)

    # first find the k-closest examples in the training set
    knnLst = self.GetNeighbors(example)

    accum = 0.0
    denom = 0.0
    for knn in knnLst:
      if knn[1] is None:
        continue
      if weightedAverage:
        dist = knn[0]
        if dist == 0.0:
          w = 1.
        else:
          w = 1. / dist
      else:
        w = 1.0
      accum += w * knn[1][-1]
      denom += w
    if denom:
      accum /= denom
    if neighborList is not None:
      neighborList.extend(knnLst)
    return accum
