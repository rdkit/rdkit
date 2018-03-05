# $Id$
#
#  Copyright (C) 2003 Rational Discovery LLC
#      All Rights Reserved
#
""" Define the class _KNNClassificationModel_, used to represent a k-nearest neighbhors
    classification model

    Inherits from _KNNModel_
"""

from rdkit.ML.KNN import KNNModel


class KNNClassificationModel(KNNModel.KNNModel):
  """ This is used to represent a k-nearest neighbor classifier

  """

  def __init__(self, k, attrs, dfunc, radius=None):
    self._setup(k, attrs, dfunc, radius)

    self._badExamples = []  # list of examples incorrectly classified

  def type(self):
    return "Classification Model"

  def SetBadExamples(self, examples):
    self._badExamples = examples

  def GetBadExamples(self):
    return self._badExamples

  def NameModel(self, varNames):
    self.SetName(self.type())

  def ClassifyExample(self, example, appendExamples=0, neighborList=None):
    """ Classify a an example by looking at its closest neighbors

    The class assigned to this example is same as the class for the mojority of its
    _k neighbors

    **Arguments**

    - examples: the example to be classified

    - appendExamples: if this is nonzero then the example will be stored on this model

    - neighborList: if provided, will be used to return the list of neighbors

    **Returns**

      - the classification of _example_
    """
    if appendExamples:
      self._examples.append(example)

    # first find the k-closest examples in the traning set
    knnLst = self.GetNeighbors(example)

    # find out how many of the neighbors belong to each of the classes
    clsCnt = {}
    for knn in knnLst:
      cls = knn[1][-1]
      if (cls in clsCnt):
        clsCnt[cls] += 1
      else:
        clsCnt[cls] = 1
    if neighborList is not None:
      neighborList.extend(knnLst)

    # now return the class with the maximum count
    mkey = -1
    mcnt = -1
    for key in clsCnt.keys():
      if (mcnt < clsCnt[key]):
        mkey = key
        mcnt = clsCnt[key]

    return mkey
