# $Id$
#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
""" code for dealing with Bayesian composite models

For a model to be useable here, it should support the following API:

  - _ClassifyExample(example)_, returns a classification

Other compatibility notes:

 1) To use _Composite.Grow_ there must be some kind of builder
    functionality which returns a 2-tuple containing (model,percent accuracy).

 2) The models should be pickleable

 3) It would be very happy if the models support the __cmp__ method so that
    membership tests used to make sure models are unique work.



"""


import numpy

from rdkit.ML.Composite import Composite


class BayesComposite(Composite.Composite):
  """a composite model using Bayesian statistics in the Decision Proxy


    **Notes**

    - typical usage:

       1) grow the composite with AddModel until happy with it

       2) call AverageErrors to calculate the average error values

       3) call SortModels to put things in order by either error or count

       4) call Train to update the Bayesian stats.

  """

  def Train(self, data, verbose=0):
    # FIX: this is wrong because it doesn't take the counts of each model into account
    nModels = len(self)
    nResults = self.nPossibleVals[-1]
    self.resultProbs = numpy.zeros(nResults, numpy.float)
    self.condProbs = [None] * nModels

    for i in range(nModels):
      self.condProbs[i] = numpy.zeros((nResults, nResults), numpy.float)
    # FIX: this is a quick hack which may slow things down a lot
    for example in data:
      act = self.QuantizeActivity(example)[-1]
      self.resultProbs[int(act)] += 1

    for example in data:
      if self._mapOrder is not None:
        example = self._RemapInput(example)
      if self.GetActivityQuantBounds():
        example = self.QuantizeActivity(example)
      if self.quantBounds is not None and 1 in self.quantizationRequirements:
        quantExample = self.QuantizeExample(example, self.quantBounds)
      else:
        quantExample = []

      trueRes = int(example[-1])

      votes = self.CollectVotes(example, quantExample)

      for i in range(nModels):
        self.condProbs[i][votes[i], trueRes] += 1

      # self.condProbs /= self.resultProbs
    for i in range(nModels):
      for j in range(nResults):
        self.condProbs[i][j] /= sum(self.condProbs[i][j])
      # self.condProbs[i] /= self.resultProbs

    self.resultProbs /= sum(self.resultProbs)

    if verbose:
      print('**** Bayesian Results')
      print('Result probabilities')
      print('\t', self.resultProbs)
      print('Model by model breakdown of conditional probs')
      for mat in self.condProbs:
        for row in mat:
          print('\t', row)
        print()

  def ClassifyExample(self, example, threshold=0, verbose=0, appendExample=0):
    """ classifies the given example using the entire composite

      **Arguments**

       - example: the data to be classified

       - threshold:  if this is a number greater than zero, then a
          classification will only be returned if the confidence is
          above _threshold_.  Anything lower is returned as -1.

      **Returns**

        a (result,confidence) tuple

    """
    if self._mapOrder is not None:
      example = self._RemapInput(example)
    if self.GetActivityQuantBounds():
      example = self.QuantizeActivity(example)
    if self.quantBounds is not None and 1 in self.quantizationRequirements:
      quantExample = self.QuantizeExample(example, self.quantBounds)
    else:
      quantExample = []
    self.modelVotes = self.CollectVotes(example, quantExample, appendExample=appendExample)

    nPossibleRes = self.nPossibleVals[-1]
    votes = [0.] * nPossibleRes
    for i in range(len(self)):
      predict = self.modelVotes[i]
      for j in range(nPossibleRes):
        votes[j] += self.condProbs[i][predict, j]

    # totVotes = sum(votes)
    res = numpy.argmax(votes)
    conf = votes[res] / len(self)
    if verbose:
      print(votes, conf, example[-1])
    if conf > threshold:
      return res, conf
    else:
      return -1, conf

  def __init__(self):
    Composite.Composite.__init__(self)
    self.resultProbs = None
    self.condProbs = None


def CompositeToBayesComposite(obj):
  """ converts a Composite to a BayesComposite

   if _obj_ is already a BayesComposite or if it is not a _Composite.Composite_ ,
    nothing will be done.

  """
  if obj.__class__ == BayesComposite:
    return
  elif obj.__class__ == Composite.Composite:
    obj.__class__ = BayesComposite
    obj.resultProbs = None
    obj.condProbs = None


def BayesCompositeToComposite(obj):
  """ converts a BayesComposite to a Composite.Composite

  """
  if obj.__class__ == Composite.Composite:
    return
  elif obj.__class__ == BayesComposite:
    obj.__class__ = Composite.Composite
    obj.resultProbs = None
    obj.condProbs = None
