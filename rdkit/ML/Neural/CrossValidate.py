#
#  Copyright (C) 2000  greg Landrum
#
""" handles doing cross validation with neural nets

This is, perhaps, a little misleading.  For the purposes of this module,
cross validation == evaluating the accuracy of a net.

"""

from rdkit.ML.Neural import Network, Trainers
from rdkit.ML.Data import SplitData
import math


def CrossValidate(net, testExamples, tolerance, appendExamples=0):
  """ Determines the classification error for the testExamples
    **Arguments**

      - tree: a decision tree (or anything supporting a _ClassifyExample()_ method)

      - testExamples: a list of examples to be used for testing

      - appendExamples: a toggle which is ignored, it's just here to maintain
         the same API as the decision tree code.

    **Returns**

      a 2-tuple consisting of:

        1) the percent error of the net

        2) a list of misclassified examples

   **Note**
     At the moment, this is specific to nets with only one output
  """
  nTest = len(testExamples)
  nBad = 0
  badExamples = []
  for i in range(nTest):
    testEx = testExamples[i]
    trueRes = testExamples[i][-1]
    res = net.ClassifyExample(testEx)
    if math.fabs(trueRes - res) > tolerance:
      badExamples.append(testEx)
      nBad = nBad + 1

  return float(nBad) / nTest, badExamples


def CrossValidationDriver(examples, attrs=[], nPossibleVals=[], holdOutFrac=.3, silent=0,
                          tolerance=0.3, calcTotalError=0, hiddenSizes=None, **kwargs):
  """
    **Arguments**

      - examples: the full set of examples

      - attrs: a list of attributes to consider in the tree building
         *This argument is ignored*

      - nPossibleVals: a list of the number of possible values each variable can adopt
         *This argument is ignored*

      - holdOutFrac: the fraction of the data which should be reserved for the hold-out set
         (used to calculate the error)

      - silent: a toggle used to control how much visual noise this makes as it goes.

      - tolerance: the tolerance for convergence of the net

      - calcTotalError: if this is true the entire data set is used to calculate
           accuracy of the net

      - hiddenSizes: a list containing the size(s) of the hidden layers in the network.
           if _hiddenSizes_ is None, one hidden layer containing the same number of nodes
           as the input layer will be used

    **Returns**

       a 2-tuple containing:

         1) the net

         2) the cross-validation error of the net

    **Note**
      At the moment, this is specific to nets with only one output

  """
  nTot = len(examples)
  if not kwargs.get('replacementSelection', 0):
    testIndices, trainIndices = SplitData.SplitIndices(nTot, holdOutFrac, silent=1, legacy=1,
                                                       replacement=0)
  else:
    testIndices, trainIndices = SplitData.SplitIndices(nTot, holdOutFrac, silent=1, legacy=0,
                                                       replacement=1)
  trainExamples = [examples[x] for x in trainIndices]
  testExamples = [examples[x] for x in testIndices]

  nTrain = len(trainExamples)
  if not silent:
    print('Training with %d examples' % (nTrain))

  nInput = len(examples[0]) - 1
  nOutput = 1
  if hiddenSizes is None:
    nHidden = nInput
    netSize = [nInput, nHidden, nOutput]
  else:
    netSize = [nInput] + hiddenSizes + [nOutput]
  net = Network.Network(netSize)
  t = Trainers.BackProp()
  t.TrainOnLine(trainExamples, net, errTol=tolerance, useAvgErr=0, silent=silent)

  nTest = len(testExamples)
  if not silent:
    print('Testing with %d examples' % nTest)
  if not calcTotalError:
    xValError, _ = CrossValidate(net, testExamples, tolerance)
  else:
    xValError, _ = CrossValidate(net, examples, tolerance)
  if not silent:
    print('Validation error was %%%4.2f' % (100 * xValError))
  net._trainIndices = trainIndices
  return net, xValError
