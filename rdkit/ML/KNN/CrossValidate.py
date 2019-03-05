#
#  Copyright (C) 2000-2008  greg Landrum
#
""" handles doing cross validation with k-nearest neighbors model

and evaluation of individual models

"""


from rdkit.ML.Data import SplitData
from rdkit.ML.KNN import DistFunctions
from rdkit.ML.KNN.KNNClassificationModel import KNNClassificationModel
from rdkit.ML.KNN.KNNRegressionModel import KNNRegressionModel


def makeClassificationModel(numNeigh, attrs, distFunc):
  return KNNClassificationModel(numNeigh, attrs, distFunc)


def makeRegressionModel(numNeigh, attrs, distFunc):
  return KNNRegressionModel(numNeigh, attrs, distFunc)


def CrossValidate(knnMod, testExamples, appendExamples=0):
  """
  Determines the classification error for the testExamples

  **Arguments**

    - tree: a decision tree (or anything supporting a _ClassifyExample()_ method)

    - testExamples: a list of examples to be used for testing

    - appendExamples: a toggle which is passed along to the tree as it does
      the classification. The trees can use this to store the examples they
      classify locally.

  **Returns**

    a 2-tuple consisting of:
      """
  nTest = len(testExamples)

  if isinstance(knnMod, KNNClassificationModel):
    badExamples = []
    nBad = 0
    for i in range(nTest):
      testEx = testExamples[i]
      trueRes = testEx[-1]
      res = knnMod.ClassifyExample(testEx, appendExamples)
      if (trueRes != res):
        badExamples.append(testEx)
        nBad += 1
    return float(nBad) / nTest, badExamples
  elif isinstance(knnMod, KNNRegressionModel):
    devSum = 0.0
    for i in range(nTest):
      testEx = testExamples[i]
      trueRes = testEx[-1]
      res = knnMod.PredictExample(testEx, appendExamples)
      devSum += abs(trueRes - res)
    return devSum / nTest, None
  raise ValueError("Unrecognized Model Type")


def CrossValidationDriver(examples, attrs, nPossibleValues, numNeigh,
                          modelBuilder=makeClassificationModel,
                          distFunc=DistFunctions.EuclideanDist, holdOutFrac=0.3, silent=0,
                          calcTotalError=0, **kwargs):
  """ Driver function for building a KNN model of a specified type

  **Arguments**

    - examples: the full set of examples

    - numNeigh: number of neighbors for the KNN model (basically k in k-NN)

    - knnModel: the type of KNN model (a classification vs regression model)

    - holdOutFrac: the fraction of the data which should be reserved for the hold-out set
      (used to calculate error)

    - silent: a toggle used to control how much visual noise this makes as it goes

    - calcTotalError: a toggle used to indicate whether the classification error
      of the tree should be calculated using the entire data set (when true) or just
      the training hold out set (when false)
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
    print("Training with %d examples" % (nTrain))

  knnMod = modelBuilder(numNeigh, attrs, distFunc)

  knnMod.SetTrainingExamples(trainExamples)
  knnMod.SetTestExamples(testExamples)

  if not calcTotalError:
    xValError, _ = CrossValidate(knnMod, testExamples, appendExamples=1)
  else:
    xValError, _ = CrossValidate(knnMod, examples, appendExamples=0)

  if not silent:
    print('Validation error was %%%4.2f' % (100 * xValError))

  knnMod._trainIndices = trainIndices
  return knnMod, xValError
