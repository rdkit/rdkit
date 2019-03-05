# $Id$
#
#  Copyright (C) 2004-2005 Rational Discovery LLC.
#   All Rights Reserved
#
""" handles doing cross validation with naive bayes models
and evaluation of individual models

"""


from rdkit.ML.Data import SplitData
from rdkit.ML.NaiveBayes.ClassificationModel import NaiveBayesClassifier

try:
  from rdkit.ML.FeatureSelect import CMIM
except ImportError:
  CMIM = None


def makeNBClassificationModel(trainExamples, attrs, nPossibleValues, nQuantBounds,
                              mEstimateVal=-1.0, useSigs=False, ensemble=None, useCMIM=0, **kwargs):
  if CMIM is not None and useCMIM > 0 and useSigs and not ensemble:
    ensemble = CMIM.SelectFeatures(trainExamples, useCMIM, bvCol=1)
  if ensemble:
    attrs = ensemble
  model = NaiveBayesClassifier(attrs, nPossibleValues, nQuantBounds, mEstimateVal=mEstimateVal,
                               useSigs=useSigs)

  model.SetTrainingExamples(trainExamples)
  model.trainModel()
  return model


def CrossValidate(NBmodel, testExamples, appendExamples=0):

  nTest = len(testExamples)
  assert nTest, 'no test examples: %s' % str(testExamples)
  badExamples = []
  nBad = 0
  preds = NBmodel.ClassifyExamples(testExamples, appendExamples)
  assert len(preds) == nTest

  for i in range(nTest):
    testEg = testExamples[i]
    trueRes = testEg[-1]
    res = preds[i]

    if (trueRes != res):
      badExamples.append(testEg)
      nBad += 1
  return float(nBad) / nTest, badExamples


def CrossValidationDriver(examples, attrs, nPossibleValues, nQuantBounds, mEstimateVal=0.0,
                          holdOutFrac=0.3, modelBuilder=makeNBClassificationModel, silent=0,
                          calcTotalError=0, **kwargs):
  nTot = len(examples)
  if not kwargs.get('replacementSelection', 0):
    testIndices, trainIndices = SplitData.SplitIndices(nTot, holdOutFrac, silent=1, legacy=1,
                                                       replacement=0)
  else:
    testIndices, trainIndices = SplitData.SplitIndices(nTot, holdOutFrac, silent=1, legacy=0,
                                                       replacement=1)

  trainExamples = [examples[x] for x in trainIndices]
  testExamples = [examples[x] for x in testIndices]

  NBmodel = modelBuilder(trainExamples, attrs, nPossibleValues, nQuantBounds, mEstimateVal,
                         **kwargs)

  if not calcTotalError:
    xValError, _ = CrossValidate(NBmodel, testExamples, appendExamples=1)
  else:
    xValError, _ = CrossValidate(NBmodel, examples, appendExamples=0)

  if not silent:
    print('Validation error was %%%4.2f' % (100 * xValError))
  NBmodel._trainIndices = trainIndices
  return NBmodel, xValError
