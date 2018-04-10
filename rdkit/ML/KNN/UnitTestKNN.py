#
#  Copyright (C) 2003  Rational Discovery LLC
#   All Rights Reserved
""" unit testing code for knn models """
import doctest
import os.path
import unittest

from rdkit import RDConfig
from rdkit import RDRandom
from rdkit.ML.Data import DataUtils
from rdkit.ML.KNN import CrossValidate, DistFunctions
from rdkit.ML.KNN import KNNModel, KNNRegressionModel


def feq(a, b, tol=1e-4):
  return abs(a - b) < tol


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(DistFunctions, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def setUp(self):
    RDRandom.seed(25)

  def test1Neighbors(self):
    fName = os.path.join(RDConfig.RDCodeDir, 'ML', 'KNN', 'test_data', 'random_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    npvals = data.GetNPossibleVals()
    nvars = data.GetNVars()
    attrs = list(range(1, nvars + 1))
    numNeigh = 11
    metric = DistFunctions.EuclideanDist
    mdl = KNNModel.KNNModel(numNeigh, attrs, metric)
    pt = examples.pop(0)
    tgt = [(metric(pt, ex, attrs), ex) for ex in examples]
    tgt.sort()
    mdl.SetTrainingExamples(examples)
    neighbors = mdl.GetNeighbors(pt)
    for i in range(numNeigh):
      assert feq(-tgt[i][0], neighbors[i][0])
      assert tgt[i][1][0] == neighbors[i][1][0]

  def test2XValClass(self):
    fName = os.path.join(RDConfig.RDCodeDir, 'ML', 'KNN', 'test_data', 'random_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    npvals = data.GetNPossibleVals()
    nvars = data.GetNVars()
    attrs = list(range(1, nvars + 1))
    numNeigh = 11
    mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, numNeigh, silent=1)
    self.assertAlmostEqual(err, 0.01075, 4)

    neighborList = []
    res = mod.ClassifyExample(examples[0], neighborList=neighborList)
    self.assertEqual(res, 1)
    self.assertEqual(neighborList[0][1], examples[0])

    self.assertEqual(mod.GetName(), '')
    mod.SetName('name')
    self.assertEqual(mod.GetName(), 'name')
    self.assertEqual(mod.type(), 'Classification Model')
    mod.NameModel('this argument is ignored')
    self.assertEqual(mod.GetName(), 'Classification Model')

  def test3Regress(self):
    # """ a carefully laid out regression data set where the results are clear: """
    fName = os.path.join(RDConfig.RDCodeDir, 'ML', 'KNN', 'test_data', 'sample_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    nvars = data.GetNVars()
    attrs = list(range(1, nvars + 1))
    numNeigh = 4
    metric = DistFunctions.EuclideanDist
    mdl = KNNRegressionModel.KNNRegressionModel(numNeigh, attrs, metric)
    mdl.SetTrainingExamples(examples)

    res = mdl.PredictExample([4, -3.5, 2.5, 0])
    assert feq(res, 1.25)
    res = mdl.PredictExample([4, 3, 2, 0])
    assert feq(res, 1.5)
    res = mdl.PredictExample([4, 3, -2.5, 0])
    assert feq(res, -1.5)
    # Use a distance dependent weight for the neighbours
    res = mdl.PredictExample([4, 3, -2.5, 0], weightedAverage=True)
    self.assertAlmostEqual(res, -1.6)
    # Check the case that the example is identical to one of the neighbours (distance = 0)
    neighborList = []
    res = mdl.PredictExample(examples[0], weightedAverage=True, neighborList=neighborList)
    self.assertAlmostEqual(res, 1.5857864)
    self.assertEqual(neighborList[0][1], examples[0])

    self.assertEqual(mdl.GetBadExamples(), [])

    self.assertEqual(mdl.GetName(), '')
    mdl.SetName('name')
    self.assertEqual(mdl.GetName(), 'name')
    self.assertEqual(mdl.type(), 'Regression Model')
    mdl.NameModel('this argument is ignored')
    self.assertEqual(mdl.GetName(), 'Regression Model')

    self.assertEqual(sorted(mdl.GetTrainingExamples() + mdl.GetTestExamples()), sorted(examples))

  def test4XValRegress(self):
    fName = os.path.join(RDConfig.RDCodeDir, 'ML', 'KNN', 'test_data', 'random_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    npvals = data.GetNPossibleVals()
    nvars = data.GetNVars()
    attrs = list(range(1, nvars + 1))
    numNeigh = 11
    _, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, numNeigh, silent=1,
                                                 modelBuilder=CrossValidate.makeRegressionModel)
    # NOTE: this number hasn't been extensively checked
    self.assertAlmostEqual(err, 0.0777, 4)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
