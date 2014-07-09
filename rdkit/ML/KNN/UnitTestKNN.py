#
#  Copyright (C) 2003  Rational Discovery LLC
#   All Rights Reserved

""" unit testing code for knn models """
from rdkit import RDConfig
import unittest
from rdkit.ML.Data import DataUtils, MLData
from rdkit.ML.KNN import CrossValidate,DistFunctions
from rdkit.ML.KNN import KNNModel,KNNClassificationModel,KNNRegressionModel
import os.path
from rdkit.six.moves import cPickle
from rdkit import RDRandom

def feq(a,b,tol=1e-4):
  return abs(a-b)<tol

class TestCase(unittest.TestCase):
  def setUp(self):
    RDRandom.seed(25)
    
  def test1Neighbors(self):
    fName = os.path.join(RDConfig.RDCodeDir,'ML','KNN','test_data','random_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    npvals = data.GetNPossibleVals()
    nvars = data.GetNVars()
    attrs = range(1,nvars+1)
    numNeigh = 11
    metric = DistFunctions.EuclideanDist
    mdl = KNNModel.KNNModel(numNeigh,attrs,metric)
    pt = examples.pop(0)
    tgt = [(metric(pt,ex,attrs),ex) for ex in examples]
    tgt.sort()
    mdl.SetTrainingExamples(examples)
    neighbors = mdl.GetNeighbors(pt)
    for i in range(numNeigh):
      assert feq(-tgt[i][0],neighbors[i][0])
      assert tgt[i][1][0]==neighbors[i][1][0]

  def test2XValClass(self):
    fName = os.path.join(RDConfig.RDCodeDir,'ML','KNN','test_data','random_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    npvals = data.GetNPossibleVals()
    nvars = data.GetNVars()
    attrs = range(1,nvars+1)
    numNeigh = 11
    mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, numNeigh,silent=1)
    self.assertAlmostEqual(err,0.01075,4)

  def test3Regress(self):
    """ a carefully laid out regression data set where the results are clear:

    """
    fName = os.path.join(RDConfig.RDCodeDir,'ML','KNN','test_data','sample_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    npvals = data.GetNPossibleVals()
    nvars = data.GetNVars()
    attrs = range(1,nvars+1)
    numNeigh = 4
    metric = DistFunctions.EuclideanDist
    mdl = KNNRegressionModel.KNNRegressionModel(numNeigh,attrs,metric)
    mdl.SetTrainingExamples(examples)

    res = mdl.PredictExample([4,-3.5,2.5,0])
    assert feq(res,1.25)
    res = mdl.PredictExample([4,3,2,0])
    assert feq(res,1.5)
    res = mdl.PredictExample([4,3,-2.5,0])
    assert feq(res,-1.5)
    
  def test4XValRegress(self):
    fName = os.path.join(RDConfig.RDCodeDir,'ML','KNN','test_data','random_pts.csv')
    data = DataUtils.TextFileToData(fName)
    examples = data.GetNamedData()
    npvals = data.GetNPossibleVals()
    nvars = data.GetNVars()
    attrs = range(1,nvars+1)
    numNeigh = 11
    mod, err = CrossValidate.CrossValidationDriver(examples, attrs, npvals, numNeigh,silent=1,
                                                   modelBuilder=CrossValidate.makeRegressionModel)
    # NOTE: this number hasn't been extensively checked
    self.assertAlmostEqual(err,0.0777,4)


      
    
    

    
if __name__ == '__main__':
  unittest.main()
  
  
