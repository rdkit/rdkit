#
#  Copyright (C) 2000  greg Landrum
#

""" unit testing code for cross validation """
import RDConfig
import unittest
from ML.DecTree import CrossValidate
import random

def feq(a,b,tol=1e-4):
  return abs(a-b)<tol

class XValTestCase(unittest.TestCase):
  def setUp(self):
    print '\n%s: '%self.shortDescription(),
    self.origTreeName = RDConfig.RDCodeDir+'/ML/DecTree/regress/XValTree.pkl'
    self.randomSeed = 23
    self.randomArraySeed = (23,42)

  def testRun(self):
    " test that the CrossValidationDriver runs "
    from ML.DecTree import randomtest
    examples,attrs,nPossibleVals = randomtest.GenRandomExamples(nExamples = 200)
    try:
      tree,frac = CrossValidate.CrossValidationDriver(examples,attrs,
                                                      nPossibleVals,silent=1)
    except:
      assert 0,'CrossValidation failed to run'

  def testResults(self):
    " test the results of CrossValidation "
    from ML.DecTree import randomtest
    import RDRandom
    RDRandom.seed(self.randomSeed)
    examples,attrs,nPossibleVals = randomtest.GenRandomExamples(nExamples = 200,
                                                                seed=self.randomArraySeed)
    tree,frac = CrossValidate.CrossValidationDriver(examples,attrs,
                                                    nPossibleVals,silent=1)
    import cPickle
    inFile = open(self.origTreeName,'r')
    oTree = cPickle.load(inFile)

    assert oTree==tree,'Random CrossValidation test failed'
    
  def testReplacementSelection(self):
    " use selection with replacement "
    from ML.DecTree import randomtest
    import RDRandom
    RDRandom.seed(self.randomSeed)
    examples,attrs,nPossibleVals = randomtest.GenRandomExamples(nExamples = 200,
                                                                seed=self.randomArraySeed)
    tree,frac = CrossValidate.CrossValidationDriver(examples,attrs,
                                                    nPossibleVals,silent=1,
                                                    replacementSelection=1)
    assert tree
    assert feq(frac,0.0333)
    
    


if __name__ == '__main__':
  unittest.main()
  
  
