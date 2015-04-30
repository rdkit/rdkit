#
#  Copyright (C) 2000  greg Landrum
#

""" unit testing code for cross validation """
from __future__ import print_function
import unittest, random
import io
from rdkit import RDConfig
from rdkit.ML.DecTree import CrossValidate

def feq(a,b,tol=1e-4):
  return abs(a-b)<tol

class XValTestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
    self.origTreeName = RDConfig.RDCodeDir+'/ML/DecTree/test_data/XValTree.pkl'
    self.randomSeed = 23
    self.randomArraySeed = (23,42)

  def testRun(self):
    " test that the CrossValidationDriver runs "
    from rdkit.ML.DecTree import randomtest
    examples,attrs,nPossibleVals = randomtest.GenRandomExamples(nExamples = 200)
    try:
      tree,frac = CrossValidate.CrossValidationDriver(examples,attrs,
                                                      nPossibleVals,silent=1)
    except:
      assert 0,'CrossValidation failed to run'

  def testResults(self):
    " test the results of CrossValidation "
    from rdkit.ML.DecTree import randomtest
    from rdkit import RDRandom
    RDRandom.seed(self.randomSeed)
    examples,attrs,nPossibleVals = randomtest.GenRandomExamples(nExamples = 200,
                                                                seed=self.randomArraySeed)
    tree,frac = CrossValidate.CrossValidationDriver(examples,attrs,
                                                    nPossibleVals,silent=1)

    from rdkit.six.moves import cPickle
    #cPickle.dump(tree,open(self.origTreeName,'w+'))
    with open(self.origTreeName,'r') as inTFile:
      buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
      inTFile.close()
    with io.BytesIO(buf) as inFile:
      oTree = cPickle.load(inFile)

    assert oTree==tree,'Random CrossValidation test failed'
    
  def testReplacementSelection(self):
    " use selection with replacement "
    from rdkit.ML.DecTree import randomtest
    from rdkit import RDRandom
    RDRandom.seed(self.randomSeed)
    examples,attrs,nPossibleVals = randomtest.GenRandomExamples(nExamples = 200,
                                                                seed=self.randomArraySeed)
    tree,frac = CrossValidate.CrossValidationDriver(examples,attrs,
                                                    nPossibleVals,silent=1,
                                                    replacementSelection=1)
    self.assertTrue(tree)
    self.assertAlmostEqual(frac,0.01666,4)
    
    


if __name__ == '__main__':
  unittest.main()
  
  
