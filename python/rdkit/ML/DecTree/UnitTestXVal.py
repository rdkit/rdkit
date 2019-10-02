#
#  Copyright (C) 2000  greg Landrum
#
""" unit testing code for cross validation """


import os
import unittest

from rdkit import RDConfig
from rdkit import RDRandom
from rdkit.ML.DecTree import CrossValidate
from rdkit.ML.DecTree import randomtest
from rdkit.TestRunner import redirect_stdout
from io import BytesIO, StringIO
import pickle


class XValTestCase(unittest.TestCase):

    def setUp(self):
        self.origTreeName = RDConfig.RDCodeDir + '/ML/DecTree/test_data/XValTree.pkl'
        self.randomSeed = 23
        self.randomArraySeed = (23, 42)

    def testRun(self):
        # " test that the CrossValidationDriver runs "
        examples, attrs, nPossibleVals = randomtest.GenRandomExamples(nExamples=200)
        f = StringIO()
        with redirect_stdout(f):
            tree, frac = CrossValidate.CrossValidationDriver(
                examples, attrs, nPossibleVals, silent=False)
        self.assertGreater(frac, 0)
        self.assertEqual('Var: 1', tree.GetName())
        self.assertIn('Validation error', f.getvalue())

        CrossValidate.CrossValidationDriver(examples, attrs, nPossibleVals, lessGreedy=True,
                                            calcTotalError=True, silent=True)

    def testResults(self):
        # " test the results of CrossValidation "
        RDRandom.seed(self.randomSeed)
        examples, attrs, nPossibleVals = randomtest.GenRandomExamples(nExamples=200,
                                                                      seed=self.randomArraySeed)
        tree, frac = CrossValidate.CrossValidationDriver(examples, attrs, nPossibleVals, silent=1)
        self.assertGreater(frac, 0)

        with open(self.origTreeName, 'r') as inTFile:
            buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
            inTFile.close()
        inFile = BytesIO(buf)
        oTree = pickle.load(inFile)

        assert oTree == tree, 'Random CrossValidation test failed'

    def testReplacementSelection(self):
        # " use selection with replacement "
        RDRandom.seed(self.randomSeed)
        examples, attrs, nPossibleVals = randomtest.GenRandomExamples(nExamples=200,
                                                                      seed=self.randomArraySeed)
        tree, frac = CrossValidate.CrossValidationDriver(examples, attrs, nPossibleVals, silent=1,
                                                         replacementSelection=1)
        self.assertTrue(tree)
        self.assertAlmostEqual(frac, 0.01666, 4)

    def test_TestRun(self):
        try:
            f = StringIO()
            with redirect_stdout(f):
                CrossValidate.TestRun()
            self.assertTrue(os.path.isfile('save.pkl'))
            s = f.getvalue()
            self.assertIn('t1 == t2 True', s)
        finally:
            if os.path.isfile('save.pkl'):
                os.remove('save.pkl')


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
