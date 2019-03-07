#
#  Copyright (C) 2001  greg Landrum
#
""" unit testing code for MLData sets

"""
import contextlib
import random
import unittest

from rdkit import RDConfig
from rdkit.ML.Data import MLData, DataUtils
from io import StringIO
import pickle


class TestCase(unittest.TestCase):

    def setUpQuantLoad(self):
        self.d = DataUtils.BuildQuantDataSet(RDConfig.RDCodeDir + '/ML/Data/test_data/test.qdat')

    def setUpGeneralLoad(self):
        self.d = DataUtils.BuildDataSet(RDConfig.RDCodeDir + '/ML/Data/test_data/test.dat')

    def testQuantLoad(self):
        # " testing QuantDataSet load"
        self.d = DataUtils.BuildQuantDataSet(RDConfig.RDCodeDir + '/ML/Data/test_data/test.qdat')

    def testQuantProps(self):
        # " testing QuantDataSet Properties"
        self.setUpQuantLoad()
        assert self.d.GetNPts() == 5, 'nPts wrong'
        assert self.d.GetNVars() == 4, 'nVars wrong'
        assert self.d.GetNResults() == 1, 'nResults wrong'
        assert self.d.GetVarNames() == ['foo1', 'foo2', 'foo3', 'foo4', 'res'], 'varNames wrong'
        assert self.d.GetPtNames() == ['p1', 'p2', 'p3', 'p4', 'p5'], 'ptNames wrong'
        assert self.d.GetNPossibleVals() == [3, 3, 2, 2, 2], 'nPossible Wrong'
        assert self.d.GetQuantBounds() == [[], [], [], [], []], 'quantBounds Wrong'
        assert self.d.GetResults() == [0, 0, 1, 1, 1], 'GetResults wrong'
        assert self.d.GetAllData()[1] == [0, 0, 0, 1, 0], 'GetAllData wrong'
        assert self.d.GetInputData()[3] == [2, 1, 0, 0], 'GetInputData wrong'
        assert self.d.GetNamedData()[2] == ['p3', 1, 0, 0, 0, 1], 'GetNamedData wrong'

    def testQuantPickle(self):
        # " testing QuantDataSet pickling "
        self.setUpQuantLoad()
        DataUtils.WritePickledData(
            RDConfig.RDCodeDir + '/ML/Data/test_data/testquant.qdat.pkl', self.d)
        with open(RDConfig.RDCodeDir + '/ML/Data/test_data/testquant.qdat.pkl', 'rb') as f:
            vNames = pickle.load(f)
            qBounds = pickle.load(f)
            ptNames = pickle.load(f)
            examples = pickle.load(f)
        d = MLData.MLQuantDataSet(examples, varNames=vNames, qBounds=qBounds, ptNames=ptNames)
        assert self.d.GetNPts() == d.GetNPts(), 'nPts wrong'
        assert self.d.GetNVars() == d.GetNVars(), 'nVars wrong'
        assert self.d.GetNResults() == d.GetNResults(), 'nResults wrong'
        assert self.d.GetVarNames() == d.GetVarNames(), 'varNames wrong'
        assert self.d.GetPtNames() == d.GetPtNames(), 'ptNames wrong'
        assert self.d.GetNPossibleVals() == d.GetNPossibleVals(), 'nPossible Wrong'
        assert self.d.GetQuantBounds() == d.GetQuantBounds(), 'quantBounds Wrong'
        assert self.d.GetResults() == d.GetResults(), 'GetResults wrong'
        assert self.d.GetAllData()[1] == d.GetAllData()[1], 'GetAllData wrong'
        assert self.d.GetInputData()[3] == d.GetInputData()[3], 'GetInputData wrong'
        assert self.d.GetNamedData()[2] == d.GetNamedData()[2], 'GetNamedData wrong'

    def testGeneralLoad(self):
        # " testing DataSet load"
        self.d = DataUtils.BuildDataSet(RDConfig.RDCodeDir + '/ML/Data/test_data/test.dat')

    def testGeneralProps(self):
        # " testing DataSet properties"
        self.setUpGeneralLoad()
        assert self.d.GetNPts() == 5, 'nPts wrong'
        assert self.d.GetNVars() == 4, 'nVars wrong'
        assert self.d.GetNResults() == 1, 'nResults wrong'
        assert self.d.GetVarNames() == ['foo1', 'foo2', 'foo3', 'foo4', 'res'], 'varNames wrong'
        assert self.d.GetPtNames() == ['p1', 'p2', 'p3', 'p4', 'p5'], 'ptNames wrong'
        assert self.d.GetNPossibleVals() == [0, 6, 0, 0, 3], 'nPossible Wrong'
        assert self.d.GetQuantBounds() == [[], [], [], [], [2, 4]], 'quantBounds Wrong'
        assert self.d.GetResults() == [1.1, 2.1, 3.1, 4.1, 5.1], 'GetResults wrong'
        assert self.d.GetAllData()[1] == ['foo', 2, 1.0, 1, 2.1], 'GetAllData wrong'
        assert self.d.GetInputData()[4] == ['foo', 5, 1.1, 1], 'GetInputData wrong'
        assert self.d.GetNamedData()[3] == ['p4', 'foo', 4, 1.0, 1, 4.1], 'GetNamedData wrong'
        assert self.d[1] == ['p2', 'foo', 2, 1.0, 1, 2.1]
        assert self.d[3] == ['p4', 'foo', 4, 1.0, 1, 4.1]
        self.d[3] = self.d[1]
        assert self.d[3] == ['p2', 'foo', 2, 1.0, 1, 2.1]

    def testGeneralPickle(self):
        # " testing DataSet pickling"
        self.setUpGeneralLoad()
        DataUtils.WritePickledData(RDConfig.RDCodeDir + '/ML/Data/test_data/testgeneral.dat.pkl',
                                   self.d)
        with open(RDConfig.RDCodeDir + '/ML/Data/test_data/testgeneral.dat.pkl', 'rb') as f:
            vNames = pickle.load(f)
            qBounds = pickle.load(f)
            ptNames = pickle.load(f)
            examples = pickle.load(f)
        d = MLData.MLDataSet(examples, varNames=vNames, qBounds=qBounds, ptNames=ptNames)
        assert self.d.GetNPts() == d.GetNPts(), 'nPts wrong'
        assert self.d.GetNVars() == d.GetNVars(), 'nVars wrong'
        assert self.d.GetNResults() == d.GetNResults(), 'nResults wrong'
        assert self.d.GetVarNames() == d.GetVarNames(), 'varNames wrong'
        assert self.d.GetPtNames() == d.GetPtNames(), 'ptNames wrong'
        assert self.d.GetNPossibleVals() == d.GetNPossibleVals(), 'nPossible Wrong'
        assert self.d.GetQuantBounds() == d.GetQuantBounds(), 'quantBounds Wrong'
        assert self.d.GetResults() == d.GetResults(), 'GetResults wrong'
        assert self.d.GetAllData()[1] == d.GetAllData()[1], 'GetAllData wrong'
        assert self.d.GetInputData()[3] == d.GetInputData()[3], 'GetInputData wrong'
        assert self.d.GetNamedData()[2] == d.GetNamedData()[2], 'GetNamedData wrong'

    def test_WriteData(self):
        self.setUpQuantLoad()
        with contextlib.closing(StringIO()) as f:
            DataUtils.WriteData(f, self.d.GetVarNames(), self.d.GetQuantBounds(), self.d.data)
            s = f.getvalue()
            self.assertIn('DataUtils', s)
            self.assertIn('foo1', s)
            self.assertIn('2 2 1 0 1', s)

    def test_CalcNPossibleUsingMap(self):
        self.setUpQuantLoad()
        order = list(range(5))
        self.assertEqual(
          DataUtils.CalcNPossibleUsingMap(self.d.data, order, self.d.GetQuantBounds()), [3, 3, 2, 2, 2])

    def test_RandomizeActivities(self):

        class RunDetails(object):
            shuffled = False
            randomized = False

        random.seed(0)
        details = RunDetails()
        self.setUpGeneralLoad()
        dataSet = self.d
        orgActivities = [d[-1] for d in dataSet]
        DataUtils.RandomizeActivities(dataSet, shuffle=True, runDetails=details)
        self.assertNotEqual(orgActivities, [d[-1] for d in dataSet])
        self.assertEqual(sorted(orgActivities), sorted([d[-1] for d in dataSet]))
        self.assertTrue(details.shuffled)
        self.assertFalse(details.randomized)

        try:
            details = RunDetails()
            self.setUpGeneralLoad()
            dataSet = self.d
            orgActivities = [d[-1] for d in dataSet]
            DataUtils.RandomizeActivities(dataSet, shuffle=False, runDetails=details)
            self.assertNotEqual(orgActivities, [d[-1] for d in dataSet])
            self.assertEqual(sorted(orgActivities), sorted([d[-1] for d in dataSet]))
            self.assertFalse(details.randomized)
            self.assertTrue(details.shuffled)
        except NameError:
            # This code branch is not working.
            pass


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
