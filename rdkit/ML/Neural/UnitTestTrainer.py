#
#  Copyright (C) 2000  greg Landrum
#
""" unit tests for the Neural network trainer implementation

   this basically works out **all** of the network code

"""


import random
import unittest

from rdkit.ML.Neural import Network, Trainers
from rdkit.ML.Neural.CrossValidate import CrossValidate, CrossValidationDriver
from rdkit.TestRunner import redirect_stdout
from io import StringIO


class TrainerTestCase(unittest.TestCase):

    def setUp(self):
        random.seed(23)
        self.trainTol = 0.3
        self.orExamples = [[0, 0, 1, 0.1], [0, 1, 1, .9], [1, 0, 1, .9], [1, 1, 1, .9]]
        self.andExamples = [[0, 0, 1, 0.1], [0, 1, 1, .1], [1, 0, 1, .1], [1, 1, 1, .9]]
        self.xorExamples = [[0, 0, 1, .1], [0, 1, 1, .9], [1, 0, 1, .9], [1, 1, 1, .1]]
        self.linExamples = [[.1, .1], [.2, .2], [.3, .3], [.4, .4], [.8, .8]]

    def _trainExamples(self, ex, arch=[3, 1], useAvgErr=False):
        net = Network.Network(arch)
        t = Trainers.BackProp()
        t.TrainOnLine(ex, net, errTol=self.trainTol, useAvgErr=useAvgErr, silent=True)
        errs = [abs(x[-1] - net.ClassifyExample(x)) for x in ex]
        return net, errs

    def testBackpropOr(self):
        # " testing backprop training on or "
        _, errs = self._trainExamples(self.orExamples)
        assert max(errs) < self.trainTol, 'net did not converge properly on or'

    def testBackpropAnd(self):
        # " testing backprop training on and "
        _, errs = self._trainExamples(self.andExamples)
        assert max(errs) < self.trainTol, 'net did not converge properly on and'

    def testBackpropLin(self):
        # " testing backprop training on a linear function "
        _, errs = self._trainExamples(self.linExamples, arch=[1, 2, 1])
        assert max(errs) < self.trainTol, 'net did not converge properly on linear fit'

        _, errs = self._trainExamples(self.linExamples, arch=[1, 2, 1], useAvgErr=True)
        assert max(errs) < 0.4, 'net did not converge properly on or'

    def test_multipleHiddenLayers(self):
        _, errs = self._trainExamples(self.linExamples, arch=[1, 1, 2, 1])
        assert max(errs) < self.trainTol, 'net did not converge properly on linear fit'

    def test_CrossValidate(self):
        # We just check here that the code works
        net, _ = self._trainExamples(self.orExamples)
        percentage, badExamples = CrossValidate(net, self.orExamples, 0.2)
        self.assertEqual(percentage, 1.0 / 4)
        self.assertEqual(len(badExamples), 1)

        percentage, badExamples = CrossValidate(net, self.orExamples, self.trainTol)
        self.assertEqual(percentage, 0.0)
        self.assertEqual(len(badExamples), 0)

        net, cvError = CrossValidationDriver(self.orExamples + self.orExamples, silent=True)
        self.assertEqual(cvError, 0.5)

        net, cvError = CrossValidationDriver(self.orExamples + self.orExamples, silent=True,
                                             replacementSelection=True)
        self.assertEqual(cvError, 0.0)

        net, cvError = CrossValidationDriver(self.orExamples + self.orExamples, silent=True,
                                             calcTotalError=True)
        self.assertEqual(cvError, 0.25)

        f = StringIO()
        with redirect_stdout(f):
            CrossValidationDriver(self.orExamples + self.orExamples)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
