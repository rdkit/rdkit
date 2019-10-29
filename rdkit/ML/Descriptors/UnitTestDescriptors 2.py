#
#  Copyright (C) 2001,2002  greg Landrum and Rational Discovery LLC
#
""" unit testing code for descriptors

"""
import unittest

from rdkit.ML.Descriptors import CompoundDescriptors
from rdkit.TestRunner import redirect_stdout
from io import StringIO


class TestCase(unittest.TestCase):

    def setUp(self):
        d = [('DED', ['NonZero', 'Mean', 'Dev']), ('M_B_electroneg', ['NonZero']),
             ('Cov_rad', ['Max', 'Min'])]
        self.desc = CompoundDescriptors.CompoundDescriptorCalculator(d)
        self.desc.BuildAtomDict()
        self.tol = 0.0001

    def testAtomDict(self):
        # " testing the atom dict "
        assert len(self.desc.atomDict.keys()) == 48, 'BuildAtomDict failed'

    def testSimpleDescriptorCalc(self):
        # " testing simple descriptor calculation "
        composList = ['Nb', 'Nb3', 'NbPt', 'Nb2Pt']
        compare = [[2.32224798203, 0.0, 1.34000003338, 1.34000003338],
                   [2.32224798203, 0.0, 1.34000003338, 1.34000003338],
                   [1.51555249095, 0.806695491076, 1.34000003338, 1.29999995232],
                   [1.78445098797, 0.717062658734, 1.34000003338, 1.29999995232]]
        for i in range(len(composList)):
            self.assertTrue(
              max(
                map(lambda x, y: abs(x - y), compare[i], self.desc.CalcSimpleDescriptorsForComposition(
                  composList[i]))) < self.tol, 'Descriptor calculation failed')

        names = self.desc.GetDescriptorNames()
        self.assertEqual(len(names), 4)
        self.assertIn('MEAN_DED', names)

    def test_exampleCode(self):
        f = StringIO()
        with redirect_stdout(f):
            CompoundDescriptors._exampleCode()


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
