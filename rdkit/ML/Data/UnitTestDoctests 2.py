#
#  Copyright (C) 2001-2008  greg Landrum
#
""" unit testing code for the Doctests

"""
import doctest
import unittest

from rdkit.ML.Data import SplitData, DataUtils
from rdkit.TestRunner import redirect_stdout
from io import StringIO


def load_tests(loader, tests, ignore):
    """ Add the Doctests from the module """
    tests.addTests(doctest.DocTestSuite(SplitData, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(DataUtils, optionflags=doctest.ELLIPSIS))
    return tests


class TestCaseSplitData(unittest.TestCase):

    def test_exceptions(self):
        self.assertRaises(ValueError, SplitData.SplitIndices, 10, -0.1)
        self.assertRaises(ValueError, SplitData.SplitIndices, 10, 1.1)

        f = StringIO()
        with redirect_stdout(f):
            SplitData.SplitIndices(10, 0.5, replacement=True, silent=False)
        s = f.getvalue()
        self.assertIn('Training', s)
        self.assertIn('hold-out', s)

    def test_SplitData(self):
        self.assertRaises(ValueError, SplitData.SplitDataSet, None, -1.1)
        self.assertRaises(ValueError, SplitData.SplitDataSet, None, 1.1)

        data = list(range(10))
        DataUtils.InitRandomNumbers((23, 42))
        f = StringIO()
        with redirect_stdout(f):
            result = SplitData.SplitDataSet(data, 0.5)
        self.assertEqual(set(result[0]).intersection(result[1]), set())
        self.assertEqual(len(result[0]), 5)
        s = f.getvalue()
        self.assertIn('Training', s)
        self.assertIn('hold-out', s)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
