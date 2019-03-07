#
#  Copyright (C) 2003  greg Landrum and Rational Discovery LLC
#
""" """
import os
import unittest

from rdkit.ML.DecTree import ID3, PruneTree, CrossValidate
from rdkit.TestRunner import redirect_stdout
from io import StringIO


def feq(a, b, tol=1e-4):
    return abs(a - b) <= tol


class TreeTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        # " testing pruning with known results "
        oPts = [
          [0, 0, 1, 0],
          [0, 1, 1, 1],
          [1, 0, 1, 1],
          [1, 1, 0, 0],
          [1, 1, 1, 1],
        ]
        tPts = oPts + [[0, 1, 1, 0], [0, 1, 1, 0]]
        tree = ID3.ID3Boot(oPts, attrs=range(3), nPossibleVals=[2] * 4)
        err, badEx = CrossValidate.CrossValidate(tree, oPts)
        assert err == 0.0, 'bad initial error'
        assert len(badEx) == 0, 'bad initial error'

        # prune with original data, shouldn't do anything
        f = StringIO()
        with redirect_stdout(f):
            PruneTree._verbose = True
            newTree, err = PruneTree.PruneTree(tree, [], oPts)
            PruneTree._verbose = False
        self.assertIn('Pruner', f.getvalue())
        assert newTree == tree, 'improper pruning'

        # prune with train data
        newTree, err = PruneTree.PruneTree(tree, [], tPts)
        assert newTree != tree, 'bad pruning'
        assert feq(err, 0.14286), 'bad error result'

    def test_exampleCode(self):
        f = StringIO()
        with redirect_stdout(f):
            try:
                PruneTree._testRandom()
                self.assertTrue(os.path.isfile('prune.pkl'))
            finally:
                if os.path.isfile('orig.pkl'):
                    os.remove('orig.pkl')
                if os.path.isfile('prune.pkl'):
                    os.remove('prune.pkl')
        self.assertIn('pruned error', f.getvalue())

        f = StringIO()
        with redirect_stdout(f):
            PruneTree._testSpecific()
        self.assertIn('pruned holdout error', f.getvalue())

        f = StringIO()
        with redirect_stdout(f):
            PruneTree._testChain()
        self.assertIn('pruned holdout error', f.getvalue())


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
