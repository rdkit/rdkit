#
#  Copyright (C) 2001,2003  greg Landrum and Rational Discovery LLC
#
""" unit tests for the QuantTree implementation """


import io
import unittest

from rdkit import RDConfig
from rdkit.ML.DecTree import BuildQuantTree
from rdkit.ML.DecTree.QuantTree import QuantTreeNode
from rdkit.TestRunner import redirect_stdout
from io import StringIO
import pickle


def cmp(t1, t2):
    return (t1 < t2) * -1 or (t1 > t2) * 1


class TestCase(unittest.TestCase):

    def setUp(self):
        self.qTree1Name = RDConfig.RDCodeDir + '/ML/DecTree/test_data/QuantTree1.pkl'
        self.qTree2Name = RDConfig.RDCodeDir + '/ML/DecTree/test_data/QuantTree2.pkl'

    def _setupTree1(self):
        examples1 = [['p1', 0, 1, 0.1, 0], ['p2', 0, 0, 0.1, 1], ['p3', 0, 0, 1.1, 2],
                     ['p4', 0, 1, 1.1, 2], ['p5', 1, 0, 0.1, 2], ['p6', 1, 0, 1.1, 2],
                     ['p7', 1, 1, 0.1, 2], ['p8', 1, 1, 1.1, 0]]
        attrs = list(range(1, len(examples1[0]) - 1))
        nPossibleVals = [0, 2, 2, 0, 3]
        boundsPerVar = [0, 0, 0, 1, 0]

        self.t1 = BuildQuantTree.QuantTreeBoot(examples1, attrs, nPossibleVals, boundsPerVar)
        self.examples1 = examples1

    def _setupTree2(self):
        examples1 = [['p1', 0.1, 1, 0.1, 0], ['p2', 0.1, 0, 0.1, 1], ['p3', 0.1, 0, 1.1, 2],
                     ['p4', 0.1, 1, 1.1, 2], ['p5', 1.1, 0, 0.1, 2], ['p6', 1.1, 0, 1.1, 2],
                     ['p7', 1.1, 1, 0.1, 2], ['p8', 1.1, 1, 1.1, 0]]
        attrs = list(range(1, len(examples1[0]) - 1))
        nPossibleVals = [0, 0, 2, 0, 3]
        boundsPerVar = [0, 1, 0, 1, 0]

        self.t2 = BuildQuantTree.QuantTreeBoot(examples1, attrs, nPossibleVals, boundsPerVar)
        self.examples2 = examples1

    def _setupTree1a(self):
        examples1 = [['p1', 0, 1, 0.1, 4.0, 0], ['p2', 0, 0, 0.1, 4.1, 1], ['p3', 0, 0, 1.1, 4.2, 2],
                     ['p4', 0, 1, 1.1, 4.2, 2], ['p5', 1, 0, 0.1, 4.2, 2], ['p6', 1, 0, 1.1, 4.2, 2],
                     ['p7', 1, 1, 0.1, 4.2, 2], ['p8', 1, 1, 1.1, 4.0, 0]]
        attrs = list(range(1, len(examples1[0]) - 1))
        nPossibleVals = [0, 2, 2, 0, 0, 3]
        boundsPerVar = [0, 0, 0, 1, -1, 0]

        self.t1 = BuildQuantTree.QuantTreeBoot(examples1, attrs, nPossibleVals, boundsPerVar)
        self.examples1 = examples1

    def test0Cmp(self):
        # " testing tree comparisons "
        self._setupTree1()
        self._setupTree2()
        assert self.t1 == self.t1, 'self equals failed'
        assert self.t2 == self.t2, 'self equals failed'
        assert self.t1 != self.t2, 'not equals failed'

    def test1Tree(self):
        # " testing tree1 "
        self._setupTree1()
        with open(self.qTree1Name, 'r') as inTFile:
            buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
            inTFile.close()
        with io.BytesIO(buf) as inFile:
            t2 = pickle.load(inFile)
        assert self.t1 == t2, 'Incorrect tree generated. '

        self.assertIn('Var: 2 []', str(self.t1))
        self.assertEqual(self.t1.GetQuantBounds(), [])

    def test2Tree(self):
        # " testing tree2 "
        self._setupTree2()
        with open(self.qTree2Name, 'r') as inTFile:
            buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
            inTFile.close()
        with io.BytesIO(buf) as inFile:
            t2 = pickle.load(inFile)
        assert self.t2 == t2, 'Incorrect tree generated.'

    def test3Classify(self):
        # " testing classification "
        self._setupTree1()
        self._setupTree2()
        for i in range(len(self.examples1)):
            self.assertEqual(
              self.t1.ClassifyExample(self.examples1[i]), self.examples1[i][-1],
              msg='examples1[%d] misclassified' % i)
        for i in range(len(self.examples2)):
            self.assertEqual(
              self.t2.ClassifyExample(self.examples2[i]), self.examples2[i][-1],
              msg='examples2[%d] misclassified' % i)

    def test4UnusedVars(self):
        # " testing unused variables "
        self._setupTree1a()
        with open(self.qTree1Name, 'r') as inTFile:
            buf = inTFile.read().replace('\r\n', '\n').encode('utf-8')
            inTFile.close()
        with io.BytesIO(buf) as inFile:
            t2 = pickle.load(inFile)
        assert self.t1 == t2, 'Incorrect tree generated.'
        for i in range(len(self.examples1)):
            self.assertEqual(
              self.t1.ClassifyExample(self.examples1[i]), self.examples1[i][-1],
              'examples1[%d] misclassified' % i)

    def test5Bug29(self):
        # """ a more extensive test of the cmp stuff using hand-built trees """
        import copy

        t1 = QuantTreeNode(None, 't1')
        t1.SetQuantBounds([1.])
        c1 = QuantTreeNode(t1, 'c1')
        c1.SetQuantBounds([2.])
        t1.AddChildNode(c1)
        c2 = QuantTreeNode(t1, 'c2')
        c2.SetQuantBounds([2.])
        t1.AddChildNode(c2)
        c11 = QuantTreeNode(c1, 'c11')
        c11.SetQuantBounds([3.])
        c1.AddChildNode(c11)
        c12 = QuantTreeNode(c1, 'c12')
        c12.SetQuantBounds([3.])
        c1.AddChildNode(c12)
        assert not cmp(t1, copy.deepcopy(t1)), 'self equality failed'

        t2 = QuantTreeNode(None, 't1')
        t2.SetQuantBounds([1.])
        c1 = QuantTreeNode(t2, 'c1')
        c1.SetQuantBounds([2.])
        t2.AddChildNode(c1)
        c2 = QuantTreeNode(t2, 'c2')
        c2.SetQuantBounds([2.])
        t2.AddChildNode(c2)
        c11 = QuantTreeNode(c1, 'c11')
        c11.SetQuantBounds([3.])
        c1.AddChildNode(c11)
        c12 = QuantTreeNode(c1, 'c12')
        c12.SetQuantBounds([3.00003])
        c1.AddChildNode(c12)
        assert cmp(t1, t2), 'inequality failed'

    def test6Bug29_2(self):
        # """ a more extensive test of the cmp stuff using pickled trees"""
        import os
        with open(os.path.join(RDConfig.RDCodeDir, 'ML', 'DecTree', 'test_data', 'CmpTree1.pkl'),
                  'r') as t1TFile:
            buf = t1TFile.read().replace('\r\n', '\n').encode('utf-8')
            t1TFile.close()
        with io.BytesIO(buf) as t1File:
            t1 = pickle.load(t1File)
        with open(os.path.join(RDConfig.RDCodeDir, 'ML', 'DecTree', 'test_data', 'CmpTree2.pkl'),
                  'r') as t2TFile:
            buf = t2TFile.read().replace('\r\n', '\n').encode('utf-8')
            t2TFile.close()
        with io.BytesIO(buf) as t2File:
            t2 = pickle.load(t2File)
        assert cmp(t1, t2), 'equality failed'

    def test7Recycle(self):
        # """ try recycling descriptors """
        examples1 = [[3, 0, 0],
                     [3, 1, 1],
                     [1, 0, 0],
                     [0, 0, 1],
                     [1, 1, 0], ]
        attrs = list(range(2))
        nPossibleVals = [2, 2, 2]
        boundsPerVar = [1, 0, 0]
        self.t1 = BuildQuantTree.QuantTreeBoot(examples1, attrs, nPossibleVals, boundsPerVar,
                                               recycleVars=1)
        assert self.t1.GetLabel() == 0, self.t1.GetLabel()
        assert self.t1.GetChildren()[0].GetLabel() == 1
        assert self.t1.GetChildren()[1].GetLabel() == 1
        assert self.t1.GetChildren()[1].GetChildren()[0].GetLabel() == 0
        assert self.t1.GetChildren()[1].GetChildren()[1].GetLabel() == 0

    def test8RandomForest(self):
        # """ try random forests descriptors """
        import random
        random.seed(23)
        nAttrs = 100
        nPts = 10
        examples = []
        for _ in range(nPts):
            descrs = [int(random.random() > 0.5) for _ in range(nAttrs)]
            act = sum(descrs) > nAttrs / 2
            examples.append(descrs + [act])
        attrs = list(range(nAttrs))
        nPossibleVals = [2] * (nAttrs + 1)
        boundsPerVar = [0] * nAttrs + [0]
        self.t1 = BuildQuantTree.QuantTreeBoot(examples, attrs, nPossibleVals, boundsPerVar, maxDepth=1,
                                               recycleVars=1, randomDescriptors=3)
        self.assertEqual(self.t1.GetLabel(), 49)
        self.assertEqual(self.t1.GetChildren()[0].GetLabel(), 3)
        self.assertEqual(self.t1.GetChildren()[1].GetLabel(), 54)

    def test_exampleCode(self):
        f = StringIO()
        with redirect_stdout(f):
            BuildQuantTree.TestTree()
        self.assertIn('Var: 2', f.getvalue())


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
