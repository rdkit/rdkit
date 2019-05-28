#
#  Copyright (C) 2000  greg Landrum
#
""" unit testing code for trees and decision trees (not learning/xvalidation) """


import copy
import os
import unittest

from rdkit import RDConfig
from rdkit.ML.DecTree import Tree
from rdkit.TestRunner import redirect_stdout
from io import StringIO
import pickle


class TreeTestCase(unittest.TestCase):

    def setUp(self):
        self.baseTree = Tree.TreeNode(None, 'root')
        self.pickleFileName = RDConfig.RDCodeDir + '/ML/DecTree/test_data/treeunit.pkl'

    def test_Tree(self):
        tree = Tree.TreeNode(None, 'root', label=0)
        self.assertEqual(tree.GetLevel(), 0)
        self.assertEqual(tree.GetName(), 'root')
        self.assertEqual(tree.GetData(), None)
        self.assertEqual(tree.GetTerminal(), False)
        self.assertEqual(tree.GetLabel(), 0)
        self.assertEqual(tree.GetParent(), None)
        self.assertEqual(tree.GetChildren(), [])

        for i in range(3):
            child = tree.AddChild('child {0}'.format(i), i + 1, data={'key': 'value'})
            self.assertEqual(child.GetLevel(), 1)
            self.assertEqual(child.GetName(), 'child {0}'.format(i))
            self.assertEqual(child.GetData(), {'key': 'value'})
            self.assertEqual(child.GetLabel(), i + 1)
            self.assertEqual(child.GetParent(), tree)
            self.assertEqual(child.GetChildren(), [])
        children = tree.GetChildren()
        self.assertEqual(len(children), 3)
        children[0].AddChild('terminal', 4, isTerminal=True)

        s = str(tree)
        self.assertIn('root', s)
        self.assertIn('    terminal', s)
        self.assertIn('  child 2', s)

        tree.NameTree(['a', 'b', 'c', 'd', 'e'])
        self.assertEqual(str(tree), 'a\n  b\n    terminal\n  c\n  d\n')

        tree.PruneChild(children[1])
        self.assertEqual(str(tree), 'a\n  b\n    terminal\n  d\n')

        f = StringIO()
        with redirect_stdout(f):
            tree.Print(showData=True)
        s = f.getvalue()
        self.assertIn('value', s)
        self.assertIn('None', s)

        f = StringIO()
        with redirect_stdout(f):
            tree.Print()
        s = f.getvalue()
        self.assertNotIn('value', s)
        self.assertNotIn('None', s)

        tree.Destroy()
        self.assertEqual(str(tree), 'a\n')

    def _readyTree(self):
        tree = self.baseTree
        tree.AddChild('child0')
        tree.AddChild('child1')

    def test5Equals(self):
        # " testing tree equals "
        nTree = Tree.TreeNode(None, 'root')
        self._readyTree()
        tTree = self.baseTree
        self.baseTree = nTree
        self._readyTree()
        assert tTree == self.baseTree, 'Equality test 1 failed. (bad Tree.__cmp__)'
        assert self.baseTree == tTree, 'Equality test 2 failed. (bad Tree.__cmp__)'
        tTree.AddChild('child2')
        assert tTree != self.baseTree, 'Inequality test 1 failed. (bad Tree.__cmp__)'
        assert self.baseTree != tTree, 'Inequality test 2 failed. (bad Tree.__cmp__)'

        self.assertTrue(tTree > self.baseTree, msg='Larger tree is greater')
        self.assertEqual(tTree.__cmp__(self.baseTree), 1)

    def test6PickleEquals(self):
        # " testing pickled tree equals "
        self._readyTree()
        pkl = pickle.dumps(self.baseTree)
        oTree = pickle.loads(pkl)

        assert oTree == self.baseTree, 'Pickle inequality test failed'
        self.assertEqual(oTree.__cmp__(self.baseTree), 0)

        self.baseTree.PruneChild(self.baseTree.GetChildren()[0])
        assert oTree != self.baseTree, 'Pickle inequality test failed (bad Tree.__cmp__)'
        self.assertEqual(abs(oTree.__cmp__(self.baseTree)), 1)

    def test7Copy(self):
        # " testing deepcopy on trees "
        self._readyTree()
        nTree = copy.deepcopy(self.baseTree)
        assert nTree == self.baseTree, 'deepcopy failed'

    def test8In(self):
        # " testing list membership "
        self._readyTree()
        nTree = copy.deepcopy(self.baseTree)
        nTree2 = copy.deepcopy(self.baseTree)
        nTree2.PruneChild(self.baseTree.GetChildren()[0])
        tList = [nTree2, nTree2, nTree]
        assert self.baseTree in tList, 'list membership (tree in list) failed'
        tList = [nTree2, nTree2]
        assert self.baseTree not in tList, 'list membership (tree not in list) failed'

    def test_exampleCode(self):
        try:
            f = StringIO()
            with redirect_stdout(f):
                Tree._exampleCode()
            self.assertTrue(os.path.isfile('save.pkl'))
            self.assertIn('tree==tree2 False', f.getvalue(), 'Example didn' 't run to end')
        finally:
            if os.path.isfile('save.pkl'):
                os.remove('save.pkl')


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
