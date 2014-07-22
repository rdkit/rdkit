#
#  Copyright (C) 2000  greg Landrum
#

""" unit testing code for trees and decision trees (not learning/xvalidation) """
from __future__ import print_function
import unittest
import copy
from rdkit import RDConfig
from rdkit.ML.DecTree import Tree,DecTree
from rdkit.six.moves import cPickle


class TreeTestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
    self.baseTree = Tree.TreeNode(None,'root')    
    self.pickleFileName = RDConfig.RDCodeDir+'/ML/DecTree/test_data/treeunit.pkl'
    
  def test0Add(self):
    " testing AddChild "
    tree = self.baseTree
    for i in range(3):
      child = tree.AddChild('child %d'%i)

  def test1GetName(self):
    " testing GetName "
    tree = self.baseTree
    # we know this works
    assert tree.GetName()=='root','Tree.GetName failed'

  def _readyTree(self):
    tree = self.baseTree
    tree.AddChild('child0')
    tree.AddChild('child1')

  def test2GetChildren(self):
    " testing GetChildren "
    self._readyTree()
    tree = self.baseTree
    child = tree.GetChildren()[1]
    assert child.GetName()=='child1','Tree.GetChildren failed'

  def test3Prune(self):
    " testing PruneChild"
    self._readyTree()
    tree = self.baseTree
    child = tree.GetChildren()[0]
    tree.PruneChild(child)
    child = tree.GetChildren()[0]
    assert child.GetName()=='child1','Tree.PruneChild failed'

  def test5Equals(self):
    " testing tree equals "
    nTree = Tree.TreeNode(None,'root')
    self._readyTree()
    tTree = self.baseTree
    self.baseTree = nTree
    self._readyTree()
    assert tTree == self.baseTree,'Equality test 1 failed. (bad Tree.__cmp__)'
    assert self.baseTree==tTree,'Equality test 2 failed. (bad Tree.__cmp__)'
    tTree.AddChild('child2')
    assert tTree != self.baseTree,'Inequality test 1 failed. (bad Tree.__cmp__)'
    assert self.baseTree != tTree,'Inequality test 2 failed. (bad Tree.__cmp__)'    
    
  def test6PickleEquals(self):
    " testing pickled tree equals "
    self._readyTree()
    pkl = cPickle.dumps(self.baseTree)
    oTree = cPickle.loads(pkl)

    assert oTree == self.baseTree,'Pickle inequality test failed'

    self.baseTree.PruneChild(self.baseTree.GetChildren()[0])
    assert oTree != self.baseTree,'Pickle inequality test failed (bad Tree.__cmp__)'    

  def test7Copy(self):
    " testing deepcopy on trees "
    self._readyTree()
    nTree = copy.deepcopy(self.baseTree)
    assert nTree == self.baseTree,'deepcopy failed'

  def test8In(self):
    " testing list membership "
    self._readyTree()
    nTree = copy.deepcopy(self.baseTree)
    nTree2 = copy.deepcopy(self.baseTree)
    nTree2.PruneChild(self.baseTree.GetChildren()[0])    
    tList = [nTree2,nTree2,nTree]
    assert self.baseTree in tList, 'list membership (tree in list) failed'
    tList = [nTree2,nTree2]    
    assert self.baseTree not in tList, 'list membership (tree not in list) failed'

if __name__ == '__main__':
  unittest.main()
  
  
