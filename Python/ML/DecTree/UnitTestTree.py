#
#  Copyright (C) 2000  greg Landrum
#

""" unit testing code for trees and decision trees (not learning/xvalidation) """
from rdkit import RDConfig
import unittest
from rdkit.ML.DecTree import Tree,DecTree
import copy
import cPickle


class TreeTestCase(unittest.TestCase):
  def setUp(self):
    print '\n%s: '%self.shortDescription(),
    self.baseTree = Tree.TreeNode(None,'root')    
    self.pickleFileName = RDConfig.RDCodeDir+'/ML/DecTree/test_data/treeunit.pkl'
    
  def testAdd(self):
    " testing AddChild "
    tree = self.baseTree
    addWorks = 1
    try:
      for i in xrange(3):
        child = tree.AddChild('child %d'%i)
    except:
      addWorks = 0
    assert addWorks,'Tree.AddChild failed'

  def testGetName(self):
    " testing GetName "
    tree = self.baseTree
    # we know this works
    assert tree.GetName()=='root','Tree.GetName failed'

  def _readyTree(self):
    tree = self.baseTree
    tree.AddChild('child0')
    tree.AddChild('child1')

  def testGetChildren(self):
    " testing GetChildren "
    self._readyTree()
    tree = self.baseTree
    child = tree.GetChildren()[1]
    assert child.GetName()=='child1','Tree.GetChildren failed'

  def testPrune(self):
    " testing PruneChild"
    self._readyTree()
    tree = self.baseTree
    child = tree.GetChildren()[0]
    tree.PruneChild(child)
    child = tree.GetChildren()[0]
    assert child.GetName()=='child1','Tree.PruneChild failed'

  def testPickle(self):
    " testing tree pickle "
    self._readyTree()
    pickleWorked=1
    try:
      self.baseTree.Pickle(self.pickleFileName)
    except:
      pickleWorked=0
    assert pickleWorked,'tree.Pickle Failed'      

  def testEquals(self):
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
    
  def testPickleEquals(self):
    " testing pickled tree equals "
    self._readyTree()
    pFile = open(self.pickleFileName,'r')
    oTree = cPickle.load(pFile)

    assert oTree == self.baseTree,'Pickle inequality test failed'

    self.baseTree.PruneChild(self.baseTree.GetChildren()[0])
    assert oTree != self.baseTree,'Pickle inequality test failed (bad Tree.__cmp__)'    

  def testCopy(self):
    " testing deepcopy on trees "
    self._readyTree()
    nTree = copy.deepcopy(self.baseTree)
    assert nTree == self.baseTree,'deepcopy failed'

  def testIn(self):
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
  
  
