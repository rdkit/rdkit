#
#  Copyright (C) 2000  greg Landrum
#
""" unit testing code for trees and decision trees (not learning/xvalidation) """


import unittest

from rdkit.ML.DecTree import TreeUtils
from rdkit.ML.DecTree.DecTree import DecTreeNode as Node


class TestTreeUtils(unittest.TestCase):

  def test_TreeUtils(self):
    # Tree is d1(d2,d3(d2,d4))
    t1 = Node(None, 'd1', 1)
    t2 = Node(None, 'd2', 2)
    t1.AddChildNode(t2)
    t2 = Node(None, 'd3', 3)
    t1.AddChildNode(t2)
    t3 = Node(None, 'd4', 4)
    t2.AddChildNode(t3)
    t3 = Node(None, 'd2', 2)
    t2.AddChildNode(t3)

    r = TreeUtils.CollectLabelLevels(t1, {})
    self.assertEqual(r, {1: 0, 2: 1, 3: 1, 4: 2})

    # Only to depth 2
    r = TreeUtils.CollectLabelLevels(t1, {}, 0, 2)
    self.assertEqual(r, {1: 0, 2: 1, 3: 1})

    # Check that we can handle subtrees:
    r = TreeUtils.CollectLabelLevels(t1, {}, 1, 2)
    self.assertEqual(r, {1: 1})

    names = TreeUtils.CollectDescriptorNames(t1, {})
    self.assertEqual(names, {1: 'd1', 2: 'd2', 3: 'd3', 4: 'd4'})

    names = TreeUtils.CollectDescriptorNames(t1, {}, 0, 2)
    self.assertEqual(names, {1: 'd1', 2: 'd2', 3: 'd3'})

    names = TreeUtils.CollectDescriptorNames(t1, {}, 1, 2)
    self.assertEqual(names, {1: 'd1'})


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
