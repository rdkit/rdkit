# $Id$
#
#  Copyright (C) 2001, 2003  greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
""" Defines the class _QuantTreeNode_, used to represent decision trees with automatic
 quantization bounds

  _QuantTreeNode_ is derived from _DecTree.DecTreeNode_

"""
from rdkit.ML.DecTree import DecTree, Tree


class QuantTreeNode(DecTree.DecTreeNode):
  """

  """

  def __init__(self, *args, **kwargs):
    DecTree.DecTreeNode.__init__(self, *args, **kwargs)
    self.qBounds = []
    self.nBounds = 0

  def ClassifyExample(self, example, appendExamples=0):
    """ Recursively classify an example by running it through the tree

      **Arguments**

        - example: the example to be classified

        - appendExamples: if this is nonzero then this node (and all children)
          will store the example

      **Returns**

        the classification of _example_

      **NOTE:**
        In the interest of speed, I don't use accessor functions
        here.  So if you subclass DecTreeNode for your own trees, you'll
        have to either include ClassifyExample or avoid changing the names
        of the instance variables this needs.

    """
    if appendExamples:
      self.examples.append(example)
    if self.terminalNode:
      return self.label
    else:
      val = example[self.label]
      if not hasattr(self, 'nBounds'):
        self.nBounds = len(self.qBounds)
      if self.nBounds:
        for i, bound in enumerate(self.qBounds):
          if val < bound:
            val = i
            break
        else:
          val = i + 1
      else:
        val = int(val)
      return self.children[val].ClassifyExample(example, appendExamples=appendExamples)

  def SetQuantBounds(self, qBounds):
    self.qBounds = qBounds[:]
    self.nBounds = len(self.qBounds)

  def GetQuantBounds(self):
    return self.qBounds

  def __cmp__(self, other):
    return (self < other) * -1 or (other < self) * 1

  def __lt__(self, other):
    if str(type(self)) < str(type(other)):
      return True
    if self.qBounds < other.qBounds:
      return True
    if Tree.TreeNode.__lt__(self, other):
      return True
    return False

  def __eq__(self, other):
    return not self < other and not other < self

  def __str__(self):
    """ returns a string representation of the tree

      **Note**

        this works recursively

    """
    here = '%s%s %s\n' % ('  ' * self.level, self.name, str(self.qBounds))
    for child in self.children:
      here = here + str(child)
    return here
