#
#  Copyright (C) 2000-2004  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#
""" Defines the class _DecTreeNode_, used to represent decision trees

  _DecTreeNode_ is derived from _Tree.TreeNode_

"""
from rdkit.ML.DecTree import Tree


class DecTreeNode(Tree.TreeNode):
  """ This is used to represent decision trees

   _DecTreeNode_s are simultaneously the roots and branches of decision trees.
   Everything is nice and recursive.

   _DecTreeNode_s can save the following pieces of internal state, accessible via
     standard setter/getter functions:

     1) _Examples_: a list of examples which have been classified

     2) _BadExamples_: a list of examples which have been misclassified

     3) _TrainingExamples_: the list of examples used to train the tree

     4) _TestExamples_: the list of examples used to test the tree

  """

  def __init__(self, *args, **kwargs):
    # apply(Tree.TreeNode.__init__,(self,)+args,kwargs)
    Tree.TreeNode.__init__(self, *args, **kwargs)
    self.examples = []
    self.badExamples = []
    self.trainingExamples = []
    self.testExamples = []

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
      return self.children[val].ClassifyExample(example, appendExamples)

  def AddChild(self, name, label=None, data=None, isTerminal=0):
    """ Constructs and adds a child with the specified data to our list

      **Arguments**

       - name: the name of the new node

       - label: the label of the new node (should be an integer)

       - data: the data to be stored in the new node

       - isTerminal: a toggle to indicate whether or not the new node is
         a terminal (leaf) node.

      **Returns*

        the _DecTreeNode_ which is constructed

    """
    child = DecTreeNode(self, name, label, data, level=self.level + 1, isTerminal=isTerminal)
    self.children.append(child)
    return child

  def GetExamples(self):
    return self.examples

  def SetExamples(self, examples):
    self.examples = examples

  def GetBadExamples(self):
    return self.badExamples

  def SetBadExamples(self, examples):
    self.badExamples = examples

  def GetTrainingExamples(self):
    return self.trainingExamples

  def SetTrainingExamples(self, examples):
    self.trainingExamples = examples

  def GetTestExamples(self):
    return self.testExamples

  def SetTestExamples(self, examples):
    self.testExamples = examples

  def ClearExamples(self):
    self.examples = []
    self.badExamples = []
    self.trainingExamples = []
    self.testExamples = []
    for child in self.GetChildren():
      child.ClearExamples()
