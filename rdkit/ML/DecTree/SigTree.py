#
#  Copyright (C) 2003-2005  Rational Discovery LLC
#   All Rights Reserved
#
""" Defines the class SigTreeNode, used to represent trees that
 use signatures (bit vectors) to represent data.  As inputs (examples),
 SigTreeNode's expect 3-sequences: (label,sig,act)

  _SigTreeNode_ is derived from _DecTree.DecTreeNode_

"""
from rdkit.ML.DecTree import DecTree
from rdkit.DataStructs.VectCollection import VectCollection
import copy


class SigTreeNode(DecTree.DecTreeNode):
  """

  """

  def NameModel(self, *args, **kwargs):
    pass

  def ClassifyExample(self, example, appendExamples=0):
    """ Recursively classify an example by running it through the tree

      **Arguments**

        - example: the example to be classified, a sequence at least
          2 long:
           ( id, sig )
          where sig is a BitVector (or something supporting __getitem__)
          additional fields will be ignored.

        - appendExamples: if this is nonzero then this node (and all children)
          will store the example

      **Returns**

        the classification of _example_

    """
    if appendExamples:
      self.examples.append(example)
    if self.terminalNode:
      return self.label
    else:
      sig = example[1]
      val = sig[self.label]
      # print 'val:',val
      if val and isinstance(sig, VectCollection):
        # we need to copy and modify the example:
        sig = copy.copy(sig)
        sig.DetachVectsNotMatchingBit(self.label)
        ex = [example[0], sig]
        if len(example) > 2:
          ex.extend(example[2:])
        example = ex
      return self.children[val].ClassifyExample(example, appendExamples=appendExamples)

  def __init__(self, *args, **kwargs):
    DecTree.DecTreeNode.__init__(self, *args, **kwargs)
