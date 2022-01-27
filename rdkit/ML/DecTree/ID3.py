#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#
""" ID3 Decision Trees

  contains an implementation of the ID3 decision tree algorithm
  as described in Tom Mitchell's book "Machine Learning"

  It relies upon the _Tree.TreeNode_ data structure (or something
    with the same API) defined locally to represent the trees

"""

import numpy

from rdkit.ML.DecTree import DecTree
from rdkit.ML.InfoTheory import entropy


def CalcTotalEntropy(examples, nPossibleVals):
  """ Calculates the total entropy of the data set (w.r.t. the results)

   **Arguments**

    - examples: a list (nInstances long) of lists of variable values + instance
              values
    - nPossibleVals: a list (nVars long) of the number of possible values each variable
      can adopt.

   **Returns**

     a float containing the informational entropy of the data set.

  """
  nRes = nPossibleVals[-1]
  resList = numpy.zeros(nRes, 'i')
  for example in examples:
    resList[int(example[-1])] += 1
  return entropy.InfoEntropy(resList)


def GenVarTable(examples, nPossibleVals, vars):
  """Generates a list of variable tables for the examples passed in.

    The table for a given variable records the number of times each possible value
    of that variable appears for each possible result of the function.

  **Arguments**

    - examples: a list (nInstances long) of lists of variable values + instance
              values

    - nPossibleVals: a list containing the number of possible values of
                   each variable + the number of values of the function.

    - vars:  a list of the variables to include in the var table


  **Returns**

      a list of variable result tables. Each table is a Numeric array
        which is varValues x nResults
  """
  nVars: int = len(vars)
  nFuncVals = nPossibleVals[-1]
  res = [numpy.zeros((nPossibleVals[vars[i]], nFuncVals), dtype=numpy.int64) for i in range(nVars)]

  for example in examples:
    value = int(example[-1])
    for i, var in enumerate(vars):
      res[i][int(example[var]), value] += 1

  return res


def ID3(examples, target, attrs, nPossibleVals, depth=0, maxDepth=-1, **kwargs):
  """ Implements the ID3 algorithm for constructing decision trees.

    From Mitchell's book, page 56

    This is *slightly* modified from Mitchell's book because it supports
      multivalued (non-binary) results.

    **Arguments**

      - examples: a list (nInstances long) of lists of variable values + instance
              values

      - target: an int

      - attrs: a list of ints indicating which variables can be used in the tree

      - nPossibleVals: a list containing the number of possible values of
                   every variable.

      - depth: (optional) the current depth in the tree

      - maxDepth: (optional) the maximum depth to which the tree
                   will be grown

    **Returns**

     a DecTree.DecTreeNode with the decision tree

    **NOTE:** This code cannot bootstrap (start from nothing...)
          use _ID3Boot_ (below) for that.
  """
  varTable = GenVarTable(examples, nPossibleVals, attrs)
  tree = DecTree.DecTreeNode(None, 'node')

  # store the total entropy... in case that is interesting
  totEntropy = CalcTotalEntropy(examples, nPossibleVals)
  tree.SetData(totEntropy)
  # tree.SetExamples(examples)

  # the matrix of results for this target:
  tMat = GenVarTable(examples, nPossibleVals, [target])[0]
  # counts of each result code:
  counts = sum(tMat)
  nzCounts = numpy.nonzero(counts)[0]

  if len(nzCounts) == 1:
    # bottomed out because there is only one result code left
    #  with any counts (i.e. there's only one type of example
    #  left... this is GOOD!).
    res = nzCounts[0]
    tree.SetLabel(res)
    tree.SetName(str(res))
    tree.SetTerminal(1)
  elif len(attrs) == 0 or (depth >= maxDepth >= 0):
    # Bottomed out: no variables left or max depth hit
    #  We don't really know what to do here, so
    #  use the heuristic of picking the most prevalent
    #  result
    v = numpy.argmax(counts)
    tree.SetLabel(v)
    tree.SetName(f'{v}?')
    tree.SetTerminal(1)
  else:
    # find the variable which gives us the largest information gain
    gains = [entropy.InfoGain(x) for x in varTable]
    best = attrs[numpy.argmax(gains)]

    # remove that variable from the lists of possible variables
    nextAttrs = list(attrs)
    if not kwargs.get('recycleVars', 0):
      nextAttrs.remove(best)

    # set some info at this node
    tree.SetName(f'Var: {best}')
    tree.SetLabel(best)
    # tree.SetExamples(examples)
    tree.SetTerminal(0)

    # loop over possible values of the new variable and
    #  build a subtree for each one
    for val in range(nPossibleVals[best]):
      nextExamples = [example for example in examples if example[best] == val]
      if len(nextExamples) == 0:
        # this particular value of the variable has no examples,
        #  so there's not much sense in recursing.
        #  This can (and does) happen.
        v = numpy.argmax(counts)
        tree.AddChild(f'{v}', label=v, data=0.0, isTerminal=1)
      else:
        # recurse
        tree.AddChildNode(
          ID3(nextExamples, best, nextAttrs, nPossibleVals, depth + 1, maxDepth, **kwargs))
  return tree


def ID3Boot(examples, attrs, nPossibleVals, initialVar=None, depth=0, maxDepth=-1, **kwargs):
  """ Bootstrapping code for the ID3 algorithm

    see ID3 for descriptions of the arguments

    If _initialVar_ is not set, the algorithm will automatically
     choose the first variable in the tree (the standard greedy
     approach).  Otherwise, _initialVar_ will be used as the first
     split.

  """
  totEntropy = CalcTotalEntropy(examples, nPossibleVals)
  varTable = GenVarTable(examples, nPossibleVals, attrs)

  tree = DecTree.DecTreeNode(None, 'node')
  # tree.SetExamples(examples)
  tree._nResultCodes = nPossibleVals[-1]

  # <perl>you've got to love any language which will let you
  # do this much work in a single line :-)</perl>
  best = attrs[numpy.argmax([entropy.InfoGain(x) for x in varTable])] if initialVar is None else initialVar

  tree.SetName(f'Var: {best}')
  tree.SetData(totEntropy)
  tree.SetLabel(best)
  tree.SetTerminal(0)
  nextAttrs = list(attrs)
  if not kwargs.get('recycleVars', 0):
    nextAttrs.remove(best)

  for val in range(nPossibleVals[best]):
    nextExamples = [example for example in examples if example[best] == val]
    tree.AddChildNode(ID3(nextExamples, best, nextAttrs, nPossibleVals, depth, maxDepth, **kwargs))
  return tree
