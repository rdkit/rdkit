#
#  Copyright (C) 2000  greg Landrum
#
""" handles doing cross validation with decision trees

This is, perhaps, a little misleading.  For the purposes of this module,
cross validation == evaluating the accuracy of a tree.


"""


import numpy

from rdkit.ML.Data import SplitData
from rdkit.ML.DecTree import ID3
from rdkit.ML.DecTree import randomtest


def ChooseOptimalRoot(examples, trainExamples, testExamples, attrs, nPossibleVals, treeBuilder,
                      nQuantBounds=[], **kwargs):
  """ loops through all possible tree roots and chooses the one which produces the best tree

  **Arguments**

    - examples: the full set of examples

    - trainExamples: the training examples

    - testExamples: the testing examples

    - attrs: a list of attributes to consider in the tree building

    - nPossibleVals: a list of the number of possible values each variable can adopt

    - treeBuilder: the function to be used to actually build the tree

    - nQuantBounds: an optional list.  If present, it's assumed that the builder
      algorithm takes this argument as well (for building QuantTrees)

  **Returns**

    The best tree found

  **Notes**

    1) Trees are built using _trainExamples_

    2) Testing of each tree (to determine which is best) is done using _CrossValidate_ and
       the entire set of data (i.e. all of _examples_)

    3) _trainExamples_ is not used at all, which immediately raises the question of
       why it's even being passed in

  """
  attrs = attrs[:]
  if nQuantBounds:
    for i in range(len(nQuantBounds)):
      if nQuantBounds[i] == -1 and i in attrs:
        attrs.remove(i)
  nAttrs = len(attrs)
  trees = [None] * nAttrs
  errs = [0] * nAttrs
  errs[0] = 1e6

  for i in range(1, nAttrs):
    argD = {'initialVar': attrs[i]}
    argD.update(kwargs)
    if nQuantBounds is None or nQuantBounds == []:
      trees[i] = treeBuilder(trainExamples, attrs, nPossibleVals, **argD)
    else:
      trees[i] = treeBuilder(trainExamples, attrs, nPossibleVals, nQuantBounds, **argD)
    if trees[i]:
      errs[i], _ = CrossValidate(trees[i], examples, appendExamples=0)
    else:
      errs[i] = 1e6
  best = numpy.argmin(errs)
  # FIX: this used to say 'trees[i]', could that possibly have been right?
  return trees[best]


def CrossValidate(tree, testExamples, appendExamples=0):
  """ Determines the classification error for the testExamples

    **Arguments**

      - tree: a decision tree (or anything supporting a _ClassifyExample()_ method)

      - testExamples: a list of examples to be used for testing

      - appendExamples: a toggle which is passed along to the tree as it does
        the classification. The trees can use this to store the examples they
        classify locally.

    **Returns**

      a 2-tuple consisting of:

        1) the percent error of the tree

        2) a list of misclassified examples

  """
  nTest = len(testExamples)
  nBad = 0
  badExamples = []
  for i in range(nTest):
    testEx = testExamples[i]
    trueRes = testEx[-1]
    res = tree.ClassifyExample(testEx, appendExamples)
    if (trueRes != res).any():
      badExamples.append(testEx)
      nBad += 1

  return float(nBad) / nTest, badExamples


def CrossValidationDriver(examples, attrs, nPossibleVals, holdOutFrac=.3, silent=0,
                          calcTotalError=0, treeBuilder=ID3.ID3Boot, lessGreedy=0, startAt=None,
                          nQuantBounds=[], maxDepth=-1, **kwargs):
  """ Driver function for building trees and doing cross validation

    **Arguments**

      - examples: the full set of examples

      - attrs: a list of attributes to consider in the tree building

      - nPossibleVals: a list of the number of possible values each variable can adopt

      - holdOutFrac: the fraction of the data which should be reserved for the hold-out set
         (used to calculate the error)

      - silent: a toggle used to control how much visual noise this makes as it goes.

      - calcTotalError: a toggle used to indicate whether the classification error
        of the tree should be calculated using the entire data set (when true) or just
        the training hold out set (when false)

      - treeBuilder: the function to call to build the tree

      - lessGreedy: toggles use of the less greedy tree growth algorithm (see
        _ChooseOptimalRoot_).

      - startAt: forces the tree to be rooted at this descriptor

      - nQuantBounds: an optional list.  If present, it's assumed that the builder
        algorithm takes this argument as well (for building QuantTrees)

      - maxDepth: an optional integer.  If present, it's assumed that the builder
        algorithm takes this argument as well

    **Returns**

       a 2-tuple containing:

         1) the tree

         2) the cross-validation error of the tree

  """
  nTot = len(examples)
  if not kwargs.get('replacementSelection', 0):
    testIndices, trainIndices = SplitData.SplitIndices(nTot, holdOutFrac, silent=1, legacy=1,
                                                       replacement=0)
  else:
    testIndices, trainIndices = SplitData.SplitIndices(nTot, holdOutFrac, silent=1, legacy=0,
                                                       replacement=1)
  trainExamples = [examples[x] for x in trainIndices]
  testExamples = [examples[x] for x in testIndices]

  nTrain = len(trainExamples)
  if not silent:
    print('Training with %d examples' % (nTrain))

  if not lessGreedy:
    if nQuantBounds is None or nQuantBounds == []:
      tree = treeBuilder(trainExamples, attrs, nPossibleVals, initialVar=startAt, maxDepth=maxDepth,
                         **kwargs)
    else:
      tree = treeBuilder(trainExamples, attrs, nPossibleVals, nQuantBounds, initialVar=startAt,
                         maxDepth=maxDepth, **kwargs)
  else:
    tree = ChooseOptimalRoot(examples, trainExamples, testExamples, attrs, nPossibleVals,
                             treeBuilder, nQuantBounds, maxDepth=maxDepth, **kwargs)

  nTest = len(testExamples)
  if not silent:
    print('Testing with %d examples' % nTest)
  if not calcTotalError:
    xValError, badExamples = CrossValidate(tree, testExamples, appendExamples=1)
  else:
    xValError, badExamples = CrossValidate(tree, examples, appendExamples=0)
  if not silent:
    print('Validation error was %%%4.2f' % (100 * xValError))
  tree.SetBadExamples(badExamples)
  tree.SetTrainingExamples(trainExamples)
  tree.SetTestExamples(testExamples)
  tree._trainIndices = trainIndices
  return tree, xValError


def TestRun():
  """ testing code """
  examples, attrs, nPossibleVals = randomtest.GenRandomExamples(nExamples=200)
  tree, _ = CrossValidationDriver(examples, attrs, nPossibleVals)

  tree.Pickle('save.pkl')

  import copy
  t2 = copy.deepcopy(tree)
  print('t1 == t2', tree == t2)
  l = [tree]
  print('t2 in [tree]', t2 in l, l.index(t2))


if __name__ == '__main__':  # pragma: nocover
  TestRun()
