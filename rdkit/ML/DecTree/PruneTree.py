#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#
""" Contains functionality for doing tree pruning

"""


import copy

import numpy

from rdkit.ML.DecTree import CrossValidate, DecTree

_verbose = 0


def MaxCount(examples):
    """ given a set of examples, returns the most common result code

     **Arguments**

        examples: a list of examples to be counted

     **Returns**

       the most common result code

    """
    resList = [x[-1] for x in examples]
    maxVal = max(resList)
    counts = [None] * (maxVal + 1)
    for i in range(maxVal + 1):
        counts[i] = sum([x == i for x in resList])

    return numpy.argmax(counts)


def _GetLocalError(node):
    nWrong = 0
    for example in node.GetExamples():
        pred = node.ClassifyExample(example, appendExamples=0)
        if pred != example[-1]:
            nWrong += 1
            # if _verbose: print('------------------>MISS:',example,pred)
    return nWrong


def _Pruner(node, level=0):
    """Recursively finds and removes the nodes whose removals improve classification

       **Arguments**

         - node: the tree to be pruned.  The pruning data should already be contained
           within node (i.e. node.GetExamples() should return the pruning data)

         - level: (optional) the level of recursion, used only in _verbose printing


       **Returns**

          the pruned version of node


       **Notes**

        - This uses a greedy algorithm which basically does a DFS traversal of the tree,
          removing nodes whenever possible.

        - If removing a node does not affect the accuracy, it *will be* removed.  We
          favor smaller trees.

    """
    if _verbose:
        print('  ' * level, '<%d>  ' % level, '>>> Pruner')
    children = node.GetChildren()[:]

    bestTree = copy.deepcopy(node)
    bestErr = 1e6
    #
    # Loop over the children of this node, removing them when doing so
    #  either improves the local error or leaves it unchanged (we're
    #  introducing a bias for simpler trees).
    #
    for i in range(len(children)):
        child = children[i]
        examples = child.GetExamples()
        if _verbose:
            print('  ' * level, '<%d>  ' % level, ' Child:', i, child.GetLabel())
            bestTree.Print()
            print()
        if len(examples):
            if _verbose:
                print('  ' * level, '<%d>  ' % level, '  Examples', len(examples))
            if child.GetTerminal():
                if _verbose:
                    print('  ' * level, '<%d>  ' % level, '    Terminal')
                continue

            if _verbose:
                print('  ' * level, '<%d>  ' % level, '    Nonterminal')

            workTree = copy.deepcopy(bestTree)
            #
            # First recurse on the child (try removing things below it)
            #
            newNode = _Pruner(child, level=level + 1)
            workTree.ReplaceChildIndex(i, newNode)
            tempErr = _GetLocalError(workTree)
            if tempErr <= bestErr:
                bestErr = tempErr
                bestTree = copy.deepcopy(workTree)
                if _verbose:
                    print('  ' * level, '<%d>  ' % level, '>->->->->->')
                    print('  ' * level, '<%d>  ' % level, 'replacing:', i, child.GetLabel())
                    child.Print()
                    print('  ' * level, '<%d>  ' % level, 'with:')
                    newNode.Print()
                    print('  ' * level, '<%d>  ' % level, '<-<-<-<-<-<')
            else:
                workTree.ReplaceChildIndex(i, child)
            #
            # Now try replacing the child entirely
            #
            bestGuess = MaxCount(child.GetExamples())
            newNode = DecTree.DecTreeNode(workTree, 'L:%d' % (
                bestGuess), label=bestGuess, isTerminal=1)
            newNode.SetExamples(child.GetExamples())
            workTree.ReplaceChildIndex(i, newNode)
            if _verbose:
                print('  ' * level, '<%d>  ' % level, 'ATTEMPT:')
                workTree.Print()
            newErr = _GetLocalError(workTree)
            if _verbose:
                print('  ' * level, '<%d>  ' % level, '---> ', newErr, bestErr)
            if newErr <= bestErr:
                bestErr = newErr
                bestTree = copy.deepcopy(workTree)
                if _verbose:
                    print('  ' * level, '<%d>  ' % level, 'PRUNING:')
                    workTree.Print()
            else:
                if _verbose:
                    print('  ' * level, '<%d>  ' % level, 'FAIL')
                # whoops... put the child back in:
                workTree.ReplaceChildIndex(i, child)
        else:
            if _verbose:
                print('  ' * level, '<%d>  ' % level, '  No Examples', len(examples))
            #
            # FIX:  we need to figure out what to do here (nodes that contain
            #   no examples in the testing set).  I can concoct arguments for
            #   leaving them in and for removing them.  At the moment they are
            #   left intact.
            #
            pass

    if _verbose:
        print('  ' * level, '<%d>  ' % level, '<<< out')
    return bestTree


def PruneTree(tree, trainExamples, testExamples, minimizeTestErrorOnly=1):
    """ implements a reduced-error pruning of decision trees

     This algorithm is described on page 69 of Mitchell's book.

     Pruning can be done using just the set of testExamples (the validation set)
     or both the testExamples and the trainExamples by setting minimizeTestErrorOnly
     to 0.

     **Arguments**

       - tree: the initial tree to be pruned

       - trainExamples: the examples used to train the tree

       - testExamples: the examples held out for testing the tree

       - minimizeTestErrorOnly: if this toggle is zero, all examples (i.e.
         _trainExamples_ + _testExamples_ will be used to evaluate the error.

     **Returns**

       a 2-tuple containing:

          1) the best tree

          2) the best error (the one which corresponds to that tree)

    """
    if minimizeTestErrorOnly:
        testSet = testExamples
    else:
        testSet = trainExamples + testExamples

    # remove any stored examples the tree may have
    tree.ClearExamples()

    #
    # screen the test data through the tree so that we end up with the
    #  appropriate points stored at each node of the tree. Results are ignored
    #
    totErr, badEx = CrossValidate.CrossValidate(tree, testSet, appendExamples=1)

    #
    # Prune
    #
    newTree = _Pruner(tree)

    #
    # And recalculate the errors
    #
    totErr, badEx = CrossValidate.CrossValidate(newTree, testSet)
    newTree.SetBadExamples(badEx)

    return newTree, totErr


# -------
#  testing code
# -------
def _testRandom():
    from rdkit.ML.DecTree import randomtest
    #   examples, attrs, nPossibleVals = randomtest.GenRandomExamples(nVars=20, randScale=0.25,
    #                                                                 nExamples=200)
    examples, attrs, nPossibleVals = randomtest.GenRandomExamples(nVars=10, randScale=0.5,
                                                                  nExamples=200)
    tree, frac = CrossValidate.CrossValidationDriver(examples, attrs, nPossibleVals)
    tree.Print()
    tree.Pickle('orig.pkl')
    print('original error is:', frac)

    print('----Pruning')
    newTree, frac2 = PruneTree(tree, tree.GetTrainingExamples(), tree.GetTestExamples())
    newTree.Print()
    print('pruned error is:', frac2)
    newTree.Pickle('prune.pkl')


def _testSpecific():
    from rdkit.ML.DecTree import ID3
    oPts = [
      [0, 0, 1, 0],
      [0, 1, 1, 1],
      [1, 0, 1, 1],
      [1, 1, 0, 0],
      [1, 1, 1, 1],
    ]
    tPts = oPts + [[0, 1, 1, 0], [0, 1, 1, 0]]

    tree = ID3.ID3Boot(oPts, attrs=range(3), nPossibleVals=[2] * 4)
    tree.Print()
    err, _ = CrossValidate.CrossValidate(tree, oPts)
    print('original error:', err)

    err, _ = CrossValidate.CrossValidate(tree, tPts)
    print('original holdout error:', err)
    newTree, frac2 = PruneTree(tree, oPts, tPts)
    newTree.Print()
    print('best error of pruned tree:', frac2)
    err, badEx = CrossValidate.CrossValidate(newTree, tPts)
    print('pruned holdout error is:', err)
    print(badEx)

#   print(len(tree), len(newTree))


def _testChain():
    from rdkit.ML.DecTree import ID3
    oPts = [
      [1, 0, 0, 0, 1],
      [1, 0, 0, 0, 1],
      [1, 0, 0, 0, 1],
      [1, 0, 0, 0, 1],
      [1, 0, 0, 0, 1],
      [1, 0, 0, 0, 1],
      [1, 0, 0, 0, 1],
      [0, 0, 1, 1, 0],
      [0, 0, 1, 1, 0],
      [0, 0, 1, 1, 1],
      [0, 1, 0, 1, 0],
      [0, 1, 0, 1, 0],
      [0, 1, 0, 0, 1],
    ]
    tPts = oPts

    tree = ID3.ID3Boot(oPts, attrs=range(len(oPts[0]) - 1), nPossibleVals=[2] * len(oPts[0]))
    tree.Print()
    err, _ = CrossValidate.CrossValidate(tree, oPts)
    print('original error:', err)

    err, _ = CrossValidate.CrossValidate(tree, tPts)
    print('original holdout error:', err)
    newTree, frac2 = PruneTree(tree, oPts, tPts)
    newTree.Print()
    print('best error of pruned tree:', frac2)
    err, badEx = CrossValidate.CrossValidate(newTree, tPts)
    print('pruned holdout error is:', err)
    print(badEx)


if __name__ == '__main__':  # pragma: nocover
    _verbose = 1
    # _testRandom()
    _testChain()
