# $Id$
#
#  Copyright (C) 2001-2008  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#
"""

"""

import numpy
import random
from rdkit.ML.DecTree import QuantTree, ID3
from rdkit.ML.InfoTheory import entropy
from rdkit.ML.Data import Quantize


def FindBest(resCodes, examples, nBoundsPerVar, nPossibleRes, nPossibleVals, attrs, exIndices=None,
             **kwargs):
    bestGain = -1e6
    best = -1
    bestBounds = []

    if exIndices is None:
        exIndices = list(range(len(examples)))

    if not len(exIndices):
        return best, bestGain, bestBounds

    nToTake = kwargs.get('randomDescriptors', 0)
    if nToTake > 0:
        nAttrs = len(attrs)
        if nToTake < nAttrs:
            ids = list(range(nAttrs))
            random.shuffle(ids, random=random.random)
            tmp = [attrs[x] for x in ids[:nToTake]]
            attrs = tmp

    for var in attrs:
        nBounds = nBoundsPerVar[var]
        if nBounds > 0:
            # vTable = map(lambda x,z=var:x[z],examples)
            try:
                vTable = [examples[x][var] for x in exIndices]
            except IndexError:
                print('index error retrieving variable: %d' % var)
                raise
            qBounds, gainHere = Quantize.FindVarMultQuantBounds(
                vTable, nBounds, resCodes, nPossibleRes)
            # print('\tvar:',var,qBounds,gainHere)
        elif nBounds == 0:
            vTable = ID3.GenVarTable((examples[x] for x in exIndices), nPossibleVals, [var])[0]
            gainHere = entropy.InfoGain(vTable)
            qBounds = []
        else:
            gainHere = -1e6
            qBounds = []
        if gainHere > bestGain:
            bestGain = gainHere
            bestBounds = qBounds
            best = var
        elif bestGain == gainHere:
            if len(qBounds) < len(bestBounds):
                best = var
                bestBounds = qBounds
    if best == -1:
        print('best unaltered')
        print('\tattrs:', attrs)
        print('\tnBounds:', numpy.take(nBoundsPerVar, attrs))
        print('\texamples:')
        for example in (examples[x] for x in exIndices):
            print('\t\t', example)

    if 0:
        print('BEST:', len(exIndices), best, bestGain, bestBounds)
        if (len(exIndices) < 10):
            print(len(exIndices), len(resCodes), len(examples))
            exs = [examples[x] for x in exIndices]
            vals = [x[best] for x in exs]
            sortIdx = numpy.argsort(vals)
            sortVals = [exs[x] for x in sortIdx]
            sortResults = [resCodes[x] for x in sortIdx]
            for i in range(len(vals)):
                print('   ', i, ['%.4f' % x for x in sortVals[i][1:-1]], sortResults[i])
    return best, bestGain, bestBounds


def BuildQuantTree(examples, target, attrs, nPossibleVals, nBoundsPerVar, depth=0, maxDepth=-1,
                   exIndices=None, **kwargs):
    """
      **Arguments**

        - examples: a list of lists (nInstances x nVariables+1) of variable
          values + instance values

        - target: an int

        - attrs: a list of ints indicating which variables can be used in the tree

        - nPossibleVals: a list containing the number of possible values of
                     every variable.

        - nBoundsPerVar: the number of bounds to include for each variable

        - depth: (optional) the current depth in the tree

        - maxDepth: (optional) the maximum depth to which the tree
                     will be grown
      **Returns**

       a QuantTree.QuantTreeNode with the decision tree

      **NOTE:** This code cannot bootstrap (start from nothing...)
            use _QuantTreeBoot_ (below) for that.
    """
    tree = QuantTree.QuantTreeNode(None, 'node')
    tree.SetData(-666)
    nPossibleRes = nPossibleVals[-1]

    if exIndices is None:
        exIndices = list(range(len(examples)))

    # counts of each result code:
    resCodes = [int(x[-1]) for x in (examples[y] for y in exIndices)]
    counts = [0] * nPossibleRes
    for res in resCodes:
        counts[res] += 1
    nzCounts = numpy.nonzero(counts)[0]

    if len(nzCounts) == 1:
        # bottomed out because there is only one result code left
        #  with any counts (i.e. there's only one type of example
        #  left... this is GOOD!).
        res = nzCounts[0]
        tree.SetLabel(res)
        tree.SetName(str(res))
        tree.SetTerminal(1)
    elif len(attrs) == 0 or (maxDepth >= 0 and depth > maxDepth):
        # Bottomed out: no variables left or max depth hit
        #  We don't really know what to do here, so
        #  use the heuristic of picking the most prevalent
        #  result
        v = numpy.argmax(counts)
        tree.SetLabel(v)
        tree.SetName('%d?' % v)
        tree.SetTerminal(1)
    else:
        # find the variable which gives us the largest information gain
        best, _, bestBounds = FindBest(resCodes, examples, nBoundsPerVar, nPossibleRes, nPossibleVals,
                                       attrs, exIndices=exIndices, **kwargs)
        # remove that variable from the lists of possible variables
        nextAttrs = attrs[:]
        if not kwargs.get('recycleVars', 0):
            nextAttrs.remove(best)

        # set some info at this node
        tree.SetName('Var: %d' % (best))
        tree.SetLabel(best)
        tree.SetQuantBounds(bestBounds)
        tree.SetTerminal(0)

        # loop over possible values of the new variable and
        #  build a subtree for each one
        indices = exIndices[:]
        if len(bestBounds) > 0:
            for bound in bestBounds:
                nextExamples = []
                for index in indices[:]:
                    ex = examples[index]
                    if ex[best] < bound:
                        nextExamples.append(index)
                        indices.remove(index)

                if len(nextExamples) == 0:
                    # this particular value of the variable has no examples,
                    #  so there's not much sense in recursing.
                    #  This can (and does) happen.
                    v = numpy.argmax(counts)
                    tree.AddChild('%d' % v, label=v, data=0.0, isTerminal=1)
                else:
                    # recurse
                    tree.AddChildNode(
                      BuildQuantTree(examples, best, nextAttrs, nPossibleVals, nBoundsPerVar, depth=depth + 1,
                                     maxDepth=maxDepth, exIndices=nextExamples, **kwargs))
            # add the last points remaining
            nextExamples = []
            for index in indices:
                nextExamples.append(index)
            if len(nextExamples) == 0:
                v = numpy.argmax(counts)
                tree.AddChild('%d' % v, label=v, data=0.0, isTerminal=1)
            else:
                tree.AddChildNode(
                  BuildQuantTree(examples, best, nextAttrs, nPossibleVals, nBoundsPerVar, depth=depth + 1,
                                 maxDepth=maxDepth, exIndices=nextExamples, **kwargs))
        else:
            for val in range(nPossibleVals[best]):
                nextExamples = []
                for idx in exIndices:
                    if examples[idx][best] == val:
                        nextExamples.append(idx)
                if len(nextExamples) == 0:
                    v = numpy.argmax(counts)
                    tree.AddChild('%d' % v, label=v, data=0.0, isTerminal=1)
                else:
                    tree.AddChildNode(
                      BuildQuantTree(examples, best, nextAttrs, nPossibleVals, nBoundsPerVar, depth=depth + 1,
                                     maxDepth=maxDepth, exIndices=nextExamples, **kwargs))
    return tree


def QuantTreeBoot(examples, attrs, nPossibleVals, nBoundsPerVar, initialVar=None, maxDepth=-1,
                  **kwargs):
    """ Bootstrapping code for the QuantTree

      If _initialVar_ is not set, the algorithm will automatically
       choose the first variable in the tree (the standard greedy
       approach).  Otherwise, _initialVar_ will be used as the first
       split.

    """
    attrs = list(attrs)
    for i in range(len(nBoundsPerVar)):
        if nBoundsPerVar[i] == -1 and i in attrs:
            attrs.remove(i)

    tree = QuantTree.QuantTreeNode(None, 'node')
    nPossibleRes = nPossibleVals[-1]
    tree._nResultCodes = nPossibleRes

    resCodes = [int(x[-1]) for x in examples]
    counts = [0] * nPossibleRes
    for res in resCodes:
        counts[res] += 1
    if initialVar is None:
        best, gainHere, qBounds = FindBest(resCodes, examples, nBoundsPerVar, nPossibleRes,
                                           nPossibleVals, attrs, **kwargs)
    else:
        best = initialVar
        if nBoundsPerVar[best] > 0:
            vTable = map(lambda x, z=best: x[z], examples)
            qBounds, gainHere = Quantize.FindVarMultQuantBounds(vTable, nBoundsPerVar[best], resCodes,
                                                                nPossibleRes)
        elif nBoundsPerVar[best] == 0:
            vTable = ID3.GenVarTable(examples, nPossibleVals, [best])[0]
            gainHere = entropy.InfoGain(vTable)
            qBounds = []
        else:
            gainHere = -1e6
            qBounds = []

    tree.SetName('Var: %d' % (best))
    tree.SetData(gainHere)
    tree.SetLabel(best)
    tree.SetTerminal(0)
    tree.SetQuantBounds(qBounds)
    nextAttrs = list(attrs)
    if not kwargs.get('recycleVars', 0):
        nextAttrs.remove(best)

    indices = list(range(len(examples)))
    if len(qBounds) > 0:
        for bound in qBounds:
            nextExamples = []
            for index in list(indices):
                ex = examples[index]
                if ex[best] < bound:
                    nextExamples.append(ex)
                    indices.remove(index)

            if len(nextExamples):
                tree.AddChildNode(
                  BuildQuantTree(nextExamples, best, nextAttrs, nPossibleVals, nBoundsPerVar, depth=1,
                                 maxDepth=maxDepth, **kwargs))
            else:
                v = numpy.argmax(counts)
                tree.AddChild('%d??' % (v), label=v, data=0.0, isTerminal=1)
        # add the last points remaining
        nextExamples = []
        for index in indices:
            nextExamples.append(examples[index])
        if len(nextExamples) != 0:
            tree.AddChildNode(
              BuildQuantTree(nextExamples, best, nextAttrs, nPossibleVals, nBoundsPerVar, depth=1,
                             maxDepth=maxDepth, **kwargs))
        else:
            v = numpy.argmax(counts)
            tree.AddChild('%d??' % (v), label=v, data=0.0, isTerminal=1)
    else:
        for val in range(nPossibleVals[best]):
            nextExamples = []
            for example in examples:
                if example[best] == val:
                    nextExamples.append(example)
            if len(nextExamples) != 0:
                tree.AddChildNode(
                  BuildQuantTree(nextExamples, best, nextAttrs, nPossibleVals, nBoundsPerVar, depth=1,
                                 maxDepth=maxDepth, **kwargs))
            else:
                v = numpy.argmax(counts)
                tree.AddChild('%d??' % (v), label=v, data=0.0, isTerminal=1)
    return tree


def TestTree():
    """ testing code for named trees

    """
    examples1 = [['p1', 0, 1, 0, 0], ['p2', 0, 0, 0, 1], ['p3', 0, 0, 1, 2], ['p4', 0, 1, 1, 2],
                 ['p5', 1, 0, 0, 2], ['p6', 1, 0, 1, 2], ['p7', 1, 1, 0, 2], ['p8', 1, 1, 1, 0]]
    attrs = list(range(1, len(examples1[0]) - 1))
    nPossibleVals = [0, 2, 2, 2, 3]
    t1 = ID3.ID3Boot(examples1, attrs, nPossibleVals, maxDepth=1)
    t1.Print()


def TestQuantTree():  # pragma: nocover
    """ Testing code for named trees

    The created pkl file is required by the unit test code.
    """
    examples1 = [['p1', 0, 1, 0.1, 0], ['p2', 0, 0, 0.1, 1], ['p3', 0, 0, 1.1, 2],
                 ['p4', 0, 1, 1.1, 2], ['p5', 1, 0, 0.1, 2], ['p6', 1, 0, 1.1, 2],
                 ['p7', 1, 1, 0.1, 2], ['p8', 1, 1, 1.1, 0]]
    attrs = list(range(1, len(examples1[0]) - 1))
    nPossibleVals = [0, 2, 2, 0, 3]
    boundsPerVar = [0, 0, 0, 1, 0]

    print('base')
    t1 = QuantTreeBoot(examples1, attrs, nPossibleVals, boundsPerVar)
    t1.Pickle('test_data/QuantTree1.pkl')
    t1.Print()

    print('depth limit')
    t1 = QuantTreeBoot(examples1, attrs, nPossibleVals, boundsPerVar, maxDepth=1)
    t1.Pickle('test_data/QuantTree1.pkl')
    t1.Print()


def TestQuantTree2():  # pragma: nocover
    """ testing code for named trees

    The created pkl file is required by the unit test code.
    """
    examples1 = [['p1', 0.1, 1, 0.1, 0], ['p2', 0.1, 0, 0.1, 1], ['p3', 0.1, 0, 1.1, 2],
                 ['p4', 0.1, 1, 1.1, 2], ['p5', 1.1, 0, 0.1, 2], ['p6', 1.1, 0, 1.1, 2],
                 ['p7', 1.1, 1, 0.1, 2], ['p8', 1.1, 1, 1.1, 0]]
    attrs = list(range(1, len(examples1[0]) - 1))
    nPossibleVals = [0, 0, 2, 0, 3]
    boundsPerVar = [0, 1, 0, 1, 0]

    t1 = QuantTreeBoot(examples1, attrs, nPossibleVals, boundsPerVar)
    t1.Print()
    t1.Pickle('test_data/QuantTree2.pkl')

    for example in examples1:
        print(example, t1.ClassifyExample(example))


if __name__ == "__main__":  # pragma: nocover
    TestTree()
    TestQuantTree()
    # TestQuantTree2()
