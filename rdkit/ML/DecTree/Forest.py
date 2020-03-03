#
#  Copyright (C) 2000-2008  greg Landrum
#
""" code for dealing with forests (collections) of decision trees

**NOTE** This code should be obsolete now that ML.Composite.Composite is up and running.

"""


import numpy

from rdkit.ML.DecTree import CrossValidate, PruneTree
import pickle


class Forest(object):
  """a forest of unique decision trees.

    adding an existing tree just results in its count field being incremented
        and the errors being averaged.

    typical usage:

      1) grow the forest with AddTree until happy with it

      2) call AverageErrors to calculate the average error values

      3) call SortTrees to put things in order by either error or count

  """

  def MakeHistogram(self):
    """ creates a histogram of error/count pairs

    """
    nExamples = len(self.treeList)
    histo = []
    i = 1
    lastErr = self.errList[0]
    countHere = self.countList[0]
    eps = 0.001
    while i < nExamples:
      if self.errList[i] - lastErr > eps:
        histo.append((lastErr, countHere))
        lastErr = self.errList[i]
        countHere = self.countList[i]
      else:
        countHere = countHere + self.countList[i]
      i = i + 1

    return histo

  def CollectVotes(self, example):
    """ collects votes across every member of the forest for the given example

      **Returns**

        a list of the results

    """
    nTrees = len(self.treeList)
    votes = [0] * nTrees
    for i in range(nTrees):
      votes[i] = self.treeList[i].ClassifyExample(example)
    return votes

  def ClassifyExample(self, example):
    """ classifies the given example using the entire forest

      **returns** a result and a measure of confidence in it.

      **FIX:** statistics sucks... I'm not seeing an obvious way to get
           the confidence intervals.  For that matter, I'm not seeing
           an unobvious way.

           For now, this is just treated as a voting problem with the confidence
           measure being the percent of trees which voted for the winning result.
    """
    self.treeVotes = self.CollectVotes(example)
    votes = [0] * len(self._nPossible)
    for i in range(len(self.treeList)):
      res = self.treeVotes[i]
      votes[res] = votes[res] + self.countList[i]

    totVotes = sum(votes)
    res = numpy.argmax(votes)
    # print 'v:',res,votes,totVotes
    return res, float(votes[res]) / float(totVotes)

  def GetVoteDetails(self):
    """ Returns the details of the last vote the forest conducted

      this will be an empty list if no voting has yet been done

    """
    return self.treeVotes

  def Grow(self, examples, attrs, nPossibleVals, nTries=10, pruneIt=0, lessGreedy=0):
    """ Grows the forest by adding trees

     **Arguments**

      - examples: the examples to be used for training

      - attrs: a list of the attributes to be used in training

      - nPossibleVals: a list with the number of possible values each variable
        (as well as the result) can take on

      - nTries: the number of new trees to add

      - pruneIt: a toggle for whether or not the tree should be pruned

      - lessGreedy: toggles the use of a less greedy construction algorithm where
        each possible tree root is used.  The best tree from each step is actually
        added to the forest.

    """
    self._nPossible = nPossibleVals
    for i in range(nTries):
      tree, frac = CrossValidate.CrossValidationDriver(examples, attrs, nPossibleVals, silent=1,
                                                       calcTotalError=1, lessGreedy=lessGreedy)
      if pruneIt:
        tree, frac2 = PruneTree.PruneTree(tree, tree.GetTrainingExamples(), tree.GetTestExamples(),
                                          minimizeTestErrorOnly=0)
        print('prune: ', frac, frac2)
        frac = frac2
      self.AddTree(tree, frac)
      if i % (nTries / 10) == 0:
        print('Cycle: % 4d' % (i))

  def Pickle(self, fileName='foo.pkl'):
    """ Writes this forest off to a file so that it can be easily loaded later

     **Arguments**

       fileName is the name of the file to be written

    """
    pFile = open(fileName, 'wb+')
    pickle.dump(self, pFile, 1)
    pFile.close()

  def AddTree(self, tree, error):
    """ Adds a tree to the forest

    If an identical tree is already present, its count is incremented

    **Arguments**

      - tree: the new tree

      - error: its error value

    **NOTE:** the errList is run as an accumulator,
        you probably want to call AverageErrors after finishing the forest

    """
    if tree in self.treeList:
      idx = self.treeList.index(tree)
      self.errList[idx] = self.errList[idx] + error
      self.countList[idx] = self.countList[idx] + 1
    else:
      self.treeList.append(tree)
      self.errList.append(error)
      self.countList.append(1)

  def AverageErrors(self):
    """ convert summed error to average error

      This does the conversion in place
    """
    self.errList = [x / y for x, y in zip(self.errList, self.countList)]

  def SortTrees(self, sortOnError=1):
    """ sorts the list of trees

      **Arguments**

        sortOnError: toggles sorting on the trees' errors rather than their counts

    """
    if sortOnError:
      order = numpy.argsort(self.errList)
    else:
      order = numpy.argsort(self.countList)

    # these elaborate contortions are required because, at the time this
    #  code was written, Numeric arrays didn't unpickle so well...
    self.treeList = [self.treeList[x] for x in order]
    self.countList = [self.countList[x] for x in order]
    self.errList = [self.errList[x] for x in order]

  def GetTree(self, i):
    return self.treeList[i]

  def SetTree(self, i, val):
    self.treeList[i] = val

  def GetCount(self, i):
    return self.countList[i]

  def SetCount(self, i, val):
    self.countList[i] = val

  def GetError(self, i):
    return self.errList[i]

  def SetError(self, i, val):
    self.errList[i] = val

  def GetDataTuple(self, i):
    """ returns all relevant data about a particular tree in the forest

      **Arguments**

        i: an integer indicating which tree should be returned

      **Returns**

        a 3-tuple consisting of:

          1) the tree

          2) its count

          3) its error
    """
    return (self.treeList[i], self.countList[i], self.errList[i])

  def SetDataTuple(self, i, tup):
    """ sets all relevant data for a particular tree in the forest

      **Arguments**

        - i: an integer indicating which tree should be returned

        - tup: a 3-tuple consisting of:

          1) the tree

          2) its count

          3) its error
    """
    self.treeList[i], self.countList[i], self.errList[i] = tup

  def GetAllData(self):
    """ Returns everything we know

    **Returns**

      a 3-tuple consisting of:

        1) our list of trees

        2) our list of tree counts

        3) our list of tree errors

    """
    return (self.treeList, self.countList, self.errList)

  def __len__(self):
    """ allows len(forest) to work

    """
    return len(self.treeList)

  def __getitem__(self, which):
    """ allows forest[i] to work.  return the data tuple

    """
    return self.GetDataTuple(which)

  def __str__(self):
    """ allows the forest to show itself as a string

    """
    outStr = 'Forest\n'
    for i in range(len(self.treeList)):
      outStr = (outStr + '  Tree % 4d:  % 5d occurrences  %%% 5.2f average error\n' %
                (i, self.countList[i], 100. * self.errList[i]))
    return outStr

  def __init__(self):
    self.treeList = []
    self.errList = []
    self.countList = []
    self.treeVotes = []


if __name__ == '__main__':
  from rdkit.ML.DecTree import DecTree
  f = Forest()
  n = DecTree.DecTreeNode(None, 'foo')
  f.AddTree(n, 0.5)
  f.AddTree(n, 0.5)
  f.AverageErrors()
  f.SortTrees()
  print(f)
