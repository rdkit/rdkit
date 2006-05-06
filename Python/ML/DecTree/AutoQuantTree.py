#
#  Copyright (C) 2001-2004  greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#
""" 

"""

from Numeric import *
from ML.InfoTheory import QuantTree


def ID3(examples,target,attrs,nPossibleVals):
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

    **Returns**
    
     a DecTree.DecTreeNode with the decision tree

    **NOTE:** This code cannot bootstrap (start from nothing...)
          use _ID3Boot_ (below) for that.
  """
  varTable = GenVarTable(examples,nPossibleVals,attrs)
  tree=DecTree.DecTreeNode(None,'node')

  # store the total entropy... in case that is interesting
  totEntropy = CalcTotalEntropy(examples,nPossibleVals)
  tree.SetData(totEntropy)
  #tree.SetExamples(examples)

  # the matrix of results for this target:
  tMat = GenVarTable(examples,nPossibleVals,[target])[0] 
  # counts of each result code:
  counts = sum(tMat)  
  nzCounts = nonzero(counts)

  if len(nzCounts) == 1:
    # bottomed out because there is only one result code left
    #  with any counts (i.e. there's only one type of example
    #  left... this is GOOD!).
    res = nzCounts[0]
    tree.SetLabel(res)
    tree.SetName(str(res))
    tree.SetTerminal(1)
  else:
    gains = map(lambda x: entropy.InfoGain(x),varTable)
    if len(attrs) == 0  or max(gains)<1e-8:
      # Bottomed out: no variables left...
      #  We don't really know what to do here, so
      #  use the heuristic of picking the most prevalent
      #  result
      v =  argmax(counts)
      tree.SetLabel(v)
      tree.SetName('%d?'%v)
      tree.SetTerminal(1)
    else:
      # find the variable which gives us the largest information gain
      best = attrs[argmax(gains)]


      # remove that variable from the lists of possible variables
      nextAttrs = attrs[:]
      nextAttrs.remove(best)

      # set some info at this node
      tree.SetName('Var: %d'%best)
      tree.SetLabel(best)
      #tree.SetExamples(examples)
      tree.SetTerminal(0)

      # loop over possible values of the new variable and
      #  build a subtree for each one
      for val in xrange(nPossibleVals[best]):
        nextExamples = []
        for example in examples:
          if example[best] == val:
            nextExamples.append(example)
        if len(nextExamples) == 0:
          # this particular value of the variable has no examples,
          #  so there's not much sense in recursing.
          #  This can (and does) happen.
          v =  argmax(counts)
          tree.AddChild('%d'%v,label=v,data=0.0,isTerminal=1)
        else:
          # recurse
          tree.AddChildNode(ID3(nextExamples,best,nextAttrs,nPossibleVals))
  return tree

def ID3Boot(examples,attrs,nPossibleVals,initialVar=None):
  """ Bootstrapping code for the ID3 algorithm

    see ID3 for descriptions of the arguments

    If _initialVar_ is not set, the algorithm will automatically
     choose the first variable in the tree (the standard greedy
     approach).  Otherwise, _initialVar_ will be used as the first
     split.
     
  """
  totEntropy = CalcTotalEntropy(examples,nPossibleVals)
  varTable = GenVarTable(examples,nPossibleVals,attrs)

  tree=DecTree.DecTreeNode(None,'node')
  #tree.SetExamples(examples)
  tree._nResultCodes = nPossibleVals[-1]

  # <perl>you've got to love any language which will let you
  # do this much work in a single line :-)</perl>
  if initialVar is None:
    best = attrs[argmax(map(lambda x: entropy.InfoGain(x),varTable))]
  else:
    best = initialVar

  tree.SetName('Var: %d'%best)
  tree.SetData(totEntropy)
  tree.SetLabel(best)
  tree.SetTerminal(0)
  nextAttrs = attrs[:]
  nextAttrs.remove(best)
  for val in xrange(nPossibleVals[best]):
    nextExamples = []
    for example in examples:
      if example[best] == val:
        nextExamples.append(example)

    tree.AddChildNode(ID3(nextExamples,best,nextAttrs,nPossibleVals))
  return tree
  
        
def TestMultiTree():
  """Testing code for generating trees with more than 2 possible results

  """
  from ML.Data import MLData
  print 'Testing MultiValue Tree Construction'
  examples = [[0,1,0,0],
              [0,0,0,1],
              [0,0,1,2],
              [0,1,1,2],
              [1,0,0,2],
              [1,0,1,2],
              [1,1,0,2],
              [1,1,1,0]
              ]
  data = MLData.MLQuantDataSet(examples)
  attrs = range(0,data.GetNVars())
  t1 = ID3Boot(data.GetAllData(),attrs,data.GetNPossibleVals())
  #t1.Print()
  t1.Pickle('multi.pkl')

  print 'Testing Pickle Load'
  import cPickle
  f = open('regress/MultiTreeRes.pkl','r')
  t2 = cPickle.load(f)
  print 'Testing Correctness'
  assert t1 == t2,'Equality Test Failed'

  print 'All Tests Passed!'
  
def TestTree():
  """Testing code for trees with a single possible result

  """
  from ML.Data import MLData

  print 'Testing Tree Construction'
  examples = [[0,0,0,0,0],
              [0,0,0,1,0],
              [1,0,0,0,1],
              [2,1,0,0,1],
              [2,2,1,0,1],
              [2,2,1,1,0],
              [1,2,1,1,1],
              [0,1,0,0,0],
              [0,2,1,0,1],
              [2,1,1,0,1],
              [0,1,1,1,1],
              [1,1,0,1,1],
              [1,0,1,0,1],
              [2,1,0,1,0]
              ]

  data = MLData.MLQuantDataSet(examples)
  attrs = range(0,data.GetNVars())
  t1 = ID3Boot(data.GetAllData(),attrs,data.GetNPossibleVals())

  print 'Testing Tree Validity'
  t2 = DecTree.DecTreeNode(None,'Var: 0',0)
  
  c = DecTree.DecTreeNode(t2,'Var: 2',2)
  t2.AddChildNode(c)
  c2 = DecTree.DecTreeNode(c,'0',0,isTerminal=1)
  c.AddChildNode(c2)
  c2 = DecTree.DecTreeNode(c,'1',1,isTerminal=1)
  c.AddChildNode(c2)
  
  c = DecTree.DecTreeNode(t2,'1',1,isTerminal=1)
  t2.AddChildNode(c)

  c = DecTree.DecTreeNode(t2,'Var: 3',3)
  t2.AddChildNode(c)
  c2 = DecTree.DecTreeNode(c,'1',1,isTerminal=1)
  c.AddChildNode(c2)
  c2 = DecTree.DecTreeNode(c,'0',0,isTerminal=1)
  c.AddChildNode(c2)

  assert t2==t1,'Trees do not match'
  #print 'Testing Printing'
  #t1.Print(showData=1)
  print 'Testing Pickle'
  t1.Pickle('save.pkl')
  print 'Classification Tests:'
  assert t1.ClassifyExample(examples[0])==examples[0][-1],'Example 0 misclassified'
  assert t1.ClassifyExample(examples[1])==examples[1][-1],'Example 1 misclassified'
  assert t1.ClassifyExample(examples[6])==examples[6][-1],'Example 6 misclassified'

  print 'Testing Copy'
  import copy
  t2 = copy.deepcopy(t1)
  assert t1==t2,'copy failed'
  print 'Testing Set Membership'
  l = [t1]
  assert t2 in l,'Set Membership failed'
  #print 't2 in [t1]', t2 in l, 'index:',l.index(t2)
  print 'All tests passed!'

def TestNamedTree():
  """ testing code for named trees

  """
  from ML.Data import MLData
  print 'Testing Named Tree Construction'
  examples = [[0,1,0,0],
              [0,0,0,1],
              [0,0,1,2],
              [0,1,1,2],
              [1,0,0,2],
              [1,0,1,2],
              [1,1,0,2],
              [1,1,1,0]
              ]
  names = ['ex1','ex2','ex3','ex4','ex5','ex6','ex7','ex8']
  data = MLData.MLQuantDataSet(examples,ptNames=names)
  attrs = range(1,data.GetNVars()+1)
  t1 = ID3Boot(data.GetNamedData(),attrs,[0]+data.GetNPossibleVals())
  print 'All tests passed!'


if __name__ == "__main__":
  TestTree()
  TestMultiTree()
  TestNamedTree()
  
