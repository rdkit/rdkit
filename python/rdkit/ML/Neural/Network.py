#
#  Copyright (C) 2000-2008  greg Landrum
#
""" Contains the class _Network_ which is used to represent neural nets

**Network Architecture:**

  A tacit assumption in all of this stuff is that we're dealing with
  feedforward networks.

  The network itself is stored as a list of _NetNode_ objects.  The list
  is ordered in the sense that nodes in earlier/later layers than a
  given node are guaranteed to come before/after that node in the list.
  This way we can easily generate the values of each node by moving
  sequentially through the list, we're guaranteed that every input for a
  node has already been filled in.

  Each node stores a list (_inputNodes_) of indices of its inputs in the
  main node list.

"""

import numpy
import random

from rdkit.ML.Neural import NetNode, ActFuncs


# FIX: this class has not been updated to new-style classes
# (RD Issue380) because that would break all of our legacy pickled
# data. Until a solution is found for this breakage, an update is
# impossible.
class Network:
  """ a neural network

  """

  def ConstructRandomWeights(self, minWeight=-1, maxWeight=1):
    """initialize all the weights in the network to random numbers

      **Arguments**

        - minWeight: the minimum value a weight can take

        - maxWeight: the maximum value a weight can take

    """
    for node in self.nodeList:
      inputs = node.GetInputs()
      if inputs:
        weights = [random.uniform(minWeight, maxWeight) for _ in range(len(inputs))]
        node.SetWeights(weights)

  def FullyConnectNodes(self):
    """ Fully connects each layer in the network to the one above it


     **Note**
       this sets the connections, but does not assign weights

    """
    nodeList = list(range(self.numInputNodes))
    nConnections = 0
    for layer in range(self.numHiddenLayers):
      for i in self.layerIndices[layer + 1]:
        self.nodeList[i].SetInputs(nodeList)
        nConnections = nConnections + len(nodeList)
      nodeList = self.layerIndices[layer + 1]

    for i in self.layerIndices[-1]:
      self.nodeList[i].SetInputs(nodeList)
      nConnections = nConnections + len(nodeList)
    self.nConnections = nConnections

  def ConstructNodes(self, nodeCounts, actFunc, actFuncParms):
    """ build an unconnected network and set node counts

      **Arguments**

        - nodeCounts: a list containing the number of nodes to be in each layer.
           the ordering is:
            (nInput,nHidden1,nHidden2, ... , nHiddenN, nOutput)

    """
    self.nodeCounts = nodeCounts
    self.numInputNodes = nodeCounts[0]
    self.numOutputNodes = nodeCounts[-1]
    self.numHiddenLayers = len(nodeCounts) - 2
    self.numInHidden = [None] * self.numHiddenLayers
    for i in range(self.numHiddenLayers):
      self.numInHidden[i] = nodeCounts[i + 1]

    numNodes = sum(self.nodeCounts)
    self.nodeList = [None] * (numNodes)
    for i in range(numNodes):
      self.nodeList[i] = NetNode.NetNode(i, self.nodeList, actFunc=actFunc,
                                         actFuncParms=actFuncParms)

    self.layerIndices = [None] * len(nodeCounts)
    start = 0
    for i in range(len(nodeCounts)):
      end = start + nodeCounts[i]
      self.layerIndices[i] = list(range(start, end))
      start = end

  def GetInputNodeList(self):
    """ returns a list of input node indices
    """
    return self.layerIndices[0]

  def GetOutputNodeList(self):
    """ returns a list of output node indices
    """
    return self.layerIndices[-1]

  def GetHiddenLayerNodeList(self, which):
    """ returns a list of hidden nodes in the specified layer
    """
    return self.layerIndices[which + 1]

  def GetNumNodes(self):
    """ returns the total number of nodes
    """
    return sum(self.nodeCounts)

  def GetNumHidden(self):
    """ returns the number of hidden layers
    """
    return self.numHiddenLayers

  def GetNode(self, which):
    """ returns a particular node
    """
    return self.nodeList[which]

  def GetAllNodes(self):
    """ returns a list of all nodes
    """
    return self.nodeList

  def ClassifyExample(self, example, appendExamples=0):
    """ classifies a given example and returns the results of the output layer.

      **Arguments**

        - example: the example to be classified

      **NOTE:**

        if the output layer is only one element long,
        a scalar (not a list) will be returned.  This is why a lot of the other
        network code claims to only support single valued outputs.

    """
    if len(example) > self.numInputNodes:
      if len(example) - self.numInputNodes > self.numOutputNodes:
        example = example[1:-self.numOutputNodes]
      else:
        example = example[:-self.numOutputNodes]
    assert len(example) == self.numInputNodes
    totNumNodes = sum(self.nodeCounts)
    results = numpy.zeros(totNumNodes, numpy.float64)
    for i in range(self.numInputNodes):
      results[i] = example[i]
    for i in range(self.numInputNodes, totNumNodes):
      self.nodeList[i].Eval(results)
    self.lastResults = results[:]
    if self.numOutputNodes == 1:
      return results[-1]
    else:
      return results

  def GetLastOutputs(self):
    """ returns the complete list of output layer values from the last time this node
    classified anything"""
    return self.lastResults

  def __str__(self):
    """ provides a string representation of the network """
    outStr = 'Network:\n'
    for i in range(len(self.nodeList)):
      outStr = outStr + '\tnode(% 3d):\n' % i
      outStr = outStr + '\t\tinputs:  %s\n' % (str(self.nodeList[i].GetInputs()))
      outStr = outStr + '\t\tweights: %s\n' % (str(self.nodeList[i].GetWeights()))

    outStr = outStr + 'Total Number of Connections: % 4d' % self.nConnections
    return outStr

  def __init__(self, nodeCounts, nodeConnections=None, actFunc=ActFuncs.Sigmoid, actFuncParms=(),
               weightBounds=1):
    """ Constructor

      This constructs and initializes the network based upon the specified
      node counts.

      A fully connected network with random weights is constructed.

      **Arguments**

        - nodeCounts: a list containing the number of nodes to be in each layer.
           the ordering is:
            (nInput,nHidden1,nHidden2, ... , nHiddenN, nOutput)

        - nodeConnections: I don't know why this is here, but it's optional.  ;-)

        - actFunc: the activation function to be used here.  Must support the API
            of _ActFuncs.ActFunc_.

        - actFuncParms: a tuple of extra arguments to be passed to the activation function
            constructor.

        - weightBounds:  a float which provides the boundary on the random initial weights

    """
    self.ConstructNodes(nodeCounts, actFunc, actFuncParms)
    self.FullyConnectNodes()
    self.ConstructRandomWeights(minWeight=-weightBounds, maxWeight=weightBounds)
    self.lastResults = []


if __name__ == '__main__':  # pragma: nocover

  print('[2,2,2]')
  net = Network([2, 2, 2])
  print(net)

  print('[2,4,1]')
  net = Network([2, 4, 1])
  print(net)

  print('[2,2]')
  net = Network([2, 2])
  print(net)
  inp = [1, 0]
  res = net.ClassifyExample(inp)
  print(inp, '->', res)
  inp = [0, 1]
  res = net.ClassifyExample(inp)
  print(inp, '->', res)
  inp = [.5, .5]
  res = net.ClassifyExample(inp)
  print(inp, '->', res)
