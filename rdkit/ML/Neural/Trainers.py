#
#  Copyright (C) 2000-2008  greg Landrum
#
""" Training algorithms for feed-forward neural nets

  Unless noted otherwise, algorithms and notation are taken from:
  "Artificial Neural Networks: Theory and Applications",
    Dan W. Patterson, Prentice Hall, 1996

"""


import numpy


class Trainer(object):
  """ "virtual base class" for network trainers

  """
  pass


class BackProp(Trainer):
  """implement back propagation (algorithm on pp 153-154 of Patterson)

   I don't *think* that I've made any assumptions about the connectivity of
     the net (i.e. full connectivity between layers is not required).

   **NOTE:** this code is currently making the assumption that the activation
     functions on the nodes in the network are capable of calculating their
     derivatives using only their values (i.e. a DerivFromVal method should
     exist).  This shouldn't be too hard to change.

  """

  def StepUpdate(self, example, net, resVect=None):
    """ does a BackProp step based upon the example

      **Arguments**

        - example: a 2-tuple:
           1) a list of variable values values
           2) a list of result values (targets)

        - net: a _Network_ (or something supporting the same API)

        - resVect: if this is nonzero, then the network is not required to
          classify the _example_

      **Returns**

        the backprop error from _network_ **before the update**

      **Note**

        In case it wasn't blindingly obvious, the weights in _network_ are modified
        in the course of taking a backprop step.

    """
    totNumNodes = net.GetNumNodes()
    if self.oldDeltaW is None:
      self.oldDeltaW = numpy.zeros(totNumNodes, numpy.float64)
    outputNodeList = net.GetOutputNodeList()
    nOutput = len(outputNodeList)
    targetVect = numpy.array(example[-nOutput:], numpy.float64)
    trainVect = example[:-nOutput]
    if resVect is None:
      # classify the example
      net.ClassifyExample(trainVect)
      resVect = net.GetLastOutputs()
    outputs = numpy.take(resVect, outputNodeList)
    errVect = targetVect - outputs

    delta = numpy.zeros(totNumNodes, numpy.float64)
    # start with the output layer
    for i in range(len(outputNodeList)):
      idx = outputNodeList[i]
      node = net.GetNode(idx)
      # the deltas here are easy
      delta[idx] = errVect[i] * node.actFunc.DerivFromVal(resVect[idx])
      # use these results to start working on the deltas of the preceding layer
      inputs = node.GetInputs()
      weights = delta[idx] * node.GetWeights()
      for j in range(len(inputs)):
        idx2 = inputs[j]
        delta[idx2] = delta[idx2] + weights[j]

    # now propagate the deltas backwards
    for layer in range(net.GetNumHidden() - 1, -1, -1):
      nodesInLayer = net.GetHiddenLayerNodeList(layer)
      for idx in nodesInLayer:
        node = net.GetNode(idx)
        # start by finishing off the error term for this guy
        delta[idx] = delta[idx] * node.actFunc.DerivFromVal(resVect[idx])

        # and then propagate our errors to the preceding layer
        if layer != 0:
          inputs = node.GetInputs()
          weights = delta[idx] * node.GetWeights()
          for i in range(len(inputs)):
            idx2 = inputs[i]
            delta[idx2] = delta[idx2] + weights[i]

        # okey dokey... we've now got the deltas for each node, use those
        #  to update the weights (whew!)
    nHidden = net.GetNumHidden()
    for layer in range(0, nHidden + 1):
      if layer == nHidden:
        idxList = net.GetOutputNodeList()
      else:
        idxList = net.GetHiddenLayerNodeList(layer)
      for idx in idxList:
        node = net.GetNode(idx)
        dW = self.speed * delta[idx] * numpy.take(resVect, node.GetInputs())
        newWeights = node.GetWeights() + dW
        node.SetWeights(newWeights)

      # return the RMS error from the OLD network
    return numpy.sqrt(errVect * errVect)[0]

  def TrainOnLine(self, examples, net, maxIts=5000, errTol=0.1, useAvgErr=1, silent=0):
    """ carries out online training of a neural net

      The definition of online training is that the network is updated after
        each example is presented.

      **Arguments**

        - examples: a list of 2-tuple:
           1) a list of variable values values
           2) a list of result values (targets)

        - net: a _Network_ (or something supporting the same API)

        - maxIts: the maximum number of *training epochs* (see below for definition) to be
          run

        - errTol: the tolerance for convergence

        - useAvgErr: if this toggle is nonzero, then the error at each step will be
          divided by the number of training examples for the purposes of checking
          convergence.

        - silent: controls the amount of visual noise produced as this runs.


      **Note**

         a *training epoch* is one complete pass through all the training examples

    """
    nExamples = len(examples)
    converged = 0
    cycle = 0

    while (not converged) and (cycle < maxIts):
      maxErr = 0
      newErr = 0
      # print('bp: ',cycle)
      for example in examples:
        localErr = self.StepUpdate(example, net)
        newErr += localErr
        if localErr > maxErr:
          maxErr = localErr
      if useAvgErr == 1:
        newErr = newErr / nExamples
      else:
        newErr = maxErr
      # print('\t',newErr,errTol)

      if newErr <= errTol:
        converged = 1

#      if cycle % 10 == 0 and not silent:
      if not silent:
        print('epoch %d, error: % 6.4f' % (cycle, newErr))

      cycle = cycle + 1
    if not silent:
      if converged:
        print('Converged after %d epochs.' % cycle)
      else:
        print('NOT Converged after %d epochs.' % cycle)
      print('final error: % 6.4f' % newErr)

  def __init__(self, speed=0.5, momentum=0.7):
    """ Constructor

      **Arguments**

        - speed: the speed parameter for back prop training

        - momentum: the momentum term for back prop training
          *Not currently used*

    """
    self.speed = speed
    self.momentum = momentum
    self.oldDeltaW = None

if __name__ == '__main__':  # pragma: nocover
  from rdkit.ML.Neural import Network

  def testAnd():
    examples = [[[0, 0, 1], [0.1]], [[0, 1, 1], [.1]], [[1, 0, 1], [.1]], [[1, 1, 1], [.9]]]
    net = Network.Network([3, 1])
    t = BackProp()
    t.TrainOnLine(examples, net)
    return net

  def testOr():
    examples = [[[0, 0, 1], [0.1]], [[0, 1, 1], [.9]], [[1, 0, 1], [.9]], [[1, 1, 1], [.9]]]
    net = Network.Network([3, 1])
    t = BackProp()
    t.TrainOnLine(examples, net, maxIts=1000, useAvgErr=0)
    print('classifications:')
    for example in examples:
      res = net.ClassifyExample(example[0])
      print('%f -> %f' % (example[1][0], res))

    return net

  def testXor():
    examples = [[[0, 0, 1], [.1]], [[0, 1, 1], [.9]], [[1, 0, 1], [.9]], [[1, 1, 1], [.1]]]
    net = Network.Network([3, 3, 1])

    t = BackProp(speed=.8)
    t.TrainOnLine(examples, net, errTol=0.2)
    return net

  def testLinear():
    examples = [
      [.1, .1],
      [.2, .2],
      [.3, .3],
      [.4, .4],
      [.8, .8],
    ]
    net = Network.Network([1, 2, 1])
    t = BackProp(speed=.8)
    t.TrainOnLine(examples, net, errTol=0.1, useAvgErr=0)
    print('classifications:')
    for example in examples:
      res = net.ClassifyExample(example[:-1])
      print('%f -> %f' % (example[-1], res))

    return net

  def runProfile(command):
    import random
    random.seed(23)
    import profile
    import pstats
    datFile = '%s.prof.dat' % (command)
    profile.run('%s()' % command, datFile)
    stats = pstats.Stats(datFile)
    stats.strip_dirs()
    stats.sort_stats('time').print_stats()

  if 0:
    net = testXor()
    print('Xor:', net)
    import pickle
    outF = open('xornet.pkl', 'wb+')
    pickle.dump(net, outF)
    outF.close()
  else:
    # runProfile('testLinear')
    net = testLinear()
    # net = testOr()
