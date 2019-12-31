# $Id$
#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#   All Rights Reserved
#
""" code for dealing with composite models

For a model to be useable here, it should support the following API:

  - _ClassifyExample(example)_, returns a classification

Other compatibility notes:

 1) To use _Composite.Grow_ there must be some kind of builder
    functionality which returns a 2-tuple containing (model,percent accuracy).

 2) The models should be pickleable

 3) It would be very happy if the models support the __cmp__ method so that
    membership tests used to make sure models are unique work.



"""

import numpy
import pickle
from rdkit.ML.Data import DataUtils


class Composite(object):
  """a composite model


    **Notes**

    - adding a model which is already present just results in its count
       field being incremented and the errors being averaged.

    - typical usage:

       1) grow the composite with AddModel until happy with it

       2) call AverageErrors to calculate the average error values

       3) call SortModels to put things in order by either error or count

    - Composites can support individual models requiring either quantized or
       nonquantized data.  This is done by keeping a set of quantization bounds
       (_QuantBounds_) in the composite and quantizing data passed in when required.
       Quantization bounds can be set and interrogated using the
       _Get/SetQuantBounds()_ methods.  When models are added to the composite,
       it can be indicated whether or not they require quantization.

    - Composites are also capable of extracting relevant variables from longer lists.
      This is accessible using _SetDescriptorNames()_ to register the descriptors about
      which the composite cares and _SetInputOrder()_ to tell the composite what the
      ordering of input vectors will be.  **Note** there is a limitation on this: each
      model needs to take the same set of descriptors as inputs.  This could be changed.

  """

  def __init__(self):
    self.modelList = []
    self.errList = []
    self.countList = []
    self.modelVotes = []
    self.quantBounds = None
    self.nPossibleVals = None
    self.quantizationRequirements = []
    self._descNames = []
    self._mapOrder = None
    self.activityQuant = []

  def SetModelFilterData(self, modelFilterFrac=0.0, modelFilterVal=0.0):
    self._modelFilterFrac = modelFilterFrac
    self._modelFilterVal = modelFilterVal

  def SetDescriptorNames(self, names):
    """ registers the names of the descriptors this composite uses

      **Arguments**

       - names: a list of descriptor names (strings).

      **NOTE**

         the _names_ list is not
         copied, so if you modify it later, the composite itself will also be modified.

    """
    self._descNames = names

  def GetDescriptorNames(self):
    """ returns the names of the descriptors this composite uses

    """
    return self._descNames

  def SetQuantBounds(self, qBounds, nPossible=None):
    """ sets the quantization bounds that the composite will use

      **Arguments**

       - qBounds:  a list of quantization bounds, each quantbound is a
             list of boundaries

       - nPossible:  a list of integers indicating how many possible values
          each descriptor can take on.

      **NOTE**

         - if the two lists are of different lengths, this will assert out

         - neither list is copied, so if you modify it later, the composite
           itself will also be modified.

    """
    if nPossible is not None:
      assert len(qBounds) == len(nPossible), 'qBounds/nPossible mismatch'
    self.quantBounds = qBounds
    self.nPossibleVals = nPossible

  def GetQuantBounds(self):
    """ returns the quantization bounds

     **Returns**

       a 2-tuple consisting of:

         1) the list of quantization bounds

         2) the nPossibleVals list

    """
    return self.quantBounds, self.nPossibleVals

  def GetActivityQuantBounds(self):
    if not hasattr(self, 'activityQuant'):
      self.activityQuant = []
    return self.activityQuant

  def SetActivityQuantBounds(self, bounds):
    self.activityQuant = bounds

  def QuantizeActivity(self, example, activityQuant=None, actCol=-1):
    if activityQuant is None:
      activityQuant = self.activityQuant
    if activityQuant:
      example = example[:]
      act = example[actCol]
      for box in range(len(activityQuant)):
        if act < activityQuant[box]:
          act = box
          break
      else:
        act = box + 1
      example[actCol] = act
    return example

  def QuantizeExample(self, example, quantBounds=None):
    """ quantizes an example

      **Arguments**

       - example: a data point (list, tuple or numpy array)

       - quantBounds:  a list of quantization bounds, each quantbound is a
             list of boundaries.  If this argument is not provided, the composite
             will use its own quantBounds

      **Returns**

        the quantized example as a list

      **Notes**

        - If _example_ is different in length from _quantBounds_, this will
           assert out.

        - This is primarily intended for internal use

    """
    if quantBounds is None:
      quantBounds = self.quantBounds
    assert len(example) == len(quantBounds), 'example/quantBounds mismatch'
    quantExample = [None] * len(example)
    for i in range(len(quantBounds)):
      bounds = quantBounds[i]
      p = example[i]
      if len(bounds):
        for box in range(len(bounds)):
          if p < bounds[box]:
            p = box
            break
        else:
          p = box + 1
      else:
        if i != 0:
          p = int(p)
      quantExample[i] = p
    return quantExample

  def MakeHistogram(self):
    """ creates a histogram of error/count pairs

     **Returns**

       the histogram as a series of (error, count) 2-tuples

    """
    nExamples = len(self.modelList)
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

  def CollectVotes(self, example, quantExample, appendExample=0, onlyModels=None):
    """ collects votes across every member of the composite for the given example

     **Arguments**

       - example: the example to be voted upon

       - quantExample: the quantized form of the example

       - appendExample: toggles saving the example on the models

       - onlyModels: if provided, this should be a sequence of model
         indices. Only the specified models will be used in the
         prediction.

     **Returns**

       a list with a vote from each member

    """
    if not onlyModels:
      onlyModels = list(range(len(self)))

    votes = [-1] * len(self)
    for i in onlyModels:
      if self.quantizationRequirements[i]:
        votes[i] = int(
          round(self.modelList[i].ClassifyExample(quantExample, appendExamples=appendExample)))
      else:
        votes[i] = int(
          round(self.modelList[i].ClassifyExample(example, appendExamples=appendExample)))

    return votes

  def ClassifyExample(self, example, threshold=0, appendExample=0, onlyModels=None):
    """ classifies the given example using the entire composite

      **Arguments**

       - example: the data to be classified

       - threshold:  if this is a number greater than zero, then a
          classification will only be returned if the confidence is
          above _threshold_.  Anything lower is returned as -1.

       - appendExample: toggles saving the example on the models

       - onlyModels: if provided, this should be a sequence of model
         indices. Only the specified models will be used in the
         prediction.

      **Returns**

        a (result,confidence) tuple


      **FIX:**
        statistics sucks... I'm not seeing an obvious way to get
           the confidence intervals.  For that matter, I'm not seeing
           an unobvious way.

        For now, this is just treated as a voting problem with the confidence
        measure being the percent of models which voted for the winning result.

    """
    if self._mapOrder is not None:
      example = self._RemapInput(example)
    if self.GetActivityQuantBounds():
      example = self.QuantizeActivity(example)
    if self.quantBounds is not None and 1 in self.quantizationRequirements:
      quantExample = self.QuantizeExample(example, self.quantBounds)
    else:
      quantExample = []

    if not onlyModels:
      onlyModels = list(range(len(self)))
    self.modelVotes = self.CollectVotes(example, quantExample, appendExample=appendExample,
                                        onlyModels=onlyModels)

    votes = [0] * self.nPossibleVals[-1]
    for i in onlyModels:
      res = self.modelVotes[i]
      votes[res] = votes[res] + self.countList[i]

    totVotes = sum(votes)
    res = numpy.argmax(votes)
    conf = float(votes[res]) / float(totVotes)
    if conf > threshold:
      return res, conf
    else:
      return -1, conf

  def GetVoteDetails(self):
    """ returns the votes from the last classification

      This will be _None_ if nothing has yet be classified
    """
    return self.modelVotes

  def _RemapInput(self, inputVect):
    """ remaps the input so that it matches the expected internal ordering

      **Arguments**

        - inputVect: the input to be reordered

      **Returns**

        - a list with the reordered (and possible shorter) data

      **Note**

        - you must call _SetDescriptorNames()_ and _SetInputOrder()_ for this to work

        - this is primarily intended for internal use

    """
    order = self._mapOrder

    if order is None:
      return inputVect
    remappedInput = [None] * len(order)

    for i in range(len(order) - 1):
      remappedInput[i] = inputVect[order[i]]
    if order[-1] == -1:
      remappedInput[-1] = 0
    else:
      remappedInput[-1] = inputVect[order[-1]]
    return remappedInput

  def GetInputOrder(self):
    """ returns the input order (used in remapping inputs)

    """
    return self._mapOrder

  def SetInputOrder(self, colNames):
    """ sets the input order

      **Arguments**

        - colNames: a list of the names of the data columns that will be passed in

      **Note**

        - you must call _SetDescriptorNames()_ first for this to work

        - if the local descriptor names do not appear in _colNames_, this will
          raise an _IndexError_ exception.
    """
    if type(colNames) != list:
      colNames = list(colNames)
    descs = [x.upper() for x in self.GetDescriptorNames()]
    self._mapOrder = [None] * len(descs)
    colNames = [x.upper() for x in colNames]

    # FIX: I believe that we're safe assuming that field 0
    #  is always the label, and therefore safe to ignore errors,
    #  but this may not be the case
    try:
      self._mapOrder[0] = colNames.index(descs[0])
    except ValueError:
      self._mapOrder[0] = 0

    for i in range(1, len(descs) - 1):
      try:
        self._mapOrder[i] = colNames.index(descs[i])
      except ValueError:
        raise ValueError('cannot find descriptor name: %s in set %s' %
                         (repr(descs[i]), repr(colNames)))
    try:
      self._mapOrder[-1] = colNames.index(descs[-1])
    except ValueError:
      # ok, there's no obvious match for the final column (activity)
      #  We'll take the last one:
      # self._mapOrder[-1] = len(descs)-1
      self._mapOrder[-1] = -1

  def Grow(self, examples, attrs, nPossibleVals, buildDriver, pruner=None, nTries=10, pruneIt=0,
           needsQuantization=1, progressCallback=None, **buildArgs):
    """ Grows the composite

      **Arguments**

       - examples: a list of examples to be used in training

       - attrs: a list of the variables to be used in training

       - nPossibleVals: this is used to provide a list of the number
          of possible values for each variable.  It is used if the
          local quantBounds have not been set (for example for when you
          are working with data which is already quantized).

       - buildDriver: the function to call to build the new models

       - pruner: a function used to "prune" (reduce the complexity of)
          the resulting model.

       - nTries: the number of new models to add

       - pruneIt: toggles whether or not pruning is done

       - needsQuantization: used to indicate whether or not this type of model
          requires quantized data

       - **buildArgs: all other keyword args are passed to _buildDriver_

      **Note**

        - new models are *added* to the existing ones

    """
    silent = buildArgs.get('silent', 0)
    buildArgs['silent'] = 1
    buildArgs['calcTotalError'] = 1

    if self._mapOrder is not None:
      examples = map(self._RemapInput, examples)
    if self.GetActivityQuantBounds():
      for i in range(len(examples)):
        examples[i] = self.QuantizeActivity(examples[i])
        nPossibleVals[-1] = len(self.GetActivityQuantBounds()) + 1
    if self.nPossibleVals is None:
      self.nPossibleVals = nPossibleVals[:]
    if needsQuantization:
      trainExamples = [None] * len(examples)
      nPossibleVals = self.nPossibleVals
      for i in range(len(examples)):
        trainExamples[i] = self.QuantizeExample(examples[i], self.quantBounds)
    else:
      trainExamples = examples

    for i in range(nTries):
      trainSet = None

      if (hasattr(self, '_modelFilterFrac')) and (self._modelFilterFrac != 0):
        trainIdx, _ = DataUtils.FilterData(trainExamples, self._modelFilterVal,
                                           self._modelFilterFrac, -1, indicesOnly=1)
        trainSet = [trainExamples[x] for x in trainIdx]

      else:
        trainSet = trainExamples

      # print("Training model %i with %i out of %i examples"%(i, len(trainSet), len(trainExamples)))
      model, frac = buildDriver(*(trainSet, attrs, nPossibleVals), **buildArgs)
      if pruneIt:
        model, frac2 = pruner(model, model.GetTrainingExamples(), model.GetTestExamples(),
                              minimizeTestErrorOnly=0)
        frac = frac2
      if (hasattr(self, '_modelFilterFrac') and self._modelFilterFrac != 0 and
          hasattr(model, '_trainIndices')):
        # correct the model's training indices:
        trainIndices = [trainIdx[x] for x in model._trainIndices]
        model._trainIndices = trainIndices

      self.AddModel(model, frac, needsQuantization)
      if not silent and (nTries < 10 or i % (nTries / 10) == 0):
        print('Cycle: % 4d' % (i))
      if progressCallback is not None:
        progressCallback(i)

  def ClearModelExamples(self):
    for i in range(len(self)):
      m = self.GetModel(i)
      try:
        m.ClearExamples()
      except AttributeError:
        pass

  def Pickle(self, fileName='foo.pkl', saveExamples=0):
    """ Writes this composite off to a file so that it can be easily loaded later

     **Arguments**

       - fileName: the name of the file to be written

       - saveExamples: if this is zero, the individual models will have
         their stored examples cleared.

    """
    if not saveExamples:
      self.ClearModelExamples()

    pFile = open(fileName, 'wb+')
    pickle.dump(self, pFile, 1)
    pFile.close()

  def AddModel(self, model, error, needsQuantization=1):
    """ Adds a model to the composite

     **Arguments**

       - model: the model to be added

       - error: the model's error

       - needsQuantization: a toggle to indicate whether or not this model
          requires quantized inputs

     **NOTE**

        - this can be used as an alternative to _Grow()_ if you already have
          some models constructed

        - the errList is run as an accumulator,
          you probably want to call _AverageErrors_ after finishing the forest

    """
    if model in self.modelList:
      try:
        idx = self.modelList.index(model)
      except ValueError:
        # FIX: we should never get here, but sometimes we do anyway
        self.modelList.append(model)
        self.errList.append(error)
        self.countList.append(1)
        self.quantizationRequirements.append(needsQuantization)
      else:
        self.errList[idx] = self.errList[idx] + error
        self.countList[idx] = self.countList[idx] + 1
    else:
      self.modelList.append(model)
      self.errList.append(error)
      self.countList.append(1)
      self.quantizationRequirements.append(needsQuantization)

  def AverageErrors(self):
    """ convert local summed error to average error

    """
    self.errList = list(map(lambda x, y: x / y, self.errList, self.countList))

  def SortModels(self, sortOnError=True):
    """ sorts the list of models

      **Arguments**

        sortOnError: toggles sorting on the models' errors rather than their counts


    """
    if sortOnError:
      order = numpy.argsort(self.errList)
    else:
      order = numpy.argsort(self.countList)

    # these elaborate contortions are required because, at the time this
    #  code was written, Numeric arrays didn't unpickle so well...
    # print(order,sortOnError,self.errList,self.countList)
    self.modelList = [self.modelList[x] for x in order]
    self.countList = [self.countList[x] for x in order]
    self.errList = [self.errList[x] for x in order]

  def GetModel(self, i):
    """ returns a particular model

    """
    return self.modelList[i]

  def SetModel(self, i, val):
    """ replaces a particular model

      **Note**

        This is included for the sake of completeness, but you need to be
        *very* careful when you use it.

    """
    self.modelList[i] = val

  def GetCount(self, i):
    """ returns the count of the _i_th model

    """
    return self.countList[i]

  def SetCount(self, i, val):
    """ sets the count of the _i_th model

    """
    self.countList[i] = val

  def GetError(self, i):
    """ returns the error of the _i_th model

    """
    return self.errList[i]

  def SetError(self, i, val):
    """ sets the error of the _i_th model

    """
    self.errList[i] = val

  def GetDataTuple(self, i):
    """ returns all relevant data about a particular model

      **Arguments**

        i: an integer indicating which model should be returned

      **Returns**

        a 3-tuple consisting of:

          1) the model

          2) its count

          3) its error
    """
    return (self.modelList[i], self.countList[i], self.errList[i])

  def SetDataTuple(self, i, tup):
    """ sets all relevant data for a particular tree in the forest

      **Arguments**

        - i: an integer indicating which model should be returned

        - tup: a 3-tuple consisting of:

          1) the model

          2) its count

          3) its error

      **Note**

        This is included for the sake of completeness, but you need to be
        *very* careful when you use it.

    """
    self.modelList[i], self.countList[i], self.errList[i] = tup

  def GetAllData(self):
    """ Returns everything we know

      **Returns**

        a 3-tuple consisting of:

          1) our list of models

          2) our list of model counts

          3) our list of model errors

    """
    return (self.modelList, self.countList, self.errList)

  def __len__(self):
    """ allows len(composite) to work

    """
    return len(self.modelList)

  def __getitem__(self, which):
    """ allows composite[i] to work, returns the data tuple

    """
    return self.GetDataTuple(which)

  def __str__(self):
    """ returns a string representation of the composite

    """
    outStr = 'Composite\n'
    for i in range(len(self.modelList)):
      outStr = (outStr + '  Model %4d:  %5d occurrences  %%%5.2f average error\n' %
                (i, self.countList[i], 100. * self.errList[i]))
    return outStr


if __name__ == '__main__':  # pragma: nocover
  if 0:
    from rdkit.ML.DecTree import DecTree
    c = Composite()
    n = DecTree.DecTreeNode(None, 'foo')
    c.AddModel(n, 0.5)
    c.AddModel(n, 0.5)
    c.AverageErrors()
    c.SortModels()
    print(c)

    qB = [[], [.5, 1, 1.5]]
    exs = [['foo', 0], ['foo', .4], ['foo', .6], ['foo', 1.1], ['foo', 2.0]]
    print('quantBounds:', qB)
    for ex in exs:
      q = c.QuantizeExample(ex, qB)
      print(ex, q)
  else:
    pass
