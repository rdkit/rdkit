#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#
""" classes to be used to help work with data sets

"""

import copy
import math

import numpy

numericTypes = (int, float)


class MLDataSet(object):
  """ A data set for holding general data (floats, ints, and strings)

     **Note**
       this is intended to be a read-only data structure
       (i.e. after calling the constructor you cannot touch it)
    """

  def __init__(self, data, nVars=None, nPts=None, nPossibleVals=None, qBounds=None, varNames=None,
               ptNames=None, nResults=1):
    """ Constructor

          **Arguments**

            - data: a list of lists containing the data. The data are copied, so don't worry
                  about us overwriting them.

            - nVars: the number of variables

            - nPts: the number of points

            - nPossibleVals: an list containing the number of possible values
                           for each variable (should contain 0 when not relevant)
                           This is _nVars_ long

            - qBounds: a list of lists containing quantization bounds for variables
                     which are to be quantized (note, this class does not quantize
                     the variables itself, it merely stores quantization bounds.
                     an empty sublist indicates no quantization for a given variable
                     This is _nVars_ long

            - varNames: a list of the names of the variables.
                     This is _nVars_ long

            - ptNames: the names (labels) of the individual data points
               This is _nPts_ long

            - nResults: the number of results columns in the data lists.  This is usually
                        1, but can be higher.
        """
    self.data = [x[:] for x in data]
    self.nResults = nResults
    if nVars is None:
      nVars = len(self.data[0]) - self.nResults
    self.nVars = nVars
    if nPts is None:
      nPts = len(data)
    self.nPts = nPts
    if qBounds is None:
      qBounds = [[]] * len(self.data[0])
    self.qBounds = qBounds
    if nPossibleVals is None:
      nPossibleVals = self._CalcNPossible(self.data)
    self.nPossibleVals = nPossibleVals
    if varNames is None:
      varNames = [''] * self.nVars
    self.varNames = varNames
    if ptNames is None:
      ptNames = [''] * self.nPts
    self.ptNames = ptNames

  def _CalcNPossible(self, data):
    """calculates the number of possible values of each variable (where possible)

          **Arguments**

             -data: a list of examples to be used

          **Returns**

             a list of nPossible values for each variable

        """
    nVars = self.GetNVars() + self.nResults
    nPossible = [-1] * nVars
    cols = list(range(nVars))
    for i, bounds in enumerate(self.qBounds):
      if len(bounds) > 0:
        nPossible[i] = len(bounds)
        cols.remove(i)

    for i, pt in enumerate(self.data):
      for col in cols[:]:
        d = pt[col]
        if type(d) in numericTypes:
          if math.floor(d) == d:
            nPossible[col] = max(math.floor(d), nPossible[col])
          else:
            nPossible[col] = -1
            cols.remove(col)
        else:
          nPossible[col] = -1
          cols.remove(col)
    return [int(x) + 1 for x in nPossible]

  def GetNResults(self):
    return self.nResults

  def GetNVars(self):
    return self.nVars

  def GetNPts(self):
    return self.nPts

  def GetNPossibleVals(self):
    return self.nPossibleVals

  def GetQuantBounds(self):
    return self.qBounds

  def __getitem__(self, idx):
    res = [self.ptNames[idx]] + self.data[idx][:]
    return res

  def __setitem__(self, idx, val):
    if len(val) != self.GetNVars() + self.GetNResults() + 1:
      raise ValueError('bad value in assignment')
    self.ptNames[idx] = val[0]
    self.data[idx] = val[1:]
    return val

  def GetNamedData(self):
    """ returns a list of named examples

         **Note**

           a named example is the result of prepending the example
            name to the data list

        """
    res = [None] * self.nPts
    for i in range(self.nPts):
      res[i] = [self.ptNames[i]] + self.data[i][:]
    return res

  def GetAllData(self):
    """ returns a *copy* of the data

        """
    return copy.deepcopy(self.data)

  def GetInputData(self):
    """ returns the input data

         **Note**

           _inputData_ means the examples without their result fields
            (the last _NResults_ entries)

        """
    v = self.GetNResults()
    return [x[:-v] for x in self.data]

  def GetResults(self):
    """ Returns the result fields from each example

        """
    if self.GetNResults() > 1:
      v = self.GetNResults()
      res = [x[-v:] for x in self.data]
    else:
      res = [x[-1] for x in self.data]
    return res

  def GetVarNames(self):
    return self.varNames

  def GetPtNames(self):
    return self.ptNames

  def AddPoint(self, pt):
    self.data.append(pt[1:])
    self.ptNames.append(pt[0])
    self.nPts += 1

  def AddPoints(self, pts, names):
    if len(pts) != len(names):
      raise ValueError("input length mismatch")
    self.data += pts
    self.ptNames += names
    self.nPts = len(self.data)


class MLQuantDataSet(MLDataSet):
  """ a data set for holding quantized data


      **Note**

        this is intended to be a read-only data structure
        (i.e. after calling the constructor you cannot touch it)

      **Big differences to MLDataSet**

        1) data are stored in a numpy array since they are homogenous

        2) results are assumed to be quantized (i.e. no qBounds entry is required)

    """

  def _CalcNPossible(self, data):
    """calculates the number of possible values of each variable

          **Arguments**

             -data: a list of examples to be used

          **Returns**

             a list of nPossible values for each variable

        """
    return [max(x) + 1 for x in numpy.transpose(data)]

  def GetNamedData(self):
    """ returns a list of named examples

         **Note**

           a named example is the result of prepending the example
            name to the data list

        """
    res = [None] * self.nPts
    for i in range(self.nPts):
      res[i] = [self.ptNames[i]] + self.data[i].tolist()
    return res

  def GetAllData(self):
    """ returns a *copy* of the data

        """
    return self.data.tolist()

  def GetInputData(self):
    """ returns the input data

         **Note**

           _inputData_ means the examples without their result fields
            (the last _NResults_ entries)

        """
    return (self.data[:, :-self.nResults]).tolist()

  def GetResults(self):
    """ Returns the result fields from each example

        """
    if self.GetNResults() > 1:
      v = self.GetNResults()
      res = [x[-v:] for x in self.data]
    else:
      res = [x[-1] for x in self.data]
    return res

  def __init__(self, data, nVars=None, nPts=None, nPossibleVals=None, qBounds=None, varNames=None,
               ptNames=None, nResults=1):
    """ Constructor

          **Arguments**

            - data: a list of lists containing the data. The data are copied, so don't worry
                  about us overwriting them.

            - nVars: the number of variables

            - nPts: the number of points

            - nPossibleVals: an list containing the number of possible values
                           for each variable (should contain 0 when not relevant)
                           This is _nVars_ long

            - qBounds: a list of lists containing quantization bounds for variables
                     which are to be quantized (note, this class does not quantize
                     the variables itself, it merely stores quantization bounds.
                     an empty sublist indicates no quantization for a given variable
                     This is _nVars_ long

            - varNames: a list of the names of the variables.
                     This is _nVars_ long

            - ptNames: the names (labels) of the individual data points
               This is _nPts_ long

            - nResults: the number of results columns in the data lists.  This is usually
                        1, but can be higher.
        """
    self.data = numpy.array(data)
    self.nResults = nResults
    if nVars is None:
      nVars = len(data[0]) - self.nResults
    self.nVars = nVars
    if nPts is None:
      nPts = len(data)
    self.nPts = nPts
    if qBounds is None:
      qBounds = [[]] * self.nVars
    self.qBounds = qBounds
    if nPossibleVals is None:
      nPossibleVals = self._CalcNPossible(data)
    self.nPossibleVals = nPossibleVals
    if varNames is None:
      varNames = [''] * self.nVars
    self.varNames = varNames
    if ptNames is None:
      ptNames = [''] * self.nPts
    self.ptNames = ptNames


if __name__ == '__main__':
  from . import DataUtils
  examples = [[0, 0, 0, 0, 0], [0, 0, 0, 1, 0], [1, 0, 0, 0, 1], [2, 1, 0, 0, 1], [2, 2, 1, 0, 1]]
  varNames = ['foo1', 'foo2', 'foo3', 'foo4', 'res']
  ptNames = ['p1', 'p2', 'p3', 'p4', 'p5']
  dataset = MLQuantDataSet(examples, varNames=varNames, ptNames=ptNames)
  DataUtils.WritePickledData('test_data/test.qdat.pkl', dataset)
  print('nVars:', dataset.GetNVars())
  print('nPts:', dataset.GetNPts())
  print('nPoss:', dataset.GetNPossibleVals())
  print('qBounds:', dataset.GetQuantBounds())
  print('data:', dataset.GetAllData())
  print('Input data:', dataset.GetInputData())
  print('results:', dataset.GetResults())

  print('nameddata:', dataset.GetNamedData())

  examples = [
    ['foo', 1, 1.0, 1, 1.1],
    ['foo', 2, 1.0, 1, 2.1],
    ['foo', 3, 1.2, 1.1, 3.1],
    ['foo', 4, 1.0, 1, 4.1],
    ['foo', 5, 1.1, 1, 5.1],
  ]
  qBounds = [[], [], [], [], [2, 4]]
  varNames = ['foo1', 'foo2', 'foo3', 'foo4', 'res']
  ptNames = ['p1', 'p2', 'p3', 'p4', 'p5']
  dataset = MLDataSet(examples, qBounds=qBounds)
  DataUtils.WritePickledData('test_data/test.dat.pkl', dataset)
  print('nVars:', dataset.GetNVars())
  print('nPts:', dataset.GetNPts())
  print('nPoss:', dataset.GetNPossibleVals())
  print('qBounds:', dataset.GetQuantBounds())
  print('data:', dataset.GetAllData())
  print('Input data:', dataset.GetInputData())
  print('results:', dataset.GetResults())

  print('nameddata:', dataset.GetNamedData())
