#
#  Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
#    All Rights Reserved
#

import random

from rdkit import RDRandom

SeqTypes = (list, tuple)


def SplitIndices(nPts, frac, silent=1, legacy=0, replacement=0):
  """ splits a set of indices into a data set into 2 pieces

    **Arguments**

     - nPts: the total number of points

     - frac: the fraction of the data to be put in the first data set

     - silent: (optional) toggles display of stats

     - legacy: (optional) use the legacy splitting approach

     - replacement: (optional) use selection with replacement

   **Returns**

     a 2-tuple containing the two sets of indices.

   **Notes**

     - the _legacy_ splitting approach uses randomly-generated floats
       and compares them to _frac_.  This is provided for
       backwards-compatibility reasons.

     - the default splitting approach uses a random permutation of
       indices which is split into two parts.

     - selection with replacement can generate duplicates.


  **Usage**:

  We'll start with a set of indices and pick from them using
  the three different approaches:
  >>> from rdkit.ML.Data import DataUtils

  The base approach always returns the same number of compounds in
  each set and has no duplicates:
  >>> DataUtils.InitRandomNumbers((23,42))
  >>> test,train = SplitIndices(10,.5)
  >>> test
  [1, 5, 6, 4, 2]
  >>> train
  [3, 0, 7, 8, 9]

  >>> test,train = SplitIndices(10,.5)
  >>> test
  [5, 2, 9, 8, 7]
  >>> train
  [6, 0, 3, 1, 4]


  The legacy approach can return varying numbers, but still has no
  duplicates.  Note the indices come back ordered:
  >>> DataUtils.InitRandomNumbers((23,42))
  >>> test,train = SplitIndices(10,.5,legacy=1)
  >>> test
  [3, 5, 7, 8, 9]
  >>> train
  [0, 1, 2, 4, 6]

  >>> test,train = SplitIndices(10,.5,legacy=1)
  >>> test
  [0, 1, 2, 3, 5, 8, 9]
  >>> train
  [4, 6, 7]

  The replacement approach returns a fixed number in the training set,
  a variable number in the test set and can contain duplicates in the
  training set.
  >>> DataUtils.InitRandomNumbers((23,42))
  >>> test,train = SplitIndices(10,.5,replacement=1)
  >>> test
  [9, 9, 8, 0, 5]
  >>> train
  [1, 2, 3, 4, 6, 7]
  >>> test,train = SplitIndices(10,.5,replacement=1)
  >>> test
  [4, 5, 1, 1, 4]
  >>> train
  [0, 2, 3, 6, 7, 8, 9]

  """
  if frac < 0. or frac > 1.:
    raise ValueError('frac must be between 0.0 and 1.0 (frac=%f)' % (frac))

  if replacement:
    nTrain = int(nPts * frac)
    resData = [None] * nTrain
    resTest = []
    for i in range(nTrain):
      val = int(RDRandom.random() * nPts)
      if val == nPts:
        val = nPts - 1
      resData[i] = val
    for i in range(nPts):
      if i not in resData:
        resTest.append(i)
  elif legacy:
    resData = []
    resTest = []
    for i in range(nPts):
      val = RDRandom.random()
      if val < frac:
        resData.append(i)
      else:
        resTest.append(i)
  else:
    perm = list(range(nPts))
    RDRandom.shuffle(perm, random=random.random)
    nTrain = int(nPts * frac)

    resData = list(perm[:nTrain])
    resTest = list(perm[nTrain:])

  if not silent:
    print('Training with %d (of %d) points.' % (len(resData), nPts))
    print('\t%d points are in the hold-out set.' % (len(resTest)))
  return resData, resTest


def SplitDataSet(data, frac, silent=0):
  """ splits a data set into two pieces

    **Arguments**

     - data: a list of examples to be split

     - frac: the fraction of the data to be put in the first data set

     - silent: controls the amount of visual noise produced.

   **Returns**

     a 2-tuple containing the two new data sets.

  """
  if frac < 0. or frac > 1.:
    raise ValueError('frac must be between 0.0 and 1.0')

  nOrig = len(data)
  train, test = SplitIndices(nOrig, frac, silent=1)
  resData = [data[x] for x in train]
  resTest = [data[x] for x in test]

  if not silent:
    print('Training with %d (of %d) points.' % (len(resData), nOrig))
    print('\t%d points are in the hold-out set.' % (len(resTest)))
  return resData, resTest


def SplitDbData(conn, fracs, table='', fields='*', where='', join='', labelCol='', useActs=0,
                nActs=2, actCol='', actBounds=[], silent=0):
  """  "splits" a data set held in a DB by returning lists of ids

  **Arguments**:

    - conn: a DbConnect object

    - frac: the split fraction. This can optionally be specified as a
      sequence with a different fraction for each activity value.

    - table,fields,where,join: (optional) SQL query parameters

    - useActs: (optional) toggles splitting based on activities
      (ensuring that a given fraction of each activity class ends
      up in the hold-out set)
      Defaults to 0

    - nActs: (optional) number of possible activity values, only
      used if _useActs_ is nonzero
      Defaults to 2

    - actCol: (optional) name of the activity column
      Defaults to use the last column returned by the query

    - actBounds: (optional) sequence of activity bounds
      (for cases where the activity isn't quantized in the db)
      Defaults to an empty sequence

    - silent: controls the amount of visual noise produced.

  **Usage**:

  Set up the db connection, the simple tables we're using have actives with even
  ids and inactives with odd ids:
  >>> from rdkit.ML.Data import DataUtils
  >>> from rdkit.Dbase.DbConnection import DbConnect
  >>> from rdkit import RDConfig
  >>> conn = DbConnect(RDConfig.RDTestDatabase)

  Pull a set of points from a simple table... take 33% of all points:
  >>> DataUtils.InitRandomNumbers((23,42))
  >>> train,test = SplitDbData(conn,1./3.,'basic_2class')
  >>> [str(x) for x in train]
  ['id-7', 'id-6', 'id-2', 'id-8']

  ...take 50% of actives and 50% of inactives:
  >>> DataUtils.InitRandomNumbers((23,42))
  >>> train,test = SplitDbData(conn,.5,'basic_2class',useActs=1)
  >>> [str(x) for x in train]
  ['id-5', 'id-3', 'id-1', 'id-4', 'id-10', 'id-8']


  Notice how the results came out sorted by activity

  We can be asymmetrical: take 33% of actives and 50% of inactives:
  >>> DataUtils.InitRandomNumbers((23,42))
  >>> train,test = SplitDbData(conn,[.5,1./3.],'basic_2class',useActs=1)
  >>> [str(x) for x in train]
  ['id-5', 'id-3', 'id-1', 'id-4', 'id-10']

  And we can pull from tables with non-quantized activities by providing
  activity quantization bounds:
  >>> DataUtils.InitRandomNumbers((23,42))
  >>> train,test = SplitDbData(conn,.5,'float_2class',useActs=1,actBounds=[1.0])
  >>> [str(x) for x in train]
  ['id-5', 'id-3', 'id-1', 'id-4', 'id-10', 'id-8']

  """
  if not table:
    table = conn.tableName
  if actBounds and len(actBounds) != nActs - 1:
    raise ValueError('activity bounds list length incorrect')
  if useActs:
    if type(fracs) not in SeqTypes:
      fracs = tuple([fracs] * nActs)
    for frac in fracs:
      if frac < 0.0 or frac > 1.0:
        raise ValueError('fractions must be between 0.0 and 1.0')
  else:
    if type(fracs) in SeqTypes:
      frac = fracs[0]
      if frac < 0.0 or frac > 1.0:
        raise ValueError('fractions must be between 0.0 and 1.0')
    else:
      frac = fracs
  # start by getting the name of the ID column:
  colNames = conn.GetColumnNames(table=table, what=fields, join=join)
  idCol = colNames[0]

  if not useActs:
    # get the IDS:
    d = conn.GetData(table=table, fields=idCol, join=join)
    ids = [x[0] for x in d]
    nRes = len(ids)
    train, test = SplitIndices(nRes, frac, silent=1)
    trainPts = [ids[x] for x in train]
    testPts = [ids[x] for x in test]
  else:
    trainPts = []
    testPts = []
    if not actCol:
      actCol = colNames[-1]
    whereBase = where.strip()
    if whereBase.find('where') != 0:
      whereBase = 'where ' + whereBase
    if where:
      whereBase += ' and '
    for act in range(nActs):
      frac = fracs[act]
      if not actBounds:
        whereTxt = whereBase + '%s=%d' % (actCol, act)
      else:
        whereTxt = whereBase
        if act != 0:
          whereTxt += '%s>=%f ' % (actCol, actBounds[act - 1])
        if act < nActs - 1:
          if act != 0:
            whereTxt += 'and '
          whereTxt += '%s<%f' % (actCol, actBounds[act])
      d = conn.GetData(table=table, fields=idCol, join=join, where=whereTxt)
      ids = [x[0] for x in d]
      nRes = len(ids)
      train, test = SplitIndices(nRes, frac, silent=1)
      trainPts.extend([ids[x] for x in train])
      testPts.extend([ids[x] for x in test])

  return trainPts, testPts


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
