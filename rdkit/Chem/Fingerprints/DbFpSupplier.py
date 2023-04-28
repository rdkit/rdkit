#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Supplies a class for working with fingerprints from databases
#DOC

"""
import pickle

from rdkit import DataStructs
from rdkit.VLib.Node import VLibNode


class DbFpSupplier(VLibNode):
  """
      new fps come back with all additional fields from the
      database set in a "_fieldsFromDb" data member

    """

  def __init__(self, dbResults, fpColName='AutoFragmentFp', usePickles=True):
    """

          DbResults should be a subclass of Dbase.DbResultSet.DbResultBase

        """
    VLibNode.__init__(self)
    self._usePickles = usePickles
    self._data = dbResults
    self._fpColName = fpColName.upper()
    self._colNames = [x.upper() for x in self._data.GetColumnNames()]
    if self._fpColName not in self._colNames:
      raise ValueError('fp column name "%s" not found in result set: %s' %
                       (self._fpColName, str(self._colNames)))
    self.fpCol = self._colNames.index(self._fpColName)
    del self._colNames[self.fpCol]
    self._colNames = tuple(self._colNames)
    self._numProcessed = 0

  def GetColumnNames(self):
    return self._colNames

  def _BuildFp(self, data):
    data = list(data)
    pkl = bytes(data[self.fpCol], encoding='Latin1')
    del data[self.fpCol]
    self._numProcessed += 1
    try:
      if self._usePickles:
        newFp = pickle.loads(pkl, encoding='bytes')
      else:
        newFp = DataStructs.ExplicitBitVect(pkl)
    except Exception:
      import traceback
      traceback.print_exc()
      newFp = None
    if newFp:
      newFp._fieldsFromDb = data
    return newFp

  def next(self):
    itm = self.NextItem()
    if itm is None:
      raise StopIteration
    return itm

  __next__ = next  # py3


class ForwardDbFpSupplier(DbFpSupplier):
  """ DbFp supplier supporting only forward iteration

    >>> from rdkit import RDConfig
    >>> from rdkit.Dbase.DbConnection import DbConnect
    >>> fName = RDConfig.RDTestDatabase
    >>> conn = DbConnect(fName,'simple_combined')
    >>> suppl = ForwardDbFpSupplier(conn.GetData())

    we can loop over the supplied fingerprints:
    
    >>> fps = []
    >>> for fp in suppl:
    ...   fps.append(fp)
    >>> len(fps)
    12

    """

  def __init__(self, *args, **kwargs):
    DbFpSupplier.__init__(self, *args, **kwargs)
    self.reset()

  def reset(self):
    DbFpSupplier.reset(self)
    self._dataIter = iter(self._data)

  def NextItem(self):
    """

          NOTE: this has side effects

        """
    d = next(self._dataIter, None)
    if d is None:
      return d
    return self._BuildFp(d)


class RandomAccessDbFpSupplier(DbFpSupplier):
  """ DbFp supplier supporting random access:

  >>> import os.path
  >>> from rdkit import RDConfig
  >>> from rdkit.Dbase.DbConnection import DbConnect
  >>> fName = RDConfig.RDTestDatabase
  >>> conn = DbConnect(fName,'simple_combined')
  >>> suppl = RandomAccessDbFpSupplier(conn.GetData())
  >>> len(suppl)
  12

  we can pull individual fingerprints:

  >>> fp = suppl[5]
  >>> fp.GetNumBits()
  128
  >>> fp.GetNumOnBits()
  54

  a standard loop over the fingerprints:

  >>> fps = []
  >>> for fp in suppl:
  ...   fps.append(fp)
  >>> len(fps)
  12

  or we can use an indexed loop:

  >>> fps = [None] * len(suppl)
  >>> for i in range(len(suppl)):
  ...   fps[i] = suppl[i]
  >>> len(fps)
  12

  """

  def __init__(self, *args, **kwargs):
    DbFpSupplier.__init__(self, *args, **kwargs)
    self.reset()

  def __len__(self):
    return len(self._data)

  def __getitem__(self, idx):
    newD = self._data[idx]
    return self._BuildFp(newD)

  def reset(self):
    self._pos = -1

  def NextItem(self):
    self._pos += 1
    res = None
    if self._pos < len(self):
      res = self[self._pos]
    return res


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
