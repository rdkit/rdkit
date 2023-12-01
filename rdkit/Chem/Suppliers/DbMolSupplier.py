#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""
Supplies a class for working with molecules from databases
"""
import sys

from rdkit import Chem
from rdkit.Chem.Suppliers.MolSupplier import MolSupplier


def warning(msg, dest=sys.stderr):
  dest.write(msg)


class DbMolSupplier(MolSupplier):
  """
    new molecules come back with all additional fields from the
    database set in a "_fieldsFromDb" data member

  """

  def __init__(self, dbResults, molColumnFormats={
    'SMILES': 'SMI',
    'SMI': 'SMI',
    'MOLPKL': 'PKL'
  }, nameCol='', transformFunc=None, **kwargs):
    """

      DbResults should be a subclass of Dbase.DbResultSet.DbResultBase

    """
    self._data = dbResults
    self._colNames = [x.upper() for x in self._data.GetColumnNames()]
    nameCol = nameCol.upper()
    self.molCol = -1
    self.transformFunc = transformFunc
    try:
      self.nameCol = self._colNames.index(nameCol)
    except ValueError:
      self.nameCol = -1
    for name in molColumnFormats:
      name = name.upper()
      try:
        idx = self._colNames.index(name)
      except ValueError:
        pass
      else:
        self.molCol = idx
        self.molFmt = molColumnFormats[name]
        break
    if self.molCol < 0:
      raise ValueError('DbResultSet has no recognizable molecule column')
    del self._colNames[self.molCol]
    self._colNames = tuple(self._colNames)
    self._numProcessed = 0

  def GetColumnNames(self):
    return self._colNames

  def _BuildMol(self, data):
    data = list(data)
    molD = data[self.molCol]
    del data[self.molCol]
    self._numProcessed += 1
    try:
      if self.molFmt == 'SMI':
        newM = Chem.MolFromSmiles(str(molD))
        if not newM:
          warning('Problems processing mol %d, smiles: %s\n' % (self._numProcessed, molD))
      elif self.molFmt == 'PKL':
        newM = Chem.Mol(str(molD))
    except Exception:
      import traceback
      traceback.print_exc()
      newM = None
    else:
      if newM and self.transformFunc:
        try:
          newM = self.transformFunc(newM, data)
        except Exception:
          import traceback
          traceback.print_exc()
          newM = None
      if newM:
        newM._fieldsFromDb = data
        nFields = len(data)
        for i in range(nFields):
          newM.SetProp(self._colNames[i], str(data[i]))
        if self.nameCol >= 0:
          newM.SetProp('_Name', str(data[self.nameCol]))
    return newM


class ForwardDbMolSupplier(DbMolSupplier):
  """ DbMol supplier supporting only forward iteration


    new molecules come back with all additional fields from the
    database set in a "_fieldsFromDb" data member

  """

  def __init__(self, dbResults, **kwargs):
    """

      DbResults should be an iterator for Dbase.DbResultSet.DbResultBase

    """
    DbMolSupplier.__init__(self, dbResults, **kwargs)
    self.Reset()

  def Reset(self):
    self._dataIter = iter(self._data)

  def NextMol(self):
    """

      NOTE: this has side effects

    """
    try:
      d = self._dataIter.next()
    except StopIteration:
      d = None
    if d is not None:
      newM = self._BuildMol(d)
    else:
      newM = None

    return newM


class RandomAccessDbMolSupplier(DbMolSupplier):

  def __init__(self, dbResults, **kwargs):
    """

      DbResults should be a Dbase.DbResultSet.RandomAccessDbResultSet

    """
    DbMolSupplier.__init__(self, dbResults, **kwargs)
    self._pos = -1

  def __len__(self):
    return len(self._data)

  def __getitem__(self, idx):
    newD = self._data[idx]
    return self._BuildMol(newD)

  def Reset(self):
    self._pos = -1

  def NextMol(self):
    self._pos += 1
    res = None
    if self._pos < len(self):
      res = self[self._pos]
    return res
