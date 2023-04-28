#  $Id$
#
#  Copyright (C) 2004 Rational Discovery LLC
#     All Rights Reserved
#

import os.path
import pickle
import sys

from rdkit import RDConfig
from rdkit.VLib.Supply import SupplyNode

if RDConfig.usePgSQL:
  from pyPgSQL import PgSQL as sql

  class _lazyDataSeq:
    """
    These classes are used to speed up (a lot) the process of
    pulling pickled objects from PostgreSQL databases.  Instead of
    having to use all of PgSQL's typechecking, we'll make a lot of
    assumptions about what's coming out of the Db and its layout.
    The results can lead to drastic improvements in performance.
    
    """

    def __init__(self, cursor, cmd, pickleCol=1, depickle=1, klass=None):
      self.cursor = cursor
      self.cmd = cmd
      self._first = 0
      self._pickleCol = pickleCol
      self._depickle = depickle
      self._klass = klass

    def _validate(self):
      curs = self.cursor
      if not curs or \
             curs.closed or \
             curs.conn is None or \
             (curs.res.resultType != sql.RESULT_DQL and curs.closed is None):
        raise ValueError('bad cursor')
      if curs.res.nfields and curs.res.nfields < 2:
        raise ValueError('invalid number of results returned (%d), must be at least 2' %
                         curs.res.nfields)
      desc1 = curs.description[self._pickleCol]
      ftv = desc1[self._pickleCol].value
      if ftv != sql.BINARY:
        raise TypeError('pickle column (%d) of bad type' % self._pickleCol)

    def __iter__(self):
      try:
        self.cursor.execute(self.cmd)
      except Exception:
        import traceback
        traceback.print_exc()
        print('COMMAND:', self.cmd)
        raise
      self._first = 1
      self._validate()
      return self

    def next(self):
      curs = self.cursor
      if not curs or \
             curs.closed or \
             curs.conn is None or \
             curs.res is None or \
             (curs.res.resultType != sql.RESULT_DQL and curs.closed is None):
        raise StopIteration
      if not self._first:
        res = curs.conn.conn.query('fetch 1 from "%s"' % self.cursor.name)

        if res.ntuples == 0:
          raise StopIteration
        else:
          if res.nfields < 2:
            raise ValueError('bad result: %s' % str(res))
          t = [res.getvalue(0, x) for x in range(res.nfields)]
          val = t[self._pickleCol]
      else:
        t = curs.fetchone()
        val = str(t[self._pickleCol])
        self._first = 0
      if self._depickle:
        if not self._klass:
          fp = pickle.loads(val)
        else:
          fp = self._klass(val)
        fields = list(t)
        del fields[self._pickleCol]
        fp._fieldsFromDb = fields
      else:
        fp = list(t)
      return fp

  class _dataSeq(_lazyDataSeq):

    def __init__(self, cursor, cmd, pickleCol=1, depickle=1):
      self.cursor = cursor
      self.cmd = cmd
      self.res = None
      self.rowCount = -1
      self.idx = 0
      self._pickleCol = pickleCol
      self._depickle = depickle

    def __iter__(self):
      self.cursor.execute(self.cmd)
      self._first = self.cursor.fetchone()
      self._validate()
      self.res = self.cursor.conn.conn.query('fetch all from "%s"' % self.cursor.name)
      self.rowCount = self.res.ntuples + 1
      self.idx = 0
      if self.res.nfields < 2:
        raise ValueError('bad query result' % str(res))

      return self

    def next(self):
      if self.idx >= self.rowCount:
        raise StopIteration

      fp = self[self.idx]
      self.idx += 1

      return fp

    def __len__(self):
      return self.rowCount

    def __getitem__(self, idx):
      if self.res is None:
        self.cursor.execute(self.cmd)
        self._first = self.cursor.fetchone()
        self._validate()
        self.res = self.cursor.conn.conn.query('fetch all from "%s"' % self.cursor.name)
        self.rowCount = self.res.ntuples + 1
        self.idx = 0
        if self.res.nfields < 2:
          raise ValueError('bad query result' % str(res))

      if idx < 0:
        idx = self.rowCount + idx
      if idx < 0 or (idx >= 0 and idx >= self.rowCount):
        raise IndexError
      if idx == 0:
        val = str(self._first[self._pickleCol])
        t = list(self._first)
      else:
        val = self.res.getvalue(self.idx - 1, self._pickleCol)
        t = [self.res.getvalue(self.idx - 1, x) for x in range(self.res.nfields)]
      if self._depickle:
        try:
          fp = pickle.loads(val)
        except Exception:
          import logging
          del t[self._pickleCol]
          logging.exception('Depickling failure in row: %s' % str(t))
          raise
        del t[self._pickleCol]
        fp._fieldsFromDb = t
      else:
        fp = t
      return fp
else:
  _dataSeq = None


class DbPickleSupplyNode(SupplyNode):
  """ Supplies pickled objects from a db result set:

  Sample Usage:
    >>> from rdkit.Dbase.DbConnection import DbConnect
  
  """

  def __init__(self, cursor, cmd, binaryCol, **kwargs):
    SupplyNode.__init__(self, **kwargs)
    self._dbResults = dbResults
    self._supplier = DbMolSupplier.RandomAccessDbMolSupplier(self._dbResults, **kwargs)

  def reset(self):
    SupplyNode.reset(self)
    self._supplier.Reset()

  def next(self):
    """

    """
    return self._supplier.next()


def GetNode(dbName, tableName):
  from rdkit.Dbase.DbConnection import DbConnect
  conn = DbConnect(dbName, tableName)
  return DbMolSupplyNode(conn.GetData())


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
