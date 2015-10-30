#
#  Copyright (C) 2003  Greg Landrum and Rational Discovery LLC
#
""" defines class _DbResultSet_ for lazy interactions with Db query results

**Note**

  this uses the Python iterator interface, so you'll need python 2.2 or above.

"""
from __future__ import print_function
import sys
from rdkit.Dbase import DbInfo

class DbResultBase(object):
  def __init__(self,cursor,conn,cmd,removeDups=-1,transform=None,extras=None):
    self.cursor = cursor
    self.removeDups = removeDups
    self.transform = transform
    self.cmd = cmd
    self.conn=conn
    self.extras = extras
    self.Reset()
    self._initColumnNamesAndTypes()
    self.Reset()

  def Reset(self):
    """ implement in subclasses

    """
    try:
      if not self.extras:
        self.cursor.execute(self.cmd)
      else:
        self.cursor.execute(self.cmd,self.extras)
    except Exception:
      sys.stderr.write('the command "%s" generated errors:\n'%(self.cmd))
      import traceback
      traceback.print_exc()


  def __iter__(self):
    self.Reset()
    return self

  def _initColumnNamesAndTypes(self):
    self.colNames = []
    self.colTypes = []
    for cName,cType in DbInfo.GetColumnInfoFromCursor(self.cursor):
      self.colNames.append(cName)
      self.colTypes.append(cType)
    self.colNames = tuple(self.colNames)
    self.colTypes = tuple(self.colTypes)

  def GetColumnNames(self):
    return self.colNames
  def GetColumnTypes(self):
    return self.colTypes
  def GetColumnNamesAndTypes(self):
    res = [None]*len(self.colNames)
    for i in range(len(self.colNames)):
      res[i] = self.colNames[i],self.colTypes[i]
    return tuple(res)
    

class DbResultSet(DbResultBase):
  """ Only supports forward iteration

  """
  def __init__(self,*args,**kwargs):
    DbResultBase.__init__(self,*args,**kwargs)
    self.seen = []
    self._stopped = 0
  def Reset(self):
    self._stopped = 0
    DbResultBase.Reset(self)

  def next(self):
    if self._stopped: 
      raise StopIteration
    r = None
    while r is None:
      r = self.cursor.fetchone()
      if not r:
        self._stopped = 1
        raise StopIteration
      if self.transform is not None:
        r = self.transform(r)
      if self.removeDups>=0:
        v = r[self.removeDups]
        if v in self.seen:
          r = None
        else:
          self.seen.append(v)
    return r

  __next__ = next # PY3
    

class RandomAccessDbResultSet(DbResultBase):
  """ Supports random access

  """
  def __init__(self,*args,**kwargs):
    DbResultBase.__init__(self,*args,**kwargs)
    self.results = []
    self.seen = []
    self._pos = -1

  def Reset(self):
    self._pos = -1
    if self.cursor is not None:
      DbResultBase.Reset(self)

  def _finish(self):
    if self.cursor:
      #sys.stderr.write('_finish:\n')
      r = self.cursor.fetchone()
      while r:
        if self.transform is not None:
          r =  self.transform(r)
        if self.removeDups >=0:
          v = r[self.removeDups]
          if v not in self.seen:
            self.seen.append(v)
            self.results.append(r)
        else:
          self.results.append(r)
        r = self.cursor.fetchone()
      self.cursor = None
  def __getitem__(self,idx):
    if idx < 0: raise IndexError("negative indices not supported")
    if self.cursor is None:
      if len(self.results):
        if idx >= len(self.results):
          raise IndexError('index %d too large (%d max)'%(idx,len(self.results)))
      else:
        raise ValueError('Invalid cursor')
      
    while idx >= len(self.results):
      r = None
      while r is None:
        r = self.cursor.fetchone()
        if not r:
          self.cursor = None
          raise IndexError('index %d too large (%d max)'%(idx,len(self.results)))

        if self.transform is not None:
          r =  self.transform(r)
        if self.removeDups>=0:
          v = r[self.removeDups]
          if v in self.seen:
            r = None
          else:
            self.results.append(r)
            self.seen.append(v)
        else:
          self.results.append(r)
          
    return self.results[idx]
      
  def __len__(self):
    if self.results is None:
      raise ValueError("len() not supported for noMemory Results Sets")
    self._finish()
    return len(self.results)

  def next(self):
    self._pos += 1
    res = None
    if self._pos < len(self):
      res = self.results[self._pos]
    else:
      raise StopIteration
    return res
  
  __next__ = next # PY3


if __name__ == '__main__':
  from rdkit.Dbase.DbConnection import DbConnect
  conn = DbConnect('TEST.GDB')
  curs = conn.GetCursor()
  print('curs:',repr(curs))
  curs.execute('select * from ten_elements')
  set = RandomAccessDbResultSet(curs)
  for i in range(12):
    try:
      val = set[i]
    except IndexError:
      assert i >= 10
    
  print('use len')
  curs = conn.GetCursor()
  curs.execute('select * from ten_elements')
  set = RandomAccessDbResultSet(curs)
  for i in range(len(set)):
    val = set[i]

  print('use iter')
  curs = conn.GetCursor()
  curs.execute('select * from ten_elements')
  set = DbResultSet(curs)
  for thing in set:
    id,val = thing

  print('dups')
  curs = conn.GetCursor()
  curs.execute('select * from ten_elements_dups')
  set = DbResultSet(curs)
  r = []
  for thing in set:
    r.append(thing)
  assert len(r)==20

  curs = conn.GetCursor()
  curs.execute('select * from ten_elements_dups')
  set = DbResultSet(curs,removeDups=0)
  r = []
  for thing in set:
    r.append(thing)
  assert len(r)==10

  curs = conn.GetCursor()
  curs.execute('select * from ten_elements_dups')
  set = RandomAccessDbResultSet(curs,removeDups=0)
  assert len(set)==10
  assert set[0] == (0,11)

  curs = conn.GetCursor()
  curs.execute('select * from ten_elements_dups')
  set = RandomAccessDbResultSet(curs,removeDups=0)
  assert set[0] == (0,11)
  assert set[1] == (2,21)
  assert set[5] == (10,61)

  curs = conn.GetCursor()
  curs.execute('select * from ten_elements_dups')
  set = RandomAccessDbResultSet(curs)
  assert set[0] == (0,11)
  assert set[1] == (0,11)


  
