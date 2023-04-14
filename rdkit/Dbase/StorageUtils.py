# $Id$
#
#  Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Various storage (molecular and otherwise) functionality

"""
from rdkit import RDConfig
from rdkit.Dbase import DbModule


def ValidateRDId(ID):
  """ returns whether or not an RDId is valid

  >>> ValidateRDId('RDCmpd-000-009-9')
  1
  >>> ValidateRDId('RDCmpd-009-000-009-8')
  1
  >>> ValidateRDId('RDCmpd-009-000-109-8')
  0
  >>> ValidateRDId('bogus')
  0

  """
  ID = ID.replace('_', '-')
  splitId = ID.split('-')
  if len(splitId) < 4:
    return 0
  accum = 0
  for entry in splitId[1:-1]:
    for char in entry:
      try:
        v = int(char)
      except ValueError:
        return 0
      accum += v
  crc = int(splitId[-1])
  return accum % 10 == crc


def RDIdToInt(ID, validate=1):
  """ Returns the integer index for a given RDId
  Throws a ValueError on error

  >>> RDIdToInt('RDCmpd-000-009-9')
  9
  >>> RDIdToInt('RDCmpd-009-000-009-8')
  9000009
  >>> RDIdToInt('RDData_000_009_9')
  9
  >>> try:
  ...   RDIdToInt('RDCmpd-009-000-109-8')
  ... except ValueError:
  ...   print('ok')
  ... else:
  ...   print('failed')
  ok
  >>> try:
  ...   RDIdToInt('bogus')
  ... except ValueError:
  ...   print('ok')
  ... else:
  ...   print('failed')
  ok

  """
  if validate and not ValidateRDId(ID):
    raise ValueError("Bad RD Id")
  ID = ID.replace('_', '-')
  terms = ID.split('-')[1:-1]
  res = 0
  factor = 1
  terms.reverse()
  for term in terms:
    res += factor * int(term)
    factor *= 1000
  return res


def IndexToRDId(idx, leadText='RDCmpd'):
  """ Converts an integer index into an RDId

  The format of the ID is:
    leadText-xxx-xxx-xxx-y
  The number blocks are zero padded and the final digit (y)
  is a checksum:

  >>> str(IndexToRDId(9))
  'RDCmpd-000-009-9'
  >>> str(IndexToRDId(9009))
  'RDCmpd-009-009-8'

  A millions block is included if it's nonzero:

  >>> str(IndexToRDId(9000009))
  'RDCmpd-009-000-009-8'

  The text at the beginning can be altered:

  >>> str(IndexToRDId(9,leadText='RDAlt'))
  'RDAlt-000-009-9'

  Negative indices are errors:

  >>> try:
  ...   IndexToRDId(-1)
  ... except ValueError:
  ...   print('ok')
  ... else:
  ...   print('failed')
  ok

  """
  if idx < 0:
    raise ValueError('indices must be >= zero')

  res = leadText + '-'
  tmpIdx = idx
  if idx >= 1e6:
    res += '%03d-' % (idx // 1e6)
    tmpIdx = idx % int(1e6)
  if tmpIdx < 1000:
    res += '000-'
  else:
    res += '%03d-' % (tmpIdx // 1000)
    tmpIdx = tmpIdx % 1000

  res += '%03d-' % (tmpIdx)
  accum = 0
  txt = str(idx)
  for char in txt:
    accum += int(char)

  res += str(accum % 10)
  return res


def GetNextId(conn, table, idColName='Id'):
  """ returns the next available Id in the database

  see RegisterItem for testing/documentation

  """
  vals = conn.GetData(table=table, fields=idColName)
  maxVal = 0
  for val in vals:
    val = RDIdToInt(val[0], validate=0)
    if val > maxVal:
      maxVal = val
  maxVal += 1
  return maxVal


def GetNextRDId(conn, table, idColName='Id', leadText=''):
  """ returns the next available RDId in the database

  see RegisterItem for testing/documentation

  """
  if not leadText:
    val = conn.GetData(table=table, fields=idColName)[0][0]
    val = val.replace('_', '-')
    leadText = val.split('-')[0]

  ID = GetNextId(conn, table, idColName=idColName)
  return IndexToRDId(ID, leadText=leadText)


def RegisterItem(conn, table, value, columnName, data=None, id='', idColName='Id',
                 leadText='RDCmpd'):
  """
  >>> from rdkit.Dbase.DbConnection import DbConnect
  >>> conn = DbConnect(tempDbName)
  >>> tblName = 'StorageTest'
  >>> conn.AddTable(tblName,'id varchar(32) not null primary key,label varchar(40),val int')
  >>> RegisterItem(conn,tblName,'label1','label',['label1',1])==(1, 'RDCmpd-000-001-1')
  True
  >>> RegisterItem(conn,tblName,'label2','label',['label2',1])==(1, 'RDCmpd-000-002-2')
  True
  >>> RegisterItem(conn,tblName,'label1','label',['label1',1])==(0, 'RDCmpd-000-001-1')
  True
  >>> str(GetNextRDId(conn,tblName))
  'RDCmpd-000-003-3'
  >>> tuple(conn.GetData(table=tblName)[0])==('RDCmpd-000-001-1', 'label1', 1)
  True

  It's also possible to provide ids by hand:

  >>> RegisterItem(conn,tblName,'label10','label',['label10',1],
  ...              id='RDCmpd-000-010-1')==(1, 'RDCmpd-000-010-1')
  True
  >>> str(GetNextRDId(conn,tblName))
  'RDCmpd-000-011-2'

  """
  curs = conn.GetCursor()
  query = 'select %s from %s where %s=%s' % (idColName, table, columnName, DbModule.placeHolder)
  curs.execute(query, (value, ))
  tmp = curs.fetchone()
  if tmp:
    return 0, tmp[0]
  ID = id or GetNextRDId(conn, table, idColName=idColName, leadText=leadText)
  if data:
    row = [ID]
    row.extend(data)
    conn.InsertData(table, row)
    conn.Commit()
  return 1, ID


def RegisterItems(conn, table, values, columnName, rows, startId='', idColName='Id',
                  leadText='RDCmpd'):
  """
  """
  if rows and len(rows) != len(values):
    raise ValueError("length mismatch between rows and values")
  nVals = len(values)
  origOrder = {}
  for i, v in enumerate(values):
    origOrder[v] = i

  curs = conn.GetCursor()
  qs = ','.join(DbModule.placeHolder * nVals)
  curs.execute("create temporary table regitemstemp (%(columnName)s)" % locals())
  curs.executemany("insert into regitemstemp values (?)", [(x, ) for x in values])
  query = ('select %(columnName)s,%(idColName)s from %(table)s ' +
           'where %(columnName)s in (select * from regitemstemp)' % locals())
  curs.execute(query)

  dbData = curs.fetchall()
  if dbData and len(dbData) == nVals:
    return 0, [x[1] for x in dbData]

  if not startId:
    startId = GetNextRDId(conn, table, idColName=idColName, leadText=leadText)
    startId = RDIdToInt(startId)
  ids = [None] * nVals
  for val, ID in dbData:
    ids[origOrder[val]] = ID

  rowsToInsert = []
  for i in range(nVals):
    if ids[i] is None:
      ID = startId
      startId += 1
      ID = IndexToRDId(ID, leadText=leadText)
      ids[i] = ID
      if rows:
        row = [ID]
        row.extend(rows[i])
        rowsToInsert.append(row)
  if rowsToInsert:
    nCols = len(rowsToInsert[0])
    qs = ','.join(DbModule.placeHolder * nCols)
    curs.executemany('insert into %(table)s values (%(qs)s)' % locals(), rowsToInsert)
    conn.Commit()
  return len(values) - len(dbData), ids


# ------------------------------------
#
#  doctest boilerplate
#

_roundtripTests = """
>>> ValidateRDId(IndexToRDId(100))
1
>>> ValidateRDId(IndexToRDId(10000,leadText='foo'))
1
>>> indices = [1,100,1000,1000000]
>>> vals = []
>>> for idx in indices:
...   vals.append(RDIdToInt(IndexToRDId(idx)))
>>> vals == indices
1

"""
__test__ = {"roundtrip": _roundtripTests}


def _test():  # pragma: nocover
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"], verbose=True)


if __name__ == '__main__':  # pragma: nocover
  import os
  import shutil
  import sys
  import tempfile
  if RDConfig.useSqlLite:
    tmpf, tempName = tempfile.mkstemp(suffix='sqlt')
    tempDbName = tempName
    shutil.copyfile(RDConfig.RDTestDatabase, tempDbName)
  else:
    tempDbName = '::RDTests'
  failed, tried = _test()
  if RDConfig.useSqlLite and os.path.exists(tempDbName):
    try:
      os.unlink(tempDbName)
    except Exception:
      import traceback
      traceback.print_exc()
  sys.exit(failed)
