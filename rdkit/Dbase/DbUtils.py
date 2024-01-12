# $Id$
#
#  Copyright (C) 2000-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" a set of functions for interacting with databases

 When possible, it's probably preferable to use a _DbConnection.DbConnect_ object

"""

import sys
from io import StringIO

from rdkit.Dbase import DbInfo, DbModule
from rdkit.Dbase.DbResultSet import DbResultSet, RandomAccessDbResultSet


def _take(fromL, what):
  """ Given a list fromL, returns an iterator of the elements specified using their
    indices in the list what """
  return [fromL[x] for x in what]


def GetColumns(dBase, table, fieldString, user='sysdba', password='masterkey', join='', cn=None):
  """ gets a set of data from a table

      **Arguments**

       - dBase: database name

       - table: table name

       - fieldString: a string with the names of the fields to be extracted,
          this should be a comma delimited list

       - user and  password:

       - join: a join clause (omit the verb 'join')


      **Returns**

       - a list of the data

    """
  if not cn:
    cn = DbModule.connect(dBase, user, password)
  c = cn.cursor()
  cmd = 'select %s from %s' % (fieldString, table)
  if join:
    if join.strip().find('join') != 0:
      join = 'join %s' % (join)
    cmd += ' ' + join
  c.execute(cmd)
  return c.fetchall()


def GetData(dBase, table, fieldString='*', whereString='', user='sysdba', password='masterkey',
            removeDups=-1, join='', forceList=0, transform=None, randomAccess=1, extras=None,
            cn=None):
  """ a more flexible method to get a set of data from a table

      **Arguments**

       - fields: a string with the names of the fields to be extracted,
            this should be a comma delimited list

       - where: the SQL where clause to be used with the DB query

       - removeDups indicates the column which should be used to screen
          out duplicates.  Only the first appearance of a duplicate will
          be left in the dataset.

      **Returns**

        - a list of the data


      **Notes**

        - EFF: this isn't particularly efficient

    """
  if forceList and (transform is not None):
    raise ValueError('forceList and transform arguments are not compatible')
  if forceList and (not randomAccess):
    raise ValueError('when forceList is set, randomAccess must also be used')
  if removeDups > -1:
    forceList = True

  if not cn:
    cn = DbModule.connect(dBase, user, password)
  c = cn.cursor()
  cmd = 'select %s from %s' % (fieldString, table)
  if join:
    if join.strip().find('join') != 0:
      join = 'join %s' % (join)
    cmd += ' ' + join
  if whereString:
    if whereString.strip().find('where') != 0:
      whereString = 'where %s' % (whereString)
    cmd += ' ' + whereString

  if forceList:
    try:
      if not extras:
        c.execute(cmd)
      else:
        c.execute(cmd, extras)
    except Exception:
      sys.stderr.write('the command "%s" generated errors:\n' % (cmd))
      import traceback
      traceback.print_exc()
      return None
    if transform is not None:
      raise ValueError('forceList and transform arguments are not compatible')
    if not randomAccess:
      raise ValueError('when forceList is set, randomAccess must also be used')
    data = c.fetchall()
    if removeDups >= 0:
      seen = set()
      for entry in data[:]:
        if entry[removeDups] in seen:
          data.remove(entry)
        else:
          seen.add(entry[removeDups])
  else:
    if randomAccess:
      klass = RandomAccessDbResultSet
    else:
      klass = DbResultSet

    data = klass(c, cn, cmd, removeDups=removeDups, transform=transform, extras=extras)

  return data


def DatabaseToText(dBase, table, fields='*', join='', where='', user='sysdba', password='masterkey',
                   delim=',', cn=None):
  """ Pulls the contents of a database and makes a deliminted text file from them

      **Arguments**
        - dBase: the name of the DB file to be used

        - table: the name of the table to query

        - fields: the fields to select with the SQL query

        - join: the join clause of the SQL query
          (e.g. 'join foo on foo.bar=base.bar')

        - where: the where clause of the SQL query
          (e.g. 'where foo = 2' or 'where bar > 17.6')

        - user: the username for DB access

        - password: the password to be used for DB access

      **Returns**

        - the CSV data (as text)

    """
  if len(where) and where.strip().find('where') == -1:
    where = 'where %s' % (where)
  if len(join) and join.strip().find('join') == -1:
    join = 'join %s' % (join)
  sqlCommand = 'select %s from %s %s %s' % (fields, table, join, where)
  if not cn:
    cn = DbModule.connect(dBase, user, password)
  c = cn.cursor()
  c.execute(sqlCommand)
  headers = []
  colsToTake = []
  # the description field of the cursor carries around info about the columns
  #  of the table
  for i in range(len(c.description)):
    item = c.description[i]
    if item[1] not in DbInfo.sqlBinTypes:
      colsToTake.append(i)
      headers.append(item[0])

  lines = []
  lines.append(delim.join(headers))

  # grab the data
  results = c.fetchall()
  for res in results:
    d = _take(res, colsToTake)
    lines.append(delim.join([str(x) for x in d]))

  return '\n'.join(lines)


def TypeFinder(data, nRows, nCols, nullMarker=None):
  """

      finds the types of the columns in _data_

      if nullMarker is not None, elements of the data table which are
        equal to nullMarker will not count towards setting the type of
        their columns.

    """
  priorities = {float: 3, int: 2, str: 1, -1: -1}
  res = [None] * nCols
  for col in range(nCols):
    typeHere = [-1, 1]
    for row in range(nRows):
      d = data[row][col]
      if d is None:
        continue
      locType = type(d)
      if locType != float and locType != int:
        locType = str
        try:
          d = str(d)
        except UnicodeError as msg:
          print('cannot convert text from row %d col %d to a string' % (row + 2, col))
          print('\t>%s' % (repr(d)))
          raise UnicodeError(msg)
      else:
        typeHere[1] = max(typeHere[1], len(str(d)))
      if isinstance(d, str):
        if nullMarker is None or d != nullMarker:
          l = max(len(d), typeHere[1])
          typeHere = [str, l]
      else:
        try:
          fD = float(int(d))
        except OverflowError:
          locType = float
        else:
          if fD == d:
            locType = int
        if not isinstance(typeHere[0], str) and priorities[locType] > priorities[typeHere[0]]:
          typeHere[0] = locType
    res[col] = typeHere
  return res


def _AdjustColHeadings(colHeadings, maxColLabelLen):
  """ *For Internal Use*

      removes illegal characters from column headings
      and truncates those which are too long.

    """
  for i in range(len(colHeadings)):
    # replace unallowed characters and strip extra white space
    colHeadings[i] = colHeadings[i].strip()
    colHeadings[i] = colHeadings[i].replace(' ', '_')
    colHeadings[i] = colHeadings[i].replace('-', '_')
    colHeadings[i] = colHeadings[i].replace('.', '_')

    if len(colHeadings[i]) > maxColLabelLen:
      # interbase (at least) has a limit on the maximum length of a column name
      newHead = colHeadings[i].replace('_', '')
      newHead = newHead[:maxColLabelLen]
      print('\tHeading %s too long, changed to %s' % (colHeadings[i], newHead))
      colHeadings[i] = newHead
  return colHeadings


def GetTypeStrings(colHeadings, colTypes, keyCol=None):
  """  returns a list of SQL type strings
    """
  typeStrs = []
  for i in range(len(colTypes)):
    typ = colTypes[i]
    if typ[0] == float:
      typeStrs.append('%s double precision' % colHeadings[i])
    elif typ[0] == int:
      typeStrs.append('%s integer' % colHeadings[i])
    else:
      typeStrs.append('%s varchar(%d)' % (colHeadings[i], typ[1]))
    if colHeadings[i] == keyCol:
      typeStrs[-1] = '%s not null primary key' % (typeStrs[-1])
  return typeStrs


def _insertBlock(conn, sqlStr, block, silent=False):
  try:
    conn.cursor().executemany(sqlStr, block)
  except Exception:
    res = 0
    conn.commit()
    for row in block:
      try:
        conn.cursor().execute(sqlStr, tuple(row))
        res += 1
      except Exception:
        if not silent:
          import traceback
          traceback.print_exc()
          print('insert failed:', sqlStr)
          print('\t', repr(row))
      else:
        conn.commit()
  else:
    res = len(block)
  return res


def _AddDataToDb(dBase, table, user, password, colDefs, colTypes, data, nullMarker=None,
                 blockSize=100, cn=None):
  """ *For Internal Use*

      (drops and) creates a table and then inserts the values

    """
  if not cn:
    cn = DbModule.connect(dBase, user, password)
  c = cn.cursor()
  try:
    c.execute('drop table %s' % (table))
  except Exception:
    print('cannot drop table %s' % (table))
  try:
    sqlStr = 'create table %s (%s)' % (table, colDefs)
    c.execute(sqlStr)
  except Exception:
    print('create table failed: ', sqlStr)
    print('here is the exception:')
    import traceback
    traceback.print_exc()
    return
  cn.commit()
  c = None

  block = []
  entryTxt = [DbModule.placeHolder] * len(data[0])
  dStr = ','.join(entryTxt)
  sqlStr = 'insert into %s values (%s)' % (table, dStr)
  nDone = 0
  for row in data:
    entries = [None] * len(row)
    for col in range(len(row)):
      if row[col] is not None and \
         (nullMarker is None or row[col] != nullMarker):
        if colTypes[col][0] == float:
          entries[col] = float(row[col])
        elif colTypes[col][0] == int:
          entries[col] = int(row[col])
        else:
          entries[col] = str(row[col])
      else:
        entries[col] = None
    block.append(tuple(entries))
    if len(block) >= blockSize:
      nDone += _insertBlock(cn, sqlStr, block)
      # note: in Python 3.12 `cn.autocommit != True`
      # is different from `not cn.autocommit` (GH #7009)
      if not hasattr(cn, 'autocommit') or cn.autocommit != True:
        cn.commit()
      block = []
  if len(block):
    nDone += _insertBlock(cn, sqlStr, block)
  # note: in Python 3.12 `cn.autocommit != True`
  # is different from `not cn.autocommit` (GH #7009)
  if not hasattr(cn, 'autocommit') or cn.autocommit != True:
    cn.commit()


def TextFileToDatabase(dBase, table, inF, delim=',', user='sysdba', password='masterkey',
                       maxColLabelLen=31, keyCol=None, nullMarker=None):
  """loads the contents of the text file into a database.

      **Arguments**

        - dBase: the name of the DB to use

        - table: the name of the table to create/overwrite

        - inF: the file like object from which the data should
          be pulled (must support readline())

        - delim: the delimiter used to separate fields

        - user: the user name to use in connecting to the DB

        - password: the password to use in connecting to the DB

        - maxColLabelLen: the maximum length a column label should be
          allowed to have (truncation otherwise)

        - keyCol: the column to be used as an index for the db

      **Notes**

        - if _table_ already exists, it is destroyed before we write
          the new data

        - we assume that the first row of the file contains the column names

    """
  table.replace('-', '_')
  table.replace(' ', '_')

  colHeadings = inF.readline().split(delim)
  _AdjustColHeadings(colHeadings, maxColLabelLen)
  nCols = len(colHeadings)
  data = []
  inL = inF.readline()
  while inL:
    inL = inL.replace('\r', '')
    inL = inL.replace('\n', '')
    splitL = inL.split(delim)
    if len(splitL) != nCols:
      print('>>>', repr(inL))
      assert len(splitL) == nCols, 'unequal length'
    tmpVect = []
    for entry in splitL:
      try:
        val = int(entry)
      except ValueError:
        try:
          val = float(entry)
        except ValueError:
          val = entry
      tmpVect.append(val)
    data.append(tmpVect)
    inL = inF.readline()
  nRows = len(data)

  # determine the types of each column
  colTypes = TypeFinder(data, nRows, nCols, nullMarker=nullMarker)
  typeStrs = GetTypeStrings(colHeadings, colTypes, keyCol=keyCol)
  colDefs = ','.join(typeStrs)

  _AddDataToDb(dBase, table, user, password, colDefs, colTypes, data, nullMarker=nullMarker)


def DatabaseToDatabase(fromDb, fromTbl, toDb, toTbl, fields='*', join='', where='', user='sysdba',
                       password='masterkey', keyCol=None, nullMarker='None'):
  """

     FIX: at the moment this is a hack

    """
  sio = StringIO()
  sio.write(
    DatabaseToText(fromDb, fromTbl, fields=fields, join=join, where=where, user=user,
                   password=password))
  sio.seek(0)
  TextFileToDatabase(toDb, toTbl, sio, user=user, password=password, keyCol=keyCol,
                     nullMarker=nullMarker)


if __name__ == '__main__':  # pragma: nocover
  sio = StringIO()
  sio.write('foo,bar,baz\n')
  sio.write('1,2,3\n')
  sio.write('1.1,4,5\n')
  sio.write('4,foo,6\n')
  sio.seek(0)
  import os

  from rdkit import RDConfig
  dirLoc = os.path.join(RDConfig.RDCodeDir, 'Dbase', 'TEST.GDB')

  TextFileToDatabase(dirLoc, 'fromtext', sio)
