# $Id$
#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import sys

from rdkit import RDConfig
from rdkit.Dbase import DbModule

sqlTextTypes = DbModule.sqlTextTypes
sqlIntTypes = DbModule.sqlIntTypes
sqlFloatTypes = DbModule.sqlFloatTypes
sqlBinTypes = DbModule.sqlBinTypes


def GetDbNames(user='sysdba', password='masterkey', dirName='.', dBase='::template1', cn=None):
  """ returns a list of databases that are available

      **Arguments**

        - user: the username for DB access

        - password: the password to be used for DB access

      **Returns**

        - a list of db names (strings)

    """
  if DbModule.getDbSql:
    if not cn:
      try:
        cn = DbModule.connect(dBase, user, password)
      except Exception:
        print('Problems opening database: %s' % (dBase))
        return []
    c = cn.cursor()
    c.execute(DbModule.getDbSql)
    if RDConfig.usePgSQL:
      names = ['::' + str(x[0]) for x in c.fetchall()]
    else:
      names = ['::' + str(x[0]) for x in c.fetchall()]
    names.remove(dBase)
  elif DbModule.fileWildcard:
    import glob
    import os.path
    names = glob.glob(os.path.join(dirName, DbModule.fileWildcard))
  else:
    names = []
  return names


def GetTableNames(dBase, user='sysdba', password='masterkey', includeViews=0, cn=None):
  """ returns a list of tables available in a database

      **Arguments**

        - dBase: the name of the DB file to be used

        - user: the username for DB access

        - password: the password to be used for DB access

        - includeViews: if this is non-null, the views in the db will
          also be returned

      **Returns**

        - a list of table names (strings)

    """
  if not cn:
    try:
      cn = DbModule.connect(dBase, user, password)
    except Exception:
      print('Problems opening database: %s' % (dBase))
      return []

  c = cn.cursor()
  if not includeViews:
    comm = DbModule.getTablesSql
  else:
    comm = DbModule.getTablesAndViewsSql
  c.execute(comm)
  names = [str(x[0]).upper() for x in c.fetchall()]
  if RDConfig.usePgSQL and 'PG_LOGDIR_LS' in names:
    names.remove('PG_LOGDIR_LS')
  return names


def GetColumnInfoFromCursor(cursor):
  if cursor is None or cursor.description is None:
    return []
  results = []
  if not RDConfig.useSqlLite:
    for item in cursor.description:
      cName = item[0]
      cType = item[1]
      if cType in sqlTextTypes:
        typeStr = 'string'
      elif cType in sqlIntTypes:
        typeStr = 'integer'
      elif cType in sqlFloatTypes:
        typeStr = 'float'
      elif cType in sqlBinTypes:
        typeStr = 'binary'
      else:
        sys.stderr.write('odd type in col %s: %s\n' % (cName, str(cType)))
      results.append((cName, typeStr))
  else:
    r = cursor.fetchone()
    if not r:
      return results
    for i, v in enumerate(r):
      cName = cursor.description[i][0]
      typ = type(v)
      if isinstance(v, str):
        typeStr = 'string'
      elif typ == int:
        typeStr = 'integer'
      elif typ == float:
        typeStr = 'float'
      elif typ in (memoryview, bytes):
        typeStr = 'binary'
      else:
        sys.stderr.write('odd type in col %s: %s\n' % (cName, typ))
      results.append((cName, typeStr))
  return results


def GetColumnNamesAndTypes(dBase, table, user='sysdba', password='masterkey', join='', what='*',
                           cn=None):
  """ gets a list of columns available in a DB table along with their types

      **Arguments**

        - dBase: the name of the DB file to be used

        - table: the name of the table to query

        - user: the username for DB access

        - password: the password to be used for DB access

        - join: an optional join clause (omit the verb 'join')

        - what: an optional clause indicating what to select

      **Returns**

        - a list of 2-tuples containing:

            1) column name

            2) column type

    """
  if not cn:
    cn = DbModule.connect(dBase, user, password)
  c = cn.cursor()
  cmd = 'select %s from %s' % (what, table)
  if join:
    cmd += ' join %s' % (join)
  c.execute(cmd)
  return GetColumnInfoFromCursor(c)


def GetColumnNames(dBase, table, user='sysdba', password='masterkey', join='', what='*', cn=None):
  """ gets a list of columns available in a DB table

      **Arguments**

        - dBase: the name of the DB file to be used

        - table: the name of the table to query

        - user: the username for DB access

        - password: the password to be used for DB access

        - join: an optional join clause  (omit the verb 'join')

        - what: an optional clause indicating what to select

      **Returns**

        -  a list of column names

    """
  if not cn:
    cn = DbModule.connect(dBase, user, password)
  c = cn.cursor()
  cmd = 'select %s from %s' % (what, table)
  if join:
    if join.strip().find('join') != 0:
      join = 'join %s' % (join)
    cmd += ' ' + join
  c.execute(cmd)
  c.fetchone()
  desc = c.description
  res = [str(x[0]) for x in desc]
  return res
