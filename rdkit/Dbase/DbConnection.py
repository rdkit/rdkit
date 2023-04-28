#
#  Copyright (C) 2000-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" defines class _DbConnect_, for abstracting connections to databases

"""

from rdkit.Dbase import DbInfo, DbModule, DbUtils


class DbError(RuntimeError):
  pass


class DbConnect(object):
  """  This class is intended to abstract away many of the details of
    interacting with databases.

    It includes some GUI functionality

  """

  def __init__(self, dbName='', tableName='', user='sysdba', password='masterkey'):
    """ Constructor

      **Arguments**  (all optional)

        - dbName: the name of the DB file to be used

        - tableName: the name of the table to be used

        - user: the username for DB access

        - password: the password to be used for DB access

    """

    self.dbName = dbName
    self.tableName = tableName
    self.user = user
    self.password = password
    self.cn = None
    self.cursor = None

  def GetTableNames(self, includeViews=0):
    """ gets a list of tables available in a database

      **Arguments**

      - includeViews: if this is non-null, the views in the db will
        also be returned

      **Returns**

          a list of table names

      **Notes**

       - this uses _DbInfo.GetTableNames_

    """
    return DbInfo.GetTableNames(self.dbName, self.user, self.password, includeViews=includeViews,
                                cn=self.cn)

  def GetColumnNames(self, table='', join='', what='*', where='', **kwargs):
    """ gets a list of columns available in the current table

      **Returns**

          a list of column names

      **Notes**

       - this uses _DbInfo.GetColumnNames_

    """
    table = table or self.tableName
    return DbInfo.GetColumnNames(self.dbName, table, self.user, self.password, join=join, what=what,
                                 cn=self.cn)

  def GetColumnNamesAndTypes(self, table='', join='', what='*', where='', **kwargs):
    """ gets a list of columns available in the current table along with their types

      **Returns**

          a list of 2-tuples containing:

            1) column name

            2) column type

      **Notes**

       - this uses _DbInfo.GetColumnNamesAndTypes_

    """
    table = table or self.tableName
    return DbInfo.GetColumnNamesAndTypes(self.dbName, table, self.user, self.password, join=join,
                                         what=what, cn=self.cn)

  def GetColumns(self, fields, table='', join='', **kwargs):
    """ gets a set of data from a table

      **Arguments**

       - fields: a string with the names of the fields to be extracted,
         this should be a comma delimited list

      **Returns**

          a list of the data

      **Notes**

        - this uses _DbUtils.GetColumns_

    """
    table = table or self.tableName
    return DbUtils.GetColumns(self.dbName, table, fields, self.user, self.password, join=join)

  def GetData(self, table=None, fields='*', where='', removeDups=-1, join='', transform=None,
              randomAccess=1, **kwargs):
    """ a more flexible method to get a set of data from a table

      **Arguments**

       - table: (optional) the table to use

       - fields: a string with the names of the fields to be extracted,
         this should be a comma delimited list

       - where: the SQL where clause to be used with the DB query

       - removeDups: indicates which column should be used to recognize
         duplicates in the data.  -1 for no duplicate removal.

      **Returns**

          a list of the data

      **Notes**

        - this uses _DbUtils.GetData_

    """
    table = table or self.tableName
    kwargs['forceList'] = kwargs.get('forceList', 0)
    return DbUtils.GetData(self.dbName, table, fieldString=fields, whereString=where,
                           user=self.user, password=self.password, removeDups=removeDups, join=join,
                           cn=self.cn, transform=transform, randomAccess=randomAccess, **kwargs)

  def GetDataCount(self, table=None, where='', join='', **kwargs):
    """ returns a count of the number of results a query will return

      **Arguments**

       - table: (optional) the table to use

       - where: the SQL where clause to be used with the DB query

       - join: the SQL join clause to be used with the DB query


      **Returns**

          an int

      **Notes**

        - this uses _DbUtils.GetData_

    """
    table = table or self.tableName
    return DbUtils.GetData(self.dbName, table, fieldString='count(*)', whereString=where,
                           cn=self.cn, user=self.user, password=self.password, join=join,
                           forceList=0)[0][0]

  def GetCursor(self):
    """ returns a cursor for direct manipulation of the DB
      only one cursor is available

    """
    if self.cursor is not None:
      return self.cursor

    self.cn = DbModule.connect(self.dbName, self.user, self.password)
    self.cursor = self.cn.cursor()
    return self.cursor

  def KillCursor(self):
    """ closes the cursor

    """
    self.cursor = None
    if self.cn is not None:
      self.cn.close()
    self.cn = None

  def AddTable(self, tableName, colString):
    """ adds a table to the database

    **Arguments**

      - tableName: the name of the table to add

      - colString: a string containing column definitions

    **Notes**

      - if a table named _tableName_ already exists, it will be dropped

      - the sqlQuery for addition is: "create table %(tableName) (%(colString))"

    """
    c = self.GetCursor()
    try:
      c.execute('drop table %s cascade' % tableName)
    except Exception:
      try:
        c.execute('drop table %s' % tableName)
      except Exception:
        pass
    self.Commit()

    addStr = 'create table %s (%s)' % (tableName, colString)
    try:
      c.execute(addStr)
    except Exception:
      import traceback
      print('command failed:', addStr)
      traceback.print_exc()
    else:
      self.Commit()

  def InsertData(self, tableName, vals):
    """ inserts data into a table

    **Arguments**

      - tableName: the name of the table to manipulate

      - vals: a sequence with the values to be inserted

    """
    c = self.GetCursor()
    if type(vals) != tuple:
      vals = tuple(vals)
    insTxt = '(' + ','.join([DbModule.placeHolder] * len(vals)) + ')'
    # insTxt = '(%s'%('%s,'*len(vals))
    # insTxt = insTxt[0:-1]+')'
    cmd = "insert into %s values %s" % (tableName, insTxt)
    try:
      c.execute(cmd, vals)
    except Exception:
      import traceback
      print('insert failed:')
      print(cmd)
      print('the error was:')
      traceback.print_exc()
      raise DbError("Insert Failed")

  def InsertColumnData(self, tableName, columnName, value, where):
    """ inserts data into a particular column of the table

    **Arguments**

      - tableName: the name of the table to manipulate

      - columnName: name of the column to update

      - value: the value to insert

      - where: a query yielding the row where the data should be inserted

    """
    c = self.GetCursor()
    cmd = "update %s set %s=%s where %s" % (tableName, columnName, DbModule.placeHolder, where)
    c.execute(cmd, (value, ))

  def AddColumn(self, tableName, colName, colType):
    """ adds a column to a table

    **Arguments**

      - tableName: the name of the table to manipulate

      - colName: name of the column to insert

      - colType: the type of the column to add

    """
    c = self.GetCursor()
    try:
      c.execute("alter table %s add %s %s" % (tableName, colName, colType))
    except Exception:
      print('AddColumn failed')

  def Commit(self):
    """ commits the current transaction

    """
    self.cn.commit()
