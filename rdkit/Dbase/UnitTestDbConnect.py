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
"""unit testing code for database connection objects

"""
import os
import shutil
import tempfile
import unittest

from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect


class TestCase(unittest.TestCase):

  def setUp(self):
    self.baseDir = os.path.join(RDConfig.RDCodeDir, 'Dbase', 'test_data')
    self.dbName = RDConfig.RDTestDatabase
    self.colHeads = ('int_col', 'floatCol', 'strCol')
    self.colTypes = ('integer', 'float', 'string')
    if RDConfig.useSqlLite:
      fd, tempName = tempfile.mkstemp(suffix='sqlt')
      self.fd = fd
      self.tempDbName = tempName
      shutil.copyfile(self.dbName, self.tempDbName)
    else:
      self.tempDbName = '::RDTests'

  def tearDown(self):
    if RDConfig.useSqlLite and os.path.exists(self.tempDbName):
      os.close(self.fd)
      try:
        os.unlink(self.tempDbName)
      except Exception:
        import traceback
        traceback.print_exc()

  def test_GetTableNames(self):
    # We can get the table names of a database with prior instantiation of a cursor
    conn = DbConnect(self.tempDbName)
    conn.GetCursor()
    names_Cursor = sorted(conn.GetTableNames())

    # and without (this tests functionality of DbInfo
    conn = DbConnect(self.tempDbName)
    names_noCursor = sorted(conn.GetTableNames())
    self.assertEqual(names_Cursor, names_noCursor)

  def testAddTable(self):
    """ tests AddTable and GetTableNames functionalities """
    newTblName = 'NEW_TABLE'
    conn = DbConnect(self.tempDbName)
    try:
      conn.GetCursor().execute('drop table %s' % (newTblName))
    except Exception:
      pass
    conn.Commit()

    self.assertNotIn(newTblName, [x.strip() for x in conn.GetTableNames()])
    conn.AddTable(newTblName, 'id int')
    self.assertIn(newTblName, [x.strip() for x in conn.GetTableNames()])

    self.assertEqual(conn.GetColumnNames(table=newTblName), ['id'])

    conn.GetCursor().execute('drop table %s' % (newTblName))

  def testCursor(self):
    """ tests GetCursor and GetTableNames functionalities """

    viewName = 'TEST_VIEW'
    conn = DbConnect(self.tempDbName)
    curs = conn.GetCursor()
    assert curs
    try:
      curs.execute('drop view %s' % (viewName))
    except Exception:
      pass
    try:
      curs.execute('create view %s as select val,id from ten_elements' % (viewName))
    except Exception:
      import traceback
      traceback.print_exc()
      raise AssertionError('create view failed')
    conn.Commit()

    self.assertNotIn(viewName, [x.strip() for x in conn.GetTableNames(includeViews=0)],
                     'improper view found')
    self.assertIn(viewName, [x.strip() for x in conn.GetTableNames(includeViews=1)],
                  'improper view not found')
    try:
      curs.execute('drop view %s' % (viewName))
    except Exception:
      raise AssertionError('drop table failed')

  def testGetData1(self):
    """ basic functionality
    """
    conn = DbConnect(self.dbName, 'ten_elements')
    d = conn.GetData(randomAccess=1)
    assert len(d) == 10
    assert tuple(d[0]) == (0, 11)
    assert tuple(d[2]) == (4, 31)
    self.assertRaises(IndexError, lambda: d[11])

    d = conn.GetColumns(fields='id,val')
    self.assertEqual(len(d), 10)
    assert tuple(d[0]) == (0, 11)
    assert tuple(d[2]) == (4, 31)

    self.assertEqual(conn.GetDataCount(), 10)

  def testGetData2(self):
    """ using removeDups
    """
    conn = DbConnect(self.dbName, 'ten_elements_dups')
    d = conn.GetData(randomAccess=1, removeDups=1)
    assert tuple(d[0]) == (0, 11)
    assert tuple(d[2]) == (4, 31)
    assert len(d) == 10
    self.assertRaises(IndexError, lambda: d[11])

  def testGetData3(self):
    """ without removeDups
    """
    conn = DbConnect(self.dbName, 'ten_elements_dups')
    d = conn.GetData(randomAccess=1, removeDups=-1)
    assert tuple(d[0]) == (0, 11)
    assert tuple(d[2]) == (2, 21)
    assert len(d) == 20
    self.assertRaises(IndexError, lambda: d[21])

    # repeat that test to make sure the table argument works
    conn = DbConnect(self.dbName, 'ten_elements')
    d = conn.GetData(table='ten_elements_dups', randomAccess=1, removeDups=-1)
    assert tuple(d[0]) == (0, 11)
    assert tuple(d[2]) == (2, 21)
    assert len(d) == 20
    self.assertRaises(IndexError, lambda: d[21])

  def testGetData4(self):
    """ non random access
    """
    conn = DbConnect(self.dbName, 'ten_elements')
    d = conn.GetData(randomAccess=0)
    self.assertRaises(TypeError, lambda: len(d))

    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs) == 10
    assert tuple(rs[0]) == (0, 11)
    assert tuple(rs[2]) == (4, 31)

  def testGetData5(self):
    """ using a RandomAccessDbResultSet with a Transform
    """

    def fn(x):
      return (x[0], x[1] * 2)

    conn = DbConnect(self.dbName, 'ten_elements')
    d = conn.GetData(randomAccess=1, transform=fn)

    assert tuple(d[0]) == (0, 22), str(d[0])
    assert tuple(d[2]) == (4, 62)
    assert len(d) == 10
    self.assertRaises(IndexError, lambda: d[11])

  def testGetData6(self):
    """ using a DbResultSet with a Transform
    """

    def fn(x):
      return (x[0], x[1] * 2)

    conn = DbConnect(self.dbName, 'ten_elements')
    d = conn.GetData(randomAccess=0, transform=fn)
    self.assertRaises(TypeError, lambda: len(d))
    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs) == 10
    assert tuple(rs[0]) == (0, 22)
    assert tuple(rs[2]) == (4, 62)

  def test_InsertData(self):
    newTblName = 'NEW_TABLE'
    conn = DbConnect(self.tempDbName)
    try:
      conn.GetCursor().execute('drop table %s' % (newTblName))
    except Exception:
      pass
    conn.Commit()
    conn.AddTable(newTblName, 'id int,val1 int, val2 int')
    for i in range(10):
      conn.InsertData(newTblName, (i, i + 1, 2 * i))
    conn.Commit()
    d = conn.GetData(table=newTblName)
    assert len(d) == 10

    self.assertEqual(len(conn.GetColumnNames(table=newTblName)), 3)
    conn.AddColumn(newTblName, 'val3', 'int')
    conn.Commit()
    self.assertEqual(len(conn.GetColumnNames(table=newTblName)), 4)
    d = conn.GetColumns('id,val3', table=newTblName)
    self.assertEqual(len(d), 10)
    self.assertTrue(all(r[1] is None for r in d))
    for r in d:
      conn.InsertColumnData(newTblName, 'val3', r[0], 'id={0}'.format(r[0]))
    conn.Commit()
    d = conn.GetColumns('id,val3', table=newTblName)
    self.assertTrue(all(r[0] == r[1] for r in d))

    d = None
    try:
      conn.GetCursor().execute('drop table %s' % (newTblName))
    except Exception:
      assert 0, 'drop table failed'

  def test_InsertColumnData(self):
    """ tests InsertData and InsertColumnData functionalities """
    newTblName = 'NEW_TABLE'
    conn = DbConnect(self.tempDbName)
    try:
      conn.GetCursor().execute('drop table %s' % (newTblName))
    except Exception:
      pass
    conn.Commit()
    conn.AddTable(newTblName, 'id int,val1 int, val2 int')
    for i in range(10):
      conn.InsertData(newTblName, (i, i + 1, 2 * i))
    conn.Commit()
    d = conn.GetData(table=newTblName)
    assert len(d) == 10

    d = None
    try:
      conn.GetCursor().execute('drop table %s' % (newTblName))
    except Exception:
      assert 0, 'drop table failed'


if __name__ == '__main__':
  unittest.main()
