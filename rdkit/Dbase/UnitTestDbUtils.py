# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the database utilities

"""
import os
import tempfile
import unittest
from io import StringIO

from rdkit import RDConfig
from rdkit.Dbase import DbUtils
from rdkit.Dbase.DbConnection import DbConnect

try:
  from contextlib import redirect_stdout
except Exception:
  from rdkit.TestRunner import redirect_stdout


class TestCase(unittest.TestCase):

  def setUp(self):
    self.baseDir = os.path.join(RDConfig.RDCodeDir, 'Dbase', 'test_data')
    self.dbName = RDConfig.RDTestDatabase
    if RDConfig.useSqlLite:
      fd, tempName = tempfile.mkstemp(suffix='sqlt')
      self.fd = fd
      self.tempDbName = tempName
    else:
      self.tempDbName = '::RDTests'
    self.colHeads = ('int_col', 'floatCol', 'strCol')
    self.colTypes = ('integer', 'float', 'string')

  def tearDown(self):
    if RDConfig.useSqlLite and os.path.exists(self.tempDbName):
      os.close(self.fd)
      try:
        os.unlink(self.tempDbName)
      except Exception:
        import traceback
        traceback.print_exc()

  def _confirm(self, tblName, dbName=None, colHeads=None, colTypes=None):
    dbName = dbName or self.dbName
    colHeads = colHeads or self.colHeads
    colTypes = colTypes or self.colTypes

    conn = DbConnect(dbName, tblName)
    res = conn.GetColumnNamesAndTypes()
    assert len(res) == len(colHeads), 'bad number of columns'
    names = [x[0] for x in res]
    for i in range(len(names)):
      assert names[i].upper() == colHeads[i].upper(), 'bad column head'
    if RDConfig.useSqlLite:
      # doesn't seem to be any column type info available
      return
    types = [x[1] for x in res]
    for i in range(len(types)):
      assert types[i] == colTypes[i], 'bad column type'

  def test1Txt(self):
    """ test reading from a text file """
    with open(os.path.join(self.baseDir, 'dbtest.csv'), 'r') as inF:
      tblName = 'fromtext'
      f = StringIO()
      with redirect_stdout(f):
        DbUtils.TextFileToDatabase(self.tempDbName, tblName, inF)
      self._confirm(tblName, dbName=self.tempDbName)

  def test3Txt(self):
    """ test reading from a text file including null markers"""
    with open(os.path.join(self.baseDir, 'dbtest.nulls.csv'), 'r') as inF:
      tblName = 'fromtext2'
      f = StringIO()
      with redirect_stdout(f):
        DbUtils.TextFileToDatabase(self.tempDbName, tblName, inF, nullMarker='NA')

      self._confirm(tblName, dbName=self.tempDbName)

  def testGetData1(self):
    """ basic functionality
    """
    d = DbUtils.GetData(self.dbName, 'ten_elements', forceList=1)
    assert len(d) == 10
    assert tuple(d[0]) == (0, 11)
    assert tuple(d[2]) == (4, 31)
    with self.assertRaisesRegexp(IndexError, ""):
      d[11]

  def testGetData2(self):
    """ using a RandomAccessDbResultSet
    """
    d = DbUtils.GetData(self.dbName, 'ten_elements', forceList=0, randomAccess=1)
    assert tuple(d[0]) == (0, 11)
    assert tuple(d[2]) == (4, 31)
    assert len(d) == 10
    with self.assertRaisesRegexp(IndexError, ""):
      d[11]

  def testGetData3(self):
    """ using a DbResultSet
    """
    d = DbUtils.GetData(self.dbName, 'ten_elements', forceList=0, randomAccess=0)
    with self.assertRaisesRegexp(TypeError, ""):
      len(d)
    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs) == 10
    assert tuple(rs[0]) == (0, 11)
    assert tuple(rs[2]) == (4, 31)

  def testGetData4(self):
    """ using a RandomAccessDbResultSet with a Transform
    """

    def fn(x):
      return (x[0], x[1] * 2)

    d = DbUtils.GetData(self.dbName, 'ten_elements', forceList=0, randomAccess=1, transform=fn)
    assert tuple(d[0]) == (0, 22)
    assert tuple(d[2]) == (4, 62)
    assert len(d) == 10
    with self.assertRaisesRegexp(IndexError, ""):
      d[11]

  def testGetData5(self):
    """ using a DbResultSet with a Transform
    """

    def fn(x):
      return (x[0], x[1] * 2)

    d = DbUtils.GetData(self.dbName, 'ten_elements', forceList=0, randomAccess=0, transform=fn)
    with self.assertRaisesRegexp(TypeError, ""):
      len(d)

    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs) == 10
    assert tuple(rs[0]) == (0, 22)
    assert tuple(rs[2]) == (4, 62)

  def test_take(self):
    self.assertEqual(list(DbUtils._take([1, 2, 3, 4], [2, 3])), [3, 4])
    self.assertEqual(list(DbUtils._take([1, 2, 3, 4], [0, 3])), [1, 4])

  def test_GetColumns(self):
    d = DbUtils.GetColumns(self.dbName, 'ten_elements', 'val')
    self.assertEqual(len(d), 10)

  def test_GetData_where(self):
    d = DbUtils.GetData(self.dbName, 'ten_elements_dups', forceList=0, randomAccess=0,
                        whereString='id<4')
    self.assertEqual(len(list(d)), 4)
    self.assertTrue(all(x[0] < 4 for x in d))

    d = DbUtils.GetData(self.dbName, 'ten_elements_dups', forceList=0, randomAccess=0,
                        whereString='id<10')
    self.assertEqual(len(list(d)), 10)
    self.assertTrue(all(x[0] < 10 for x in d))

    d = DbUtils.GetData(self.dbName, 'ten_elements_dups', removeDups=1, forceList=True)
    self.assertEqual(len(list(d)), 10)

  def test_DatabaseToText(self):
    txt = DbUtils.DatabaseToText(self.dbName, 'ten_elements')
    self.assertIn('id,val', txt)
    self.assertIn('0,11', txt)
    self.assertIn('18,101', txt)
    self.assertEqual(len(txt.split('\n')), 11)

    txt = DbUtils.DatabaseToText(self.dbName, 'ten_elements', fields='val')
    self.assertNotIn('id', txt)
    self.assertNotIn(',', txt)

    txt = DbUtils.DatabaseToText(self.dbName, 'ten_elements', where='id<4')
    self.assertIn('id,val', txt)
    self.assertEqual(len(txt.split('\n')), 3)

  def test_TypeFinder(self):
    data = [('-', 1.45, 'abc', None), (20, 3, 'defgh', None)]
    self.assertEqual(
      DbUtils.TypeFinder(data, 2, 4, nullMarker='-'), [[int, 2], [float, 4], [str, 5], [-1, 1]])

  def test_AdjustColHeadings(self):
    headers = ['abc def', ' abc def', 'abc-def ', 'abc.def']
    self.assertEqual(DbUtils._AdjustColHeadings(headers, 7), ['abc_def'] * 4)
    f = StringIO()
    with redirect_stdout(f):
      headers = ['abc def', ' abc def', 'abc-def ', 'abc.def']
      self.assertEqual(DbUtils._AdjustColHeadings(headers, 3), ['abc'] * 4)
    self.assertIn('Heading', f.getvalue())

  def test_GetTypeStrings(self):
    headers = ['pk', 'a', 'b', 'c']
    colTypes = [(int, 2), (int, 3), (float, 5), (str, 10)]
    self.assertEqual(
      DbUtils.GetTypeStrings(headers, colTypes),
      ['pk integer', 'a integer', 'b double precision', 'c varchar(10)'])
    self.assertEqual(
      DbUtils.GetTypeStrings(headers, colTypes, keyCol='pk'),
      ['pk integer not null primary key', 'a integer', 'b double precision', 'c varchar(10)'])

  def test_DatabaseToDatabase(self):
    tblName = 'db2db'
    f = StringIO()
    with redirect_stdout(f):
      DbUtils.DatabaseToDatabase(self.dbName, 'ten_elements', self.tempDbName, tblName)
    self._confirm(tblName, dbName=self.tempDbName, colHeads=['id', 'val'])


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
