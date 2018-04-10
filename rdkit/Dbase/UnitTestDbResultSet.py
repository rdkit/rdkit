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
"""unit testing code for the DbResultSet object

"""
import unittest

from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.Dbase.DbResultSet import DbResultSet, RandomAccessDbResultSet


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dbName = RDConfig.RDTestDatabase
    self.conn = DbConnect(self.dbName)
    self.curs = self.conn.GetCursor()

  def test1(self):
    """ test indexing in, ensure acceptable error conditions """
    cmd = 'select * from ten_elements'
    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd)
    self.assertRaises(IndexError, resultSet.__getitem__, -1)
    for i in range(12):
      try:
        _ = resultSet[i]
      except IndexError:
        assert i >= 10
    self.assertEqual(resultSet.GetColumnNames(), ('id', 'val'))
    self.assertEqual(resultSet.GetColumnTypes(), ('integer', 'integer'))
    self.assertEqual(resultSet.GetColumnNamesAndTypes(), (('id', 'integer'), ('val', 'integer')))

    cmd = 'select * from ten_elements_dups'
    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd)
    for i in range(22):
      try:
        _ = resultSet[i]
      except IndexError:
        assert i >= 20

    cmd = 'select * from ten_elements_dups'
    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd, removeDups=0)
    for i in range(22):
      try:
        _ = resultSet[i]
      except IndexError:
        assert i >= 10

    # Test iterator
    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd, removeDups=0)
    self.assertEqual(next(resultSet), resultSet[0])
    self.assertEqual(len(list(resultSet)), 10)

  def test2(self):
    cmd = 'select * from ten_elements'
    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd)
    assert len(resultSet) == 10
    for i in range(len(resultSet)):
      _ = resultSet[i]

  def test3(self):
    cmd = 'select * from ten_elements'
    resultSet = DbResultSet(self.curs, self.conn, cmd)
    r = [obj for obj in resultSet]
    self.assertEqual(len(r), 10)

    # Test iterator
    resultSet = DbResultSet(self.curs, self.conn, cmd)
    self.assertEqual(next(resultSet), (0, 11))
    self.assertEqual(len(list(resultSet)), 10)

  def test4(self):
    cmd = 'select * from ten_elements_dups'
    resultSet = DbResultSet(self.curs, self.conn, cmd)
    r = [obj for obj in resultSet]
    self.assertEqual(len(r), 20)

    resultSet = DbResultSet(self.curs, self.conn, cmd, removeDups=0)
    r = [obj for obj in resultSet]
    self.assertEqual(len(r), 10)

  def test5(self):
    cmd = 'select * from ten_elements_dups'
    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd)
    self.assertEqual(len(resultSet), 20)
    r = [obj for obj in resultSet]
    self.assertEqual(len(r), 20)

    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd)
    r = [obj for obj in resultSet]
    self.assertEqual(len(r), 20)
    self.assertEqual(len(resultSet), 20)

    resultSet = RandomAccessDbResultSet(self.curs, self.conn, cmd, removeDups=0)
    self.assertEqual(len(resultSet), 10)
    r = [obj for obj in resultSet]
    self.assertEqual(len(r), 10)

  def test6(self):
    cmd = 'select * from ten_elements_dups'
    resultSet = DbResultSet(self.curs, self.conn, cmd, removeDups=0)
    r = [obj for obj in resultSet]
    self.assertEqual(len(r), 10)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
