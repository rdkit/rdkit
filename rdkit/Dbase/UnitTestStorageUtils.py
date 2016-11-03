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
import traceback
import unittest

from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.Dbase.StorageUtils import RDIdToInt, ValidateRDId, IndexToRDId, RegisterItem, GetNextRDId


class TestCase(unittest.TestCase):

  def test_ValidateRDId(self):
    self.assertEqual(ValidateRDId('RDCmpd-000-009-9'), True)
    self.assertEqual(ValidateRDId('RDCmpd-009-000-009-8'), True)
    self.assertEqual(ValidateRDId('RDCmpd-009-000-109-8'), False)
    self.assertEqual(ValidateRDId('bogus'), False)
    self.assertEqual(ValidateRDId('RDCmpd-bogus-000-9'), False)

  def test_RDIdToInt(self):
    self.assertEqual(RDIdToInt('RDCmpd-000-009-9'), 9)
    self.assertEqual(RDIdToInt('RDCmpd-009-000-009-8'), 9000009)
    self.assertEqual(RDIdToInt('RDData_000_009_9'), 9)
    self.assertRaises(ValueError, RDIdToInt, 'RDCmpd-009-000-109-8')
    self.assertRaises(ValueError, RDIdToInt, 'bogus')

  def test_roundTrip(self):
    self.assertTrue(ValidateRDId(IndexToRDId(100)))
    self.assertTrue(ValidateRDId(IndexToRDId(10000, leadText='foo')))
    indices = [1, 100, 1000, 1000000]
    vals = [RDIdToInt(IndexToRDId(idx)) for idx in indices]
    self.assertEqual(indices, vals)

  def test_IndexToRDId(self):
    self.assertEqual(str(IndexToRDId(9)), 'RDCmpd-000-009-9')
    self.assertEqual(str(IndexToRDId(9009)), 'RDCmpd-009-009-8')
    self.assertEqual(str(IndexToRDId(9000009)), 'RDCmpd-009-000-009-8')
    self.assertEqual(str(IndexToRDId(9, leadText='RDAlt')), 'RDAlt-000-009-9')
    self.assertRaises(ValueError, IndexToRDId, -1)

  def test_RegisterItem(self):
    if RDConfig.useSqlLite:
      _, tempName = tempfile.mkstemp(suffix='sqlt')
      tempDbName = tempName
      shutil.copyfile(RDConfig.RDTestDatabase, tempDbName)
    else:  # pragma: nocover
      tempDbName = '::RDTests'
    conn = DbConnect(tempDbName)
    tblName = 'StorageTest'
    conn.AddTable(tblName, 'id varchar(32) not null primary key,label varchar(40),val int')
    self.assertEqual(
      RegisterItem(conn, tblName, 'label1', 'label', data=['label1', 1]), (1, 'RDCmpd-000-001-1'))
    self.assertEqual(
      RegisterItem(conn, tblName, 'label2', 'label', data=['label2', 1]), (1, 'RDCmpd-000-002-2'))
    self.assertEqual(
      RegisterItem(conn, tblName, 'label1', 'label', data=['label1', 1]), (0, 'RDCmpd-000-001-1'))
    self.assertEqual(str(GetNextRDId(conn, tblName)), 'RDCmpd-000-003-3')
    self.assertEqual(tuple(conn.GetData(table=tblName)[0]), ('RDCmpd-000-001-1', 'label1', 1))
    self.assertEqual(
      RegisterItem(conn, tblName, 'label10', 'label', data=['label10', 1], id='RDCmpd-000-010-1'),
      (1, 'RDCmpd-000-010-1'))
    self.assertEqual(str(GetNextRDId(conn, tblName)), 'RDCmpd-000-011-2')
    if RDConfig.useSqlLite and os.path.exists(tempDbName):
      try:
        os.unlink(tempDbName)
      except:  # pragma: nocover
        traceback.print_exc()


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
