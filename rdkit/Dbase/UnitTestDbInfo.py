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
import unittest

from rdkit.Dbase import DbInfo, DbModule


class TestCase(unittest.TestCase):

  def setUp(self):
    self._oldFileWildcard = DbModule.fileWildcard

  def tearDown(self):
    DbModule.fileWildcard = self._oldFileWildcard

  def test_GetDbNames(self):
    DbModule.fileWildcard = '*.sqlite'
    self.assertEqual(
      len(DbInfo.GetDbNames(dirName=os.path.join(os.path.dirname(__file__), 'test_data'))), 2)

    DbModule.fileWildcard = '*.notexisting'
    self.assertEqual(DbInfo.GetDbNames(), [])

    DbModule.fileWildcard = None
    self.assertEqual(DbInfo.GetDbNames(), [])


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
