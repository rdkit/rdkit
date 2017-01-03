#
#  Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import unittest
from rdkit.Chem import pyPeriodicTable


class TestCase(unittest.TestCase):
  def test_Issue1194(self):
    assert isinstance(pyPeriodicTable.metalNames, list)
    self.assertIn('Al', pyPeriodicTable.metalNames)
    self.assertIn('Sc', pyPeriodicTable.metalNames)
    self.assertIn('Lr', pyPeriodicTable.metalNames)
    self.assertEqual(len(pyPeriodicTable.metalNames), len(pyPeriodicTable.metalNumList))
    self.assertEqual(len(pyPeriodicTable.metalNames), 69)

  def test_KierHall(self):
    self.assertEqual(len(pyPeriodicTable.hallKierAlphas), 10)

if __name__ == '__main__':  # pragma: nocover
  unittest.main()
