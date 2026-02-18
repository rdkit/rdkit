#
#  Copyright (C) 2025 Niels Maeder and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Testing Flat matrix functionality from rdkit.DataStructs.__init__.py."""

import unittest

from rdkit.DataStructs import getElementFromFlatMatrix, getNForFlatMatrix


class TestCase(unittest.TestCase):
  def testGetNFromFlatMatrix(self):
    self.assertEqual(getNForFlatMatrix([0] * 1), 2)
    self.assertEqual(getNForFlatMatrix([0] * 3), 3)
    self.assertEqual(getNForFlatMatrix([0] * 6), 4)

  def testGetElementFromFlatMatrix(self):
    # (0,1)=a, (0,2)=b, (1,2)=c, (0,3)=d, ...
    flat = ["a", "b", "c", "d", "e", "f"]
    self.assertEqual(getElementFromFlatMatrix(flat, 1, 0), "a")
    self.assertEqual(getElementFromFlatMatrix(flat, 0, 1), "a")
    self.assertEqual(getElementFromFlatMatrix(flat, 2, 0), "b")
    self.assertEqual(getElementFromFlatMatrix(flat, 0, 2), "b")
    self.assertEqual(getElementFromFlatMatrix(flat, 2, 1), "c")
    self.assertEqual(getElementFromFlatMatrix(flat, 1, 2), "c")
    self.assertEqual(getElementFromFlatMatrix(flat, 0, 3), "d")
    self.assertEqual(getElementFromFlatMatrix(flat, 3, 0), "d")
    self.assertEqual(getElementFromFlatMatrix(flat, 1, 3), "e")
    self.assertEqual(getElementFromFlatMatrix(flat, 3, 1), "e")
    self.assertEqual(getElementFromFlatMatrix(flat, 2, 3), "f")
    self.assertEqual(getElementFromFlatMatrix(flat, 3, 2), "f")

    self.assertEqual(getElementFromFlatMatrix(flat, 0, 0), 0.0)


