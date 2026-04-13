#
#  Copyright (C) 2026  Emily Rhodes
#         All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import unittest

from rdkit import Chem


class TestCase(unittest.TestCase):

  def testGetPropDefault(self):
    m = Chem.MolFromSmiles("CC")

    self.assertTrue(m.GetProp("missing_prop") is None)
    self.assertTrue(m.GetProp("missing_prop", defaultValue=None) is None)
    self.assertEqual(m.GetProp("missing_prop", defaultValue="fallback"), "fallback")

    m.SetProp("present_prop", "value")
    self.assertEqual(m.GetProp("present_prop", defaultValue="fallback"), "value")


if __name__ == '__main__':
  print("Testing GetProp default value")
  unittest.main()
