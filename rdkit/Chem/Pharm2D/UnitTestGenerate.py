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
"""unit testing code for the signature utils

"""
import doctest
import itertools
import unittest

from rdkit.Chem.Pharm2D import Generate


class TestCase(unittest.TestCase):

  def test_Gen2DFingerprint(self):
    # The code is run from most of the other tests. This is only 
    # to test the assertions
    self.assertRaises(ValueError, Generate.Gen2DFingerprint, 'mol', 'incorrectType')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
