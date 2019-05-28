#
#  Copyright (C) 2004-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#


import unittest

from rdkit.Chem.Pharm3D.ExcludedVolume import ExcludedVolume


class TestCase(unittest.TestCase):

  def test_ExcludedVolume(self):
    # featInfo must have a length
    self.assertRaises(ValueError, ExcludedVolume, 123)
    self.assertRaises(ValueError, ExcludedVolume, [])
    self.assertRaises(ValueError, ExcludedVolume, [123, ])
    self.assertRaises(ValueError, ExcludedVolume, [[], ])

    featInfo = ([(0, ), 0.5, 1.0], )
    excludedVolume = ExcludedVolume(featInfo)
    self.assertEqual(excludedVolume.featInfo, featInfo)
    self.assertAlmostEqual(excludedVolume.exclusionDist, 3.0)
    excludedVolume = ExcludedVolume(featInfo, exclusionDist=3.14)
    self.assertEqual(excludedVolume.featInfo, featInfo)
    self.assertAlmostEqual(excludedVolume.exclusionDist, 3.14)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
