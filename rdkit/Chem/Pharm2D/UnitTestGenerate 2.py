"""
unit testing code rdkit.Chem.Pharm2D.Generate

"""
import unittest

from rdkit.Chem.Pharm2D import Generate


class TestCase(unittest.TestCase):

  def test_Gen2DFingerprint(self):
    # The code is run from most of the other tests. This is only
    # to test the assertions
    self.assertRaises(ValueError, Generate.Gen2DFingerprint, 'mol', 'incorrectType')


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
