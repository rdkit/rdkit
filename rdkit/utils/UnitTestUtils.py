'''
Created on 26 Oct 2016

@author: peter
'''
import unittest
import doctest
from rdkit.utils import fileutils, chemutils, listutils


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(listutils, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def test_fileutils(self):
    filename = fileutils.__file__
    with open(filename) as inFile:
      line = fileutils.MoveToMatchingLine(inFile, 'pass')
      self.assertEqual(line, '  pass\n')
      self.assertRaises(fileutils.NoMatchFoundError, fileutils.MoveToMatchingLine, inFile, 'pass')
    with open(filename) as inFile:
      self.assertRaises(fileutils.NoMatchFoundError, fileutils.MoveToMatchingLine, inFile, 'pass',
                        fullMatch=True)

  def test_chemutils(self):
    self.assertEqual(chemutils.SplitComposition('Fe'), [('Fe', 1)])
    self.assertEqual(chemutils.SplitComposition('Fe3Al'), [('Fe', 3.0), ('Al', 1)])
    self.assertEqual(chemutils.SplitComposition('Fe99PdAl'), [('Fe', 99.0), ('Pd', 1), ('Al', 1)])
    self.assertEqual(
      chemutils.SplitComposition('TiNiSiSO12P'),
      [('Ti', 1), ('Ni', 1), ('Si', 1), ('S', 1), ('O', 12.0), ('P', 1)])
    temp = ['[Xe] 4f^12 6s^2', '[Xe] 4f^14 5d^6 6s^2', '[Xe] 4f^14 5d^10 6s^2',
            '[Xe] 4f^14 5d^10 6s^2 6p^1', '[Xe] 5d^10']
    for entry, expected in zip(temp, (14, 8, 2, 3, 10)):
      self.assertEqual(
        chemutils.ConfigToNumElectrons(entry, ignoreFullD=True, ignoreFullF=True), expected)

    for entry, expected in zip(temp, (14, 8, 12, 13, 10)):
      self.assertEqual(
        chemutils.ConfigToNumElectrons(entry, ignoreFullD=False, ignoreFullF=True), expected)

    for entry, expected in zip(temp, (14, 22, 16, 17, 10)):
      self.assertEqual(
        chemutils.ConfigToNumElectrons(entry, ignoreFullD=True, ignoreFullF=False), expected)

    for entry, expected in zip(temp, (14, 22, 26, 27, 10)):
      self.assertEqual(
        chemutils.ConfigToNumElectrons(entry, ignoreFullD=False, ignoreFullF=False), expected)


if __name__ == "__main__":  # pragma: nocover
  unittest.main()
