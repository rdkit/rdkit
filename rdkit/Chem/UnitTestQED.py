from __future__ import print_function

from collections import namedtuple
import doctest
import os.path
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import QED


doLong = False
TestData = namedtuple('TestData', 'lineNo,smiles,mol,expected')
dataNCI200 = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'QED', 'NCI_200_qed.csv')
dataRegression = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'QED', 'Regression_qed.csv')


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(QED, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def testNCI200(self):
    for d in readTestData(dataNCI200):
      self.assertAlmostEqual(QED.qed(d.mol), d.expected,
                             msg='QED not equal to expected in line {}'.format(d.lineNo))
      # Check that adding hydrogens will not change the result
      # This is currently not the case. Hydrogens change the number of rotatable bonds and the
      # number of alerts.
      if False:
        mol = Chem.AddHs(d.mol)
        self.assertAlmostEqual(QED.qed(mol), d.expected,
                               msg='QED not equal to expected in line {}'.format(d.lineNo))

  def testRegression(self):
    if not doLong:
      raise unittest.SkipTest('long test')
    for d in readTestData(dataRegression):
      self.assertAlmostEqual(QED.qed(d.mol), d.expected,
                             msg='QED not equal to expected in line {}'.format(d.lineNo))

  def test_properties(self):
    m = Chem.MolFromSmiles('N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O')
    p = QED.properties(m)
    self.assertAlmostEqual(p.MW, 337.456)
    self.assertAlmostEqual(p.ALOGP, -0.55833)
    self.assertAlmostEqual(p.HBA, 6)
    self.assertAlmostEqual(p.HBD, 5)
    self.assertAlmostEqual(p.PSA, 173.33)
    self.assertAlmostEqual(p.ROTB, 7)
    self.assertAlmostEqual(p.AROM, 1)
    self.assertAlmostEqual(p.ALERTS, 3)

    if False:
      p = QED.properties(Chem.AddHs(m))
      self.assertAlmostEqual(p.MW, 337.456)
      self.assertAlmostEqual(p.ALOGP, -0.55833)
      self.assertAlmostEqual(p.HBA, 6)
      self.assertAlmostEqual(p.HBD, 5)
      self.assertAlmostEqual(p.PSA, 173.33)
      self.assertAlmostEqual(p.ROTB, 7)
      self.assertAlmostEqual(p.AROM, 1)
      self.assertAlmostEqual(p.ALERTS, 3)


def readTestData(filename):
  """ Read test data from file """
  with open(filename, 'r') as f:
    for lineNo, line in enumerate(f, 1):
      if line[0] == '#':
        continue
      smiles, expected = line.strip().split(',')
      mol = Chem.MolFromSmiles(smiles)
      if not mol:
        raise AssertionError('molecule construction failed on line %d' % lineNo)
      yield TestData(lineNo, smiles, mol, float(expected))


def updateTestData():
  """ Update the test data. This should only be done if the method changes! """
  for filename in (dataNCI200, dataRegression,):
    data = list(readTestData(filename))
    with open(filename, 'w') as f:
      print('# Test data for QED descriptor', file=f)
      for d in data:
        expected = QED.qed(d.mol)
        print('{0.smiles},{1}'.format(d, expected), file=f)


if __name__ == '__main__':  # pragma: nocover
  import argparse
  import sys
  parser = argparse.ArgumentParser()
  parser.add_argument('-l', default=False, action='store_true', dest='doLong')
  parser.add_argument('-u', default=False, action='store_true', dest='updateTestData')
  args = parser.parse_args()

  # Handle possible arguments
  doLong = args.doLong
  if args.doLong:
    sys.argv.remove('-l')

  if args.updateTestData:
    updateTestData()
    sys.argv.remove('-u')

  unittest.main()
