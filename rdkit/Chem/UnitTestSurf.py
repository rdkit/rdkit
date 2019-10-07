#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""
unit testing code for calculations in rdkit.Chem.MolSurf
"""


from collections import namedtuple
import os.path
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import MolSurf

doLong = False
TestData = namedtuple('TestData', 'lineNo,smiles,mol,expected')


class TestCase(unittest.TestCase):
  dataNCI200 = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.tpsa.csv')
  dataNCI5000 = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'NCI_5K_TPSA.csv')
  dataTPSAregr = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'tpsa_regr.csv')

  @classmethod
  def readNCI_200(cls):
    for data in readPSAtestData(cls.dataNCI200):
      yield data

  @classmethod
  def readNCI_5000(cls):
    for data in readPSAtestData(cls.dataNCI5000):
      yield data

  @classmethod
  def readTPSAregres(cls):
    for data in readPSAtestData(cls.dataTPSAregr):
      yield data

  def testTPSAShort(self):
    for data in self.readNCI_200():
      calc = MolSurf.TPSA(data.mol)
      self.assertAlmostEqual(
        calc, data.expected, delta=1e-4,
        msg='bad TPSA for SMILES {0.smiles} ({1:.2f} != {0.expected:.2f})'.format(data, calc))

  def testTPSALong(self):
    if not doLong:
      raise unittest.SkipTest('long test')
    for data in self.readNCI_5000():
      try:
        calc = MolSurf.TPSA(data.mol)
      except Exception:
        raise AssertionError(
          'Line {0.lineNo}: TPSA Calculation failed for SMILES {0.smiles}'.format(data))
      self.assertAlmostEqual(
        calc, data.expected, delta=1e-4,
        msg='bad TPSA for SMILES {0.smiles} ({1:.2f} != {0.expected:.2f})'.format(data, calc))

  def testTPSALongNCI(self):
    if not doLong:
      raise unittest.SkipTest('long test')
    for data in self.readTPSAregres():
      try:
        calc = MolSurf.TPSA(data.mol)
      except Exception:
        raise AssertionError(
          'Line {0.lineNo}: TPSA Calculation failed for SMILES {0.smiles}'.format(data))
      self.assertAlmostEqual(
        calc, data.expected, delta=1e-4,
        msg='bad TPSA for SMILES {0.smiles} ({1:.2f} != {0.expected:.2f})'.format(data, calc))

  def testHsAndTPSA(self):
    """
     testing the impact of Hs in the graph on PSA calculations
     This was sf.net issue 1969745
    """
    mol = Chem.MolFromSmiles('c1c[nH]cc1')
    molH = Chem.AddHs(mol)
    psa = MolSurf.TPSA(mol)
    psaH = MolSurf.TPSA(molH)

    if psa != psaH:
      psac = MolSurf.rdMolDescriptors._CalcTPSAContribs(mol)
      psaHc = MolSurf.rdMolDescriptors._CalcTPSAContribs(molH)
      for i, v in enumerate(psac):
        print('\t', i, '\t', v, '\t', psaHc[i])
      while i < len(psaHc):
        print('\t\t\t', psaHc[i])
        i += 1
    self.assertEqual(psa, psaH)

    for data in self.readNCI_200():
      mol = Chem.AddHs(data.mol)
      calc = MolSurf.TPSA(mol)
      self.assertAlmostEqual(
        calc, data.expected, delta=1e-4,
        msg='bad TPSA for SMILES {0.smiles} ({1:.2f} != {0.expected:.2f})'.format(data, calc))

    if doLong:
      for data in self.readNCI_5000():
        mol = Chem.AddHs(data.mol)
        calc = MolSurf.TPSA(mol)
        self.assertAlmostEqual(
          calc, data.expected, delta=1e-4,
          msg='bad TPSA for SMILES {0.smiles} ({1:.2f} != {0.expected:.2f})'.format(data, calc))


class TestCase_descriptorRegression(unittest.TestCase):

  def __testDesc(self, fileN, col, func):
    for data in readRegressionData(fileN, col):
      if abs(data.expected - 666.0) < 1e-4:
        print(data)
        raise AssertionError('check this data entry')
      try:
        val = func(data.mol)
      except Exception:
        val = 666
      self.assertAlmostEqual(
        val, data.expected, delta=1e-4,
        msg='line {0.lineNo}, mol {0.smiles} (calc = {1}) should have val = {0.expected}'.format(
          data, val))

  def testLabuteASALong(self):
    if not doLong:
      raise unittest.SkipTest('long test')
    col = 6
    self.__testDesc('PP_descrs_regress.csv', col, lambda x: MolSurf.LabuteASA(x, includeHs=1))
    self.__testDesc('PP_descrs_regress.2.csv', col, lambda x: MolSurf.LabuteASA(x, includeHs=1))

  def testTPSALong(self):
    col = 28
    self.__testDesc('PP_descrs_regress.csv', col, MolSurf.TPSA)
    if doLong:
      self.__testDesc('PP_descrs_regress.2.csv', col, MolSurf.TPSA)

  def testMOELong(self):
    if not doLong:
      raise unittest.SkipTest('long test')
    fName = 'PP_descrs_regress.VSA.csv'
    self.__testDesc(fName, 1, MolSurf.SMR_VSA1)
    self.__testDesc(fName, 2, MolSurf.SMR_VSA10)
    self.__testDesc(fName, 3, MolSurf.SMR_VSA2)
    self.__testDesc(fName, 4, MolSurf.SMR_VSA3)
    self.__testDesc(fName, 5, MolSurf.SMR_VSA4)
    self.__testDesc(fName, 6, MolSurf.SMR_VSA5)
    self.__testDesc(fName, 7, MolSurf.SMR_VSA6)
    self.__testDesc(fName, 8, MolSurf.SMR_VSA7)
    self.__testDesc(fName, 9, MolSurf.SMR_VSA8)
    self.__testDesc(fName, 10, MolSurf.SMR_VSA9)
    self.__testDesc(fName, 11, MolSurf.SlogP_VSA1)
    self.__testDesc(fName, 12, MolSurf.SlogP_VSA10)
    self.__testDesc(fName, 13, MolSurf.SlogP_VSA11)
    self.__testDesc(fName, 14, MolSurf.SlogP_VSA12)

    fName = 'PP_descrs_regress.VSA.2.csv'
    self.__testDesc(fName, 1, MolSurf.SMR_VSA1)
    self.__testDesc(fName, 2, MolSurf.SMR_VSA10)
    self.__testDesc(fName, 11, MolSurf.SlogP_VSA1)
    self.__testDesc(fName, 12, MolSurf.SlogP_VSA10)
    self.__testDesc(fName, 13, MolSurf.SlogP_VSA11)
    self.__testDesc(fName, 14, MolSurf.SlogP_VSA12)


class TestCase_python(unittest.TestCase):
  """ Test the Python implementation of the various descriptors """

  def test_pyTPSA(self):
    for data in TestCase.readNCI_200():
      molPy = Chem.MolFromSmiles(data.smiles)
      self.assertAlmostEqual(MolSurf.TPSA(data.mol), MolSurf._pyTPSA(molPy))

  def test_pyLabuteASA(self):
    for data in TestCase.readNCI_200():
      molPy = Chem.MolFromSmiles(data.smiles)
      self.assertAlmostEqual(MolSurf.LabuteASA(data.mol), MolSurf.pyLabuteASA(molPy))

  def test_pyPEOE_VSA_(self):
    for data in TestCase.readNCI_200():
      molPy = Chem.MolFromSmiles(data.smiles)
      for calcC, calcPy in zip(MolSurf.PEOE_VSA_(data.mol), MolSurf.pyPEOE_VSA_(molPy,
                                                                                force=False)):
        self.assertAlmostEqual(calcC, calcPy)

  def test_pySlogP_VSA_(self):
    for data in TestCase.readNCI_200():
      molPy = Chem.MolFromSmiles(data.smiles)
      for calcC, calcPy in zip(
          MolSurf.SlogP_VSA_(data.mol), MolSurf.pySlogP_VSA_(molPy, force=False)):
        self.assertAlmostEqual(calcC, calcPy)

  def test_pySMR_VSA_(self):
    for data in TestCase.readNCI_200():
      molPy = Chem.MolFromSmiles(data.smiles)
      for calcC, calcPy in zip(MolSurf.SMR_VSA_(data.mol), MolSurf.pySMR_VSA_(molPy, force=False)):
        self.assertAlmostEqual(calcC, calcPy)

  def test_pyLabuteHelper(self):
    for data in TestCase.readNCI_200():
      molPy = Chem.MolFromSmiles(data.smiles)
      for calcC, calcPy in zip(MolSurf._LabuteHelper(data.mol), MolSurf._pyLabuteHelper(molPy)):
        self.assertAlmostEqual(calcC, calcPy)


def readPSAtestData(filename):
  """ Read test data for PSA method from file """
  with open(filename, 'r') as f:
    for lineNo, line in enumerate(f, 1):
      if line[0] == '#':
        continue
      smiles, expected = line.strip().split(',')
      mol = Chem.MolFromSmiles(smiles)
      if not mol:
        raise AssertionError('molecule construction failed on line %d' % lineNo)
      yield TestData(lineNo, smiles, mol, float(expected))


def readRegressionData(filename, col):
  """ Return entries form regression dataset.
  Returns the line number, smiles, molecule, and the value found in column col
  """
  with open(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', filename), 'r') as inF:
    for lineNum, line in enumerate(inF, 1):
      if line[0] == '#':
        continue
      splitL = line.split(',')
      smi = splitL[0]
      mol = Chem.MolFromSmiles(smi)
      if mol is None:
        raise AssertionError('line %d, smiles: %s' % (lineNum, smi))
      expected = float(splitL[col])
      yield TestData(lineNum, smi, mol, expected)


if __name__ == '__main__':  # pragma: nocover
  import argparse
  import sys
  parser = argparse.ArgumentParser()
  parser.add_argument('-l', default=False, action='store_true', dest='doLong')
  args = parser.parse_args()
  doLong = args.doLong

  # Remove the -l flag if present so that it isn't interpreted by unittest.main()
  if 'l' in sys.argv:
    sys.argv.remove('-l')
  unittest.main()
