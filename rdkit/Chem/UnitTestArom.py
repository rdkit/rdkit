#
#  Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""
unit testing code for the aromaticity detection
"""
from __future__ import print_function

import os
import unittest
from collections import namedtuple

from rdkit import Chem
from rdkit import RDConfig


TestData = namedtuple('TestData', 'lineNo,smiles,smiles2,tgtCount,tgtAtoms')


class TestCase(unittest.TestCase):

  def setUp(self):
    testdir = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data')
    self.fName1 = os.path.join(testdir, 'aromat_regress.txt')
    self.fName2 = os.path.join(testdir, 'NCI_aromat_regress.txt')

  def test1(self, maxFailures=50):
    self.checkFile(self.fName1, maxFailures=maxFailures)

  def test2(self, maxFailures=50):
    self.checkFile(self.fName2, maxFailures=maxFailures)

  def _readData(self, fName):
    with open(fName, 'r') as fdata:
      for lineNo, line in enumerate(fdata, 1):
        if len(line.strip()) == 0 or line[0] == '#':
          continue
        splitL = line.split('\t')
        if len(splitL) == 4:
          smi1, smi2, count, ats = splitL
          yield TestData(lineNo, smi1, smi2, int(count), eval(ats))

  def checkFile(self, filename, maxFailures=50):
    nFailed = 0
    for data in self._readData(filename):
      mol = Chem.MolFromSmiles(data.smiles)
      if mol is None:
        print('failure({0.lineNo}): {0.smiles}'.format(data))
        print('-----------------------------')
      else:
        aroms = [at.GetIdx() for at in mol.GetAtoms() if at.GetIsAromatic()]
        count = len(aroms)
        if count != data.tgtCount:
          print('Fail({0.lineNo}): {0.smiles}, {1}'.format(data, Chem.MolToSmiles(mol)))
          print('\t {0} != {1.tgtCount}'.format(count, data))
          print('\t {0}'.format(aroms))
          print('\t {0.tgtAtoms}'.format(data))
          print('-----------------------------')
          nFailed += 1
          if nFailed >= maxFailures:
            raise AssertionError('Too many failures')


if __name__ == '__main__':
  unittest.main()
