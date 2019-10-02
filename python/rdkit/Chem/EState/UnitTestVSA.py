# $Id$
#
#  Copyright (C) 2004-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for MOE-type descriptors with EStates

"""


import os
import unittest

from rdkit import Chem
from rdkit.Chem.EState import EState_VSA


class TestCase(unittest.TestCase):

  @staticmethod
  def referenceData():
    filename = os.sep.join(
      [os.path.dirname(os.path.abspath(__file__)), 'test_data', 'EState_VSA.csv'])
    with open(filename) as fin:
      header = fin.readline()
      header = [s.strip() for s in header.split(',')][1:]

      funcEstates = dict((k, getattr(EState_VSA, k)) for k in header)
      yield funcEstates

      for line in fin:
        line = [s.strip() for s in line.split(',')]
        smiles = line.pop(0)
        mol = Chem.MolFromSmiles(smiles)
        data = dict((k, float(v)) for k, v in zip(header, line))
        yield smiles, mol, data

  def test1(self):
    referenceData = self.referenceData()
    funcEstates = next(referenceData)
    for smiles, mol, data in referenceData:
      for name in funcEstates:
        calc = funcEstates[name](mol)
        exp = data[name]
        self.assertAlmostEqual(calc, exp, delta=1e-4,
                               msg='{0}: {1:.4f}!={2:.4f}'.format(smiles, calc, exp))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
