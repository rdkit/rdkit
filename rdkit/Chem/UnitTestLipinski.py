#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for Lipinski parameter calculation

  This provides a workout for the SMARTS matcher

"""


import os
import unittest

from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import Lipinski, rdMolDescriptors, Crippen

NonStrict = "NUM_ROTATABLEBONDS_O"
Strict = "NUM_ROTATABLEBONDS"

doLong = False


class TestCase(unittest.TestCase):

  def setUp(self):
    self.inFileName = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_200.props.sdf')

  def test1(self):
    " testing first 200 mols from NCI "
    # figure out which rotor version we are using
    m = Chem.MolFromSmiles("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C")
    if Lipinski.NumRotatableBonds(m) == 2:
      rot_prop = NonStrict
    else:
      rot_prop = Strict

    suppl = Chem.SDMolSupplier(self.inFileName)
    idx = 1
    for m in suppl:
      if m:
        calc = Lipinski.NHOHCount(m)
        orig = int(m.GetProp('NUM_LIPINSKIHDONORS'))
        assert calc == orig, 'bad num h donors for mol %d (%s): %d != %d' % (
          idx, m.GetProp('SMILES'), calc, orig)

        calc = Lipinski.NOCount(m)
        orig = int(m.GetProp('NUM_LIPINSKIHACCEPTORS'))
        assert calc == orig, 'bad num h acceptors for mol %d (%s): %d != %d' % (
          idx, m.GetProp('SMILES'), calc, orig)

        calc = Lipinski.NumHDonors(m)
        orig = int(m.GetProp('NUM_HDONORS'))
        assert calc == orig, 'bad num h donors for mol %d (%s): %d != %d' % (
          idx, m.GetProp('SMILES'), calc, orig)

        calc = Lipinski.NumHAcceptors(m)
        orig = int(m.GetProp('NUM_HACCEPTORS'))
        assert calc == orig, 'bad num h acceptors for mol %d (%s): %d != %d' % (
          idx, m.GetProp('SMILES'), calc, orig)

        calc = Lipinski.NumHeteroatoms(m)
        orig = int(m.GetProp('NUM_HETEROATOMS'))
        assert calc == orig, 'bad num heteroatoms for mol %d (%s): %d != %d' % (
          idx, m.GetProp('SMILES'), calc, orig)

        calc = Lipinski.NumRotatableBonds(m)
        orig = int(m.GetProp(rot_prop))
        assert calc == orig, 'bad num rotors for mol %d (%s): %d != %d' % (idx, m.GetProp('SMILES'),
                                                                           calc, orig)

        # test the underlying numrotatable bonds
        calc = rdMolDescriptors.CalcNumRotatableBonds(
          m, rdMolDescriptors.NumRotatableBondsOptions.NonStrict)
        orig = int(m.GetProp(NonStrict))
        assert calc == orig, 'bad num rotors for mol %d (%s): %d != %d' % (idx, m.GetProp('SMILES'),
                                                                           calc, orig)

        calc = rdMolDescriptors.CalcNumRotatableBonds(
          m, rdMolDescriptors.NumRotatableBondsOptions.Strict)
        orig = int(m.GetProp(Strict))
        assert calc == orig, 'bad num rotors for mol %d (%s): %d != %d' % (idx, m.GetProp('SMILES'),
                                                                           calc, orig)

      idx += 1

  def testIssue2183420(self):
    " testing a problem with the acceptor definition "
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('NC')) == 1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CNC')) == 1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CN(C)C')) == 1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('NC(=O)')) == 1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('NC(=O)C')) == 1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CNC(=O)')) == 1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('CNC(=O)C')) == 1)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('O=CNC(=O)C')) == 2)
    self.assertTrue(Lipinski.NumHAcceptors(Chem.MolFromSmiles('O=C(C)NC(=O)C')) == 2)


class TestCase_Regression(unittest.TestCase):

  @staticmethod
  def _regressionData(filename, col):
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
        yield lineNum, smi, mol, expected

  def __testDesc(self, fileN, col, func):
    for lineNum, smi, mol, expected in self._regressionData(fileN, col):
      if abs(expected - 666.0) < 1e-4:  # Unavailable data point
        continue
      try:
        val = func(mol)
      except Exception:
        val = 666
      assert abs(val - expected) < 1e-4, 'line %d, mol %s (calc = %f) should have val = %f' % (
        lineNum, smi, val, expected)

  def testLipinskiLong(self):
    """ Lipinski parameter """
    if not doLong:
      raise unittest.SkipTest('long test')
    fName = 'PP_descrs_regress.csv'
    self.__testDesc(fName, 30, Lipinski.NumHDonors)
    self.__testDesc(fName, 31, Lipinski.NumHeteroatoms)
    self.__testDesc(fName, 32, Lipinski.NumRotatableBonds)
    self.__testDesc(fName, 33, lambda x: Crippen.MolLogP(x, includeHs=1))

    fName = 'Block_regress.Lip.csv'
    self.__testDesc(fName, 1, Lipinski.NumHAcceptors)
    self.__testDesc(fName, 2, Lipinski.NumHDonors)
    self.__testDesc(fName, 3, Lipinski.NumHeteroatoms)
    self.__testDesc(fName, 4, Lipinski.NumRotatableBonds)

    fName = 'PP_descrs_regress.2.csv'
    self.__testDesc(fName, 33, lambda x: Crippen.MolLogP(x, includeHs=1))


if __name__ == '__main__':
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
