# $Id$
#
#  Copyright (C) 2003-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for graph-theoretical descriptors

"""

import os.path
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import GraphDescriptors

doLong = False
_THREE_RING = Chem.MolFromSmarts('*1~*~*~1')


def feq(n1, n2, tol=1e-4):
  return abs(n1 - n2) <= tol


def _skip3rings(mol):
  """ The chi3v and chi3n descriptors changed for molecules with 3-rings """
  return mol.HasSubstructMatch(_THREE_RING)


def _hasAromaticAtoms(mol):
  """ The BertzCT descriptor changed for molecules with aromatic rings """
  return any(atom.GetIsAromatic() for atom in mol.GetAtoms())


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


class TestCase(unittest.TestCase):

  def __testDesc(self, fileN, col, func, molFilter=None):
    """ Regression tests for descriptor calculator.

    The optional argument molFilter can be used to skip examples. """
    for lineNum, smi, mol, expected in _regressionData(fileN, col):
      if molFilter and molFilter(mol):
        continue
      if feq(expected, 666.0):  # Unavailable data point
        continue
      try:
        val = func(mol)
      except Exception:
        val = 666
      assert feq(
        val, expected,
        1e-4), 'line %d, mol %s (calc = %f) should have val = %f' % (lineNum, smi, val, expected)

  def test_FullRegression(self):
    if not doLong:
      raise unittest.SkipTest('Full regression tests of descriptors skipped.')

  def testBertzCT(self):
    # test calculation of Bertz 'C(T)' index """
    data = [('C=CC=C', 21.01955), ('O=CC=O', 25.01955), ('FCC(=O)CF', 46.7548875),
            ('O=C1C=CC(=O)C=C1', 148.705216), ('C12C(F)=C(O)C(F)C1C(F)=C(O)C(F)2', 315.250442),
            ('C12CC=CCC1C(=O)C3CC=CCC3C(=O)2', 321.539522)]

    for smi, expected in data:
      m = Chem.MolFromSmiles(smi)
      newCT = GraphDescriptors.BertzCT(m, forceDMat=1)
      self.assertAlmostEqual(
        newCT, expected, delta=1e-3,
        msg='mol %s (CT calc = %f) should have CT = %f' % (smi, newCT, expected))

    if doLong:
      # We need to skip molecules with aromatic rings, due to changes in the
      # treatment of aromatic atoms. (Tests pass actually even without the filter!)
      self.__testDesc('PP_descrs_regress.2.csv', 1, GraphDescriptors.BertzCT,
                      molFilter=_hasAromaticAtoms)

  def testHallKierAlpha(self):
    self.__testDesc('PP_descrs_regress.csv', 3, GraphDescriptors.HallKierAlpha)
    if doLong:
      self.__testDesc('PP_descrs_regress.2.csv', 3, GraphDescriptors.HallKierAlpha)

  def testIpc(self):
    data = [('CCCCC', 1.40564, 11.24511), ('CCC(C)C', 1.37878, 9.65148),
            ('CC(C)(C)C', 0.72193, 3.60964), ('CN(CC)CCC', 1.67982, 31.91664),
            ('C1CCCCC1', 1.71997, 34.39946), ('CC1CCCCC1', 1.68562, 47.19725),
            ('Cc1ccccc1', 1.68562, 47.19725), ('CC(C)=C(C)C', 1.36096, 13.60964),
            ('C#N', 1.00000, 2.00000), ('OC#N', 0.91830, 2.75489)]
    for smi, res1, res2 in data:
      m = Chem.MolFromSmiles(smi)
      Ipc = GraphDescriptors.Ipc(m, forceDMat=1)
      Ipc_avg = GraphDescriptors.Ipc(m, avg=1, forceDMat=1)
      self.assertAlmostEqual(
        Ipc_avg, res1, delta=1e-3,
        msg='mol %s (Ipc_avg=%f) should have Ipc_avg=%f' % (smi, Ipc_avg, res1))
      avgIpc = GraphDescriptors.AvgIpc(m, forceDMat=1)
      self.assertEqual(Ipc_avg, avgIpc)

      self.assertAlmostEqual(Ipc, res2, delta=1e-3,
                             msg='mol %s (Ipc=%f) should have Ipc=%f' % (smi, Ipc, res2))

      Ipc = GraphDescriptors.Ipc(m)
      Ipc_avg = GraphDescriptors.Ipc(m, avg=1)
      self.assertAlmostEqual(
        Ipc_avg, res1, delta=1e-3,
        msg='2nd pass: mol %s (Ipc_avg=%f) should have Ipc_avg=%f' % (smi, Ipc_avg, res1))
      self.assertAlmostEqual(Ipc, res2, delta=1e-3,
                             msg='2nd pass: mol %s (Ipc=%f) should have Ipc=%f' % (smi, Ipc, res2))

      if doLong:
        self.__testDesc('PP_descrs_regress.csv', 4, GraphDescriptors.Ipc)
        self.__testDesc('PP_descrs_regress.2.csv', 4, GraphDescriptors.Ipc)

  def testKappa1(self):
    """ test calculation of the Hall-Kier kappa1 value

     corrected data from Tables 3 and 6 of Rev. Comp. Chem. vol 2, 367-422, (1991)
    """
    data = [('C12CC2C3CC13', 2.344), ('C1CCC12CC2', 3.061), ('C1CCCCC1', 4.167), ('CCCCCC', 6.000),
            ('CCC(C)C1CCC(C)CC1', 9.091), ('CC(C)CC1CCC(C)CC1', 9.091),
            ('CC(C)C1CCC(C)CCC1', 9.091)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      kappa = GraphDescriptors.Kappa1(m)
      assert feq(kappa, res, 1e-3), 'mol %s (kappa1=%f) should have kappa1=%f' % (smi, kappa, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 31, GraphDescriptors.Kappa1)

  def testKappa2(self):
    """ test calculation of the Hall-Kier kappa2 value

     corrected data from Tables 5 and 6 of Rev. Comp. Chem. vol 2, 367-422, (1991)

    """
    data = [('[C+2](C)(C)(C)(C)(C)C', 0.667), ('[C+](C)(C)(C)(C)(CC)', 1.240),
            ('C(C)(C)(C)(CCC)', 2.3444), ('CC(C)CCCC', 4.167), ('CCCCCCC', 6.000),
            ('CCCCCC', 5.000), ('CCCCCCC', 6.000), ('C1CCCC1', 1.440), ('C1CCCC1C', 1.633),
            ('C1CCCCC1', 2.222), ('C1CCCCCC1', 3.061), ('CCCCC', 4.00), ('CC=CCCC', 4.740),
            ('C1=CN=CN1', 0.884), ('c1ccccc1', 1.606), ('c1cnccc1', 1.552), ('n1ccncc1', 1.500),
            ('CCCCF', 3.930), ('CCCCCl', 4.290), ('CCCCBr', 4.480), ('CCC(C)C1CCC(C)CC1', 4.133),
            ('CC(C)CC1CCC(C)CC1', 4.133), ('CC(C)C1CCC(C)CCC1', 4.133)]
    for smi, res in data:
      # we have some molecules in there with bogus valences, so we need
      # to do custom sanitization
      m = Chem.MolFromSmiles(smi, sanitize=False)
      m.UpdatePropertyCache(strict=False)
      Chem.SanitizeMol(
        m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
      kappa = GraphDescriptors.Kappa2(m)
      assert feq(kappa, res, 1e-3), 'mol %s (kappa2=%f) should have kappa2=%f' % (smi, kappa, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 32, GraphDescriptors.Kappa2)

  def testKappa3(self):
    """ test calculation of the Hall-Kier kappa3 value

     corrected data from Tables 3 and 6 of Rev. Comp. Chem. vol 2, 367-422, (1991)

    """
    data = [('C[C+](C)(C)(C)C(C)(C)C', 2.000),
            ('CCC(C)C(C)(C)(CC)', 2.380), ('CCC(C)CC(C)CC', 4.500), ('CC(C)CCC(C)CC', 5.878),
            ('CC(C)CCCC(C)C', 8.000), ('CCC(C)C1CCC(C)CC1', 2.500), ('CC(C)CC1CCC(C)CC1', 3.265),
            ('CC(C)C1CCC(C)CCC1', 2.844)]
    for smi, res in data:
      # we have some molecules in there with bogus valences, so we need
      # to do custom sanitization
      m = Chem.MolFromSmiles(smi, sanitize=False)
      m.UpdatePropertyCache(strict=False)
      Chem.SanitizeMol(
        m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
      kappa = GraphDescriptors.Kappa3(m)
      assert feq(kappa, res, 1e-3), 'mol %s (kappa3=%f) should have kappa3=%f' % (smi, kappa, res)

    self.__testDesc('PP_descrs_regress.csv', 5, GraphDescriptors.Kappa3)
    if doLong:
      self.__testDesc('PP_descrs_regress.2.csv', 5, GraphDescriptors.Kappa3)

  def testBalabanJ(self):
    """ test calculation of the Balaban J value

      J values are from Balaban's paper and have had roundoff
      errors and typos corrected.
    """
    data = [  # alkanes
      ('CC', 1.0),
      ('CCC', 1.6330),
      ('CCCC', 1.9747),
      ('CC(C)C', 2.3238),
      ('CCCCC', 2.1906),
      ('CC(C)CC', 2.5396),
      ('CC(C)(C)C', 3.0237),
      ('CCCCCC', 2.3391),
      ('CC(C)CCC', 2.6272),
      ('CCC(C)CC', 2.7542),
      ('CC(C)(C)CC', 3.1685),
      ('CC(C)C(C)C', 2.9935),

      # cycloalkanes
      ('C1CCCCC1', 2.0000),
      ('C1C(C)CCCC1', 2.1229),
      ('C1C(CC)CCCC1', 2.1250),
      ('C1C(C)C(C)CCC1', 2.2794),
      ('C1C(C)CC(C)CC1', 2.2307),
      ('C1C(C)CCC(C)C1', 2.1924),
      ('C1C(CCC)CCCC1', 2.0779),
      ('C1C(C(C)C)CCCC1', 2.2284),
      ('C1C(CC)C(C)CCC1', 2.2973),
      ('C1C(CC)CC(C)CC1', 2.2317),
      ('C1C(CC)CCC(C)C1', 2.1804),
      ('C1C(C)C(C)C(C)CC1', 2.4133),
      ('C1C(C)C(C)CC(C)C1', 2.3462),
      ('C1C(C)CC(C)CC1(C)', 2.3409),
      # aromatics
      ('c1ccccc1', 3.0000),
      ('c1c(C)cccc1', 3.0215),
      ('c1c(CC)cccc1', 2.8321),
      ('c1c(C)c(C)ccc1', 3.1349),
      ('c1c(C)cc(C)cc1', 3.0777),
      ('c1c(C)ccc(C)c1', 3.0325),
      ('c1c(CCC)cccc1', 2.6149),
      ('c1c(C(C)C)cccc1', 2.8483),
      ('c1c(CC)c(C)ccc1', 3.0065),
      ('c1c(CC)cc(C)cc1', 2.9369),
      ('c1c(CC)ccc(C)c1', 2.8816),
      ('c1c(C)c(C)c(C)cc1', 3.2478),
      ('c1c(C)c(C)cc(C)c1', 3.1717),
      ('c1c(C)cc(C)cc1(C)', 3.1657)
    ]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      j = GraphDescriptors.BalabanJ(m, forceDMat=1)
      assert feq(j, res), 'mol %s (J=%f) should have J=%f' % (smi, j, res)
      j = GraphDescriptors.BalabanJ(m)
      assert feq(j, res), 'second pass: mol %s (J=%f) should have J=%f' % (smi, j, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 1, GraphDescriptors.BalabanJ)

  def testChi0Long(self):
    self.__testDesc('PP_descrs_regress.csv', 2, GraphDescriptors.Chi0)
    if doLong:
      self.__testDesc('PP_descrs_regress.2.csv', 2, GraphDescriptors.Chi0)
      self.__testDesc('PP_descrs_regress.rest.2.csv', 5, GraphDescriptors.Chi0)

  def testChi1Long(self):
    if not doLong:
      raise unittest.SkipTest('long test')
    self.__testDesc('PP_descrs_regress.rest.2.csv', 8, GraphDescriptors.Chi1)

  def testChi0v(self):
    data = [('CCCCCC', 4.828), ('CCC(C)CC', 4.992), ('CC(C)CCC', 4.992), ('CC(C)C(C)C', 5.155),
            ('CC(C)(C)CC', 5.207), ('CCCCCO', 4.276), ('CCC(O)CC', 4.439), ('CC(O)(C)CC', 4.654),
            ('c1ccccc1O', 3.834), ('CCCl', 2.841), ('CCBr', 3.671), ('CCI', 4.242)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi0v(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi0v=%f) should have Chi0V=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 7, GraphDescriptors.Chi0v)

  def testChi1v(self):
    data = [('CCCCCC', 2.914), ('CCC(C)CC', 2.808), ('CC(C)CCC', 2.770), ('CC(C)C(C)C', 2.643),
            ('CC(C)(C)CC', 2.561), ('CCCCCO', 2.523), ('CCC(O)CC', 2.489), ('CC(O)(C)CC', 2.284),
            ('c1ccccc1O', 2.134)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi1v(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi1v=%f) should have Chi1V=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 10, GraphDescriptors.Chi1v)

  def testChi2v(self):
    data = [
      ('CCCCCC', 1.707),
      ('CCC(C)CC', 1.922),
      ('CC(C)CCC', 2.183),
      ('CC(C)C(C)C', 2.488),
      ('CC(C)(C)CC', 2.914),
      ('CCCCCO', 1.431),
      ('CCC(O)CC', 1.470),
      ('CC(O)(C)CC', 2.166),
      ('c1ccccc1O', 1.336),
    ]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi2v(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi2v=%f) should have Chi2V=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 12, GraphDescriptors.Chi2v)

  def testChi3v(self):
    data = [('CCCCCC', 0.957), ('CCC(C)CC', 1.394), ('CC(C)CCC', 0.866), ('CC(C)C(C)C', 1.333),
            ('CC(C)(C)CC', 1.061), ('CCCCCO', 0.762), ('CCC(O)CC', 0.943), ('CC(O)(C)CC', 0.865),
            ('c1ccccc1O', 0.756)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi3v(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi3v=%f) should have Chi3V=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 14, GraphDescriptors.Chi3v,
                      molFilter=_skip3rings)

  def testChi4v(self):
    data = [('CCCCCC', 0.500), ('CCC(C)CC', 0.289), ('CC(C)CCC', 0.577), ('CC(C)C(C)C', 0.000),
            ('CC(C)(C)CC', 0.000), ('CCCCCO', 0.362), ('CCC(O)CC', 0.289), ('CC(O)(C)CC', 0.000),
            ('c1ccccc1O', 0.428)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi4v(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi4v=%f) should have Chi4V=%f' % (smi, chi, res)

  def testChi5v(self):
    data = [('CCCCCC', 0.250), ('CCC(C)CC', 0.000), ('CC(C)CCC', 0.000), ('CC(C)C(C)C', 0.000),
            ('CC(C)(C)CC', 0.000), ('CCCCCO', 0.112), ('CCC(O)CC', 0.000), ('CC(O)(C)CC', 0.000),
            ('c1ccccc1O', 0.242)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.ChiNv_(m, 5)
      assert feq(chi, res, 1e-3), 'mol %s (Chi5v=%f) should have Chi5V=%f' % (smi, chi, res)

  def testChi0n(self):
    data = [
      ('CCCCCC', 4.828),
      ('CCC(C)CC', 4.992),
      ('CC(C)CCC', 4.992),
      ('CC(C)C(C)C', 5.155),
      ('CC(C)(C)CC', 5.207),
      ('CCCCCO', 4.276),
      ('CCC(O)CC', 4.439),
      ('CC(O)(C)CC', 4.654),
      ('c1ccccc1O', 3.834),
      ('CCCl', 2.085),
      ('CCBr', 2.085),
      ('CCI', 2.085),
    ]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi0n(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi0n=%f) should have Chi0n=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 6, GraphDescriptors.Chi0n)

  def testChi1n(self):
    """ test calculation of Chi1n

    """
    data = [('CCCCCC', 2.914), ('CCC(C)CC', 2.808), ('CC(C)CCC', 2.770), ('CC(C)C(C)C', 2.643),
            ('CC(C)(C)CC', 2.561), ('CCCCCO', 2.523), ('CCC(O)CC', 2.489), ('CC(O)(C)CC', 2.284),
            ('c1ccccc1O', 2.134)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi1n(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi1n=%f) should have Chi1N=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 9, GraphDescriptors.Chi1n)

  def testChi2n(self):
    data = [('CCCCCC', 1.707), ('CCC(C)CC', 1.922), ('CC(C)CCC', 2.183), ('CC(C)C(C)C', 2.488),
            ('CC(C)(C)CC', 2.914), ('CCCCCO', 1.431), ('CCC(O)CC', 1.470), ('CC(O)(C)CC', 2.166),
            ('c1ccccc1O', 1.336)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi2n(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi2n=%f) should have Chi2N=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 11, GraphDescriptors.Chi2n)

  def testChi3n(self):
    data = [('CCCCCC', 0.957), ('CCC(C)CC', 1.394), ('CC(C)CCC', 0.866), ('CC(C)C(C)C', 1.333),
            ('CC(C)(C)CC', 1.061), ('CCCCCO', 0.762), ('CCC(O)CC', 0.943), ('CC(O)(C)CC', 0.865),
            ('c1ccccc1O', 0.756)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi3n(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi3n=%f) should have Chi3N=%f' % (smi, chi, res)

    if doLong:
      self.__testDesc('PP_descrs_regress.rest.2.csv', 13, GraphDescriptors.Chi3n,
                      molFilter=_skip3rings)

  def testChi4n(self):
    data = [('CCCCCC', 0.500), ('CCC(C)CC', 0.289), ('CC(C)CCC', 0.577), ('CC(C)C(C)C', 0.000),
            ('CC(C)(C)CC', 0.000), ('CCCCCO', 0.362), ('CCC(O)CC', 0.289), ('CC(O)(C)CC', 0.000),
            ('c1ccccc1O', 0.428)]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      chi = GraphDescriptors.Chi4n(m)
      assert feq(chi, res, 1e-3), 'mol %s (Chi4n=%f) should have Chi4N=%f' % (smi, chi, res)

  def testIssue125(self):
    # test an issue with calculating BalabanJ
    smi = 'O=C(OC)C1=C(C)NC(C)=C(C(OC)=O)C1C2=CC=CC=C2[N+]([O-])=O'
    m1 = Chem.MolFromSmiles(smi)
    m2 = Chem.MolFromSmiles(smi)
    Chem.MolToSmiles(m1)
    j1 = GraphDescriptors.BalabanJ(m1)
    j2 = GraphDescriptors.BalabanJ(m2)
    assert feq(j1, j2)

  def testOrderDepend(self):
    data = [('C=CC=C', 21.01955, 2.73205), ('O=CC=O', 25.01955, 2.73205),
            ('FCC(=O)CF', 46.7548875, 2.98816), ('O=C1C=CC(=O)C=C1', 148.705216, 2.8265),
            ('C12C(F)=C(O)C(F)C1C(F)=C(O)C(F)2', 315.250442, 2.4509),
            ('C12CC=CCC1C(=O)C3CC=CCC3C(=O)2', 321.539522, 1.95986)]

    for smi, CT, bal in data:
      m = Chem.MolFromSmiles(smi)
      newBal = GraphDescriptors.BalabanJ(m, forceDMat=1)
      assert feq(newBal, bal, 1e-4), 'mol %s %f!=%f' % (smi, newBal, bal)
      m = Chem.MolFromSmiles(smi)
      newCT = GraphDescriptors.BertzCT(m, forceDMat=1)
      assert feq(newCT, CT, 1e-4), 'mol %s (CT calc = %f) should have CT = %f' % (smi, newCT, CT)
      m = Chem.MolFromSmiles(smi)
      newCT = GraphDescriptors.BertzCT(m, forceDMat=1)
      assert feq(newCT, CT, 1e-4), 'mol %s (CT calc = %f) should have CT = %f' % (smi, newCT, CT)
      newBal = GraphDescriptors.BalabanJ(m, forceDMat=1)
      assert feq(newBal, bal, 1e-4), 'mol %s %f!=%f' % (smi, newBal, bal)

      m = Chem.MolFromSmiles(smi)
      newBal = GraphDescriptors.BalabanJ(m, forceDMat=1)
      assert feq(newBal, bal, 1e-4), 'mol %s %f!=%f' % (smi, newBal, bal)
      newCT = GraphDescriptors.BertzCT(m, forceDMat=1)
      assert feq(newCT, CT, 1e-4), 'mol %s (CT calc = %f) should have CT = %f' % (smi, newCT, CT)

  def testPathCounts(self):
    """ FIX: this should be in some other file

    """
    data = [
      ('CCCCCC', (6, 5, 4, 3, 2, 1)),
      ('CCC(C)CC', (6, 5, 5, 4, 1, 0)),
      ('CC(C)CCC', (6, 5, 5, 3, 2, 0)),
      ('CC(C)C(C)C', (6, 5, 6, 4, 0, 0)),
      ('CC(C)(C)CC', (6, 5, 7, 3, 0, 0)),
      ('CCCCCO', (6, 5, 4, 3, 2, 1)),
      ('CCC(O)CC', (6, 5, 5, 4, 1, 0)),
      ('CC(O)(C)CC', (6, 5, 7, 3, 0, 0)),
      ('c1ccccc1O', (7, 7, 8, 8, 8, 8)),
    ]
    for smi, res in data:
      m = Chem.MolFromSmiles(smi)
      for i in range(1, 6):
        cnt = len(Chem.FindAllPathsOfLengthN(m, i, useBonds=1))
        assert cnt == res[i], (smi, i, cnt, res[i], Chem.FindAllPathsOfLengthN(m, i, useBonds=1))
        cnt = len(Chem.FindAllPathsOfLengthN(m, i + 1, useBonds=0))
        assert cnt == res[i], (smi, i, cnt, res[i],
                               Chem.FindAllPathsOfLengthN(m, i + 1, useBonds=1))


class TestCase_python(unittest.TestCase):
  """ Test the Python implementation of the various descriptors """

  def test_equivalence(self):
    if not doLong:
      raise unittest.SkipTest('long test')
    self._compareImplementations('PP_descrs_regress.csv', 3, GraphDescriptors.HallKierAlpha,
                                 GraphDescriptors._pyHallKierAlpha)
    self._compareImplementations('PP_descrs_regress.2.csv', 3, GraphDescriptors.HallKierAlpha,
                                 GraphDescriptors._pyHallKierAlpha)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 31, GraphDescriptors.Kappa1,
                                 GraphDescriptors._pyKappa1)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 32, GraphDescriptors.Kappa2,
                                 GraphDescriptors._pyKappa2)
    self._compareImplementations('PP_descrs_regress.csv', 5, GraphDescriptors.Kappa3,
                                 GraphDescriptors._pyKappa3)
    self._compareImplementations('PP_descrs_regress.2.csv', 5, GraphDescriptors.Kappa3,
                                 GraphDescriptors._pyKappa3)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 7, GraphDescriptors.Chi0v,
                                 GraphDescriptors._pyChi0v)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 10, GraphDescriptors.Chi1v,
                                 GraphDescriptors._pyChi1v)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 12, GraphDescriptors.Chi2v,
                                 GraphDescriptors._pyChi2v)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 14, GraphDescriptors.Chi3v,
                                 GraphDescriptors._pyChi3v, molFilter=_skip3rings)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 6, GraphDescriptors.Chi0n,
                                 GraphDescriptors._pyChi0n)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 9, GraphDescriptors.Chi1n,
                                 GraphDescriptors._pyChi1n)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 11, GraphDescriptors.Chi2n,
                                 GraphDescriptors._pyChi2n)
    self._compareImplementations('PP_descrs_regress.rest.2.csv', 13, GraphDescriptors.Chi3n,
                                 GraphDescriptors._pyChi3n, molFilter=_skip3rings)

  def _compareImplementations(self, filename, col, cFunc, pyFunc, molFilter=None):
    """ Comparison of two implementations.
    The optional argument molFilter can be used to skip examples. """
    for lineNum, smi, mol, expected in _regressionData(filename, col):
      if molFilter and molFilter(mol):
        continue
      if feq(expected, 666.0):  # Unavailable data point
        continue
      cVal = cFunc(mol)
      pyVal = pyFunc(mol)
      assert feq(cVal, pyVal,
                 1e-4), 'line %d, mol %s (c = %f, py = %f)' % (lineNum, smi, cVal, pyVal)


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
