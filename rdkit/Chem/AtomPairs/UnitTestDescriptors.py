# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import doctest
import gzip
import os
import pickle
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem.AtomPairs import Pairs, Sheridan, Torsions, Utils


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(Pairs, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(Sheridan, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(Torsions, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(Utils, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def setUp(self):
    self.testDataPath = os.path.join(RDConfig.RDCodeDir, 'Chem', 'AtomPairs', 'test_data')
    inF = gzip.open(os.path.join(self.testDataPath, 'mols1000.pkl.gz'), 'rb')
    self.mols = pickle.load(inF, encoding='bytes')

  def testPairsRegression(self):
    inF = gzip.open(os.path.join(self.testDataPath, 'mols1000.aps.pkl.gz'), 'rb')
    atomPairs = pickle.load(inF, encoding='bytes')
    for i, m in enumerate(self.mols):
      ap = Pairs.GetAtomPairFingerprint(m)
      if ap != atomPairs[i]:  # pragma: nocover
        debugFingerprint(m, ap, atomPairs[i])
      self.assertEqual(ap, atomPairs[i])
      self.assertNotEqual(ap, atomPairs[i - 1])

  def testTorsionsRegression(self):
    inF = gzip.open(os.path.join(self.testDataPath, 'mols1000.tts.pkl.gz'), 'rb')
    torsions = pickle.load(inF, encoding='bytes')
    for i, m in enumerate(self.mols):
      tt = Torsions.GetTopologicalTorsionFingerprintAsIntVect(m)
      if tt != torsions[i]:  # pragma: nocover
        debugFingerprint(m, tt, torsions[i])
      self.assertEqual(tt, torsions[i])
      self.assertNotEqual(tt, torsions[i - 1])

  def testGetTopologicalTorsionFingerprintAsIds(self):
    mol = Chem.MolFromSmiles('C1CCCCN1')
    tt = Torsions.GetTopologicalTorsionFingerprint(mol)
    self.assertEqual(tt.GetNonzeroElements(), {4437590049: 2, 8732557345: 2, 4445978657: 2})
    tt = Torsions.GetTopologicalTorsionFingerprintAsIds(mol)
    self.assertEqual(sorted(tt),
                     [4437590049, 4437590049, 4445978657, 4445978657, 8732557345, 8732557345])
    tt = Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)
    self.assertEqual(tt.GetNonzeroElements(), {4437590049: 2, 8732557345: 2, 4445978657: 2})

  def testGithub334(self):
    m1 = Chem.MolFromSmiles('N#C')
    self.assertEqual(Chem.GetNumPiElectrons(m1.GetAtomWithIdx(0)), 2)
    self.assertEqual(Chem.GetNumPiElectrons(m1.GetAtomWithIdx(1)), 2)

    m1 = Chem.MolFromSmiles('N#[CH]')
    self.assertEqual(Chem.GetNumPiElectrons(m1.GetAtomWithIdx(0)), 2)
    self.assertEqual(Chem.GetNumPiElectrons(m1.GetAtomWithIdx(1)), 2)

    m = Chem.MolFromSmiles('C=C')
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(0)), 1)

    m = Chem.MolFromSmiles('C#CC')
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(0)), 2)
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(1)), 2)

    m = Chem.MolFromSmiles('O=C=CC')
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(0)), 1)
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(1)), 2)
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(2)), 1)
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(3)), 0)

    m = Chem.MolFromSmiles('c1ccccc1')
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(0)), 1)

    # FIX: this behaves oddly in these cases:

    m = Chem.MolFromSmiles('S(=O)(=O)')
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(0)), 2)

    m = Chem.MolFromSmiles('S(=O)(=O)(O)O')
    self.assertEqual(Chem.GetNumPiElectrons(m.GetAtomWithIdx(0)), 0)


def debugFingerprint(mol, fpCalc, fpExpected):  # pragma: nocover
  print(Chem.MolToSmiles(mol))
  pd = fpCalc.GetNonzeroElements()
  rd = fpExpected.GetNonzeroElements()
  for k, v in pd.items():
    if k in rd:
      if rd[k] != v:
        print('>>>1', k, v, rd[k])
    else:
      print('>>>2', k, v)
  for k, v in rd.items():
    if k in pd:
      if pd[k] != v:
        print('>>>3', k, v, pd[k])
    else:
      print('>>>4', k, v)


if __name__ == '__main__':
  unittest.main()
