#
#  Copyright (C) 2007-2017 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" General descriptor testing code

"""

import doctest
import os.path
import pickle
import unittest

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Descriptors3D, Lipinski, rdMolDescriptors


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(Descriptors, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def testGithub1287(self):
    smis = ('CCC', )
    for smi in smis:
      m = Chem.MolFromSmiles(smi)
      self.assertTrue(m)
      for nm, fn in Descriptors._descList:
        try:
          _ = fn(m)
        except Exception:
          import traceback
          traceback.print_exc()
          raise AssertionError('SMILES: %s; Descriptor: %s' % (smi, nm))

  def testBadAtomHandling(self):
    smis = ('CC[Pu]', 'CC[*]')
    for smi in smis:
      m = Chem.MolFromSmiles(smi)
      self.assertTrue(m)
      for nm, fn in Descriptors._descList:
        try:
          v = fn(m)
        except RuntimeError:
          # 3D descriptors fail since the mol has no conformers
          pass
        except Exception:
          import traceback
          traceback.print_exc()
          raise AssertionError('SMILES: %s; Descriptor: %s' % (smi, nm))

  def testMolFormula(self):
    for (smiles, expected) in (
      ("[NH4+]", "H4N+"),
      ("c1ccccc1", "C6H6"),
      ("C1CCCCC1", "C6H12"),
      ("c1ccccc1O", "C6H6O"),
      ("C1CCCCC1O", "C6H12O"),
      ("C1CCCCC1=O", "C6H10O"),
      ("N[Na]", "H2NNa"),
      ("[C-][C-]", "C2-2"),
      ("[H]", "H"),
      ("[H-1]", "H-"),
      ("[H-1]", "H-"),
      ("[CH2]", "CH2"),
      ("[He-2]", "He-2"),
      ("[U+3]", "U+3"),
    ):
      mol = Chem.MolFromSmiles(smiles)
      actual = AllChem.CalcMolFormula(mol)
      self.assertEqual(actual, expected)

  def testMQNDetails(self):
    refFile = os.path.join(os.path.dirname(__file__), 'test_data', 'MQNs_regress.pkl')
    refFile2 = os.path.join(os.path.dirname(__file__), 'test_data', 'MQNs_non_strict_regress.pkl')
    # figure out which definition we are currently using
    m = Chem.MolFromSmiles("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C")
    if Lipinski.NumRotatableBonds(m) == 2:
      refFile = refFile2

    with open(refFile, 'rb') as intf:
      refData = pickle.load(intf)
    fn = os.path.join(os.path.dirname(__file__), 'test_data', 'aromat_regress.txt')
    ms = [x for x in Chem.SmilesMolSupplier(fn, delimiter='\t')]
    for i, m in enumerate(ms):
      mqns = rdMolDescriptors.MQNs_(m)
      if mqns != refData[i][1]:
        indices = [(j, x, y) for j, x, y in zip(range(len(mqns)), mqns, refData[i][1]) if x != y]
        print(i, Chem.MolToSmiles(m), indices)
      self.assertEqual(mqns, refData[i][1])

  def testMQN(self):
    m = Chem.MolFromSmiles("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C")
    if Lipinski.NumRotatableBonds(m) == 2:
      tgt = [
        42917, 274, 870, 621, 135, 1582, 29, 3147, 5463, 6999, 470, 62588, 19055, 4424, 309, 24061,
        17820, 1, 9303, 24146, 16076, 5560, 4262, 646, 746, 13725, 5430, 2629, 362, 24211, 15939,
        292, 41, 20, 1852, 5642, 31, 9, 1, 2, 3060, 1750
      ]
    else:
      tgt = [
        42917, 274, 870, 621, 135, 1582, 29, 3147, 5463, 6999, 470, 62588, 19055, 4424, 309, 24061,
        17820, 1, 8314, 24146, 16076, 5560, 4262, 646, 746, 13725, 5430, 2629, 362, 24211, 15939,
        292, 41, 20, 1852, 5642, 31, 9, 1, 2, 3060, 1750
      ]
      tgt = [
        42917, 274, 870, 621, 135, 1582, 29, 3147, 5463, 6999, 470, 62588, 19055, 4424, 309, 24059,
        17822, 1, 8314, 24146, 16076, 5560, 4262, 646, 746, 13725, 5430, 2629, 362, 24211, 15939,
        292, 41, 20, 1852, 5642, 31, 9, 1, 2, 3060, 1750
      ]
    fn = os.path.join(os.path.dirname(__file__), 'test_data', 'aromat_regress.txt')
    ms = [x for x in Chem.SmilesMolSupplier(fn, delimiter='\t')]
    vs = np.zeros((42, ), np.int32)

    for m in ms:
      vs += rdMolDescriptors.MQNs_(m)
    self.assertEqual(list(vs), tgt)

  def test_FpDensityMorgan(self):
    self.assertEqual(Descriptors.FpDensityMorgan1.version, '1.0.0')
    self.assertEqual(Descriptors.FpDensityMorgan2.version, '1.0.0')
    self.assertEqual(Descriptors.FpDensityMorgan3.version, '1.0.0')

    m = Chem.MolFromSmiles('C')
    self.assertAlmostEqual(Descriptors.FpDensityMorgan2(m), 1)
    m = Chem.MolFromSmiles('CC')
    self.assertAlmostEqual(Descriptors.FpDensityMorgan2(m), 1)
    m = Chem.MolFromSmiles('CCC')
    self.assertAlmostEqual(Descriptors.FpDensityMorgan2(m), 4.0 / 3)
    m = Chem.MolFromSmiles('C' * 10)
    self.assertAlmostEqual(Descriptors.FpDensityMorgan2(m), 8.0 / 10)
    m = Chem.MolFromSmiles('C' * 100)
    self.assertAlmostEqual(Descriptors.FpDensityMorgan2(m), 8.0 / 100)

    m = Chem.MolFromSmiles('CCCc1ccccc1')
    fpd1 = Descriptors.FpDensityMorgan1(m)
    fpd2 = Descriptors.FpDensityMorgan2(m)
    fpd3 = Descriptors.FpDensityMorgan3(m)
    self.assertAlmostEqual(fpd1, 10.0 / 9)
    self.assertLess(fpd1, fpd2)
    self.assertLess(fpd2, fpd3)

  @unittest.skipIf(not hasattr(rdMolDescriptors, 'BCUT2D'), "BCUT descriptor not available")
  def testVectorDescriptors(self):
    m = Chem.MolFromSmiles('CCCc1ccccc1')
    results = rdMolDescriptors.BCUT2D(m)
    names = [
      "BCUT2D_%s" % s
      for s in ('MWHI', "MWLOW", "CHGHI", "CHGLO", "LOGPHI", "LOGPLOW", "MRHI", "MRLOW")
    ]
    for i, n in enumerate(names):
      f = getattr(Descriptors, n)
      self.assertEqual(results[i], f(m))

    results = rdMolDescriptors.CalcAUTOCORR2D(m)
    names = ["AUTOCORR2D_%s" % str(i + 1) for i in range(192)]
    for i, n in enumerate(names):
      f = getattr(Descriptors, n)
      self.assertEqual(results[i], f(m))

  @unittest.skipIf(
    not hasattr(rdMolDescriptors, 'BCUT2D') or not hasattr(rdMolDescriptors, 'CalcAUTOCORR2D'),
    "BCUT or AUTOCORR descriptors not available")
  def testVectorDescriptorsInDescList(self):
    # First try only bcuts should exist
    descriptors = set([n for n, _ in Descriptors.descList])
    names = set([
      "BCUT2D_%s" % i
      for i in ('MWHI', "MWLOW", "CHGHI", "CHGLO", "LOGPHI", "LOGPLOW", "MRHI", "MRLOW")
    ])
    self.assertEqual(descriptors.intersection(names), names)

    Descriptors.setupAUTOCorrDescriptors()
    descriptors2 = set([n for n, _ in Descriptors.descList])
    self.assertEqual(descriptors2.intersection(descriptors), descriptors)

    names = set(["AUTOCORR2D_%s" % str(i + 1) for i in range(192)])
    self.assertEqual(descriptors2.intersection(names), names)

  def test_issue_4567(self):
    """
      Github issue #4567:
      Requesting "density" fingerprint Hydrogen molecule fails with exception
      """
    mol = AllChem.MolFromSmiles('[HH]')
    self.assertEqual(Descriptors.FpDensityMorgan1(mol), 0)

  def testGetMolDescriptors(self):
    mol = Chem.MolFromSmiles('CCCO')
    descs = Descriptors.CalcMolDescriptors(mol)
    self.assertTrue('MolLogP' in descs)
    self.assertEqual(descs['NumHDonors'], 1)

  def testGet3DMolDescriptors(self):
    mol = Chem.MolFromSmiles('CCCO')

    # check ValueError raised when no 3D coordinates supplied
    with self.assertRaises(ValueError):
      Descriptors3D.CalcMolDescriptors3D(mol)

    # test function returns expected outputs
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    descs = Descriptors3D.CalcMolDescriptors3D(mol)
    self.assertTrue('InertialShapeFactor' in descs)
    self.assertAlmostEqual(descs['PMI1'], 20.9582649071385, delta=1e-4)


if __name__ == '__main__':
  unittest.main()
