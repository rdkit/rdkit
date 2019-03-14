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


import io
import os.path
import unittest
import doctest

import numpy as np
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors
import pickle


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(Descriptors, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def testGithub1287(self):
    smis = ('CCC',)
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
    for (smiles, expected) in (("[NH4+]", "H4N+"),
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
                               ("[U+3]", "U+3"),):
      mol = Chem.MolFromSmiles(smiles)
      actual = AllChem.CalcMolFormula(mol)
      self.assertEqual(actual, expected)

  def testMQNDetails(self):
    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'MQNs_regress.pkl')
    refFile2 = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'MQNs_non_strict_regress.pkl')
    # figure out which definition we are currently using
    m = Chem.MolFromSmiles("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C")
    if Lipinski.NumRotatableBonds(m) == 2:
      refFile = refFile2

    with open(refFile, 'r') as intf:
      buf = intf.read().replace('\r\n', '\n').encode('utf-8')
      intf.close()
    with io.BytesIO(buf) as inf:
      pkl = inf.read()
    refData = pickle.loads(pkl, encoding='bytes')
    fn = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'aromat_regress.txt')
    ms = [x for x in Chem.SmilesMolSupplier(fn, delimiter='\t')]
    refData2 = []
    for i, m in enumerate(ms):
      mqns = rdMolDescriptors.MQNs_(m)
      refData2.append((m, mqns))
      if mqns != refData[i][1]:
        indices = [(j, x, y) for j, x, y in zip(range(len(mqns)), mqns, refData[i][1]) if x != y]
        print(i, Chem.MolToSmiles(m), indices)
      self.assertEqual(mqns, refData[i][1])

  def testMQN(self):
    m = Chem.MolFromSmiles("CC(C)(C)c1cc(O)c(cc1O)C(C)(C)C")
    if Lipinski.NumRotatableBonds(m) == 2:
      tgt = np.array(
        [42917, 274, 870, 621, 135, 1582, 29, 3147, 5463, 6999, 470, 62588, 19055, 4424, 309, 24061,
         17820, 1, 9303, 24146, 16076, 5560, 4262, 646, 746, 13725, 5430, 2629, 362, 24211, 15939,
         292, 41, 20, 1852, 5642, 31, 9, 1, 2, 3060, 1750])
    else:
      tgt = np.array(
        [42917, 274, 870, 621, 135, 1582, 29, 3147, 5463, 6999, 470, 62588, 19055, 4424, 309, 24061,
         17820, 1, 8314, 24146, 16076, 5560, 4262, 646, 746, 13725, 5430, 2629, 362, 24211, 15939,
         292, 41, 20, 1852, 5642, 31, 9, 1, 2, 3060, 1750])
    fn = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'aromat_regress.txt')
    ms = [x for x in Chem.SmilesMolSupplier(fn, delimiter='\t')]
    vs = np.zeros((42,), np.int32)
    for m in ms:
      vs += rdMolDescriptors.MQNs_(m)
    self.assertFalse(False in (vs == tgt))

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


if __name__ == '__main__':
  unittest.main()
