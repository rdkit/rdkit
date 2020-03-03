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

import unittest
import doctest
from rdkit.DataStructs import BitUtils, VectCollection, LazySignature, FingerprintSimilarity
from rdkit import DataStructs, Chem
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprinterDetails


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(BitUtils, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(LazySignature, optionflags=doctest.ELLIPSIS))
  tests.addTests(doctest.DocTestSuite(VectCollection, optionflags=doctest.ELLIPSIS))
  return tests


class TestCaseAdditional(unittest.TestCase):

  def test_VectCollection(self):
    # We mainly test the use of Reset
    bv1 = DataStructs.ExplicitBitVect(10)
    bv1.SetBitsFromList((1, 3, 5))
    bv2 = DataStructs.ExplicitBitVect(10)
    bv2.SetBitsFromList((6, 8))

    vc = VectCollection.VectCollection()
    self.assertEqual(vc.GetOrVect(), None)
    vc.AddVect(1, bv1)
    vc.AddVect(2, bv2)

    onBits = set([1, 3, 5, 6, 8])
    for i, onOff in enumerate(vc.GetOrVect()):
      self.assertEqual(i in onBits, onOff == 1)
    vc.Reset()
    self.assertEqual(onBits, set(vc.GetOnBits()))

    vc = VectCollection.VectCollection()
    self.assertEqual(vc.GetOrVect(), None)
    vc.AddVect(1, bv1)
    vc.AddVect(2, bv2)
    for i in onBits:
      self.assertEqual(vc[i], 1)

  def test_LazySig(self):
    self.assertRaises(ValueError, LazySignature.LazySig, lambda x: 1, 0)

    # Check cache works
    obj = LazySignature.LazySig(lambda x: x, 10)
    self.assertNotIn(8, obj._cache)
    self.assertEqual(obj[8], 8)
    self.assertIn(8, obj._cache)
    obj._cache[8] = 'cached'
    self.assertEqual(obj[8], 'cached')

  def test__init__(self):
    from rdkit.Chem.Fingerprints import FingerprintMols
    ms = [Chem.MolFromSmiles('CCOC'), Chem.MolFromSmiles('CCO'), Chem.MolFromSmiles('COC')]
    fps = [FingerprintMols.FingerprintMol(x) for x in ms]
    self.assertAlmostEqual(FingerprintSimilarity(fps[0], fps[1]), 0.6, places=2)

    details = FingerprinterDetails()
    fpArgs = details.__dict__
    fps = []
    for i, x in enumerate(ms, 1):
      fpArgs['fpSize'] = 16 * i
      fps.append(FingerprintMols.FingerprintMol(x, **fpArgs))
    self.assertAlmostEqual(FingerprintSimilarity(fps[0], fps[1]), 0.555, places=2)
    self.assertAlmostEqual(FingerprintSimilarity(fps[1], fps[0]), 0.555, places=2)

    fpArgs['fpSize'] = 1024
    fpArgs['tgtDensity'] = 0.8
    fp = FingerprintMols.FingerprintMol(ms[0], **fpArgs)
    self.assertEqual(len(fp), 64)
    fp = DataStructs.FoldToTargetDensity(fp, density=0.1, minLength=2)
    self.assertEqual(len(fp), 4)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
