#
#  Copyright (C) 2002-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the signatures

"""
import os.path
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D import Generate, SigFactory


class TestCase(unittest.TestCase):

  def setUp(self):
    fdefFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'Pharm2D', 'test_data', 'BaseFeatures.fdef')
    featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
    self.factory = SigFactory.SigFactory(featFactory, minPointCount=2, maxPointCount=3,
                                         trianglePruneBins=False)
    self.factory.SetBins([(0, 2), (2, 5), (5, 8)])
    self.factory.Init()
    SigFactory._verbose = False

  def test1Sizes(self):
    self.factory.maxPointCount = 2
    self.factory.Init()
    sig = self.factory.GetSignature()
    self.assertEqual(len(sig), 45)

    self.factory.maxPointCount = 3
    self.factory.Init()
    sig = self.factory.GetSignature()
    self.assertEqual(len(sig), 990)

    self.factory.maxPointCount = 4
    self.factory.Init()
    sig = self.factory.GetSignature()
    self.assertEqual(len(sig), 18000)

  def test2BitIdx(self):
    data = [
      ((0, 0), [0], 0),
      ((0, 0), [2], 1),
      ((0, 0), [5], 2),
      ((0, 1), [5], 5),
      ((1, 1), [4], 16),
      ((1, 1), [7], 17),
      ((0, 0, 0), [1, 1, 1], 45),
      ((0, 0, 1), [1, 1, 1], 72),
      ((0, 0, 1), [1, 1, 3], 75),
      ((0, 0, 1), [3, 1, 1], 81),
      ((0, 0, 1), [3, 3, 1], 84),
    ]
    for tpl in data:
      patts, dists, bit = tpl
      idx = self.factory.GetBitIdx(patts, dists)
      self.assertEqual(bit, idx)

      cnt, feats, _ = self.factory.GetBitInfo(bit)
      self.assertEqual(cnt, len(patts))
      self.assertEqual(feats, patts)

  def test3BitIdx(self):
    # test 3 point p'cophore ids, you can never have too much of this stuff
    self.factory.SetBins(((0, 2), (2, 4), (4, 8)))
    self.factory.Init()
    self.assertEqual(self.factory.GetSigSize(), 990)
    probes = [((0, 0, 0), (1, 3, 1), 54),
              ((0, 0, 0), (3, 1, 1), 54),
              ((0, 0, 0), (1, 1, 3), 54),
              ((0, 0, 0), (1, 3, 3), 57),
              ((0, 0, 1), (1, 3, 1), 75), ]
    for patts, dists, ans in probes:
      idx = self.factory.GetBitIdx(patts, dists)
      self.assertEqual(idx, ans)

      cnt, feats, _ = self.factory.GetBitInfo(ans)
      self.assertEqual(cnt, len(patts))
      self.assertEqual(feats, patts)

  def test4BitIdx(self):
    self.factory.trianglePruneBins = True
    self.factory.Init()
    sig = self.factory.GetSignature()
    self.assertEqual(len(sig), 885)

    probes = [((0, 0, 0), (1, 3, 1), 52),
              ((0, 0, 0), (1, 1, 3), 52),
              ((0, 0, 0), (3, 1, 1), 52),
              ((0, 0, 0), (1, 3, 3), 55),
              ((0, 0, 1), (1, 3, 1), 71), ]
    for patts, dists, ans in probes:
      idx = self.factory.GetBitIdx(patts, dists)
      self.assertEqual(idx, ans)
      cnt, feats, _ = self.factory.GetBitInfo(ans)
      self.assertEqual(cnt, len(patts))
      self.assertEqual(feats, patts)

    # Test exception handlers
    self.assertRaises(NotImplementedError, self.factory.GetBitIdx, (0, 0, 0, 0), (1, 3, 1, 2))
    self.assertRaises(IndexError, self.factory.GetBitIdx, (0, ), (1, ))
    self.factory.maxPointCount = 2
    self.assertRaises(IndexError, self.factory.GetBitIdx, (0, 0, 0), (1, 3, 1))
    self.factory.maxPointCount = 3
    self.assertRaises(IndexError, self.factory.GetBitIdx, (0, -1, 0), (1, 3, 1))
    self.assertRaises(IndexError, self.factory.GetBitIdx, (0, self.factory._nFeats, 0), (1, 3, 1))

  def test5SimpleSig(self):
    factory = self.factory
    factory.SetBins([(1, 3), (3, 7), (7, 10)])
    factory.minPointCount = 2
    factory.maxPointCount = 3
    factory.Init()

    mol = Chem.MolFromSmiles('O=CCC=O')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(len(sig), 990)
    bs = tuple(sig.GetOnBits())
    self.assertEqual(bs, (1, ))

    mol = Chem.MolFromSmiles('O=CC(CC=O)CCC=O')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(len(sig), 990)
    bs = tuple(sig.GetOnBits())
    self.assertEqual(bs, (1, 2, 67))

  def test6SimpleSigCounts(self):
    factory = self.factory
    factory.SetBins([(1, 3), (3, 7), (7, 10)])
    factory.minPointCount = 2
    factory.maxPointCount = 3
    factory.useCounts = True
    factory.Init()

    mol = Chem.MolFromSmiles('O=CCC=O')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(sig.GetLength(), 990)
    cs = tuple(sig.GetNonzeroElements().items())
    self.assertEqual(cs, ((1, 1), ))

    mol = Chem.MolFromSmiles('O=CC(CC=O)CCC=O')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(sig.GetLength(), 990)
    elems = sig.GetNonzeroElements()
    bs = list(elems.keys())
    bs.sort()
    cs = [(x, elems[x]) for x in bs]
    self.assertEqual(tuple(cs), ((1, 2), (2, 1), (67, 1)))

  def test7SimpleSigSkip(self):
    factory = self.factory
    factory.SetBins([(1, 3), (3, 7), (7, 10)])
    factory.minPointCount = 2
    factory.maxPointCount = 3
    factory.skipFeats = 'Acceptor'
    factory.Init()

    mol = Chem.MolFromSmiles('O=CCC=O')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(len(sig), 570)
    bs = tuple(sig.GetOnBits())
    self.assertEqual(bs, ())

  def test8MultiPointMatches(self):
    factory = self.factory
    factory.SetBins([(1, 3), (3, 7), (7, 10)])
    factory.minPointCount = 2
    factory.maxPointCount = 3
    factory.Init()

    mol = Chem.MolFromSmiles('O=Cc1ccccc1')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(len(sig), 990)
    bs = tuple(sig.GetOnBits())
    self.assertEqual(bs, (3, ))

    mol = Chem.MolFromSmiles('O=CCCCCCCCCc1ccccc1')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(len(sig), 990)
    bs = tuple(sig.GetOnBits())
    self.assertEqual(bs, ())

  # FIX: add test for perms argument to Gen2DFingerprint

  def test9BondOrderSigs(self):
    # test sigs where bond order is used
    factory = self.factory
    factory.SetBins([(1, 4), (4, 7), (7, 10)])
    factory.minPointCount = 2
    factory.maxPointCount = 3
    factory.Init()

    mol = Chem.MolFromSmiles('[O-]CCC(=O)')
    sig = Generate.Gen2DFingerprint(mol, self.factory)
    self.assertEqual(len(sig), 990)
    bs = tuple(sig.GetOnBits())
    self.assertEqual(bs, (1, ))

    self.factory.includeBondOrder = True
    sig = Generate.Gen2DFingerprint(mol, self.factory)
    self.assertEqual(len(sig), 990)
    bs = tuple(sig.GetOnBits())
    self.assertEqual(bs, (0, ))

  def testDefaultFactory(self):
    from rdkit.Chem import Pharm2D
    factory = Pharm2D.DefaultSigFactory()
    self.assertEqual(factory.GetNumBins(), 7)
    
    # Generate._verbose=True
    mol = Chem.MolFromSmiles('OCCC(=O)')
    sig = Generate.Gen2DFingerprint(mol, factory)
    self.assertEqual(len(sig), 19355)
    self.assertEqual(tuple(sig.GetOnBits()), (2,
                                              16,
                                              21,
                                              84,
                                              1274,
                                              4361, ))
    nPts, combo, scaffold, labels, dMat = factory._GetBitSummaryData(21)
    self.assertEqual(nPts, 2)
    self.assertEqual(labels, ['Acceptor', 'Hydrophobe'])
    self.assertEqual(list(dMat[0]), [0, 0])
    self.assertEqual(list(dMat[1]), [0, 0])

    txt = factory.GetBitDescription(21)
    self.assertEqual(txt, 'Acceptor Hydrophobe |0 0|0 0|')

    txt = factory.GetBitDescription(21)
    self.assertEqual(txt, 'Acceptor Hydrophobe |0 0|0 0|')

    nPts, combo, scaffold, labels, dMat = factory._GetBitSummaryData(2)
    self.assertEqual(nPts, 2)
    self.assertEqual(labels, ['Acceptor', 'Acceptor'])
    self.assertEqual(list(dMat[0]), [0, 2])
    self.assertEqual(list(dMat[1]), [2, 0])

    nPts, combo, scaffold, labels, dMat = factory._GetBitSummaryData(4361)
    self.assertEqual(nPts, 3)
    self.assertEqual(labels, ['Acceptor', 'Donor', 'Hydrophobe'])
    self.assertEqual(list(dMat[0]), [0, 2, 0])
    self.assertEqual(list(dMat[1]), [2, 0, 0])
    self.assertEqual(list(dMat[2]), [0, 0, 0])
    self.assertEqual(
      factory.GetBitDescription(4361), 'Acceptor Donor Hydrophobe |0 2 0|2 0 0|0 0 0|')

    self.assertRaises(NotImplementedError, factory.GetBitDescriptionAsText, 21)

if __name__ == '__main__':  # pragma: nocover
  unittest.main()
