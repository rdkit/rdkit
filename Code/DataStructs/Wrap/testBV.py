import pickle
import random
import unittest

import numpy

from rdkit import DataStructs, RDConfig


def feq(a, b, tol=1e-4):
  return abs(a - b) < tol


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test0FromList(self):
    bv1 = DataStructs.SparseBitVect(1000)
    bv2 = DataStructs.SparseBitVect(1000)
    obits = range(0, 1000, 3)

    for bit in obits:
      bv1.SetBit(bit)

    bv2.SetBitsFromList(obits)

    for i in range(1000):
      assert bv1.GetBit(i) == bv2.GetBit(i)

    self.assertTrue(bv1 == bv2)
    bv2.SetBit(1)
    self.assertTrue(bv1 != bv2)
    bv2.UnSetBit(1)
    self.assertTrue(bv1 == bv2)

    bv2.UnSetBitsFromList(obits)
    for i in range(1000):
      assert bv2.GetBit(i) == 0

    bv1 = DataStructs.ExplicitBitVect(1000)
    bv2 = DataStructs.ExplicitBitVect(1000)
    obits = range(0, 1000, 3)

    for bit in obits:
      bv1.SetBit(bit)

    bv2.SetBitsFromList(obits)

    for i in range(1000):
      assert bv1.GetBit(i) == bv2.GetBit(i)

    bv2.UnSetBitsFromList(obits)
    for i in range(1000):
      assert bv2.GetBit(i) == 0

  def test01BVWithAllOnes(self):
    bv1 = DataStructs.ExplicitBitVect(10, True)
    for i in range(10):
      assert bv1.GetBit(i) == 1

  def test1SparsePickle(self):
    nbits = 10000
    bv1 = DataStructs.SparseBitVect(nbits)
    for i in range(1000):
      x = random.randrange(0, nbits)
      bv1.SetBit(x)

    pkl = pickle.dumps(bv1, 1)
    bv2 = pickle.loads(pkl)
    for i in range(nbits):
      assert bv1[i] == bv2[i]

  def test2ExplicitPickle(self):
    nbits = 10000
    bv1 = DataStructs.ExplicitBitVect(nbits)
    for i in range(1000):
      x = random.randrange(0, nbits)
      bv1.SetBit(x)

    pkl = pickle.dumps(bv1, 1)
    bv2 = pickle.loads(pkl)
    for i in range(nbits):
      assert bv1[i] == bv2[i]

  def test3Bounds(self):
    nbits = 10
    bv1 = DataStructs.ExplicitBitVect(nbits)
    bv1[0]
    with self.assertRaisesRegex(IndexError, ""):
      bv1[11]

  def test4OnBitsInCommon(self):
    sz = 100
    bv1 = DataStructs.ExplicitBitVect(sz)
    bv2 = DataStructs.ExplicitBitVect(sz)
    for i in range(0, sz, 2):
      bv1.SetBit(i)
      if i < 3 * sz / 4:
        bv2.SetBit(i)
    self.assertTrue(DataStructs.AllProbeBitsMatch(bv1, bv1.ToBinary()))
    self.assertTrue(DataStructs.AllProbeBitsMatch(bv2, bv1.ToBinary()))
    self.assertFalse(DataStructs.AllProbeBitsMatch(bv1, bv2.ToBinary()))
    self.assertTrue(DataStructs.AllProbeBitsMatch(bv2, bv2.ToBinary()))

  def test5FromBitString(self):
    s1 = '1010'
    bv = DataStructs.CreateFromBitString(s1)
    self.assertTrue(len(bv) == 4)
    self.assertTrue(list(bv.GetOnBits()) == [0, 2])

  def _bulkTest(self, bvs):
    for metric in 'Tanimoto', 'Dice', 'AllBit', 'OnBit', 'RogotGoldberg':
      bulk = getattr(DataStructs, f'Bulk{metric}Similarity')
      single = getattr(DataStructs, f'{metric}Similarity')
    sims = bulk(bvs[0], bvs)
    for i, bbvs in enumerate(bvs):
      sim = single(bvs[0], bbvs)
      self.assertEqual(sim, sims[i])
      self.assertEqual(sim, single(bvs[0], bbvs.ToBinary()))
    dists = bulk(bvs[0], bvs, returnDistance=True)
    for i, bbvs in enumerate(bvs):
      dist = single(bvs[0], bbvs, returnDistance=True)
      self.assertEqual(dist, dists[i])
      self.assertEqual(dist, single(bvs[0], bbvs.ToBinary(), returnDistance=True))

    sims = DataStructs.BulkTverskySimilarity(bvs[0], bvs, 1, 1)
    for i, bbvs in enumerate(bvs):
      sim = DataStructs.TverskySimilarity(bvs[0], bbvs, 1, 1)
      self.assertEqual(sim, sims[i])
      sim = DataStructs.TanimotoSimilarity(bvs[0], bbvs)
      self.assertEqual(sim, sims[i])

    sims = DataStructs.BulkTverskySimilarity(bvs[0], bvs, 1, 1, returnDistance=True)
    for i, bbvs in enumerate(bvs):
      sim = DataStructs.TverskySimilarity(bvs[0], bbvs, 1, 1, returnDistance=True)
      self.assertEqual(sim, sims[i])
      sim = DataStructs.TanimotoSimilarity(bvs[0], bbvs, returnDistance=True)
      self.assertEqual(sim, sims[i])

  def test6BulkOps(self):
    nbits = 10000
    bvs = []
    for _ in range(10):
      bv = DataStructs.ExplicitBitVect(nbits)
      for _ in range(nbits):
        x = random.randrange(0, nbits)
        bv.SetBit(x)
      bvs.append(bv)
    self._bulkTest(bvs)

  def test10BulkOps2(self):
    nbits = 10000
    bvs = []
    for _ in range(10):
      bv = DataStructs.ExplicitBitVect(nbits)
      for _ in range(nbits):
        x = random.randrange(0, nbits)
        bv.SetBit(x)
      bvs.append(bv)
    bvs = tuple(bvs)
    self._bulkTest(bvs)

  def test10BulkOps3(self):
    nbits = 10000
    bvs = numpy.empty((10, ), DataStructs.ExplicitBitVect)
    for bvi in range(10):
      bv = DataStructs.ExplicitBitVect(nbits)
      for _ in range(nbits):
        x = random.randrange(0, nbits)
        bv.SetBit(x)
      bvs[bvi] = bv
    self._bulkTest(bvs)

  def test7FPS(self):
    bv = DataStructs.ExplicitBitVect(32)
    bv.SetBit(0)
    bv.SetBit(1)
    bv.SetBit(17)
    bv.SetBit(23)
    bv.SetBit(31)

    self.assertEqual(DataStructs.BitVectToFPSText(bv), "03008280")
    bv2 = DataStructs.CreateFromFPSText("03008280")
    self.assertEqual(bv, bv2)

    self.assertRaises(ValueError, lambda: DataStructs.CreateFromFPSText("030082801"))

    bv2 = DataStructs.CreateFromFPSText("")
    self.assertEqual(bv2.GetNumBits(), 0)

  def test8BinText(self):
    bv = DataStructs.ExplicitBitVect(32)
    bv.SetBit(0)
    bv.SetBit(1)
    bv.SetBit(17)
    bv.SetBit(23)
    bv.SetBit(31)

    bv2 = DataStructs.CreateFromBinaryText(DataStructs.BitVectToBinaryText(bv))
    self.assertEqual(bv, bv2)

    bv2 = DataStructs.CreateFromBinaryText("")
    self.assertEqual(bv2.GetNumBits(), 0)

  def test9ToNumpy(self):
    import numpy
    for typ in (DataStructs.ExplicitBitVect, ):
      bv = typ(32)
      bv.SetBit(0)
      bv.SetBit(1)
      bv.SetBit(17)
      bv.SetBit(23)
      bv.SetBit(31)
      arr = numpy.zeros((32, ), 'i')
      DataStructs.ConvertToNumpyArray(bv, arr)
      for i in range(bv.GetNumBits()):
        self.assertEqual(bv[i], arr[i])

    for typ in (DataStructs.IntSparseIntVect, DataStructs.LongSparseIntVect,
                DataStructs.UIntSparseIntVect, DataStructs.ULongSparseIntVect):
      iv = typ(32)
      iv[0] = 1
      iv[1] = 1
      iv[17] = 1
      iv[23] = 1
      iv[31] = 1
      arr = numpy.zeros((32, ), 'i')
      DataStructs.ConvertToNumpyArray(iv, arr)
      for i in range(iv.GetLength()):
        self.assertEqual(iv[i], arr[i])

  def test11BulkNeighbors(self):
    nbits = 2048
    bvs = []
    for _ in range(1000):
      bv = DataStructs.ExplicitBitVect(nbits)
      for _ in range(nbits):
        x = random.randrange(0, nbits)
        bv.SetBit(x)
      bvs.append(bv)
    qs = bvs[:10]
    db = bvs[10:]
    for metric in [
        'Tanimoto', 'Cosine', 'Kulczynski', 'Dice', 'Sokal', 'McConnaughey', 'Asymmetric',
        'BraunBlanquet', 'Russel', 'RogotGoldberg'
    ]:
      bulkSim = getattr(DataStructs, f'Bulk{metric}Similarity')
      nbrSim = getattr(DataStructs, f'{metric}SimilarityNeighbors')
      tgts = []
      for q in qs:
        sims = bulkSim(q, db)
        sim, idx = max((sim, -idx) for idx, sim in enumerate(sims))
        tgts.append((-idx, sim))
      nbrs = nbrSim(qs, db)
      self.assertEqual(tgts, nbrs)

  def test12ToList(self):
    nbits = 2048
    for cls in [DataStructs.ExplicitBitVect, DataStructs.SparseBitVect]:
      bv = cls(nbits)
      l = [0] * 2048

      # test no bits set
      l2 = list(bv)
      l3 = bv.ToList()
      self.assertEqual(l, l2)
      self.assertEqual(l, l3)

      for _ in range(nbits):
        x = random.randrange(0, nbits)
        l[x] = 1
        bv.SetBit(x)

      l2 = list(bv)
      l3 = bv.ToList()
      self.assertEqual(l, l2)
      self.assertEqual(l, l3)

  def test13Base64(self):
    nbits = 2048
    for cls in [DataStructs.ExplicitBitVect, DataStructs.SparseBitVect]:
      bv = cls(nbits)
      for _ in range(nbits):
        x = random.randrange(0, nbits)
        bv.SetBit(x)

      bv2 = cls(nbits)
      bv2.FromBase64(bv.ToBase64())
      self.assertEqual(bv, bv2)

  def test14NegativeIndices(self):
    nbits = 2048
    for cls in [DataStructs.ExplicitBitVect, DataStructs.SparseBitVect]:
      bv = cls(nbits)
      bv2 = cls(nbits)
      for _ in range(nbits):
        x = random.randrange(0, nbits)
        bv.SetBit(x)
        bv2[-(nbits - x)] = 1

      self.assertEqual(bv, bv2)
      for j in range(nbits):
        self.assertEqual(bv[j], bv[-(nbits - j)])
      with self.assertRaises(IndexError):
        bv[-(nbits + 1)]
      with self.assertRaises(IndexError):
        bv2[-(nbits + 1)] = 1

  def test15FoldFingerprint(self):
    for cls in [DataStructs.ExplicitBitVect, DataStructs.SparseBitVect]:
      fp = cls(8)
      fp[0] = 1
      fp[1] = 1
      fp[6] = 1
      ffp = DataStructs.FoldFingerprint(fp)
      self.assertTrue(ffp[0])
      self.assertTrue(ffp[1])
      self.assertTrue(ffp[2])
      self.assertFalse(ffp[3])


if __name__ == '__main__':
  unittest.main()
