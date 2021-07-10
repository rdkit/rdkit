from rdkit import RDConfig
import unittest, os
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.DataManip.Metric import rdMetricMatrixCalc as rdmmc
from rdkit import DataStructs
import numpy
import random


class TestCase(unittest.TestCase):

  def setUp(self):
    self.n = 1000
    self.m = 80
    self.d = 2
    self.dataPts = numpy.zeros((self.n, self.d), 'd')
    for i in range(self.n):
      for j in range(self.d):
        self.dataPts[i, j] = random.random()
    self.dMat = rdmmc.GetEuclideanDistMat(self.dataPts)

  def test0MaxMin(self):
    pkr = rdSimDivPickers.MaxMinPicker()
    maxmin = pkr.Pick(self.dMat, self.n, self.m, (886, 112))
    self.assertEqual(maxmin[0], 886)
    self.assertEqual(maxmin[1], 112)

    def func(i, j):
      if i == j:
        return 0.0
      if i < j:
        j, i = i, j
      return self.dMat[i * (i - 1) // 2 + j]

    lmaxmin = pkr.LazyPick(func, self.n, self.m, (886, 112))
    self.assertEqual(list(lmaxmin), list(maxmin))

    lmaxmin = pkr.LazyPick(func, self.n, self.m, (886, 112), useCache=False)
    self.assertEqual(list(lmaxmin), list(maxmin))

    self.assertRaises(ValueError, lambda: pkr.Pick(self.dMat, self.n, self.m, (1012, )))
    self.assertRaises(ValueError, lambda: pkr.Pick(self.dMat, self.n, self.m, (-1, )))

    maxmin = pkr.Pick(self.dMat, self.n, self.m)
    self.assertTrue(maxmin)
    lmaxmin = pkr.LazyPick(func, self.n, self.m)
    self.assertTrue(lmaxmin)

  def test1HierarchPick(self):
    fname = os.path.join(RDConfig.RDBaseDir, 'Code', 'SimDivPickers', 'Wrap', 'test_data',
                         'points.csv')
    with open(fname) as infil:
      lines = infil.readlines()
    self.dataPts = numpy.zeros((len(lines), 2), 'd')
    labels = []
    i = 0
    for line in lines:
      tlst = line.strip().split(',')
      self.dataPts[i, 0] = float(tlst[1])
      self.dataPts[i, 1] = float(tlst[2])
      labels.append(int(tlst[3]))
      i += 1
    self.dMat = rdmmc.GetEuclideanDistMat(self.dataPts)
    pkr = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    clusters = pkr.Cluster(self.dMat, i, 2)
    # check that each of the clusters have the same label
    for cl in clusters:
      clbl = labels[cl[0]]
      for id in cl:
        assert clbl == labels[id]
    hierarch = pkr.Pick(self.dMat, i, 2)
    self.assertEqual(tuple(hierarch), (1, 30))

  def testIssue208(self):
    sz = 10
    N = 3
    m = []
    for i in range(sz):
      for j in range(i + 1, sz):
        m.append(random.random())
    m = numpy.array(m)
    picker = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    p1 = list(picker.Pick(m, sz, N))
    p1.sort()
    p2 = list(picker.Pick(m, sz, N))
    p2.sort()
    self.assertEqual(p1, p2)

  def testInts(self):
    """ make sure we can handle ints too """
    sz = 10
    N = 3
    m = []
    for i in range(sz):
      for j in range(i + 1, sz):
        m.append(int(100 * random.random()))
    m = numpy.array(m)
    picker = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    p1 = list(picker.Pick(m, sz, N))
    p1.sort()
    p2 = list(picker.Pick(m, sz, N))
    p2.sort()
    self.assertEqual(p1, p2)

  def testNonUniqueCrash(self):
    from rdkit import DataStructs
    sz = 300
    nbits = 40
    nBitsToSet = int(nbits * .3)
    N = 8
    vs = []
    for i in range(sz):
      bv = DataStructs.ExplicitBitVect(nbits)
      for j in range(nBitsToSet):
        val = int(nbits * random.random())
        bv.SetBit(val)
      vs.append(bv)
      vs.append(bv)

    def taniFunc(i, j, bvs=vs):
      d = 1 - DataStructs.FingerprintSimilarity(bvs[i], bvs[j])
      return d

    picker = rdSimDivPickers.MaxMinPicker()
    mm1 = picker.LazyPick(taniFunc, len(vs), N)
    self.assertEqual(len(mm1), N)
    picker = None

    picker = rdSimDivPickers.MaxMinPicker()
    mm2 = picker.LazyBitVectorPick(vs, len(vs), N)
    self.assertEqual(len(mm2), N)

    picker = rdSimDivPickers.MaxMinPicker()
    mm3 = picker.LazyBitVectorPick(vs, len(vs), N)
    self.assertEqual(len(mm3), N)

    # we get the occasional dupe randomly,
    # make sure we don't get three dupes in a row
    self.assertTrue(tuple(mm2) != tuple(mm1)) or (tuple(mm3) != tuple(mm1))
    picker = None

    ds = []
    nvs = len(vs)
    for i in range(nvs):
      for j in range(i + 1, nvs):
        d = taniFunc(i, j)
        ds.append(d)
    m = numpy.array(ds)
    picker = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    p1 = list(picker.Pick(m, nvs, N))

  def testBitVectorMaxMin(self):
    from rdkit import DataStructs
    sz = 100
    nbits = 200
    nBitsToSet = int(nbits * .1)
    N = 10
    vs = []
    for i in range(sz):
      bv = DataStructs.ExplicitBitVect(nbits)
      for j in range(nBitsToSet):
        val = int(nbits * random.random())
        bv.SetBit(val)
      vs.append(bv)

    def func(i, j, bvs=vs):
      d = DataStructs.TanimotoSimilarity(bvs[i], bvs[j], returnDistance=True)
      return d

    picker = rdSimDivPickers.MaxMinPicker()
    mm1 = picker.LazyPick(func, len(vs), N, seed=42)
    self.assertEqual(len(mm1), N)

    mm2 = picker.LazyPick(func, len(vs), N, useCache=False, seed=42)
    self.assertEqual(len(mm2), N)
    self.assertEqual(list(mm1), list(mm2))

    mm2 = picker.LazyBitVectorPick(vs, len(vs), N, seed=42)
    self.assertEqual(len(mm2), N)
    self.assertEqual(list(mm1), list(mm2))

    mm2 = picker.LazyBitVectorPick(vs, len(vs), N, useCache=False, seed=42)
    self.assertEqual(len(mm2), N)
    self.assertEqual(list(mm1), list(mm2))

  def testBitVectorMaxMin2(self):
    fps = [
      "11110010101000000000", "00000000000010010000", "11001010000000000001",
      "00100110101000001000", "01010110000100011001", "11000110101001000011",
      "00000000001100001111", "00011110110000001101", "00000011011110100010",
      "11000010110001000000", "00000100010000010000", "10000001000010110010",
      "00010010000000010100", "00011100100110101000", "10001001100110100000",
      "10000110100110010000", "00101110000101000000", "11011101100011100000",
      "10000110000100101000", "00101000100000010001", "01000001000010000000",
      "00101101010100000110", "10001000100110110001", "00011000010100000001",
      "00101000001000100011", "00010000100010011001", "01100001000100010001",
      "10000101000001101101", "00001000011001011000", "11110000100100100000",
      "10100110000000011010", "00110100010110010010", "00000000000001010010",
      "00100000000010100001", "11110011000010001000", "10110001010100001000",
      "00001100100110011011", "00010010100100001110", "10100101100010100010",
      "01100100010100000001", "10101110011100000000", "01011000000001000001",
      "00000011100110100010", "01100001010001001001", "00001000000001001100",
      "10011001110000000100", "10110000001001100100", "00011000000001001011",
      "11001011010001100010", "10010000000001001011", "00010000100111100000",
      "00001000001110001000", "11010000010001100110", "01101001100000111000",
      "01001000001110111000", "10000000000100010010", "11001000010010000000",
      "01010010000100110001", "00010001010100100001", "01110010000000010000",
      "10001010000011000001", "00000110000000100100", "00010000010001000000",
      "11101100011010000011", "00000010100001010001", "00010000110010000101",
      "00010001001000111001", "01000010001100100110", "00110110000000100001",
      "00100010010110110010", "01000000110011001111", "00011000001000110010",
      "01111010101000110100", "00001010000010110110", "00110011000011011010",
      "00111010111010000110", "00010011101010000011", "00000001011000010000",
      "00011011101110110000", "00010001101000000001", "00010000001010011010",
      "00000010100100100010", "00000010001011000100", "11010000000001011100",
      "00001000110101000001", "00000010000000110010", "10000000010011000001",
      "11110110100100010000", "10001111000110001001", "00100110000110000100",
      "00000100100000100100", "00110000101100010100", "00001010100000100000",
      "01011000000011000111", "00010000100001010001", "10000010100000010000",
      "00001000000000110010", "00001000101011010001", "00011110000100100000", "11001001010001010100"
    ]
    N = 5
    fps = [DataStructs.CreateFromBitString(x) for x in fps]
    picker = rdSimDivPickers.MaxMinPicker()
    mm1 = picker.LazyBitVectorPick(fps, len(fps), N, seed=42)
    self.assertEqual(len(mm1), N)
    self.assertEqual(list(mm1), [37, 1, 43, 38, 16])

    mm2 = picker.LazyBitVectorPick(fps, len(fps), N, useCache=False, seed=42)
    self.assertEqual(len(mm2), N)
    self.assertEqual(list(mm1), list(mm2))

  def testBitVectorMaxMin3(self):
    fname = os.path.join(RDConfig.RDBaseDir, 'Code', 'SimDivPickers', 'Wrap', 'test_data',
                         'chembl_cyps.head.fps')
    fps = []
    with open(fname) as infil:
      for line in infil:
        fp = DataStructs.CreateFromFPSText(line.strip())
        fps.append(fp)
    mmp = rdSimDivPickers.MaxMinPicker()
    ids = list(mmp.LazyBitVectorPick(fps, len(fps), 20, seed=42))
    self.assertEqual(ids,[374,720,690,339,875,842,404,725,120,385,115,868,630,
                          881,516,497,412,718,869,407])

    ids = list(
      mmp.LazyBitVectorPick(fps, len(fps), 20, firstPicks=[374, 720, 690, 339, 875], seed=42))
    self.assertEqual(ids,[374,720,690,339,875,842,404,725,120,385,115,868,630,
                          881,516,497,412,718,869,407])

  def testBitVectorMaxMin4(self):
    # threshold tests
    fname = os.path.join(RDConfig.RDBaseDir, 'Code', 'SimDivPickers', 'Wrap', 'test_data',
                         'chembl_cyps.head.fps')
    fps = []
    with open(fname) as infil:
      for line in infil:
        fp = DataStructs.CreateFromFPSText(line.strip())
        fps.append(fp)
    mmp = rdSimDivPickers.MaxMinPicker()
    ids, threshold = mmp.LazyBitVectorPickWithThreshold(fps, len(fps), 20, -1.0, seed=42)
    self.assertEqual(list(ids),[374,720,690,339,875,842,404,725,120,385,115,868,630,
                          881,516,497,412,718,869,407])

    self.assertAlmostEqual(threshold, 0.8977, 4)

    ids, threshold = mmp.LazyBitVectorPickWithThreshold(fps, len(fps), 20, 0.91, seed=42)
    self.assertEqual(list(ids), [374, 720, 690, 339, 875, 842, 404, 725, 120, 385, 115, 868, 630])
    self.assertTrue(threshold >= 0.91)

  def testBitVectorLeader1(self):
    # threshold tests
    fname = os.path.join(RDConfig.RDBaseDir, 'Code', 'SimDivPickers', 'Wrap', 'test_data',
                         'chembl_cyps.head.fps')
    fps = []
    with open(fname) as infil:
      for line in infil:
        fp = DataStructs.CreateFromFPSText(line.strip())
        fps.append(fp)
    mmp = rdSimDivPickers.LeaderPicker()
    thresh = 0.8
    ids = mmp.LazyBitVectorPick(fps, len(fps), thresh)
    self.assertEqual(len(ids), 146)
    for i in range(len(ids)):
      for j in range(i):
        self.assertGreaterEqual(1 - DataStructs.TanimotoSimilarity(fps[ids[i]], fps[ids[j]]),
                                thresh)
    thresh = 0.9
    ids = mmp.LazyBitVectorPick(fps, len(fps), thresh)
    self.assertEqual(len(ids), 14)
    for i in range(len(ids)):
      for j in range(i):
        self.assertGreaterEqual(1 - DataStructs.TanimotoSimilarity(fps[ids[i]], fps[ids[j]]),
                                thresh)

    ids = mmp.LazyBitVectorPick(fps, len(fps), thresh, pickSize=10)
    self.assertEqual(len(ids), 10)
    for i in range(len(ids)):
      for j in range(i):
        self.assertGreaterEqual(1 - DataStructs.TanimotoSimilarity(fps[ids[i]], fps[ids[j]]),
                                thresh)

  def testLazyLeader(self):
    pkr = rdSimDivPickers.LeaderPicker()

    def func(i, j):
      if i == j:
        return 0.0
      if i < j:
        j, i = i, j
      return i - j

    lres = pkr.LazyPick(func, 100, 20)
    self.assertEqual(list(lres), [0, 21, 42, 63, 84])


if __name__ == '__main__':
  unittest.main()
