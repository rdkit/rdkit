# testRealValueVect.py
#
# 14.04.2014 by David Hahn, hahnda6

from __future__ import print_function, unicode_literals
from rdkit import RDConfig
import os, sys, pickle
import unittest
from rdkit import DataStructs as ds


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1Real(self):
    v1 = ds.RealValueVect(30)
    for i in range(30):
      v1[i] = i / 10.0

    self.assertTrue(len(v1) == 30)
    self.assertAlmostEqual(v1.GetTotalVal(), 43.5)

    for i in range(len(v1)):
      self.assertAlmostEqual(v1[i], 0.1 * i)

  def test2RealVectDistances(self):
    v1 = ds.RealValueVect(30)
    v2 = ds.RealValueVect(30)
    for i in range(15):
      v1[2 * i] = 1.3
      v2[2 * i] = 1.3
    self.assertAlmostEqual(ds.ComputeL1Norm(v1, v2), 0)
    for i in range(30):
      if (i % 3 == 0):
        v2[i] = 1.3
      else:
        v2[i] = 0
    self.assertAlmostEqual(ds.ComputeL1Norm(v1, v2), 19.5)

  def test3Pickles(self):
    #outF = file('../testData/rvvs.pkl','wb+')
    with open(os.path.join(RDConfig.RDBaseDir, 'Code/DataStructs/Wrap/testData/rvvs.pkl'),
              'rb') as inF:
      v1 = ds.RealValueVect(30)
      for i in range(15):
        v1[2 * i] = 1.3
      v2 = pickle.loads(pickle.dumps(v1))
      self.assertAlmostEqual(ds.ComputeL1Norm(v1, v2), 0)
      #pickle.dump(v1,outF)
      v2 = pickle.load(inF, encoding='bytes')
      self.assertAlmostEqual(ds.ComputeL1Norm(v1, v2), 0)
      self.assertAlmostEqual(v1.GetTotalVal(), v2.GetTotalVal())
      self.assertTrue(v2.GetTotalVal() != 0)
    #outF.close()

  def test4RealVectOps(self):
    v1 = ds.RealValueVect(8)
    for i in range(4):
      v1[2 * i] = 2.1
    self.assertAlmostEqual(v1.GetTotalVal(), 8.4)
    v2 = ds.RealValueVect(8)
    for i in range(4):
      v2[2 * i + 1] = 2.1
      v2[2 * i] = 1.1
    self.assertAlmostEqual(v2.GetTotalVal(), 12.8)

    v3 = v1 | v2
    self.assertTrue(len(v3) == len(v2))
    self.assertAlmostEqual(v3.GetTotalVal(), 16.8)

    v3 = v1 & v2
    self.assertTrue(len(v3) == len(v2))
    self.assertAlmostEqual(v3.GetTotalVal(), 4.4)

    v4 = v1 + v2
    self.assertTrue(len(v4) == len(v2))
    self.assertAlmostEqual(v4.GetTotalVal(), 21.2)

    v4 = v1 - v2
    self.assertAlmostEqual(v4.GetTotalVal(), -4.4)
    v4 = v2 - v1
    self.assertAlmostEqual(v4.GetTotalVal(), 4.4)

    v4 = v2
    v4 -= v1
    self.assertAlmostEqual(v4.GetTotalVal(), 4.4)
    v4 -= v4
    self.assertAlmostEqual(v4.GetTotalVal(), 0)

  def testIterator(self):
    v1 = ds.RealValueVect(30)
    for i in range(15):
      v1[2 * i] = 1.1
    l1 = list(v1)
    self.assertTrue(len(l1) == len(v1))
    for i, v in enumerate(v1):
      self.assertAlmostEqual(l1[i], v)
    with self.assertRaises(IndexError):
      v1[40]

  def test9ToNumpy(self):
    import numpy
    bv = ds.RealValueVect(32)
    bv[0] = 1.1
    bv[1] = 4.4
    bv[17] = 1.2
    bv[23] = 8.3
    bv[31] = 12.2
    arr = numpy.zeros((32, ), numpy.double)
    ds.ConvertToNumpyArray(bv, arr)
    for i in range(len(bv)):
      self.assertAlmostEqual(bv[i], arr[i])


if __name__ == '__main__':
  unittest.main()
