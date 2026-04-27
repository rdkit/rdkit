# $Id$
#
# Copyright (C) 2007,2008 Greg Landrum
#
#  @@ All Rights Reserved @@
#
import io
import os
import pickle
import random
import sys
import unittest

from rdkit import DataStructs as ds
from rdkit import RDConfig


def feq(v1, v2, tol=1e-4):
  return abs(v1 - v2) < tol


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1Int(self):
    """

    """
    v1 = ds.IntSparseIntVect(5)
    self.assertRaises(IndexError, lambda: v1[5])
    v1[0] = 1
    v1[2] = 2
    v1[3] = 3
    self.assertTrue(v1 == v1)
    self.assertTrue(v1.GetLength() == 5)

    v2 = ds.IntSparseIntVect(5)
    self.assertTrue(v1 != v2)
    v2 |= v1
    self.assertTrue(v2 == v1)

    v3 = v2 | v1
    self.assertTrue(v3 == v1)

    onVs = v1.GetNonzeroElements()
    self.assertTrue(onVs == {0: 1, 2: 2, 3: 3})

  def test2Long(self):
    """

    """
    l = 1 << 42
    v1 = ds.LongSparseIntVect(l)
    self.assertRaises(IndexError, lambda: v1[l])
    v1[0] = 1
    v1[2] = 2
    v1[1 << 35] = 3
    self.assertTrue(v1 == v1)
    self.assertTrue(v1.GetLength() == l)

    v2 = ds.LongSparseIntVect(l)
    self.assertTrue(v1 != v2)
    v2 |= v1
    self.assertTrue(v2 == v1)

    v3 = v2 | v1
    self.assertTrue(v3 == v1)

    onVs = v1.GetNonzeroElements()
    self.assertTrue(onVs == {0: 1, 2: 2, 1 << 35: 3})

  def test3Pickle1(self):
    """

    """
    l = 1 << 42
    v1 = ds.LongSparseIntVect(l)
    self.assertRaises(IndexError, lambda: v1[l + 1])
    v1[0] = 1
    v1[2] = 2
    v1[1 << 35] = 3
    self.assertTrue(v1 == v1)

    v2 = pickle.loads(pickle.dumps(v1))
    self.assertTrue(v2 == v1)

    v3 = ds.LongSparseIntVect(v2.ToBinary())
    self.assertTrue(v2 == v3)
    self.assertTrue(v1 == v3)

    #pickle.dump(v1,file('lsiv.pkl','wb+'))
    with open(os.path.join(RDConfig.RDBaseDir, 'Code/DataStructs/Wrap/testData/lsiv.pkl'),
              'r') as tf:
      buf = tf.read().replace('\r\n', '\n').encode('utf-8')
      tf.close()
    with io.BytesIO(buf) as f:
      v3 = pickle.load(f)
      self.assertTrue(v3 == v1)

  def test3Pickle2(self):
    """

    """
    l = 1 << 21
    v1 = ds.IntSparseIntVect(l)
    self.assertRaises(IndexError, lambda: v1[l + 1])
    v1[0] = 1
    v1[2] = 2
    v1[1 << 12] = 3
    self.assertTrue(v1 == v1)

    v2 = pickle.loads(pickle.dumps(v1))
    self.assertTrue(v2 == v1)

    v3 = ds.IntSparseIntVect(v2.ToBinary())
    self.assertTrue(v2 == v3)
    self.assertTrue(v1 == v3)

    #pickle.dump(v1,file('isiv.pkl','wb+'))
    with open(os.path.join(RDConfig.RDBaseDir, 'Code/DataStructs/Wrap/testData/isiv.pkl'),
              'r') as tf:
      buf = tf.read().replace('\r\n', '\n').encode('utf-8')
      tf.close()
    with io.BytesIO(buf) as f:
      v3 = pickle.load(f)
      self.assertTrue(v3 == v1)

  def test4Update(self):
    """

    """
    v1 = ds.IntSparseIntVect(5)
    self.assertRaises(IndexError, lambda: v1[6])
    v1[0] = 1
    v1[2] = 2
    v1[3] = 3
    self.assertTrue(v1 == v1)

    v2 = ds.IntSparseIntVect(5)
    v2.UpdateFromSequence((0, 2, 3, 3, 2, 3))
    self.assertTrue(v1 == v2)

  def test5Dice(self):
    """

    """
    v1 = ds.IntSparseIntVect(5)
    v1[4] = 4
    v1[0] = 2
    v1[3] = 1
    self.assertTrue(feq(ds.DiceSimilarity(v1, v1), 1.0))

    v1 = ds.IntSparseIntVect(5)
    v1[0] = 2
    v1[2] = 1
    v1[3] = 4
    v1[4] = 6
    v2 = ds.IntSparseIntVect(5)
    v2[1] = 2
    v2[2] = 3
    v2[3] = 4
    v2[4] = 4
    self.assertTrue(feq(ds.DiceSimilarity(v1, v2), 18.0 / 26.))
    self.assertTrue(feq(ds.DiceSimilarity(v2, v1), 18.0 / 26.))

  def test6BulkDice(self):
    """

    """
    sz = 10
    nToSet = 5
    nVs = 6
    import random
    vs = []
    for _ in range(nVs):
      v = ds.IntSparseIntVect(sz)
      for _ in range(nToSet):
        v[random.randint(0, sz - 1)] = random.randint(1, 10)
      vs.append(v)

    baseDs = [ds.DiceSimilarity(vs[0], vs[x]) for x in range(1, nVs)]
    bulkDs = ds.BulkDiceSimilarity(vs[0], vs[1:])
    for bbaseDs, bbulkDs in zip(baseDs, bulkDs):
      self.assertTrue(feq(bbaseDs, bbulkDs))

  def test6BulkTversky(self):
    """

    """
    sz = 10
    nToSet = 5
    nVs = 6
    import random
    vs = []
    for _ in range(nVs):
      v = ds.IntSparseIntVect(sz)
      for _ in range(nToSet):
        v[random.randint(0, sz - 1)] = random.randint(1, 10)
      vs.append(v)

    baseDs = [ds.TverskySimilarity(vs[0], vs[x], .5, .5) for x in range(1, nVs)]
    bulkDs = ds.BulkTverskySimilarity(vs[0], vs[1:], 0.5, 0.5)
    diceDs = [ds.DiceSimilarity(vs[0], vs[x]) for x in range(1, nVs)]
    for bbaseDs, bbulkDs, ddiceDs in zip(baseDs, bulkDs, diceDs):
      self.assertTrue(feq(bbaseDs, bbulkDs))
      self.assertTrue(feq(bbaseDs, ddiceDs))
      self.assertTrue(feq(bbulkDs, ddiceDs))

    bulkDs = ds.BulkTverskySimilarity(vs[0], vs[1:], 1.0, 1.0)
    taniDs = [ds.TanimotoSimilarity(vs[0], vs[x]) for x in range(1, nVs)]
    for bbulkDs, ttaniDs in zip(bulkDs, taniDs):
      self.assertTrue(feq(bbulkDs, ttaniDs))
    taniDs = ds.BulkTanimotoSimilarity(vs[0], vs[1:])
    for bbulkDs, ttaniDs in zip(bulkDs, taniDs):
      self.assertTrue(feq(bbulkDs, ttaniDs))

  def test7ToList(self):
    l = [0] * 2048
    nbits = 2048
    bv = ds.IntSparseIntVect(nbits)
    for _ in range(nbits):
      x = random.randrange(0, nbits)
      l[x] = x
      bv[x] = x

    l2 = list(bv)
    l3 = bv.ToList()
    self.assertEqual(l, l2)
    self.assertEqual(l, l3)


if __name__ == '__main__':
  unittest.main()
