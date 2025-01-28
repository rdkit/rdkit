import random
import unittest

import numpy

from rdkit import DataStructs, RDConfig
from rdkit.DataManip.Metric import rdMetricMatrixCalc as rdmmc


def feq(v1, v2, tol2=1e-4):
  return abs(v1 - v2) <= tol2


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test0DistsArray(self):
    exp = numpy.array([1., 1.414213, 1.0], 'd')

    # initialize a double array and check if get back the expected distances
    desc = numpy.zeros((3, 2), 'd')
    desc[1, 0] = 1.0
    desc[2, 0] = 1.0
    desc[2, 1] = 1.0

    dmat = rdmmc.GetEuclideanDistMat(desc)
    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

    # repeat with an float array
    desc = numpy.zeros((3, 2), 'f')
    desc[1, 0] = 1.0
    desc[2, 0] = 1.0
    desc[2, 1] = 1.0

    dmat = rdmmc.GetEuclideanDistMat(desc)
    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

    # finally with an integer array
    desc = numpy.zeros((3, 2), 'i')
    desc[1, 0] = 1
    desc[2, 0] = 1
    desc[2, 1] = 1

    dmat = rdmmc.GetEuclideanDistMat(desc)
    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

  def ctest1DistsListArray(self):
    exp = numpy.array([1., 1.414213, 1.0], 'd')

    desc = [
      numpy.array([0.0, 0.0], 'd'),
      numpy.array([1.0, 0.0], 'd'),
      numpy.array([1.0, 1.0], 'd')
    ]
    dmat = rdmmc.GetEuclideanDistMat(desc)

    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

    # repeat the test with a list of numpy.arrays of floats
    desc = [
      numpy.array([0.0, 0.0], 'f'),
      numpy.array([1.0, 0.0], 'f'),
      numpy.array([1.0, 1.0], 'f')
    ]
    dmat = rdmmc.GetEuclideanDistMat(desc)
    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

    # repeat the test with a list of numpy.arrays of ints
    desc = [numpy.array([0, 0], 'i'), numpy.array([1, 0], 'i'), numpy.array([1, 1], 'i')]
    dmat = rdmmc.GetEuclideanDistMat(desc)
    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

  def test2DistListList(self):
    exp = numpy.array([1., 1.414213, 1.0], 'd')

    desc = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]
    dmat = rdmmc.GetEuclideanDistMat(desc)
    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

    # test with ints
    desc = [[0, 0], [1, 0], [1, 1]]
    dmat = rdmmc.GetEuclideanDistMat(desc)
    for i in range(numpy.shape(dmat)[0]):
      assert feq(dmat[i], exp[i])

  def test3Compare(self):
    n = 30
    m = 5

    dscArr = numpy.random.default_rng().random((n,m))
    dmatArr = rdmmc.GetEuclideanDistMat(dscArr)

    dscLL = []
    for i in range(n):
      row = []
      for j in range(m):
        row.append(dscArr[i, j])
      dscLL.append(row)
    dmatLL = rdmmc.GetEuclideanDistMat(dscLL)

    assert numpy.shape(dmatArr) == numpy.shape(dmatLL)

    for i in range(n * (n - 1) // 2):
      assert feq(dmatArr[i], dmatLL[i])

  def test4ebv(self):

    n = 30
    m = 2048
    dm = 800
    lst = []
    for i in range(n):
      v = DataStructs.ExplicitBitVect(m)
      for j in range(dm):
        v.SetBit(random.randrange(0, m))
      lst.append(v)

    dMat = rdmmc.GetTanimotoDistMat(lst)

    sMat = rdmmc.GetTanimotoSimMat(lst)

    for i in range(n * (n - 1) // 2):
      assert feq(sMat[i] + dMat[i], 1.0)

  def test5sbv(self):

    n = 30
    m = 2048
    dm = 800
    lst = []
    for i in range(n):
      v = DataStructs.SparseBitVect(m)
      for j in range(dm):
        v.SetBit(random.randrange(0, m))
      lst.append(v)

    dMat = rdmmc.GetTanimotoDistMat(lst)

    sMat = rdmmc.GetTanimotoSimMat(lst)

    for i in range(n * (n - 1) // 2):
      assert feq(sMat[i] + dMat[i], 1.0)


if __name__ == '__main__':
  unittest.main()
