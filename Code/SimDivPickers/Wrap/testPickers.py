from rdkit import RDConfig
import unittest
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.DataManip.Metric import rdMetricMatrixCalc as rdmmc
import numpy
import random

class TestCase(unittest.TestCase):
  def setUp(self) :
      self.n = 1000
      self.m = 80
      self.d = 2
      self.dataPts = numpy.zeros((self.n, self.d), 'd')
      for i in range(self.n):
        for j in range(self.d):
          self.dataPts[i,j] = random.random()
      self.dMat = rdmmc.GetEuclideanDistMat(self.dataPts)
       
  def test0MaxMin(self):
    pkr = rdSimDivPickers.MaxMinPicker()
    maxmin = pkr.Pick(self.dMat, self.n, self.m,(886,112))
    self.failUnless(maxmin[0]==886)
    self.failUnless(maxmin[1]==112)

    def func(i,j):
      if i==j:
        return 0.0
      if i<j:
        j,i=i,j
      return self.dMat[i*(i-1)/2+j]
    lmaxmin = pkr.LazyPick(func, self.n, self.m,(886,112))
    self.failUnless(list(lmaxmin)==list(maxmin))

    self.failUnlessRaises(ValueError,lambda:pkr.Pick(self.dMat, self.n, self.m,(1012,)))
    self.failUnlessRaises(ValueError,lambda:pkr.Pick(self.dMat, self.n, self.m,(-1,)))

    maxmin = pkr.Pick(self.dMat, self.n, self.m)
    self.failUnless(maxmin)
    lmaxmin = pkr.LazyPick(func, self.n, self.m)
    self.failUnless(lmaxmin)

  def test1HierarchPick(self) :
    infil = open("test_data/points.csv", 'r')
    lines = infil.readlines()
    infil.close()
    self.dataPts = numpy.zeros((len(lines), 2), 'd')
    labels = []
    i = 0
    for line in lines :
      tlst = line.strip().split(',')
      self.dataPts[i, 0] = float(tlst[1])
      self.dataPts[i, 1] = float(tlst[2])
      labels.append(int(tlst[3]))
      i += 1
    self.dMat = rdmmc.GetEuclideanDistMat(self.dataPts)
    pkr = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    clusters = pkr.Cluster(self.dMat, i, 2)
    # check that each of the clusters have the same label
    for cl in clusters :
      clbl = labels[cl[0]]
      for id in cl:
        assert clbl == labels[id]
    hierarch = pkr.Pick(self.dMat, i, 2)
    assert tuple(hierarch) == (1,30)


  def testIssue208(self) :
    sz = 10
    N=3
    m = []
    for i in range(sz):
      for j in range(i+1,sz):
        m.append(random.random())
    m = numpy.array(m)
    picker = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    p1 = list(picker.Pick(m,sz,N))
    p1.sort()
    p2 = list(picker.Pick(m,sz,N))
    p2.sort()
    self.failUnless(p1==p2)

  def testInts(self) :
    """ make sure we can handle ints too """
    sz = 10
    N=3
    m = []
    for i in range(sz):
      for j in range(i+1,sz):
        m.append(int(100*random.random()))
    m = numpy.array(m)
    picker = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    p1 = list(picker.Pick(m,sz,N))
    p1.sort()
    p2 = list(picker.Pick(m,sz,N))
    p2.sort()
    self.failUnless(p1==p2)

            
  def testNonUniqueCrash(self) :
    from rdkit import DataStructs
    sz = 10
    nbits=20
    nBitsToSet=int(nbits*.3)
    N=12
    vs = []
    for i in range(sz):
      bv = DataStructs.ExplicitBitVect(nbits)
      for j in range(nBitsToSet):
        val= int(nbits*random.random())
        bv.SetBit(val)
      vs.append(bv)
      vs.append(bv)
    def taniFunc(i,j,bvs = vs):
      d = 1-DataStructs.FingerprintSimilarity(bvs[i],bvs[j])
      return d
    picker = rdSimDivPickers.MaxMinPicker()
    try:
      mm = picker.LazyPick(taniFunc,len(vs),N)
    except:
      ok=False
    else:
      ok=True
    self.failUnless(ok)
    self.failUnless(len(mm)==N)
    picker = None

    ds = []
    nvs = len(vs)
    for i in range(nvs):
      for j in range(i+1,nvs):
        d = taniFunc(i,j)
        ds.append(d)
    m = numpy.array(ds)
    picker = rdSimDivPickers.HierarchicalClusterPicker(rdSimDivPickers.ClusterMethod.WARD)
    p1 = list(picker.Pick(m,nvs,N))

            
if __name__ == '__main__':
    unittest.main()


