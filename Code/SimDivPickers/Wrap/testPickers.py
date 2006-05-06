import RDConfig
import unittest
from SimDivFilters import rdSimDivPickers as rdsimdiv
from DataManip.Metric import rdMetricMatrixCalc as rdmmc
from Numeric import *
import random

class TestCase(unittest.TestCase):
  def setUp(self) :
      self.n = 1000
      self.m = 80
      self.d = 2
      self.dataPts = zeros((self.n, self.d), 'd')
      for i in range(self.n):
        for j in range(self.d):
          self.dataPts[i,j] = random.random()
      self.dMat = rdmmc.GetEuclideanDistMat(self.dataPts)
       
  def test0MaxMin(self):
    pkr = rdsimdiv.MaxMinPicker()
    maxmin = pkr.Pick(self.dMat, self.n, self.m)

  def test1HierarchPick(self) :
    infil = open("test_data/points.csv", 'r')
    lines = infil.readlines()
    infil.close()
    self.dataPts = zeros((len(lines), 2), 'd')
    labels = []
    i = 0
    for line in lines :
      tlst = line.strip().split(',')
      self.dataPts[i, 0] = float(tlst[1])
      self.dataPts[i, 1] = float(tlst[2])
      labels.append(int(tlst[3]))
      i += 1
    self.dMat = rdmmc.GetEuclideanDistMat(self.dataPts)
    pkr = rdsimdiv.HierarchicalClusterPicker(rdsimdiv.ClusterMethod.WARD)
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
    m = array(m)
    picker = rdsimdiv.HierarchicalClusterPicker(rdsimdiv.ClusterMethod.WARD)
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
    m = array(m)
    picker = rdsimdiv.HierarchicalClusterPicker(rdsimdiv.ClusterMethod.WARD)
    p1 = list(picker.Pick(m,sz,N))
    p1.sort()
    p2 = list(picker.Pick(m,sz,N))
    p2.sort()
    self.failUnless(p1==p2)

            
if __name__ == '__main__':
    unittest.main()


