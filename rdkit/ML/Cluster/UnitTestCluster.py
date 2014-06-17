# $Id$
#
#  Copyright (C) 2001-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for clustering

"""
from rdkit import RDConfig
import unittest
from rdkit.ML.Cluster import Standardize,ClusterUtils
from rdkit.ML.Cluster import Clusters
try:
  from rdkit.ML.Cluster import Murtagh
except ImportError:
  Murtagh=None
import numpy
from rdkit.six.moves import cPickle


class TestCase(unittest.TestCase):
  def setUp(self):
    # this is the data set used by Romesburg in "Cluster Analysis for Researchers"
    #  to demonstrate the different clustering methods
    #print '\n%s: '%self.shortDescription(),
    self.d = numpy.array([[10.,5.],[20.,20.],[30.,10.],[30.,15.],[5.,10.]])
    self.names = ['p1','p2','p3','p4','p5']
  def testDivide(self):
    " tests the cluster division algorithms "
    ca = Clusters.Cluster(index=1)
    cb = Clusters.Cluster(index=2)
    cc = Clusters.Cluster(index=3)
    cd = Clusters.Cluster(index=4)
    ce = Clusters.Cluster(index=5)
    cf = Clusters.Cluster(index=6)

    c1 = Clusters.Cluster(metric=10,children=[ca,cb],index=7)
    c2 = Clusters.Cluster(metric=15,children=[cc,cd],index=8)
    c3 = Clusters.Cluster(metric=20,children=[ce,cf],index=9)
    c4 = Clusters.Cluster(metric=25,children=[c2,c3],index=10)
    c5 = Clusters.Cluster(metric=30,children=[c4,c1],index=11)

    cs = ClusterUtils.SplitIntoNClusters(c5,4,breadthFirst=1)
    assert len(cs)==4,'bad split length'
    indices = [x.GetIndex() for x in cs]
    for index in [9,8,1,2]:
      assert index in indices,'index %d not found in %s'%(index,str(indices))
    # we may not want to preserve order, but test it for now
    assert indices==[9,8,1,2],'bad index order'
    
    cs2 = ClusterUtils.SplitIntoNClusters(c5,4,breadthFirst=0)
    indices = [x.GetIndex() for x in cs2]
    for index in [8,7,5,6]:
      assert index in indices,'index %d not found in %s'%(index,str(indices))
    # we may not want to preserve order, but test it for now
    assert indices==[8,7,5,6],'bad index order'

  def testMurtaghUPGMA(self):
    if Murtagh is None: return
    nPts = 5
    sz = 5
    dataP = numpy.random.random((nPts,sz))
    newClust = Murtagh.ClusterData(dataP,nPts,Murtagh.UPGMA)[0]
    ds = []
    for i in range(nPts):
      for j in range(i):
        d = dataP[i]-dataP[j]
        ds.append(sum(d*d))
    ds = numpy.array(ds)    
    newClust2 = Murtagh.ClusterData(ds,nPts,Murtagh.UPGMA,isDistData=1)[0]

    assert len(newClust)==len(newClust2),'length mismatch2'

    assert not newClust.Compare(newClust2,ignoreExtras=0),'equality failed3'

    newClust2 = Murtagh.ClusterData(dataP,nPts,Murtagh.UPGMA,isDistData=0)[0]
    assert len(newClust)==len(newClust2),'length mismatch2'

    assert not newClust.Compare(newClust2,ignoreExtras=0),'equality failed3'

  def testGetNodeList(self):
    """ tests the GetNodeList functionality """
    root = Clusters.Cluster(index=1)
    c1 = Clusters.Cluster(index=10)
    c1.AddChild(Clusters.Cluster(index=30))
    c1.AddChild(Clusters.Cluster(index=31))
    c1.AddChild(Clusters.Cluster(index=32))
    
    c2 = Clusters.Cluster(index=11)
    c2.AddChild(Clusters.Cluster(index=40))
    c2.AddChild(Clusters.Cluster(index=41))

    root.AddChild(c1)
    root.AddChild(c2)
    nodes = ClusterUtils.GetNodeList(root)
    
    indices = [x.GetIndex() for x in nodes]
    assert indices==[30,31,32,10,40,41,11,1],'bad indices'

profileTest=0


if __name__ == '__main__':
  unittest.main()

