# $Id$
#
# Copyright (C) 2007-2008 Greg Landrum
#  All Rights Reserved
#
import unittest

from rdkit.six.moves import cPickle
from rdkit import RDConfig
from rdkit.ML.Cluster import Butina


class TestCase(unittest.TestCase):
  def test1(self):
    dists = [1,
             2,1,
             4,3,2,
             6,5,4,2,
             7,6,5,3,1
             ]
    nPts = 6
    cs = Butina.ClusterData(dists,nPts,1.1,isDistData=1)
    self.assertTrue(len(cs)==3)

    self.assertTrue(cs[0]==(1,0,2))
    self.assertTrue(cs[1]==(5,4))
    self.assertTrue(cs[2]==(3,))

  def test2(self):
    dists = [.5,
             1,.5,
             2,1.5,1,
             3,2.5,2,1,
             5,4.5,4,3,2,
             8,7.5,7,6,5,3,
             9,8.5,8,7,6,4,1,
             ]
    nPts = 8
    cs = Butina.ClusterData(dists,nPts,2.1,isDistData=1)

    self.assertTrue(len(cs)==3)

    self.assertTrue(cs[0]==(3,0,1,2,4))
    self.assertTrue(cs[1]==(7,6))
    self.assertTrue(cs[2]==(5,))

  def test3(self):
    " edge case: everything a singleton "
    dists = [1,
             2,1,
             ]
    nPts = 3
    cs = Butina.ClusterData(dists,nPts,0.9,isDistData=1)
    self.assertTrue(len(cs)==3)

    self.assertTrue(cs[0]==(2,))
    self.assertTrue(cs[1]==(1,))
    self.assertTrue(cs[2]==(0,))

  def test4(self):
    " edge case: everything in one cluster "
    dists = [1,
             2,1,
             3,2,1,
             ]
    nPts = 4
    cs = Butina.ClusterData(dists,nPts,2,isDistData=1)
    self.assertTrue(len(cs)==1)
    self.assertTrue(cs[0]==(3,0,1,2))

  def test4(self):
    " edge case: one in the middle leaves the edges lonely "
    dists = [1.5,
             2.5,1,
             3.5,2,1,
             5,3.5,2.5,1.5,
             ]
    nPts = 5
    cs = Butina.ClusterData(dists,nPts,1.1,isDistData=1)
    self.assertTrue(len(cs)==3)
    self.assertTrue(cs[0]==(2,1,3))
    self.assertTrue(cs[1]==(4,))
    self.assertTrue(cs[2]==(0,))

  def test6(self):
    " edge case: zero distances: "
    dists = [1,
             2,0,
             2,0,0,
             4,2,2,2,
             ]
    nPts = 5
    cs = Butina.ClusterData(dists,nPts,0.9,isDistData=1)
    self.assertTrue(len(cs)==3)
    self.assertTrue(cs[0]==(3,1,2))
    self.assertTrue(cs[1]==(4,))
    self.assertTrue(cs[2]==(0,))

  def test7(self):
    " reordering: no changes "
    dists = [1,
             2,1,
             4,3,2,
             6,5,4,2,
             7,6,5,3,1
             ]
    nPts = 6
    cs = Butina.ClusterData(dists,nPts,1.1,isDistData=1,reordering=True)
    self.assertTrue(len(cs)==3)

    self.assertTrue(cs[0]==(1,0,2))
    self.assertTrue(cs[1]==(5,4))
    self.assertTrue(cs[2]==(3,))

  def test8(self):
    " reordering: changes"
    dists = [2,
             3.5,1.5,
             5,3,1.5,
             7,5,3.5,2,
             8,6,4.5,3,1,
             9,7,5.5,4,2,1,
             ]
    nPts = 7
    # without reordering
    cs = Butina.ClusterData(dists,nPts,2.1,isDistData=1)
    self.assertTrue(len(cs)==3)
    self.assertTrue(cs[0]==(4,3,5,6))
    self.assertTrue(cs[1]==(2,1))
    self.assertTrue(cs[2]==(0,))

    # with reordering
    cs = Butina.ClusterData(dists,nPts,2.1,isDistData=1,reordering=True)
    self.assertTrue(len(cs)==2)
    self.assertTrue(cs[0]==(4,3,5,6))
    self.assertTrue(cs[1]==(1,0,2))
    

profileTest=0


if __name__ == '__main__':
  unittest.main()

