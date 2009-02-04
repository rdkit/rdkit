#
#  Copyright (C) 2001,2003  greg Landrum and Rational Discovery LLC
#

""" unit testing code for variable quantization

"""
import unittest
from rdkit import RDConfig
from rdkit.ML.Data import DataUtils

class TestCase(unittest.TestCase):
  def setUp(self):
    self.d1 = [[0,0,0],[1,0,1],[2,0,0],[3,0,1],[4,0,0],[5,0,0]]

  def test1(self):
    """ basics """
    probes = [
      (.5,4,2),
      (.7,3,3),
      (.75,3,3),
      (.333,6,0),
      (.25,4,2),
      ]
    for frac,nKeep,nRej in probes:
      k,r = DataUtils.FilterData(self.d1,1,frac)
      assert len(k)==nKeep,'bad nKeep (%d != %d)'%(len(k),nKeep)
      assert len(r)==nRej,'bad nRej (%d != %d)'%(len(r),nRej)

  def test3(self):
    """ indicesOnly """
    probes = [
      (.5,4,2),
      (.7,3,3),
      (.75,3,3),
      (.333,6,0),
      (.25,4,2),
      ]
    for frac,nKeep,nRej in probes:
      DataUtils.InitRandomNumbers((23,42))
      k,r = DataUtils.FilterData(self.d1,1,frac,indicesOnly=1)
      assert len(k)==nKeep,'bad nKeep (%d != %d)'%(len(k),nKeep)
      assert len(r)==nRej,'bad nRej (%d != %d)'%(len(r),nRej)
      # make sure the indices are actually correct
      keep = [self.d1[x] for x in k]
      rej = [self.d1[x] for x in r]
      DataUtils.InitRandomNumbers((23,42))
      tgtKeep,tgtRej = DataUtils.FilterData(self.d1,1,frac)
      assert keep==tgtKeep,'%.2f: %s!=%s'%(frac,str(keep),str(tgtKeep))
      assert rej==tgtRej,'%.2f: %s!=%s'%(frac,str(rej),str(tgtRej))


  def test4(self):
    """ indicesOnly with indicesToUse """
    probes = [
      (.5,4,2),
      (.7,3,3),
      (.75,3,3),
      (.333,6,0),
      (.25,4,2),
      ]
    nPts = len(self.d1)
    for frac,nKeep,nRej in probes:
      DataUtils.InitRandomNumbers((23,42))
      k,r = DataUtils.FilterData(self.d1,1,frac,indicesToUse=range(nPts),
                                 indicesOnly=1)
      assert len(k)==nKeep,'bad nKeep (%d != %d)'%(len(k),nKeep)
      assert len(r)==nRej,'bad nRej (%d != %d)'%(len(r),nRej)
      # make sure the indices are actually correct
      keep = [self.d1[x] for x in k]
      rej = [self.d1[x] for x in r]
      DataUtils.InitRandomNumbers((23,42))
      tgtKeep,tgtRej = DataUtils.FilterData(self.d1,1,frac)
      assert keep==tgtKeep,'%.2f: %s!=%s'%(frac,str(keep),str(tgtKeep))
      assert rej==tgtRej,'%.2f: %s!=%s'%(frac,str(rej),str(tgtRej))

  def test5(self):
    """ indicesToUse """
    probes = [
      (.5,4,2),
      (.7,3,3),
      (.75,3,3),
      (.333,6,0),
      (.25,4,2),
      ]
    nPts = len(self.d1)
    for frac,nKeep,nRej in probes:
      DataUtils.InitRandomNumbers((23,42))
      k,r = DataUtils.FilterData(self.d1,1,frac,indicesToUse=range(nPts))
      assert len(k)==nKeep,'bad nKeep (%d != %d)'%(len(k),nKeep)
      assert len(r)==nRej,'bad nRej (%d != %d)'%(len(r),nRej)
      keep,rej = k,r
      # make sure the examples are actually correct
      DataUtils.InitRandomNumbers((23,42))
      tgtKeep,tgtRej = DataUtils.FilterData(self.d1,1,frac)
      assert keep==tgtKeep,'%.2f: %s!=%s'%(frac,str(keep),str(tgtKeep))
      assert rej==tgtRej,'%.2f: %s!=%s'%(frac,str(rej),str(tgtRej))


if __name__ == '__main__':
  unittest.main()
