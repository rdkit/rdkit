# $Id: UnitTestBitEnsemble.py 5009 2006-02-22 15:19:49Z glandrum $
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" unit testing code for BitEnsembles
"""
import RDConfig
import os
import unittest
from DataStructs.BitEnsemble import BitEnsemble
from DataStructs import BitEnsembleDb
from DataStructs import SparseBitVect


class TestCase(unittest.TestCase):
  def test1(self):
    ensemble = BitEnsemble()
    ensemble.SetBits([1,11,21,31])
    bv = SparseBitVect(100)
    bv.SetBit(1)
    bv.SetBit(11)
    bv.SetBit(13)

    score = ensemble.ScoreWithOnBits(bv) 
    assert score==2,'bad score: %d'%(score)
    score = ensemble.ScoreWithIndex(bv) 
    assert score==2,'bad score: %d'%(score)
    
  def test2(self):
    ensemble = BitEnsemble([1,11,21,31])
    bv = SparseBitVect(100)
    bv.SetBit(1)
    bv.SetBit(11)
    bv.SetBit(13)

    score = ensemble.ScoreWithOnBits(bv) 
    assert score==2,'bad score: %d'%(score)
    score = ensemble.ScoreWithIndex(bv) 
    assert score==2,'bad score: %d'%(score)
    
  def test3(self):
    ensemble = BitEnsemble()
    for bit in [1,11,21,31]:
      ensemble.AddBit(bit)
    bv = SparseBitVect(100)
    bv.SetBit(1)
    bv.SetBit(11)
    bv.SetBit(13)

    score = ensemble.ScoreWithOnBits(bv) 
    assert score==2,'bad score: %d'%(score)
    score = ensemble.ScoreWithIndex(bv) 
    assert score==2,'bad score: %d'%(score)

  def _setupDb(self):
    from Dbase.DbConnection import DbConnect
    if not RDConfig.usePgSQL:
      fName = os.path.join(RDConfig.RDCodeDir,'Dbase','testData','TEST.GDB')
    else:
      fName = os.path.join('::RDTests')
    self.conn  = DbConnect(fName)
    self.dbTblName = 'bit_ensemble_test'
    return self.conn
    
  def testdb1(self):
    """ test the sig - db functionality """
    conn = self._setupDb()
    ensemble = BitEnsemble()
    for bit in [1,3,4]:
      ensemble.AddBit(bit)

    sigBs = [([0,0,0,0,0,0],(0,0,0)),
             ([0,1,0,1,0,0],(1,1,0)),
             ([0,1,0,0,1,0],(1,0,1)),
             ([0,1,0,0,1,1],(1,0,1)),
             ]
    ensemble.InitScoreTable(conn,self.dbTblName)
    for bs,tgt in sigBs:
      ensemble.ScoreToDb(bs,conn)
    conn.Commit()

    d = conn.GetData(table=self.dbTblName)
    assert len(d) == len(sigBs),'bad number of results returned'
    for i in range(len(sigBs)):
      bs,tgt = tuple(sigBs[i])
      dbRes = tuple(d[i])
      assert dbRes==tgt,'bad bits returned: %s != %s'%(str(dbRes),str(tgt))
    d = None
    self.conn = None
  def testdb2(self):
    """ test the sig - db functionality """
    conn = self._setupDb()
    ensemble = BitEnsemble()
    for bit in [1,3,4]:
      ensemble.AddBit(bit)

    sigBs = [([0,0,0,0,0,0],(0,0,0)),
             ([0,1,0,1,0,0],(1,1,0)),
             ([0,1,0,0,1,0],(1,0,1)),
             ([0,1,0,0,1,1],(1,0,1)),
             ]
    ensemble.InitScoreTable(conn,self.dbTblName,idInfo='id varchar(10)',actInfo='act int')
    for bs,tgt in sigBs:
      ensemble.ScoreToDb(bs,conn,id='foo',act=1)
    conn.Commit()

    d = conn.GetData(table=self.dbTblName)
    assert len(d) == len(sigBs),'bad number of results returned'
    for i in range(len(sigBs)):
      bs,tgt = tuple(sigBs[i])
      dbRes = tuple(d[i])
      assert dbRes[1:-1]==tgt,'bad bits returned: %s != %s'%(str(dbRes[1:-1]),str(tgt))
      
    d = None
    self.conn = None
    
    
if __name__ == '__main__':
  unittest.main()
