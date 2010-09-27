# $Id$
#
#  Copyright (C) 2003-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the ScreenComposite functionality

"""
from rdkit import RDConfig
import unittest,os
from rdkit.ML import BuildComposite
from rdkit.ML import ScreenComposite
import cPickle as pickle

def feq(a,b,tol=1e-4):
  if abs(a-b)>tol: return 0
  else: return 1
  
class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    self.baseDir = os.path.join(RDConfig.RDCodeDir,'ML','test_data')
    self.dbName = RDConfig.RDTestDatabase
    self.details = ScreenComposite.SetDefaults()
    self.details.dbName = self.dbName
    self.details.dbUser = RDConfig.defaultDBUser
    self.details.dbPassword = RDConfig.defaultDBPassword

  def test1(self):
    """ basics """
    self.details.tableName = 'ferro_quant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 7
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==93
    assert misCount==2
    assert nSkipped==0
    assert feq(avgGood,.9849),avgGood
    assert feq(avgBad,.8500),avgBad
    assert tbl[0,0] == 54,tbl
    assert tbl[1,1] == 39
    assert tbl[0,1] == 2
    assert tbl[1,0] == 0
    
  def test2(self):
    """ include holdout data only """
    self.details.tableName = 'ferro_quant'
    self.details.doHoldout=1
    self.details.doTraining=0
    
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 7
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==28
    assert misCount==1
    assert nSkipped==0
    assert feq(avgGood,.9857),avgGood
    assert feq(avgBad,1.000),avgBad
    assert tbl[0,0] == 16,tbl
    assert tbl[1,1] == 12
    assert tbl[0,1] == 1
    assert tbl[1,0] == 0
    
  def test3(self):
    """ include training data only """
    self.details.tableName = 'ferro_quant'
    self.details.doHoldout=0
    self.details.doTraining=1

    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 7
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==65
    assert misCount==1
    assert nSkipped==0
    assert feq(avgGood,.9846),avgGood
    assert feq(avgBad,.7000),avgBad
    assert tbl[0,0] == 38,tbl
    assert tbl[1,1] == 27
    assert tbl[0,1] == 1
    assert tbl[1,0] == 0
    
  def test4(self):
    """ include thresholding """
    self.details.tableName = 'ferro_quant'
    self.details.threshold = 0.80
    self.details.doHoldout=0
    self.details.doTraining=0

    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 7
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==87,str(nGood)
    assert misCount==1
    assert nSkipped==7,nSkipped
    assert feq(avgGood,1.0),avgGood
    assert feq(avgBad,1.000),avgBad
    assert feq(avgSkip,.7571),avgSkip
    assert tbl[0,0] == 50
    assert tbl[1,1] == 37
    assert tbl[0,1] == 1
    assert tbl[1,0] == 0
    
  def test5(self):
    """ basics """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_auto_10_3.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)

    tpl = ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl

    assert nGood==93,nGood
    assert misCount==10
    assert nSkipped==0
    assert feq(avgGood,.9699),avgGood
    assert feq(avgBad,.8100),avgBad
    assert tbl[0,0] == 48,tbl
    assert tbl[1,1] == 45
    assert tbl[0,1] == 7
    assert tbl[1,0] == 3
    
  def test6(self):
    """ multiple models """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_auto_10_3.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)
    composites = [compos,compos]
    tpl = ScreenComposite.ScreenFromDetails(composites,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    assert feq(nGood[0],93),nGood
    assert feq(misCount[0],10)
    assert feq(nSkipped[0],0)
    assert feq(avgGood[0],.9699),avgGood
    assert feq(avgBad[0],.8100),avgBad
    assert feq(nGood[1],0)
    assert feq(misCount[1],0)
    assert feq(nSkipped[1],0)
    assert feq(avgGood[1],0)
    assert feq(avgBad[1],0)
    assert feq(tbl[0,0],48),tbl
    assert feq(tbl[1,1],45)
    assert feq(tbl[0,1],7)
    assert feq(tbl[1,0],3)
    
  def test7(self):
    """ shuffle """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)
    self.details.shuffleActivities=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==50,nGood
    assert misCount==53
    assert nSkipped==0
    assert feq(avgGood,.7380),avgGood
    assert feq(avgBad,.7660),avgBad
    assert tbl[0,0] == 30,tbl
    assert tbl[1,1] == 20
    assert tbl[0,1] == 25
    assert tbl[1,0] == 28
    
  def test8(self):
    """ shuffle with segmentation """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)
    self.details.shuffleActivities=1
    self.details.doHoldout=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==19,nGood
    assert misCount==12
    assert nSkipped==0
    assert feq(avgGood,.7737),avgGood
    assert feq(avgBad,.7500),avgBad
    assert tbl[0,0] == 12,tbl
    assert tbl[1,1] == 7
    assert tbl[0,1] == 6
    assert tbl[1,0] == 6
    
  def test9(self):
    """ shuffle with segmentation2 """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)
    self.details.shuffleActivities=1
    self.details.doTraining=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==31,nGood
    assert misCount==41
    assert nSkipped==0
    assert feq(avgGood,.7161),avgGood
    assert feq(avgBad,.7707),avgBad
    assert tbl[0,0] == 18,tbl
    assert tbl[1,1] == 13
    assert tbl[0,1] == 19
    assert tbl[1,0] == 22
    
  def test10(self):
    """ filtering """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_filt_10_3.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)
    self.details.filterVal=1
    self.details.filterFrac=.33

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood==90
    assert misCount==13
    assert nSkipped==0
    assert feq(avgGood,.9578)
    assert feq(avgBad,.8538)
    assert tbl[0,0] == 54
    assert tbl[1,1] == 36
    assert tbl[0,1] == 1
    assert tbl[1,0] == 12
    
  def test11(self):
    """ filtering with segmentation """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_filt_10_3.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)
    self.details.doHoldout=1
    self.details.filterVal=1
    self.details.filterFrac=.33

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)

    assert nGood==37,nGood
    assert misCount==6
    assert nSkipped==0
    assert feq(avgGood,.9594)
    assert feq(avgBad,.85)
    assert tbl[0,0] == 14,tbl
    assert tbl[1,1] == 23
    assert tbl[0,1] == 1
    assert tbl[1,0] == 5

  def test12(self):
    """ test the naive bayes composite"""
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_NaiveBayes.pkl'),
                              'rb'))
    tgt = 10
    assert len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt)
    self.details.doHoldout=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    assert nGood == 27,nGood
    assert misCount == 4,misCount
    assert nSkipped == 0,nSkipped
    assert feq(avgGood, 0.9407),avgGood
    assert feq(avgBad, 0.875),avgBad
    assert tbl[0,0] == 11,tbl
    assert tbl[0,1] == 4
    assert tbl[1,0] == 0
    assert tbl[1,1] == 16
    

if __name__ == '__main__':
  unittest.main()

