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
from rdkit.six.moves import cPickle as pickle

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
    with open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 5
    self.failUnlessEqual(len(compos),tgt)

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,93)
    self.failUnlessEqual(misCount,2)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.9871,4)
    self.failUnlessAlmostEqual(avgBad,.8000,4)
    self.failUnlessEqual(tbl[0,0] , 54)
    self.failUnlessEqual(tbl[1,1] , 39)
    self.failUnlessEqual(tbl[0,1] , 2)
    self.failUnlessEqual(tbl[1,0] , 0)
    
  def test2(self):
    """ include holdout data only """
    self.details.tableName = 'ferro_quant'
    self.details.doHoldout=1
    self.details.doTraining=0
    
    with open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 5
    self.failUnlessEqual(len(compos),tgt)

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,28)
    self.failUnlessEqual(misCount,1)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.9964,4)
    self.failUnlessAlmostEqual(avgBad,1.000,4)
    self.failUnlessEqual(tbl[0,0] , 16)
    self.failUnlessEqual(tbl[1,1] , 12)
    self.failUnlessEqual(tbl[0,1] , 1)
    self.failUnlessEqual(tbl[1,0] , 0)
    
  def test3(self):
    """ include training data only """
    self.details.tableName = 'ferro_quant'
    self.details.doHoldout=0
    self.details.doTraining=1

    with open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 5
    self.failUnlessEqual(len(compos),tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,65)
    self.failUnlessEqual(misCount,1)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.98307,4)
    self.failUnlessAlmostEqual(avgBad,0.600,4)
    self.failUnlessEqual(tbl[0,0] , 38,tbl)
    self.failUnlessEqual(tbl[1,1] , 27)
    self.failUnlessEqual(tbl[0,1] , 1)
    self.failUnlessEqual(tbl[1,0] , 0)
    
  def test4(self):
    """ include thresholding """
    self.details.tableName = 'ferro_quant'
    self.details.threshold = 0.80
    self.details.doHoldout=0
    self.details.doTraining=0

    with open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 5
    self.failUnlessEqual(len(compos),tgt)

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,91)
    self.failUnlessEqual(misCount,1)
    self.failUnlessEqual(nSkipped,3)
    self.failUnlessAlmostEqual(avgGood,0.9956,4)
    self.failUnlessAlmostEqual(avgBad,1.000,4)
    self.failUnlessAlmostEqual(avgSkip,0.6000,4)
    self.failUnlessEqual(tbl[0,0] , 54)
    self.failUnlessEqual(tbl[1,1] , 37)
    self.failUnlessEqual(tbl[0,1] , 1)
    self.failUnlessEqual(tbl[1,0] , 0)
    
  def test5(self):
    """ basics """
    self.details.tableName = 'ferro_noquant'

    with open(os.path.join(self.baseDir,'ferromag_auto_10_3.pkl'),'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)

    tpl = ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl

    self.failUnlessEqual(nGood,95)
    self.failUnlessEqual(misCount,8)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.9684,4)
    self.failUnlessAlmostEqual(avgBad,.8375,4)
    self.failUnlessEqual(tbl[0,0] , 50)
    self.failUnlessEqual(tbl[1,1] , 45)
    self.failUnlessEqual(tbl[0,1] , 5)
    self.failUnlessEqual(tbl[1,0] , 3)
    
  def test6(self):
    """ multiple models """
    self.details.tableName = 'ferro_noquant'
    with open(os.path.join(self.baseDir,'ferromag_auto_10_3.pkl'),'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)
    composites = [compos,compos]
    tpl = ScreenComposite.ScreenFromDetails(composites,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnlessEqual(nGood[0],95)
    self.failUnlessEqual(misCount[0],8)
    self.failUnlessEqual(nSkipped[0],0)
    self.failUnlessAlmostEqual(avgGood[0],.9684,4)
    self.failUnlessAlmostEqual(avgBad[0],.8375,4)
    self.failUnlessEqual(nGood[1],0)
    self.failUnlessEqual(misCount[1],0)
    self.failUnlessEqual(nSkipped[1],0)
    self.failUnlessEqual(avgGood[1],0)
    self.failUnlessEqual(avgBad[1],0)
    self.failUnlessEqual(tbl[0,0],50)
    self.failUnlessEqual(tbl[1,1],45)
    self.failUnlessEqual(tbl[0,1],5)
    self.failUnlessEqual(tbl[1,0],3)
    
  def test7(self):
    """ shuffle """
    self.details.tableName = 'ferro_noquant'
    with open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)
    self.details.shuffleActivities=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,50)
    self.failUnlessEqual(misCount,53)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.7380,4)
    self.failUnlessAlmostEqual(avgBad,.7660,4)
    self.failUnlessEqual(tbl[0,0] , 30)
    self.failUnlessEqual(tbl[1,1] , 20)
    self.failUnlessEqual(tbl[0,1] , 25)
    self.failUnlessEqual(tbl[1,0] , 28)
    
  def test8(self):
    """ shuffle with segmentation """
    self.details.tableName = 'ferro_noquant'
    with open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
              'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)
    self.details.shuffleActivities=1
    self.details.doHoldout=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,19)
    self.failUnlessEqual(misCount,12)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.7737,4)
    self.failUnlessAlmostEqual(avgBad,.7500,4)
    self.failUnlessEqual(tbl[0,0] , 12)
    self.failUnlessEqual(tbl[1,1] , 7)
    self.failUnlessEqual(tbl[0,1] , 6)
    self.failUnlessEqual(tbl[1,0] , 6)
    
  def test9(self):
    """ shuffle with segmentation2 """
    self.details.tableName = 'ferro_noquant'
    with open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
              'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)
    self.details.shuffleActivities=1
    self.details.doTraining=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,31)
    self.failUnlessEqual(misCount,41)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.7161,4)
    self.failUnlessAlmostEqual(avgBad,.7707,4)
    self.failUnlessEqual(tbl[0,0] , 18)
    self.failUnlessEqual(tbl[1,1] , 13)
    self.failUnlessEqual(tbl[0,1] , 19)
    self.failUnlessEqual(tbl[1,0] , 22)
    
  def test10(self):
    """ filtering """
    self.details.tableName = 'ferro_noquant'
    with open(os.path.join(self.baseDir,'ferromag_filt_10_3.pkl'),
              'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)
    self.details.filterVal=1
    self.details.filterFrac=.33

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood,90)
    self.failUnlessEqual(misCount,13)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.9578,4)
    self.failUnlessAlmostEqual(avgBad,.8538,4)
    self.failUnlessEqual(tbl[0,0] , 54)
    self.failUnlessEqual(tbl[1,1] , 36)
    self.failUnlessEqual(tbl[0,1] , 1)
    self.failUnlessEqual(tbl[1,0] , 12)
    
  def test11(self):
    """ filtering with segmentation """
    self.details.tableName = 'ferro_noquant'
    with open(os.path.join(self.baseDir,'ferromag_filt_10_3.pkl'),
              'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)
    self.details.doHoldout=1
    self.details.filterVal=1
    self.details.filterFrac=.33

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)

    self.failUnlessEqual(nGood,37)
    self.failUnlessEqual(misCount,6)
    self.failUnlessEqual(nSkipped,0)
    self.failUnlessAlmostEqual(avgGood,.95946,4)    
    self.failUnlessAlmostEqual(avgBad,.85,4)
    self.failUnlessEqual(tbl[0,0] , 14)
    self.failUnlessEqual(tbl[1,1] , 23)
    self.failUnlessEqual(tbl[0,1] , 1)
    self.failUnlessEqual(tbl[1,0] , 5)

  def test12(self):
    """ test the naive bayes composite"""
    self.details.tableName = 'ferro_noquant'
    with open(os.path.join(self.baseDir,'ferromag_NaiveBayes.pkl'),
              'rb') as pklF:
      compos = pickle.load(pklF)
    tgt = 10
    self.failUnlessEqual(len(compos),tgt)
    self.details.doHoldout=1
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnlessEqual(nGood , 25)
    self.failUnlessEqual(misCount , 6)
    self.failUnlessEqual(nSkipped , 0)
    self.failUnlessAlmostEqual(avgGood, 0.9800,4)
    self.failUnlessAlmostEqual(avgBad, 0.86667,4)
    self.failUnlessEqual(tbl[0,0] , 9)
    self.failUnlessEqual(tbl[0,1] , 6)
    self.failUnlessEqual(tbl[1,0] , 0)
    self.failUnlessEqual(tbl[1,1] , 16)
    

if __name__ == '__main__':
  unittest.main()

