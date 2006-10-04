# $Id$
#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the ScreenComposite functionality

"""
import RDConfig
import unittest,os
from ML import BuildComposite
from ML import ScreenComposite
import cPickle as pickle

def feq(a,b,tol=1e-4):
  if abs(a-b)>tol: return 0
  else: return 1
  
class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    self.baseDir = os.path.join(RDConfig.RDCodeDir,'ML','test_data')
    if not RDConfig.usePgSQL:
      self.dbName = os.path.join(self.baseDir,'test.gdb')
    else:
      self.dbName = "::RDTests"
    self.details = ScreenComposite.SetDefaults()
    self.details.dbName = self.dbName
    self.details.dbUser = RDConfig.defaultDBUser
    self.details.dbPassword = RDConfig.defaultDBPassword

  def test1(self):
    """ basics """
    self.details.tableName = 'ferro_quant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 6
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))

    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==90)
    self.failUnless(misCount==5)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.9878))
    self.failUnless(feq(avgBad,.7400))
    self.failUnless(tbl[0,0] == 56)
    self.failUnless(tbl[1,1] == 34)
    self.failUnless(tbl[0,1] == 0)
    self.failUnless(tbl[1,0] == 5)
    
  def test2(self):
    """ include holdout data only """
    self.details.tableName = 'ferro_quant'
    self.details.doHoldout=1
    self.details.doTraining=0
    
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 6
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))

    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==25)
    self.failUnless(misCount==4)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.9880))
    self.failUnless(feq(avgBad,.7500))
    self.failUnless(tbl[0,0] == 16)
    self.failUnless(tbl[1,1] == 9)
    self.failUnless(tbl[0,1] == 0)
    self.failUnless(tbl[1,0] == 4)
    
  def test3(self):
    """ include training data only """
    self.details.tableName = 'ferro_quant'
    self.details.doHoldout=0
    self.details.doTraining=1

    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 6
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))

    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==65)
    self.failUnless(misCount==1)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.9877))
    self.failUnless(feq(avgBad,.7000))
    self.failUnless(tbl[0,0] == 40)
    self.failUnless(tbl[1,1] == 25)
    self.failUnless(tbl[0,1] == 0)
    self.failUnless(tbl[1,0] == 1)
    
  def test4(self):
    """ include thresholding """
    self.details.tableName = 'ferro_quant'
    self.details.threshold = 0.80
    self.details.doHoldout=0
    self.details.doTraining=0

    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_quant_10.pkl'),
                              'rb'))
    tgt = 6
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))

    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==88,str(nGood))
    self.failUnless(misCount==0)
    self.failUnless(nSkipped==7)
    self.failUnless(feq(avgGood,.9932))
    self.failUnless(feq(avgBad,.000))
    self.failUnless(feq(avgSkip,.7428))
    self.failUnless(tbl[0,0] == 54)
    self.failUnless(tbl[1,1] == 34)
    self.failUnless(tbl[0,1] == 0)
    self.failUnless(tbl[1,0] == 0)
    
  def test5(self):
    """ basics """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_auto_10_3.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))

    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==94)
    self.failUnless(misCount==9)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.9298))
    self.failUnless(feq(avgBad,.6333))
    self.failUnless(tbl[0,0] == 55)
    self.failUnless(tbl[1,1] == 39)
    self.failUnless(tbl[0,1] == 0)
    self.failUnless(tbl[1,0] == 9)
    
  def test6(self):
    """ multiple models """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_auto_10_3.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))
    composites = [compos,compos]
    tpl=ScreenComposite.ScreenFromDetails(composites,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(feq(nGood[0],94))
    self.failUnless(feq(misCount[0],9))
    self.failUnless(feq(nSkipped[0],0))
    self.failUnless(feq(avgGood[0],.9298))
    self.failUnless(feq(avgBad[0],.6333))
    self.failUnless(feq(nGood[1],0))
    self.failUnless(feq(misCount[1],0))
    self.failUnless(feq(nSkipped[1],0))
    self.failUnless(feq(avgGood[1],0))
    self.failUnless(feq(avgBad[1],0))
    self.failUnless(feq(tbl[0,0],55))
    self.failUnless(feq(tbl[1,1],39))
    self.failUnless(feq(tbl[0,1],0))
    self.failUnless(feq(tbl[1,0],9))
    
  def test7(self):
    """ shuffle """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))
    self.details.shuffleActivities=1
    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==46)
    self.failUnless(misCount==57)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.7565))
    self.failUnless(feq(avgBad,.7491))
    self.failUnless(tbl[0,0] == 28)
    self.failUnless(tbl[1,1] == 18)
    self.failUnless(tbl[0,1] == 27)
    self.failUnless(tbl[1,0] == 30)
    
  def test8(self):
    """ shuffle with segmentation """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))
    self.details.shuffleActivities=1
    self.details.doHoldout=1
    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==13)
    self.failUnless(misCount==18)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.7923))
    self.failUnless(feq(avgBad,.6889))
    self.failUnless(tbl[0,0] == 7)
    self.failUnless(tbl[1,1] == 6)
    self.failUnless(tbl[0,1] == 8)
    self.failUnless(tbl[1,0] == 10)
    
  def test9(self):
    """ shuffle with segmentation2 """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_shuffle_10_3.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))
    self.details.shuffleActivities=1
    self.details.doTraining=1
    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood==33)
    self.failUnless(misCount==39)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.7424))
    self.failUnless(feq(avgBad,.7769))
    self.failUnless(tbl[0,0] == 21)
    self.failUnless(tbl[1,1] == 12)
    self.failUnless(tbl[0,1] == 19)
    self.failUnless(tbl[1,0] == 20)
    
  def test10(self):
    """ filtering """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_filt_10_3.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))
    self.details.filterVal=1
    self.details.filterFrac=.33

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    self.failUnless(nGood==90)
    self.failUnless(misCount==13)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.9578))
    self.failUnless(feq(avgBad,.8538))
    self.failUnless(tbl[0,0] == 54)
    self.failUnless(tbl[1,1] == 36)
    self.failUnless(tbl[0,1] == 1)
    self.failUnless(tbl[1,0] == 12)
    
  def _test11(self):
    """ filtering with segmentation """
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_filt_10_3.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))
    self.details.doHoldout=1
    self.details.filterVal=1
    self.details.filterFrac=.33

    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = ScreenComposite.ScreenFromDetails(compos,self.details)
    print nGood,misCount,nSkipped,avgGood,avgBad,avgSkip
    self.failUnless(nGood==46,nGood)
    self.failUnless(misCount==6)
    self.failUnless(nSkipped==0)
    self.failUnless(feq(avgGood,.9500))
    self.failUnless(feq(avgBad,.8333))
    self.failUnless(tbl[0,0] == 37)
    self.failUnless(tbl[1,1] == 9)
    self.failUnless(tbl[0,1] == 1)
    self.failUnless(tbl[1,0] == 5)

  def test12(self):
    """ test the naive bayes composite"""
    self.details.tableName = 'ferro_noquant'
    compos = pickle.load(open(os.path.join(self.baseDir,'ferromag_NaiveBayes.pkl'),
                              'rb'))
    tgt = 10
    self.failUnless(len(compos)==tgt,'bad composite loaded: %d != %d'%(len(compos),tgt))
    self.details.doHoldout=1
    tpl=ScreenComposite.ScreenFromDetails(compos,self.details)
    nGood,misCount,nSkipped,avgGood,avgBad,avgSkip,tbl = tpl
    self.failUnless(nGood == 25)
    self.failUnless(misCount == 6)
    self.failUnless(nSkipped == 0)
    self.failUnless(feq(avgGood, 0.9800))
    self.failUnless(feq(avgBad, 0.9500))
    self.failUnless(tbl[0,0] == 13)
    self.failUnless(tbl[0,1] == 6)
    self.failUnless(tbl[1,0] == 0)
    self.failUnless(tbl[1,1] == 12)
    

if __name__ == '__main__':
  unittest.main()

