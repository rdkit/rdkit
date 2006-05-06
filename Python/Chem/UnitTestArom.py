# $Id: UnitTestArom.py 5083 2006-03-11 17:56:19Z NightlyBuild $
#
#  Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the Crippen clogp and MR calculators

"""
import RDConfig
import unittest,sys,os,cPickle
import Chem
from Chem import Crippen
from Numeric import *

def feq(n1,n2,tol=1e-5):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    self.fName = os.path.join(RDConfig.RDCodeDir,'Chem/tests','aromat_regress.txt')
    self.fName2 = os.path.join(RDConfig.RDCodeDir,'Chem/tests','NCI_aromat_regress.txt')

  def _readData(self,fName):
    d = []
    lineNo=0
    for line in open(fName,'r').xreadlines():
      lineNo+=1
      if len(line) and line[0] != '#':
        splitL = line.split('\t')
        if len(splitL)==4:
          smi1,smi2,count,ats = splitL
          d.append((lineNo,smi1,smi2,int(count),eval(ats)))
    self.data = d
          
  def test1(self,maxFailures=50):
    self._readData(self.fName)
    nMols = len(self.data)
    nFailed = 0
    for i in range(nMols):
      lineNo,smi,smi2,tgtCount,tgtAts = self.data[i]
      try:
        mol = Chem.MolFromSmiles(smi)
        if not mol: raise ValueError
      except:
        mol = None
        print 'failure(%d): '%lineNo,smi
        print '-----------------------------'
      else:
        count = 0
        aroms = []
        for at in mol.GetAtoms():
          if at.GetIsAromatic():
            aroms.append(at.GetIdx())
            count+=1
        if count != tgtCount:
          print 'Fail(%d): %s, %s'%(lineNo,smi,Chem.MolToSmiles(mol))
          print '\t %d != %d'%(count,tgtCount)
          print '\t ',repr(aroms)
          print '\t ',repr(tgtAts)
          print '-----------------------------'
          nFailed += 1
          if nFailed >= maxFailures:
            assert 0

        
  def test2(self,maxFailures=50):
    self._readData(self.fName2)
    nMols = len(self.data)
    nFailed = 0
    for i in range(nMols):
      lineNo,smi,smi2,tgtCount,tgtAts = self.data[i]
      try:
        mol = Chem.MolFromSmiles(smi)
        if not mol: raise ValueError
      except:
        mol = None
        print 'failure(%d): '%lineNo,smi
        print '-----------------------------'
      else:
        count = 0
        aroms = []
        for at in mol.GetAtoms():
          if at.GetIsAromatic():
            aroms.append(at.GetIdx())
            count+=1
        if count != tgtCount:
          print 'Fail(%d): %s, %s'%(lineNo,smi,Chem.MolToSmiles(mol))
          print '\t %d != %d'%(count,tgtCount)
          print '\t ',repr(aroms)
          print '\t ',repr(tgtAts)
          print '-----------------------------'
          nFailed += 1
          if nFailed >= maxFailures:
            assert 0

        
      


          
          
      
  
    
if __name__ == '__main__':
  unittest.main()

