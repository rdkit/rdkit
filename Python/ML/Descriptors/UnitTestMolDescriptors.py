#
#  Copyright (C) 2002  greg Landrum and Rational Discovery LLC
#

""" unit testing code for molecular descriptor calculators

"""
import unittest,cPickle,os.path
import RDConfig
from Numeric import *
from ML.Descriptors import MoleculeDescriptors
from Chem import *

class TestCase(unittest.TestCase):
  def setUp(self):
    self.descs = ['MolLogP','Chi1v']
    self.calc = MoleculeDescriptors.MolecularDescriptorCalculator(self.descs)
    self.testD = [
      ('CCOC',    (0.6497, 1.40403)),
      ('CC=O',    (0.2022, 0.81305)),
      ('CCC(=O)O',(0.475, 1.48839))]

  def testGetNames(self):
    assert self.calc.GetDescriptorNames()==tuple(self.descs),'bad descriptor names: %s'%(self.calc.GetDescriptorNames())
    
  def _testVals(self,calc,testD):
    for smi,vals in testD:
      mol = MolFromSmiles(smi)
      ans = array(vals)
      res = array(calc.CalcDescriptors(mol))
      assert max(abs(res-ans))<1e-4,'bad descriptor values for SMILES %s (%s)'%(smi,str(res))
    
  def testCalcVals(self):
    self._testVals(self.calc,self.testD)

  def testSaveState(self):
    fName = os.path.join(RDConfig.RDCodeDir,'ML/Descriptors/test_data','molcalc.dsc')
    ok = 1
    try:
      inF = open(fName,'rb')
    except:
      ok = 0
    assert ok,'problems opening saved file %s'%(fName)
    try:
      calc = cPickle.load(inF)
    except:
      ok = 0
    assert ok,'problems reading saved file %s'%(fName)
      

    assert calc.GetDescriptorNames()==tuple(self.descs),'bad descriptor names'
    self._testVals(calc,self.testD)
    
    
if __name__ == '__main__':
  unittest.main()
