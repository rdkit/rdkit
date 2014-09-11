#
#  Copyright (C) 2001  greg Landrum
#

""" unit testing code for the descriptor COM server

"""
from __future__ import print_function
from rdkit import RDConfig
import unittest
import Parser
from win32com.client import Dispatch
from Numeric import *

class TestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
  def testConnect(self):
    " testing connection "
    ok = 1
    try:
      c = Dispatch('RD.DescCalc')
    except:
      ok = 0
    assert ok and c is not None, 'connection to COM server failed'
  def testLoad(self):
    " testing load "
    c = Dispatch('RD.DescCalc')
    ok = 1
    try:
      c.LoadCalculator(RDConfig.RDCodeDir+'/ml/descriptors/test_data/ferro.dsc')
    except:
      ok = 0
    assert ok, 'LoadCalculator failed'
  def testNames(self):
    " testing GetDescriptorNames "
    c = Dispatch('RD.DescCalc')
    c.LoadCalculator(RDConfig.RDCodeDir+'/ml/descriptors/test_data/ferro.dsc')
    names = c.GetDescriptorNames()
    expectedNames = ('MAX_DED','has3d','has4d','has5d','elconc','atvol')
    assert names==expectedNames, 'GetDescriptorNames failed (%s != %s)'%(repr(names),
                                                                         repr(expectedNames))

  def testCalc(self):
    " testing descriptor calculation "
    argV = ['CrPt3','fcc','AuCu3',58.09549962,1,4,0.228898,8.876,1]
    nameV = ['Compound','Structure','Structure_Type','Volume',
             'Z','Atoms_per_Formula_Unit','Hardness','RawDOS_Ef',
             'IsFerromagnetic']
    c = Dispatch('RD.DescCalc')
    c.LoadCalculator(RDConfig.RDCodeDir+'/ml/descriptors/test_data/ferro.dsc')
    ok = 1
    descVect = array(c.CalcDescriptors(argV,nameV))
    expected = array((3.67481803894, 1, 0, 1, 0.619669341609, 14.523874905))
    diffV = abs(descVect-expected)
    assert ok, 'CalcDescriptors failed'
    assert max(diffV)<0.0001,'bad descriptors: %s, %s'%(str(expected),str(descVect))
    
def TestSuite():
  suite = unittest.TestSuite()
  suite.addTest(TestCase('testConnect'))
  suite.addTest(TestCase('testLoad'))
  suite.addTest(TestCase('testNames'))
  suite.addTest(TestCase('testCalc'))
  return suite


if __name__ == '__main__':
  suite = TestSuite()
  unittest.TextTestRunner().run(suite)
