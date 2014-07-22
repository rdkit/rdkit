#
#  Copyright (C) 2001  greg Landrum
#

""" unit testing code for compound descriptors

"""
from __future__ import print_function
import unittest
import Parser

from rdkit.six.moves import xrange    

class TestCase(unittest.TestCase):
  def setUp(self):
    print('\n%s: '%self.shortDescription(),end='')
    self.piece1 = [['d1','d2'],['d1','d2']]
    self.aDict = {'Fe':{'d1':1,'d2':2},'Pt':{'d1':10,'d2':20}}
    self.pDict = {'d1':100.,'d2':200.}
    self.compos = [('Fe',1),('Pt',1)]
    self.cExprs = ["SUM($1)","SUM($1)+SUM($2)","MEAN($1)","DEV($2)","MAX($1)","MIN($2)","SUM($1)/$a"]
    self.results = [11.,33.,5.5,9.,10.,2.,0.11]
    self.tol = 0.0001
  def testSingleCalcs(self):
    " testing calculation of a single descriptor "
    for i in xrange(len(self.cExprs)):
      cExpr= self.cExprs[i]
      argVect = self.piece1 + [cExpr]
      res = Parser.CalcSingleCompoundDescriptor(self.compos,argVect,self.aDict,self.pDict)
      self.assertAlmostEqual(res,self.results[i],2)
  def testMultipleCalcs(self):
    " testing calculation of multiple descriptors "
    for i in xrange(len(self.cExprs)):
      cExpr= self.cExprs[i]
      argVect = self.piece1 + [cExpr]
      res = Parser.CalcMultipleCompoundsDescriptor([self.compos,self.compos],argVect,
                                                  self.aDict,[self.pDict,self.pDict])
      self.assertAlmostEqual(res[0],self.results[i],2)
      self.assertAlmostEqual(res[1],self.results[i],2)
      #self.assertTrue(abs(res[0]-self.results[i])<self.tol,'Expression %s failed'%(cExpr))
      #self.assertTrue((res[1]-self.results[i])<self.tol,'Expression %s failed'%(cExpr))

def TestSuite():
  suite = unittest.TestSuite()
  suite.addTest(TestCase('testSingleCalcs'))
  suite.addTest(TestCase('testMultipleCalcs'))
  return suite


if __name__ == '__main__':
  suite = TestSuite()
  unittest.TextTestRunner().run(suite)
