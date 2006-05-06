#
#  Copyright (C) 2000  greg Landrum
#

""" unit testing code for the virtual library enumerator

"""
import unittest
import random

from EnumNodes import *

class TestCase(unittest.TestCase):
  def setUp(self):
    random.seed(23)

  def testSupply(self):
    s = SupplyNode(['foo','bar','baz'])
    s.AddSupply('grn')
    res = [s.GetNext(),s.GetNext()]
    assert res==['foo','bar'], 'testSupply Failed'

  def testSupplySlice(self):
    s = SupplyNode(['foo','bar','baz'])
    s.AddSupply('grn')
    res = s[:]
    assert res==['foo','bar','baz','grn'], 'testSupplySlice Failed'

  def testRandomSupply(self):
    s = SupplyNode(['foo','bar','baz'])
    s.AddSupply('grn')
    res = [s.GetRandom(),s.GetRandom(),s.GetRandom()]
    expect = ['foo','foo','baz']
    assert res==expect, 'testRandomSupply Failed: %s != %s'%(str(res),str(expect))

  def testCombine(self):
    s1 = SupplyNode([1,2,3],'s1')
    s2 = SupplyNode(['a','b','c'],'s2')
    s3 = SupplyNode([100,200],'s3')
    c = CombineNode()
    c.AddInput(s1)
    c.AddInput(s2)
    c.AddInput(s3)
    res = [c.GetNext(),c.GetNext()]
    assert res==[(1, 'a', 100), (1, 'a', 200)], 'testCombine Failed'

  def testCombine2(self):
    s1 = SupplyNode(name='s1')
    s2 = SupplyNode(['a','b','c'],'s2')
    c = CombineNode()
    c.AddInput(s1)
    c.AddInput(s2)
    s1.AddSupply(1)
    s1.AddSupply(2)
    assert len(c) == 6, 'testCombine2 Add failed'
    s1.RemoveSupply(2)
    assert len(c) == 3, 'testCombine2 Remove failed'
    
  def testRandomCombine(self):
    s1 = SupplyNode([1,2,3],'s1')
    s2 = SupplyNode(['a','b','c'],'s2')
    s3 = SupplyNode([100,200],'s3')
    c = CombineNode()
    c.AddInput(s1)
    c.AddInput(s2)
    c.AddInput(s3)
    res = [c.GetRandom(),c.GetRandom(),c.GetRandom()]
    expect = [(1, 'b', 100), (1, 'b', 100),(2,'b',200)]
    assert res==expect, 'testRandomCombine Failed: %s != %s'%(str(res),str(expect))

  def _composeCompare(self,compos1,compos2):
    if compos1[0] != compos2[0]:
      return 0
    if compos1[1] != compos2[1]:
      return 0
    if abs(compos1[2][0]-compos2[2][0])>0.001:
      return 0
    if abs(compos1[2][1]-compos2[2][1])>0.001:
      return 0
    return 1

  def testAlloy(self):
    s1 = SupplyNode(['a','b'],'s1')
    s2 = SupplyNode(['d','e'],'s2')
    s3 = SupplyNode(['g','h'],'s3')
    a = AlloyNode()
    a.AddInput(s1)
    a.AddInput(s2,compRange=(0.3,0.6))
    a.AddInput(s3,compRange=(0.5,0.8))

    a.RemoveInput(s2)
    res = a[0:3]
    cVect = map(self._composeCompare,res,[('a', 'g', [0.2, 0.8]),
                                          ('a', 'g', [0.3, 0.7]),
                                          ('a', 'g', [0.4, 0.6])])
    assert cVect == [1,1,1],'testAlloy Failed'

  def testRandomAlloy(self):
    s1 = SupplyNode(['a','b'],'s1')
    s3 = SupplyNode(['g','h'],'s3')
    a = AlloyNode()
    a.AddInput(s1)
    a.AddInput(s3,compRange=(0.5,0.8))
    res = [a.GetRandom(),a.GetRandom()]
    expect = [('a', 'g', [0.4, 0.6]),
              ('a', 'g', [0.3, 0.7])]
    cVect = map(self._composeCompare,res,expect)
    assert cVect == [1,1],'testRandomAlloy Failed: %s != %s'%(res,expect)
    
  def testTransform(self):
    s1 = SupplyNode(['foo','bar'],'s1')
    s2 = SupplyNode(['baz','grn'],'s2')
    t = TransformNode(lambda x,y:x+y)
    t.AddInput(s1)
    t.AddInput(s2)
    res=t[0:]
    assert res==['foobaz','bargrn'],'testTransform Failed'

  def testFilter(self):
    s1 = SupplyNode([1,2,3,4],'s1')
    s2 = SupplyNode([1,2,3,4],'s2')
    f = FilterNode(lambda x,y:x+y>3)
    f.AddInput(s1)
    f.AddInput(s2)
    res = [f.GetNext(),f.GetNext()]
    assert res == [(2,2),(3,3)],'testFilter Failed'

  def testRandomFilter(self):
    s1 = SupplyNode([1,2,3,4],'s1')
    s2 = SupplyNode([1,2,3,4],'s2')
    f = FilterNode(lambda x,y:x+y>3)
    f.AddInput(s1)
    f.AddInput(s2)
    res = [f.GetRandom(),f.GetRandom()]
    expect = [(3,3),(4,4)]
    assert res==expect,'testRandomFilter Failed: %s != %s'%(res,expect)
    

def TestSuite():
  suite = unittest.TestSuite()
  suite.addTest(TestCase('testSupply'))
  suite.addTest(TestCase('testSupplySlice'))
  suite.addTest(TestCase('testRandomSupply'))
  suite.addTest(TestCase('testCombine'))
  suite.addTest(TestCase('testCombine2'))
  suite.addTest(TestCase('testRandomCombine'))
  suite.addTest(TestCase('testAlloy'))
  suite.addTest(TestCase('testRandomAlloy'))
  suite.addTest(TestCase('testTransform'))
  suite.addTest(TestCase('testFilter'))
  suite.addTest(TestCase('testRandomFilter'))
  return suite

if __name__ == '__main__':
  suite = TestSuite()
  unittest.TextTestRunner().run(suite)
