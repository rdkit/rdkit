import RDConfig
import os,sys,cPickle
import unittest
import DataStructs as ds

class TestCase(unittest.TestCase):
  def setUp(self) :
    pass

  def test1Int(self):
    """

    """
    v1 = ds.IntSparseIntVect(5)
    self.failUnlessRaises(IndexError,lambda:v1[6])
    v1[0]=1
    v1[2]=2
    v1[3]=3
    self.failUnless(v1==v1)

    v2= ds.IntSparseIntVect(5)
    self.failUnless(v1!=v2)
    v2|=v1
    self.failUnless(v2==v1)

    v3=v2|v1
    self.failUnless(v3==v1)

  def test2Long(self):
    """

    """
    l=1L<<42
    v1 = ds.LongSparseIntVect(l)
    self.failUnlessRaises(IndexError,lambda:v1[l+1])
    v1[0]=1
    v1[2]=2
    v1[1L<<35]=3
    self.failUnless(v1==v1)

    v2= ds.LongSparseIntVect(l)
    self.failUnless(v1!=v2)
    v2|=v1
    self.failUnless(v2==v1)

    v3=v2|v1
    self.failUnless(v3==v1)

  def test3Pickle(self):
    """

    """
    l=1L<<42
    v1 = ds.LongSparseIntVect(l)
    self.failUnlessRaises(IndexError,lambda:v1[l+1])
    v1[0]=1
    v1[2]=2
    v1[1L<<35]=3
    self.failUnless(v1==v1)

    v2=  cPickle.loads(cPickle.dumps(v1))
    self.failUnless(v2==v1)
    
    v3=  ds.LongSparseIntVect(v2.ToBinary())
    self.failUnless(v2==v3)
    self.failUnless(v1==v3)
    
  def test4Update(self):
    """

    """
    v1 = ds.IntSparseIntVect(5)
    self.failUnlessRaises(IndexError,lambda:v1[6])
    v1[0]=1
    v1[2]=2
    v1[3]=3
    self.failUnless(v1==v1)

    v2 = ds.IntSparseIntVect(5)
    v2.UpdateFromSequence((0,2,3,3,2,3))
    self.failUnless(v1==v2)
    

    
    
    
if __name__ == '__main__':
    unittest.main()
