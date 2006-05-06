#  Copyright (C) 2004-2005  Rational Discovery LLC
#    All Rights Reserved
#
"""unit testing code for the SVD wrapper

"""
import RDConfig
import unittest,os.path
import PySVD
from PySVD import cSVD

from Numeric import *

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):
  def setUp(self):
    if doLong:
      print '\n%s: '%self.shortDescription(),

  def _checkOrthog(self,u):
    for i in range(u.shape[0]):
      for j in range(i,u.shape[0]):
        dotProd = dot(u[i],u[j])
        if i!=j: assert feq(dotProd,0)
        if i==j: assert feq(dotProd,1)

  def test1(self):
    " simple SVD case, using defaults "
    indices = ((0,3),(1,3),(2,3))
    vals = (1.0,1.0,2.0,1.5,3.0,1.7)
    u,s,v = cSVD.SparseSVD(indices,vals,3,4)
    assert u.shape==(3,3)
    assert s.shape==(3,)
    assert v is None
    self._checkOrthog(u)

  def test2(self):
    " simple SVD case, changing defaults "
    indices = ((0,3),(1,3),(2,3))
    vals = (1.0,1.0,2.0,1.5,3.0,1.7)
    u,s,v = cSVD.SparseSVD(indices,vals,3,4,3)
    assert u.shape==(3,3)
    assert s.shape==(3,)
    assert v is None
    self._checkOrthog(u)

    u,s,v = cSVD.SparseSVD(indices,vals,3,4,2)
    assert u.shape==(2,3)
    assert s.shape==(2,)
    assert v is None
    self._checkOrthog(u)

    u,s,v = cSVD.SparseSVD(indices,vals,3,4,2,1)
    assert u.shape==(2,3)
    assert s.shape==(2,)
    assert v is not None
    assert v.shape==(2,4)
    self._checkOrthog(u)
    self._checkOrthog(v)

    try:
      u,s,v = cSVD.SparseSVD(indices,vals,3,4,5)
    except ValueError:
      ok=1
    else:
      ok=0
    assert ok
     
      
  def test3(self):
    " simple SVD bitmatrix cases "
    indices = ((0,3),(1,3),(2,3))
    u,s,v = cSVD.SparseSVDBitMatrix(indices,3,4)
    assert u.shape==(3,3)
    assert s.shape==(3,)
    assert v is None
    self._checkOrthog(u)

    u,s,v = cSVD.SparseSVDBitMatrix(indices,3,4,2,1)
    assert u.shape==(2,3)
    assert s.shape==(2,)
    assert v is not None
    assert v.shape==(2,4)
    self._checkOrthog(u)
    self._checkOrthog(v)

  def test4(self):
    " SVD with dupes "
    indices = ((0,3),(1,3),(0,3))
    vals = (1.0,1.0,2.0,1.5,1.0,2.0)
    u,s,v = cSVD.SparseSVD(indices,vals,3,4)
    assert u.shape==(3,3)
    assert s.shape==(3,)
    assert v is None
    self._checkOrthog(u)




if __name__ == '__main__':
  import sys,getopt,re
  doLong=0
  if len(sys.argv) >1:
    args,extras=getopt.getopt(sys.argv[1:],'l')
    for arg,val in args:
      if arg=='-l':
        doLong=1
      sys.argv.remove('-l')
  if doLong:
    for methName in dir(TestCase):
      if re.match('_test',methName):
        newName = re.sub('_test','test',methName)
        exec('TestCase.%s = TestCase.%s'%(newName,methName))
        
  unittest.main()

