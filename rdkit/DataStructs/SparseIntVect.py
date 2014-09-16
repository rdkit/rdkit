# $Id$
#
#  Copyright (C) 2007 Greg Landrum
#   All Rights Reserved
#
import copy

class pySparseIntVect(object):
  """ this class is pretty much obsolete (it's in C++ now) """
  size=0
  container={}
  def __init__(self,size):
    self.size=size
    self.container={}

  def UpdateFromSequence(self,seq):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1.UpdateFromSequence((0,1,1,5))
    >>> [x for x in c1]
    [(0, 1), (1, 2), (5, 1)]
    >>> c1.UpdateFromSequence((0,3))
    >>> [x for x in c1]
    [(0, 2), (1, 2), (3, 1), (5, 1)]

    """
    for v in seq:
      self[v] += 1
  def InitFromSequence(self,seq):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1.InitFromSequence((0,1,1,5))
    >>> [x for x in c1]
    [(0, 1), (1, 2), (5, 1)]

    """
    self.container={}
    self.UpdateFromSequence(seq)
      
  def Sum(self,useAbs=False):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c1.Sum()
    10

    >>> c1[2] = -2
    >>> c1.Sum()
    6
    >>> c1.Sum(useAbs=True)
    10
    """
    res=0
    for v in self.container.values():
      if not useAbs:
        res+=v
      else:
        res+=abs(v)
    return res
  def GetTotalVal(self,useAbs=False):
    return self.Sum(useAbs=useAbs)
  

  def __eq__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 3
    >>> c2[2] = 2
    >>> c1 == c2
    False
    >>> c1 == c1
    True
    """
    if not isinstance(other,pySparseIntVect):
      raise TypeError
    if self.size != other.size:
      return 0
    return self.container==other.container
    

  def __iand__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[2] = -2
    >>> c2[5] = 6
    >>> c1 &= c2
    >>> [x for x in c1]
    [(0, 2), (2, -2)]

    """
    if not isinstance(other,pySparseIntVect):
      raise TypeError
    if self.size != other.size:
      raise ValueError

    newC = {}
    for idx,v in self.container.iteritems():
      ov = other.container.get(idx,0)
      if ov:
        if v<ov:
          newC[idx]=v
        else:
          newC[idx]=ov
    self.container=newC
    return self
  def __ior__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[2] = -2
    >>> c2[5] = 6
    >>> c1 |= c2
    >>> [x for x in c1]
    [(0, 3), (2, 2), (4, 5), (5, 6)]

    """
    if not isinstance(other,pySparseIntVect):
      raise TypeError
    if self.size != other.size:
      raise ValueError

    newC = {}
    for idx,v in self.container.iteritems():
      ov = other.container.get(idx,0)
      if v<ov:
        newC[idx]=ov
      else:
        newC[idx]=v
    for k,v in other.container.iteritems():
      if not newC.has_key(k):
        newC[k]=v
    self.container=newC
    return self
  
  def __iadd__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[2] = -2
    >>> c2[5] = 6
    >>> c1 += c2
    >>> [x for x in c1]
    [(0, 5), (4, 5), (5, 6)]

    """
    if not isinstance(other,pySparseIntVect):
      raise TypeError
    if self.size != other.size:
      raise ValueError
    seen={}
    for idx in self.container.keys():
      seen[idx]=1
      v = self.container[idx]+other[idx]
      if v:
        self.container[idx]=v
      else:
        del self.container[idx]
    for idx,v in other:
      if not seen.has_key(idx):
        self.container[idx]=v
    return self
  
  def __isub__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[2] = 2
    >>> c2[5] = 6
    >>> c1 -= c2
    >>> [x for x in c1]
    [(0, 1), (4, 5), (5, -6)]

    """
    if not isinstance(other,pySparseIntVect):
      raise TypeError
    if self.size != other.size:
      raise ValueError
    seen={}
    for idx in self.container.keys():
      seen[idx]=1
      v = self.container[idx]-other[idx]
      if v:
        self.container[idx]=v
      else:
        del self.container[idx]
    for idx,v in other:
      if not seen.has_key(idx):
        self.container[idx]=-v
    return self

  def __imul__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[5] = 6
    >>> c1 *= c2
    >>> [x for x in c1]
    [(0, 6)]

    """
    if not isinstance(other,pySparseIntVect):
      raise TypeError
    if self.size != other.size:
      raise ValueError
    for idx in self.container.keys():
      v = self.container[idx]*other[idx]
      if v:
        self.container[idx]=v
      else:
        del self.container[idx]
    return self
  
  def __add__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[5] = 6
    >>> c3 = c2+c1
    >>> [x for x in c3]
    [(0, 5), (4, 5), (5, 6)]

    """
    res = pySparseIntVect(self.size)
    res.container = copy.deepcopy(self.container)
    res += other
    return res
  def __sub__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[2] = 2
    >>> c2[5] = 6
    >>> c3 = c1-c2
    >>> [x for x in c3]
    [(0, 1), (4, 5), (5, -6)]
    >>> [x for x in c1]
    [(0, 3), (2, 2), (4, 5)]

    """
    res = pySparseIntVect(self.size)
    res.container = copy.deepcopy(self.container)
    res -= other
    return res
  def __mul__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[5] = 6
    >>> c3 = c1*c2
    >>> [x for x in c3]
    [(0, 6)]
    >>> [x for x in c1]
    [(0, 3), (4, 5)]

    """
    res = pySparseIntVect(self.size)
    res.container = copy.deepcopy(self.container)
    res *= other
    return res
  def __and__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[2] = -2
    >>> c2[5] = 6
    >>> c3 = c1 & c2
    >>> [x for x in c3]
    [(0, 2), (2, -2)]
    >>> [x for x in c1]
    [(0, 3), (2, 2), (4, 5)]

    """
    res = pySparseIntVect(self.size)
    res.container = copy.deepcopy(self.container)
    res &= other
    return res
  def __or__(self,other):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[2] = 2
    >>> c1[4] = 5
    >>> c2=pySparseIntVect(10)
    >>> c2[0] = 2
    >>> c2[2] = -2
    >>> c2[5] = 6
    >>> c3 = c1 | c2
    >>> [x for x in c3]
    [(0, 3), (2, 2), (4, 5), (5, 6)]
    >>> [x for x in c1]
    [(0, 3), (2, 2), (4, 5)]

    """
    res = pySparseIntVect(self.size)
    res.container = copy.deepcopy(self.container)
    res |= other
    return res
  
  def __len__(self):
    return self.size
  def __getitem__(self,which):
    """
    >>> c1=pySparseIntVect(10)
    >>> c1[0] = 3
    >>> c1[4] = 5
    >>> c1[0]
    3
    >>> c1[1]
    0
    
    """
    if abs(which)>=self.size:
      raise IndexError(which)
    if which<0:
      which = self.size-which
    return self.container.get(which,0)
  def __setitem__(self,which,val):
    if abs(which)>=self.size:
      raise IndexError(which)
    val = int(val)
    if which<0:
      which = self.size-which
    self.container[which]=val
  def __iter__(self):
    """
    >>> c=pySparseIntVect(10)
    >>> c[0] = 3
    >>> c[4] = 5
    >>> c[7] = -1
    >>> for idx,v in c:
    ...  print idx,v
    0 3
    4 5
    7 -1

    """
    return self.container.iteritems()

    
from rdkit import DataStructs
DiceSimilarity=DataStructs.DiceSimilarity
def pyDiceSimilarity(v1,v2,bounds=None,useAbs=False):
  """ Implements the DICE similarity metric.

  >>> v1 = DataStructs.IntSparseIntVect(10)
  >>> v2 = DataStructs.IntSparseIntVect(10)
  >>> v1.UpdateFromSequence((1,2,3))
  >>> v2.UpdateFromSequence((1,2,3))
  >>> DiceSimilarity(v1,v2)
  1.0

  >>> v2 = DataStructs.IntSparseIntVect(10)
  >>> v2.UpdateFromSequence((5,6))
  >>> DiceSimilarity(v1,v2)
  0.0

  >>> v1 = DataStructs.IntSparseIntVect(10)
  >>> v2 = DataStructs.IntSparseIntVect(10)
  >>> v1.UpdateFromSequence((1,2,3,4))
  >>> v2.UpdateFromSequence((1,3,5,7))
  >>> DiceSimilarity(v1,v2)
  0.5

  >>> v1 = DataStructs.IntSparseIntVect(10)
  >>> v2 = DataStructs.IntSparseIntVect(10)
  >>> v1.UpdateFromSequence((1,2,3,4,5,6))
  >>> v2.UpdateFromSequence((1,3))
  >>> DiceSimilarity(v1,v2)
  0.5

  """
  denom = 1.0*(v1.GetTotalVal(useAbs=useAbs)+v2.GetTotalVal(useAbs=useAbs))
  if not denom:
    res = 0.0
  else:
    if bounds and (min(len(v1),len(v2))/denom) < bounds:
      numer = 0.0
    else:
      numer=0.0
      v3=v1&v2
      numer=v3.GetTotalVal(useAbs=useAbs)
    res = 2.*numer/denom

  return res

def DotProduct(bv1,bv2):
  res = 0.0
  for k,v in bv1.GetNonzeroElements().iteritems():
    res += v*bv2[k]
  return res

#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
  
