# $Id$
#
#  Copyright (C) 2001-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" This has all be re-implemented in the C++ code
"""
from __future__ import print_function
import math

class BitVect:
  def __init__(self,nBits):
    self.nBits = nBits
    self.bits = [0]*nBits

  def NumOnBits(self):
    return len(self.GetOnBits())

  def GetOnBits(self,sort=1,reverse=0):
    l = [idx for idx in xrange(self.nBits) if self.bits[idx] == 1]
    if reverse:
      l.reverse()
    return l
  
  def TanimotoSimilarity(self,other):
    if not isinstance(other,BitVect):
      raise TypeError("Tanimoto similarities can only be calculated between two BitVects")
    if len(self)!=len(other):
      raise ValueError("BitVects must be the same length")
    bc = len(self & other)
    b1 = self.NumOnBits()
    b2 = other.NumOnBits()
    return float(bc) / float(b1 + b2 - bc)
  TanimotoSimilarity = TanimotoSimilarity

  def EuclideanDistance(self,other):
    if not isinstance(other,BitVect):
      raise TypeError("Tanimoto similarities can only be calculated between two BitVects")
    bt = len(self)
    bi = len(self ^ (~ other))
    return math.sqrt(bt-bi)/bt
  
  def __getitem__(self,which):
    if which >= self.nBits or which < 0:
      raise ValueError('bad index')
    return self.bits[which]

  def __setitem__(self,which,val):
    if which >= self.nBits or which < 0:
      raise ValueError('bad index')
    if val not in [0,1]:
      raise ValueError('val must be 0 or 1')

    self.bits[which] = val
      
  def __len__(self):
    return self.nBits

  def __and__(self,other):
    if not isinstance(other,BitVect):
      raise TypeError("BitVects can only be &'ed with other BitVects")
    if len(self) != len(other):
      raise ValueError("BitVects must be of the same length")

    l1 = self.GetOnBits()
    l2 = other.GetOnBits()
    r = [bit for bit in l1 if bit in l2]
    return r

  def __or__(self,other):
    if not isinstance(other,BitVect):
      raise TypeError("BitVects can only be |'ed with other BitVects")
    if len(self) != len(other):
      raise ValueError("BitVects must be of the same length")
    l1 = self.GetOnBits()
    l2 = other.GetOnBits()
    r = l1 + [bit for bit in l2 if bit not in l1]
    r.sort()
    return r

  def __xor__(self,other):
    if not isinstance(other,BitVect):
      raise TypeError("BitVects can only be ^'ed with other BitVects")
    if len(self) != len(other):
      raise ValueError("BitVects must be of the same length")
    
    l1 = self.GetOnBits()
    l2 = other.GetOnBits()
    r = [bit for bit in l1 if bit not in l2] + [bit for bit in l2 if bit not in l1]
    r.sort()
    return r

  def __invert__(self):
    res = BitVect(len(self))
    for i in xrange(len(self)):
      res[i] = not self[i]
    return res
  
class SparseBitVect(BitVect):
  def __init__(self,nBits):
    self.nBits = nBits
    self.bits = []

  def NumOnBits(self):
    return len(self.bits)

  def GetOnBits(self,sort=1,reverse=0):
    l = self.bits[:]
    if sort:
      l.sort()
    if reverse:
      l.reverse()
    return l
  
  def __getitem__(self,which):
    if which >= self.nBits or which < 0:
      raise ValueError('bad index')
    if which in self.bits:
      return 1
    else:
      return 0

  def __setitem__(self,which,val):
    if which >= self.nBits or which < 0:
      raise ValueError('bad index')
    if val == 0:
      if which in self.bits:
        self.bits.remove(which)
    else:
      self.bits.append(which)
      
  def __len__(self):
    return self.nBits

if __name__ == '__main__':
  b1 = BitVect(10)
  b2 = SparseBitVect(10)
  b1[0] = 1
  b2[0] = 1
  b1[3] = 1
  b2[4] = 1
  b2[5] = 1
  b2[5] = 0
  print('b1:',b1.GetOnBits())
  print('b2:',b2.GetOnBits())
  print('&:', b1 & b2)
  print('|:', b1 | b2)
  print('^:', b1 ^ b2)
  print('b1.Tanimoto(b2):',b1.TanimotoSimilarity(b2))
  print('b1.Tanimoto(b1):',b1.TanimotoSimilarity(b1))
  print('b2.Tanimoto(b2):',b2.TanimotoSimilarity(b2))
  print('b2.Tanimoto(b1):',b2.TanimotoSimilarity(b1))
  
  
