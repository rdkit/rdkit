# $Id$
#
#  Copyright (C) 2005 Rational Discovery LLC
#   All Rights Reserved
#


class LazySig:
  def __init__(self,computeFunc,sigSize):
    """
    computeFunc should take a single argument, the integer bit id
    to compute

    """
    if sigSize<=0:
      raise ValueError('zero size')
    self.computeFunc=computeFunc
    self.size=sigSize
    self._cache={}

  def __len__(self):
    """

     >>> obj = LazySig(lambda x:1,10)
     >>> len(obj)
     10

    """
    return self.size

  def __getitem__(self,which):
    """

     >>> obj = LazySig(lambda x:x,10)
     >>> obj[1]
     1
     >>> obj[-1]
     9
     >>> try:
     ...   obj[10]
     ... except IndexError:
     ...   1
     ... else:
     ...   0
     1
     >>> try:
     ...   obj[-10]
     ... except IndexError:
     ...   1
     ... else:
     ...   0
     1
     
    """
    if which<0:
      # handle negative indices
      which = self.size+which
    
    if which<=0 or which>=self.size:
      raise IndexError('bad index')
    
    if which in self._cache:
      v= self._cache[which]
    else:
      v = self.computeFunc(which)
      self._cache[which]=v
    return v
  
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

