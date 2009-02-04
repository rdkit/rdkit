# $Id$
#
#  Copyright (C) 2005-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
def ConstructEnsembleBV(bv,bitsToKeep):
  """

  >>> import DataStructs
  >>> bv = DataStructs.ExplicitBitVect(128)
  >>> bv.SetBitsFromList((1,5,47,99,120))
  >>> r = ConstructEnsembleBV(bv,(0,1,2,3,45,46,47,48,49))
  >>> r.GetNumBits()
  9
  >>> r.GetBit(0)
  0
  >>> r.GetBit(1)
  1
  >>> r.GetBit(5)
  0
  >>> r.GetBit(6)  # old bit 47
  1
  
  

  """
  finalSize=len(bitsToKeep)
  res = bv.__class__(finalSize)


  for i,bit in enumerate(bitsToKeep):
    if bv.GetBit(bit):
      res.SetBit(i)
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
