# $Id$
#
#  Copyright (C) 2005-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#


def ConstructEnsembleBV(bv, bitsToKeep):
  """

  >>> from rdkit import DataStructs
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
  finalSize = len(bitsToKeep)
  res = bv.__class__(finalSize)

  for i, bit in enumerate(bitsToKeep):
    if bv.GetBit(bit):
      res.SetBit(i)
  return res


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
