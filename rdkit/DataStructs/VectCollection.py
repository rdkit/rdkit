# $Id$
#
#  Copyright (C) 2005-2006 greg landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import copy
import struct

from rdkit import DataStructs


class VectCollection(object):
  """

    >>> vc = VectCollection()
    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((1,3,5))
    >>> vc.AddVect(1,bv1)
    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((6,8))
    >>> vc.AddVect(2,bv1)
    >>> len(vc)
    10
    >>> vc.GetNumBits()
    10
    >>> vc[0]
    0
    >>> vc[1]
    1
    >>> vc[9]
    0
    >>> vc[6]
    1
    >>> vc.GetBit(6)
    1
    >>> list(vc.GetOnBits())
    [1, 3, 5, 6, 8]

    keys must be unique, so adding a duplicate replaces the
    previous values:

    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((7,9))
    >>> vc.AddVect(1,bv1)
    >>> len(vc)
    10
    >>> vc[1]
    0
    >>> vc[9]
    1
    >>> vc[6]
    1

    we can also query the children:

    >>> vc.NumChildren()
    2
    >>> cs = vc.GetChildren()
    >>> id,fp = cs[0]
    >>> id
    1
    >>> list(fp.GetOnBits())
    [7, 9]
    >>> id,fp = cs[1]
    >>> id
    2
    >>> list(fp.GetOnBits())
    [6, 8]

    attach/detach operations:

    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((5,6))
    >>> vc.AddVect(3,bv1)
    >>> vc.NumChildren()
    3
    >>> list(vc.GetOnBits())
    [5, 6, 7, 8, 9]
    >>> vc.DetachVectsNotMatchingBit(6)
    >>> vc.NumChildren()
    2
    >>> list(vc.GetOnBits())
    [5, 6, 8]


    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((7,9))
    >>> vc.AddVect(1,bv1)
    >>> vc.NumChildren()
    3
    >>> list(vc.GetOnBits())
    [5, 6, 7, 8, 9]
    >>> vc.DetachVectsMatchingBit(6)
    >>> vc.NumChildren()
    1
    >>> list(vc.GetOnBits())
    [7, 9]


    to copy VectCollections, use the copy module:

    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((5,6))
    >>> vc.AddVect(3,bv1)
    >>> list(vc.GetOnBits())
    [5, 6, 7, 9]
    >>> vc2 = copy.copy(vc)
    >>> vc.DetachVectsNotMatchingBit(6)
    >>> list(vc.GetOnBits())
    [5, 6]
    >>> list(vc2.GetOnBits())
    [5, 6, 7, 9]

    The Uniquify() method can be used to remove duplicate vectors:

    >>> vc = VectCollection()
    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((7,9))
    >>> vc.AddVect(1,bv1)
    >>> vc.AddVect(2,bv1)
    >>> bv1 = DataStructs.ExplicitBitVect(10)
    >>> bv1.SetBitsFromList((2,3,5))
    >>> vc.AddVect(3,bv1)
    >>> vc.NumChildren()
    3
    >>> vc.Uniquify()
    >>> vc.NumChildren()
    2


    """

  def __init__(self):
    self.__vects = {}
    self.__orVect = None
    self.__numBits = -1
    self.__needReset = True

  def GetOrVect(self):
    if self.__needReset:
      self.Reset()
    return self.__orVect

  orVect = property(GetOrVect)

  def AddVect(self, idx, vect):
    self.__vects[idx] = vect
    self.__needReset = True

  def Reset(self):
    if not self.__needReset:
      return
    self.__orVect = None
    if not self.__vects:
      return
    ks = list(iter(self.__vects))
    self.__orVect = copy.copy(self.__vects[ks[0]])
    self.__numBits = self.__orVect.GetNumBits()
    for i in range(1, len(ks)):
      self.__orVect |= self.__vects[ks[i]]
    self.__needReset = False

  def NumChildren(self):
    return len(self.__vects.keys())

  def GetChildren(self):
    return tuple(self.__vects.items())

  def __getitem__(self, idx):
    if self.__needReset:
      self.Reset()
    return self.__orVect.GetBit(idx)

  GetBit = __getitem__

  def __len__(self):
    if self.__needReset:
      self.Reset()
    return self.__numBits

  GetNumBits = __len__

  def GetOnBits(self):
    if self.__needReset:
      self.Reset()
    return self.__orVect.GetOnBits()

  def DetachVectsNotMatchingBit(self, bit):
    items = list(self.__vects.items())
    for k, v in items:
      if not v.GetBit(bit):
        del (self.__vects[k])
        self.__needReset = True

  def DetachVectsMatchingBit(self, bit):
    items = list(self.__vects.items())
    for k, v in items:
      if v.GetBit(bit):
        del (self.__vects[k])
        self.__needReset = True

  def Uniquify(self, verbose=False):
    obls = {}
    for k, v in self.__vects.items():
      obls[k] = list(v.GetOnBits())

    keys = list(self.__vects.keys())
    nKeys = len(keys)
    keep = list(self.__vects.keys())
    for i in range(nKeys):
      k1 = keys[i]
      if k1 in keep:
        obl1 = obls[k1]
        idx = keys.index(k1)
        for j in range(idx + 1, nKeys):
          k2 = keys[j]
          if k2 in keep:
            obl2 = obls[k2]
            if obl1 == obl2:
              keep.remove(k2)

    self.__needsReset = True
    tmp = {}
    for k in keep:
      tmp[k] = self.__vects[k]
    if verbose:
      print('uniquify:', len(self.__vects), '->', len(tmp))
    self.__vects = tmp

  #
  # set up our support for pickling:
  #
  def __getstate__(self):
    pkl = struct.pack('<I', len(self.__vects))
    for k, v in self.__vects.items():
      pkl += struct.pack('<I', k)
      p = v.ToBinary()
      l = len(p)
      pkl += struct.pack('<I', l)
      pkl += struct.pack('%ds' % (l), p)
    return pkl

  def __setstate__(self, pkl):
    if isinstance(pkl, str):
      pkl = bytes(pkl, encoding='Latin1')

    self.__vects = {}
    self.__orVect = None
    self.__numBits = -1
    self.__needReset = True
    szI = struct.calcsize('I')
    offset = 0
    nToRead = struct.unpack('<I', pkl[offset:offset + szI])[0]
    offset += szI
    for _ in range(nToRead):
      k = struct.unpack('<I', pkl[offset:offset + szI])[0]
      offset += szI
      l = struct.unpack('<I', pkl[offset:offset + szI])[0]
      offset += szI
      sz = struct.calcsize('%ds' % l)
      bv = DataStructs.ExplicitBitVect(struct.unpack('%ds' % l, pkl[offset:offset + sz])[0])
      offset += sz
      self.AddVect(k, bv)


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
