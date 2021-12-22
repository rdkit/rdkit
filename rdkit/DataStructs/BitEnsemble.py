# $Id$
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" #DOC


"""


class BitEnsemble(object):
  """  used to store a collection of bits and score
  BitVects (or signatures) against them.

  """

  def __init__(self, bits=None):
    self._bits = list(bits) if bits is not None else []

  def SetBits(self, bits):
    self._bits = list(bits) if bits is not None else []

  def AddBit(self, bit):
    self._bits.append(bit)

  def GetBits(self):
    return tuple(self._bits)

  def GetNumBits(self):
    return len(self._bits)

  def ScoreWithOnBits(self, other):
    """ other must support GetOnBits() """
    obl = other.GetOnBits()
    thisBits = self.GetBits()
    return len([None for bit in thisBits if bit in obl])

  def ScoreWithIndex(self, other):
    """ other must support __getitem__() """
    thisBits = self.GetBits()
    return len([None for bit in thisBits if other[bit]])


if __name__ == '__main__':

  pass
