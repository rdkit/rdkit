# $Id$
#
#  Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""basic unit testing code for the Bond wrapper

"""
import unittest
from rdkit import Chem


class TestCase(unittest.TestCase):

  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    self.m = Chem.MolFromSmiles('CCCC1=CC=C1')

  def test1Get(self):
    " testing GetBond "
    ok = 1
    try:
      b = self.m.GetBondBetweenAtoms(0, 1)
    except Exception:
      ok = 0
    assert ok, 'GetBond failed'

  def test2Setters(self):
    " testing setting bond props "
    b = self.m.GetBondBetweenAtoms(0, 1)
    assert b.GetBondType() == Chem.BondType.SINGLE
    b.SetBondDir(Chem.BondDir.BEGINWEDGE)
    assert self.m.GetBondBetweenAtoms(0, 1).GetBondDir() == Chem.BondDir.BEGINWEDGE
    b = self.m.GetBondBetweenAtoms(0, 1)

  def test3Props(self):
    " testing bond props "
    b = self.m.GetBondBetweenAtoms(0, 1)
    assert b.GetBondType() == Chem.BondType.SINGLE
    assert b.GetBeginAtom().GetIdx() == self.m.GetAtomWithIdx(0).GetIdx()
    assert b.GetBeginAtomIdx() == 0
    assert b.GetEndAtom().GetIdx() == self.m.GetAtomWithIdx(1).GetIdx()
    assert b.GetEndAtomIdx() == 1

  def test4Props2(self):
    " testing more bond props "
    b = self.m.GetBondBetweenAtoms(3, 4)
    assert b.GetBondType() == Chem.BondType.DOUBLE
    b2 = self.m.GetBondBetweenAtoms(1, 2)
    assert b2.GetBondType() == Chem.BondType.SINGLE
    assert b.GetIsConjugated()
    assert not b2.GetIsConjugated()


if __name__ == '__main__':
  unittest.main()
