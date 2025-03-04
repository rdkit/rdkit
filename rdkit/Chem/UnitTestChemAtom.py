#
#  Copyright (C) 2003-2024 Greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""basic unit testing code for atoms

"""
import unittest

from rdkit import Chem


class TestCase(unittest.TestCase):

  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    self.m = Chem.MolFromSmiles('CC(=O)CCSC')

  def test1Implicit(self):
    " testing ImplicitValence "
    a = self.m.GetAtoms()[0]
    iV = a.GetValence(which=Chem.ValenceType.IMPLICIT)
    assert iV == 3
    assert self.m.GetAtomWithIdx(1).GetValence(which=Chem.ValenceType.IMPLICIT) == 0
    assert self.m.GetAtomWithIdx(2).GetValence(which=Chem.ValenceType.IMPLICIT) == 0
    assert self.m.GetAtomWithIdx(3).GetValence(which=Chem.ValenceType.IMPLICIT) == 2

  def test2BondIter(self):
    " testing bond iteration "
    a = self.m.GetAtomWithIdx(1)
    bs = a.GetBonds()
    r = []
    for b in bs:
      r.append(b)
    assert len(r) == 3

  def test3GetBond(self):
    " testing GetBondBetweenAtoms(idx,idx) "
    b = self.m.GetBondBetweenAtoms(1, 2)
    assert b.GetBondType() == Chem.BondType.DOUBLE, 'GetBond failed'

  def test4Props(self):
    " testing atomic props "
    a = self.m.GetAtomWithIdx(1)
    assert a.GetSymbol() == 'C'
    assert a.GetAtomicNum() == 6
    assert a.GetFormalCharge() == 0
    assert a.GetDegree() == 3
    assert a.GetValence(which=Chem.ValenceType.IMPLICIT) == 0
    assert a.GetValence(which=Chem.ValenceType.EXPLICIT) == 4

  def test5Setters(self):
    " testing setting atomic props "
    a = Chem.Atom(6)
    assert a.GetSymbol() == 'C'
    assert a.GetAtomicNum() == 6
    a.SetFormalCharge(1)
    assert a.GetFormalCharge() == 1
    assert a.GetValence(which=Chem.ValenceType.IMPLICIT) == 0


if __name__ == '__main__':
  unittest.main()
