# $Id: UnitTestFragments.py 5007 2006-02-22 15:14:41Z glandrum $
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Unit Test code for cFragments module

"""
import unittest
from Numeric import *
import LinearAlgebra
import Chem

from cFragments import *

import copy

_verbose = 0
_zeroTol = 1e-6

class TestCase(unittest.TestCase):
  def setUp(self) :
    m = Mol(io_types.SMI,io_types.SMI)
    a = Atom()
    a.SetAtomicNum(6)
    for i in range(6): m.AddAtom(a)
    m.GetAtom(1).SetAtomicNum(8)

    for thing in [(1,2),(1,5),(1,6),(2,3),(3,4),(4,5)]:
      b = m.CreateBond()
      b.SetBegin(m.GetAtom(thing[0]))
      b.SetEnd(m.GetAtom(thing[1]))
      m.AddBond(b)
    self.mol1 = m

    # following 3 molecules are from Rucker paper, JCICS, 41, 2, 2001, 314-318
    self.mol2 = MolFromSmi("CCC(O)C(c1ccccc1)CC(C)N(C)C")

    self.mol3 = MolFromSmi("O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))")

    m = Mol(io_types.SMI,io_types.SMI)
    a = Atom()
    a.SetAtomicNum(6)
    for i in range(8): m.AddAtom(a)
    m.GetAtom(1).SetAtomicNum(7)
    
    for thing in [(1,2),(1,4),(1,6),(2,3),(2,7),(3,4),(3,8),(4,5),(5,6),(5,8),(6,7),(7,8)]:
      b = m.CreateBond()
      b.SetBegin(m.GetAtom(thing[0]))
      b.SetEnd(m.GetAtom(thing[1]))
      m.AddBond(b)
    self.mol4 = m
    
  def test1(self):
    
    m = self.mol1
    tgts = {1:6,2:7,3:8,4:9,5:6,6:1}
    for i in range(1,7):
      #print Test_frag(m)
      r = FindAllSubgraphs(m,i,1, [], [])
      #print '%d: %d walks'%(i,len(r))
      assert len(r) == tgts[i],'bad length'

  def test2(self) :
    m = self.mol1
    tgts = {3:6, 4:8, 5:6, 6:1}
    for i in range(3,7):
      r = FindAllSubgraphs(m,i,0,[[1], [2]], [], 1, 0)
      #print sort(pt)
      #print '%d: %d walks'%(i,len(r))
      assert len(r) == tgts[i],'bad length'

  def test3(self) :
    m = self.mol1
    tgts = {3:3, 4:4, 5:5, 6:1}
    for i in range(3,7):
      r = FindAllSubgraphs(m,i,1,[[1,3,3], [2]], [], 1, 0)
      #print sort(pt)
      #print '%d: %d walks'%(i,len(r))
      assert len(r) == tgts[i],'bad length'

  def test4(self) :
    m = self.mol1
    tgts = {4:4}
    for i in range(4,5) :
      r = FindUniqueSubgraphs(m, i, 0, 0)
      assert len(r) == tgts[i], 'wrong number of unique subgraphs'
      

  def test5(self) :
    at = 1990
    ut = 894
    nt = 0
    ns = 0
    mol = self.mol2
    for i in range(1, 18):
      rall = FindAllSubgraphs(mol, i, 0, [], [])
      nt += len(rall)
      rui = FindUniqueSubgraphs(mol, i, 0, 0)
      ns += len(rui)
      
    assert nt == at, "Wrong number of all subgraphs"
    assert ns == ut, "wrong number of unique subgraphs"

  def test6(self) :
    mol = self.mol3
    at = 6435
    ut = 5618
    nt = 0
    ns = 0
    for i in range(1, 26):
      rall = FindAllSubgraphs(mol, i, 0, [], [])
      nt += len(rall)
      rui = FindUniqueSubgraphs(mol, i, 0, 0)
      ns += len(rui)
      
    assert nt == at, "Wrong number of all subgraphs"
    assert ns == ut, "wrong number of unique subgraphs" 

  def test7(self) :
    m = self.mol4
    a = Atom()
    at = 2433
    ut = 300
    nt = 0;
    ns = 0;
    for i in range(1,13) :
      rall = FindAllSubgraphs(m, i, 0, [], [])
      nt += len(rall)
      rui = FindUniqueSubgraphs(m, i, 0, 0)
      ns += len(rui)
    assert nt == at, "Wrong number of all subgraphs"
    assert ns == ut, "wrong number of unique subgraphs" 
    
if __name__ == '__main__':
  from Chem import *
  unittest.main()
