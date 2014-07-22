# $Id$
#
#  Copyright (C) 2003-2013  Greg Landrum and Rational Discovery LLC
#         All Rights Reserved
#
""" This is a rough coverage test of the python wrapper

it's intended to be shallow, but broad

"""
from __future__ import print_function
import os,sys,tempfile
import unittest
from rdkit import RDConfig,rdBase
from rdkit import DataStructs
from rdkit import Chem
from rdkit import six

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

class TestCase(unittest.TestCase):
  def setUp(self):
    pass

  def test0Except(self):

    try:
      Chem.tossit()
    except IndexError:
      ok=1
    else:
      ok=0
    assert ok


  def test1Table(self):

    tbl = Chem.GetPeriodicTable()
    self.assertTrue(tbl)

    self.assertTrue(feq(tbl.GetAtomicWeight(6),12.011))
    self.assertTrue(feq(tbl.GetAtomicWeight("C"),12.011))
    self.assertTrue(tbl.GetAtomicNumber('C')==6)
    self.assertTrue(feq(tbl.GetRvdw(6),1.950))
    self.assertTrue(feq(tbl.GetRvdw("C"),1.950))
    self.assertTrue(feq(tbl.GetRcovalent(6),0.680))
    self.assertTrue(feq(tbl.GetRcovalent("C"),0.680))
    self.assertTrue(tbl.GetDefaultValence(6)==4)
    self.assertTrue(tbl.GetDefaultValence("C")==4)
    self.assertTrue(tuple(tbl.GetValenceList(6))==(4,))
    self.assertTrue(tuple(tbl.GetValenceList("C"))==(4,))
    self.assertTrue(tuple(tbl.GetValenceList(16))==(2,4,6))
    self.assertTrue(tuple(tbl.GetValenceList("S"))==(2,4,6))
    self.assertTrue(tbl.GetNOuterElecs(6)==4)
    self.assertTrue(tbl.GetNOuterElecs("C")==4)

  def test2Atom(self):
    atom = Chem.Atom(6)
    self.assertTrue(atom)
    self.assertTrue(atom.GetAtomicNum()==6)
    atom.SetAtomicNum(8)
    self.assertTrue(atom.GetAtomicNum()==8)

    atom = Chem.Atom("C")
    self.assertTrue(atom)
    self.assertTrue(atom.GetAtomicNum()==6)


  def test3Bond(self):
    # No longer relevant, bonds are not constructible from Python
    pass

  def test4Mol(self):
    mol = Chem.Mol()
    self.assertTrue(mol)

  def test5Smiles(self):
    mol = Chem.MolFromSmiles('n1ccccc1')
    self.assertTrue(mol)
    self.assertTrue(mol.GetNumAtoms()==6)
    self.assertTrue(mol.GetNumAtoms(1)==6)
    self.assertTrue(mol.GetNumAtoms(0)==11)
    at = mol.GetAtomWithIdx(2)
    self.assertTrue(at.GetAtomicNum()==6)
    at = mol.GetAtomWithIdx(0)
    self.assertTrue(at.GetAtomicNum()==7)
    
  def _test6Bookmarks(self):
    mol = Chem.MolFromSmiles('n1ccccc1')
    self.assertTrue(mol)

    self.assertTrue(not mol.HasAtomBookmark(0))
    mol.SetAtomBookmark(mol.GetAtomWithIdx(0),0)
    mol.SetAtomBookmark(mol.GetAtomWithIdx(1),1)
    self.assertTrue(mol.HasAtomBookmark(0))
    self.assertTrue(mol.HasAtomBookmark(1))

    if 1:
      self.assertTrue(not mol.HasBondBookmark(0))
      self.assertTrue(not mol.HasBondBookmark(1))
      mol.SetBondBookmark(mol.GetBondWithIdx(0),0)
      mol.SetBondBookmark(mol.GetBondWithIdx(1),1)
      self.assertTrue(mol.HasBondBookmark(0))
      self.assertTrue(mol.HasBondBookmark(1))


    at = mol.GetAtomWithBookmark(0)
    self.assertTrue(at)
    self.assertTrue(at.GetAtomicNum()==7)
    mol.ClearAtomBookmark(0)
    self.assertTrue(not mol.HasAtomBookmark(0))
    self.assertTrue(mol.HasAtomBookmark(1))
    mol.ClearAllAtomBookmarks()
    self.assertTrue(not mol.HasAtomBookmark(0))
    self.assertTrue(not mol.HasAtomBookmark(1))
    
    mol.SetAtomBookmark(mol.GetAtomWithIdx(1),1)

    if 1:
      self.assertTrue(mol.HasBondBookmark(0))
      self.assertTrue(mol.HasBondBookmark(1))
      bond = mol.GetBondWithBookmark(0)
      self.assertTrue(bond)
      mol.ClearBondBookmark(0)
      self.assertTrue(not mol.HasBondBookmark(0))
      self.assertTrue(mol.HasBondBookmark(1))
      mol.ClearAllBondBookmarks()
      self.assertTrue(not mol.HasBondBookmark(0))
      self.assertTrue(not mol.HasBondBookmark(1))

      self.assertTrue(mol.HasAtomBookmark(1))
    
  def test7Atom(self):
    mol = Chem.MolFromSmiles('n1ccccc1C[CH2-]')
    self.assertTrue(mol)
    Chem.SanitizeMol(mol)
    a0 = mol.GetAtomWithIdx(0)
    a1 = mol.GetAtomWithIdx(1)
    a6 = mol.GetAtomWithIdx(6)
    a7 = mol.GetAtomWithIdx(7)
    
    self.assertTrue(a0.GetAtomicNum()==7)
    self.assertTrue(a0.GetSymbol()=='N')
    self.assertTrue(a0.GetIdx()==0)

    aList = [a0,a1,a6,a7]
    self.assertTrue(a0.GetDegree()==2)
    self.assertTrue(a1.GetDegree()==2)
    self.assertTrue(a6.GetDegree()==2)
    self.assertTrue(a7.GetDegree()==1)
    self.assertTrue([x.GetDegree() for x in aList]==[2,2,2,1])

    self.assertTrue([x.GetTotalNumHs() for x in aList]==[0,1,2,2])
    self.assertTrue([x.GetNumImplicitHs() for x in aList]==[0,1,2,0])
    self.assertTrue([x.GetExplicitValence() for x in aList]==[3,3,2,3])
    self.assertTrue([x.GetImplicitValence() for x in aList]==[0,1,2,0])
    self.assertTrue([x.GetFormalCharge() for x in aList]==[0,0,0,-1])
    self.assertTrue([x.GetNoImplicit() for x in aList]==[0,0,0,1])
    self.assertTrue([x.GetNumExplicitHs() for x in aList]==[0,0,0,2])
    self.assertTrue([x.GetIsAromatic() for x in aList]==[1,1,0,0])
    self.assertTrue([x.GetHybridization() for x in aList]==[Chem.HybridizationType.SP2,Chem.HybridizationType.SP2,
                                                   Chem.HybridizationType.SP3,Chem.HybridizationType.SP3],\
                                                   [x.GetHybridization() for x in aList])
    


  def test8Bond(self):
    mol = Chem.MolFromSmiles('n1ccccc1CC(=O)O')
    self.assertTrue(mol)
    Chem.SanitizeMol(mol)
    # note bond numbering is funny because of ring closure
    b0 = mol.GetBondWithIdx(0)
    b6 = mol.GetBondWithIdx(6)
    b7 = mol.GetBondWithIdx(7)
    b8 = mol.GetBondWithIdx(8)

    bList = [b0,b6,b7,b8]
    self.assertTrue([x.GetBondType() for x in bList] ==
           [Chem.BondType.AROMATIC,Chem.BondType.SINGLE,
            Chem.BondType.DOUBLE,Chem.BondType.SINGLE])
    self.assertTrue([x.GetIsAromatic() for x in bList] == 
           [1,0,0,0])
    self.assertEqual(bList[0].GetBondTypeAsDouble(),1.5)
    self.assertEqual(bList[1].GetBondTypeAsDouble(),1.0)
    self.assertEqual(bList[2].GetBondTypeAsDouble(),2.0)
    
    self.assertTrue([x.GetIsConjugated()!=0 for x in bList] ==
           [1,0,1,1],[x.GetIsConjugated()!=0 for x in bList])
    self.assertTrue([x.GetBeginAtomIdx() for x in bList] ==
           [0,6,7,7],[x.GetBeginAtomIdx() for x in bList])
    self.assertTrue([x.GetBeginAtom().GetIdx() for x in bList] ==
           [0,6,7,7])
    self.assertTrue([x.GetEndAtomIdx() for x in bList] ==
           [1,7,8,9])
    self.assertTrue([x.GetEndAtom().GetIdx() for x in bList] ==
           [1,7,8,9])


  def test9Smarts(self):
    query1 = Chem.MolFromSmarts('C(=O)O')
    self.assertTrue(query1)
    query2 = Chem.MolFromSmarts('C(=O)[O,N]')
    self.assertTrue(query2)
    query3 = Chem.MolFromSmarts('[$(C(=O)O)]')
    self.assertTrue(query3)


    mol = Chem.MolFromSmiles('CCC(=O)O')
    self.assertTrue(mol)

    self.assertTrue(mol.HasSubstructMatch(query1))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.HasSubstructMatch(query3))

    mol = Chem.MolFromSmiles('CCC(=O)N')
    self.assertTrue(mol)

    self.assertTrue(not mol.HasSubstructMatch(query1))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(not mol.HasSubstructMatch(query3))

  def test10Iterators(self):
    mol = Chem.MolFromSmiles('CCOC')
    self.assertTrue(mol)

    for atom in mol.GetAtoms():
      self.assertTrue(atom)
    ats = mol.GetAtoms()
    try:
      ats[1]
    except:
      ok = 0
    else:
      ok = 1
    self.assertTrue(ok)
    try:
      ats[12]
    except IndexError:
      ok = 1
    else:
      ok = 0
    self.assertTrue(ok)

    if 0:
      for atom in mol.GetHeteros():
        self.assertTrue(atom)
      ats = mol.GetHeteros()
      try:
        ats[0]
      except:
        ok = 0
      else:
        ok = 1
      self.assertTrue(ok)
      self.assertTrue(ats[0].GetIdx()==2)
      try:
        ats[12]
      except IndexError:
        ok = 1
      else:
        ok = 0
      self.assertTrue(ok)


    for bond in mol.GetBonds():
      self.assertTrue(bond)
    bonds = mol.GetBonds()
    try:
      bonds[1]
    except:
      ok = 0
    else:
      ok = 1
    self.assertTrue(ok)
    try:
      bonds[12]
    except IndexError:
      ok = 1
    else:
      ok = 0
    self.assertTrue(ok)
      
    if 0:
      mol = Chem.MolFromSmiles('c1ccccc1C')
      for atom in mol.GetAromaticAtoms():
        self.assertTrue(atom)
      ats = mol.GetAromaticAtoms()
      try:
        ats[0]
      except:
        ok = 0
      else:
        ok = 1
      self.assertTrue(ok)
      self.assertTrue(ats[0].GetIdx()==0)
      self.assertTrue(ats[1].GetIdx()==1)
      self.assertTrue(ats[2].GetIdx()==2)
      try:
        ats[12]
      except IndexError:
        ok = 1
      else:
        ok = 0
      self.assertTrue(ok)




  def test11MolOps(self) :
    mol = Chem.MolFromSmiles('C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3')
    self.assertTrue(mol)
    smi = Chem.MolToSmiles(mol)
    Chem.SanitizeMol(mol)
    nr = Chem.GetSymmSSSR(mol)

    self.assertTrue((len(nr) == 3))
    
  def test12Smarts(self):
    query1 = Chem.MolFromSmarts('C(=O)O')
    self.assertTrue(query1)
    query2 = Chem.MolFromSmarts('C(=O)[O,N]')
    self.assertTrue(query2)
    query3 = Chem.MolFromSmarts('[$(C(=O)O)]')
    self.assertTrue(query3)


    mol = Chem.MolFromSmiles('CCC(=O)O')
    self.assertTrue(mol)

    self.assertTrue(mol.HasSubstructMatch(query1))
    self.assertTrue(mol.GetSubstructMatch(query1)==(2,3,4))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.GetSubstructMatch(query2)==(2,3,4))
    self.assertTrue(mol.HasSubstructMatch(query3))
    self.assertTrue(mol.GetSubstructMatch(query3)==(2,))

    mol = Chem.MolFromSmiles('CCC(=O)N')
    self.assertTrue(mol)

    self.assertTrue(not mol.HasSubstructMatch(query1))
    self.assertTrue(not mol.GetSubstructMatch(query1))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.GetSubstructMatch(query2)==(2,3,4))
    self.assertTrue(not mol.HasSubstructMatch(query3))

    mol = Chem.MolFromSmiles('OC(=O)CC(=O)O')
    self.assertTrue(mol)
    self.assertTrue(mol.HasSubstructMatch(query1))
    self.assertTrue(mol.GetSubstructMatch(query1)==(1,2,0))
    self.assertTrue(mol.GetSubstructMatches(query1)==((1,2,0),(4,5,6)))
    self.assertTrue(mol.HasSubstructMatch(query2))
    self.assertTrue(mol.GetSubstructMatch(query2)==(1,2,0))
    self.assertTrue(mol.GetSubstructMatches(query2)==((1,2,0),(4,5,6)))
    self.assertTrue(mol.HasSubstructMatch(query3))
    self.assertTrue(mol.GetSubstructMatches(query3)==((1,),(4,)))
    
  def test13Smarts(self):
    # previous smarts problems:
    query = Chem.MolFromSmarts('N(=,-C)')
    self.assertTrue(query)
    mol = Chem.MolFromSmiles('N#C')
    self.assertTrue(not mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('N=C')
    self.assertTrue(mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('NC')
    self.assertTrue(mol.HasSubstructMatch(query))

    query = Chem.MolFromSmarts('[Cl,$(O)]')
    mol = Chem.MolFromSmiles('C(=O)O')
    self.assertTrue(len(mol.GetSubstructMatches(query))==2)
    mol = Chem.MolFromSmiles('C(=N)N')
    self.assertTrue(len(mol.GetSubstructMatches(query))==0)
    
    query = Chem.MolFromSmarts('[$([O,S]-[!$(*=O)])]')
    mol = Chem.MolFromSmiles('CC(S)C(=O)O')
    self.assertTrue(len(mol.GetSubstructMatches(query))==1)
    mol = Chem.MolFromSmiles('C(=O)O')
    self.assertTrue(len(mol.GetSubstructMatches(query))==0)
    
  def test14Hs(self):
    m = Chem.MolFromSmiles('CC(=O)[OH]')
    self.assertEqual(m.GetNumAtoms(),4)
    m2 = Chem.AddHs(m)
    self.assertEqual(m2.GetNumAtoms(),8)
    m2 = Chem.RemoveHs(m2)
    self.assertEqual(m2.GetNumAtoms(),4)

    m = Chem.MolFromSmiles('CC[H]',False)
    self.assertEqual(m.GetNumAtoms(),3)
    m2 = Chem.MergeQueryHs(m)
    self.assertEqual(m2.GetNumAtoms(),2)
    self.assertTrue(m2.GetAtomWithIdx(1).HasQuery())
    
  def test15Neighbors(self):
    m = Chem.MolFromSmiles('CC(=O)[OH]')
    self.assertTrue(m.GetNumAtoms()==4)
    
    a = m.GetAtomWithIdx(1)
    ns = a.GetNeighbors()
    self.assertTrue(len(ns)==3)

    bs = a.GetBonds()
    self.assertTrue(len(bs)==3)

    for b in bs:
      try:
        a2 = b.GetOtherAtom(a)
      except:
        a2=None
      self.assertTrue(a2)
    self.assertTrue(len(bs)==3)

  def test16Pickle(self):
    from rdkit.six.moves import cPickle
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    pkl = cPickle.dumps(m)
    m2 = cPickle.loads(pkl)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1==smi2)
    
  def test16Props(self):
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    self.assertTrue(not m.HasProp('prop1'))
    self.assertTrue(not m.HasProp('prop2'))
    self.assertTrue(not m.HasProp('prop2'))
    m.SetProp('prop1','foob')
    self.assertTrue(not m.HasProp('prop2'))
    self.assertTrue(m.HasProp('prop1'))
    self.assertTrue(m.GetProp('prop1')=='foob')
    self.assertTrue(not m.HasProp('propo'))
    try:
      m.GetProp('prop2')
    except KeyError:
      ok=1
    else:
      ok=0
    self.assertTrue(ok)

    # test computed properties
    m.SetProp('cprop1', 'foo', 1)
    m.SetProp('cprop2', 'foo2', 1)

    m.ClearComputedProps()
    self.assertTrue(not m.HasProp('cprop1'))
    self.assertTrue(not m.HasProp('cprop2'))

  def test17Kekulize(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi=='c1ccccc1')

    Chem.Kekulize(m)
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi=='c1ccccc1')

    m = Chem.MolFromSmiles('c1ccccc1')
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi=='c1ccccc1')

    Chem.Kekulize(m,1)
    smi = Chem.MolToSmiles(m)
    self.assertTrue(smi=='C1=CC=CC=C1', smi)

  def test18Paths(self):


    m = Chem.MolFromSmiles("C1CC2C1CC2")
    #self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7)
    #print(Chem.FindAllPathsOfLengthN(m,3,useBonds=0))
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==10,
           Chem.FindAllPathsOfLengthN(m,2,useBonds=1))
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==14)
    

    m = Chem.MolFromSmiles('C1CC1C')
    self.assertTrue(m)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==4)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3,Chem.FindAllPathsOfLengthN(m,3,useBonds=1))
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,4,useBonds=1))==1,Chem.FindAllPathsOfLengthN(m,4,useBonds=1))
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,5,useBonds=1))==0,Chem.FindAllPathsOfLengthN(m,5,useBonds=1))

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==4)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3)
    
    m = Chem.MolFromSmiles("CCC(C)CC")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==4,Chem.FindAllPathsOfLengthN(m,3,useBonds=1))
    
    m = Chem.MolFromSmiles("CCCC(C)C")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3)
    
    m = Chem.MolFromSmiles("CC(C)C(C)C")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==6)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==4)
    
    m = Chem.MolFromSmiles("CC(C)(C)CC")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==7)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3,Chem.FindAllPathsOfLengthN(m,3,useBonds=1))
    
    m = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==6)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==6)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==6)
    
    m = Chem.MolFromSmiles("C1CC2C1CC2")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==10,
           Chem.FindAllPathsOfLengthN(m,2,useBonds=1))
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==14)
    
    
    m = Chem.MolFromSmiles("CC2C1CCC12")
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7)
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==11)
    # FIX: this result disagrees with the paper (which says 13),
    #   but it seems right
    self.assertTrue(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==15,
           Chem.FindAllPathsOfLengthN(m,3,useBonds=1))
    
    
    
  def test19Subgraphs(self):
    m = Chem.MolFromSmiles('C1CC1C')
    self.assertTrue(m)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1,0))==4)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==4)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,4))==1)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,5))==0)

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==4)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==3)
    
    l = Chem.FindAllSubgraphsOfLengthMToN(m,1,3)
    self.assertEqual(len(l),3)
    self.assertEqual(len(l[0]),5)
    self.assertEqual(len(l[1]),4)
    self.assertEqual(len(l[2]),3)
    self.assertRaises(ValueError,lambda :Chem.FindAllSubgraphsOfLengthMToN(m,4,3))

    
    m = Chem.MolFromSmiles("CCC(C)CC")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==5)
    
    m = Chem.MolFromSmiles("CCCC(C)C")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==4)
    
    m = Chem.MolFromSmiles("CC(C)C(C)C")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==6)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==6)
    
    m = Chem.MolFromSmiles("CC(C)(C)CC")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==7)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==7,Chem.FindAllSubgraphsOfLengthN(m,3))
    
    m = Chem.MolFromSmiles("C1CCCCC1")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==6)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==6)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==6)
    #self.assertTrue(len(Chem.FindUniqueSubgraphsOfLengthN(m,1))==1)
    self.assertTrue(len(Chem.FindUniqueSubgraphsOfLengthN(m,2))==1)
    self.assertTrue(len(Chem.FindUniqueSubgraphsOfLengthN(m,3))==1)
    
    m = Chem.MolFromSmiles("C1CC2C1CC2")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==7)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==10)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==16)
    

    m = Chem.MolFromSmiles("CC2C1CCC12")
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,1))==7)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,2))==11)
    self.assertTrue(len(Chem.FindAllSubgraphsOfLengthN(m,3))==18,
           len(Chem.FindAllSubgraphsOfLengthN(m,3)))
    
    
  def test20IsInRing(self):
    m = Chem.MolFromSmiles('C1CCC1C')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).IsInRingSize(4))
    self.assertTrue(m.GetAtomWithIdx(1).IsInRingSize(4))
    self.assertTrue(m.GetAtomWithIdx(2).IsInRingSize(4))
    self.assertTrue(m.GetAtomWithIdx(3).IsInRingSize(4))
    self.assertTrue(not m.GetAtomWithIdx(4).IsInRingSize(4))
    
    self.assertTrue(not m.GetAtomWithIdx(0).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(1).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(2).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(3).IsInRingSize(3))
    self.assertTrue(not m.GetAtomWithIdx(4).IsInRingSize(3))
    
    self.assertTrue(m.GetBondWithIdx(0).IsInRingSize(4))
    self.assertTrue(not m.GetBondWithIdx(3).IsInRingSize(4))
    self.assertTrue(not m.GetBondWithIdx(0).IsInRingSize(3))
    self.assertTrue(not m.GetBondWithIdx(3).IsInRingSize(3))
    
  def test21Robustification(self):
    ok = False
    # FIX: at the moment I can't figure out how to catch the
    # actual exception that BPL is throwinng when it gets
    # invalid arguments (Boost.Python.ArgumentError)
    try:
      Chem.MolFromSmiles('C=O').HasSubstructMatch(Chem.MolFromSmarts('fiib'))
    #except ValueError:
    #  ok=True
    except:
      ok=True
    self.assertTrue(ok  )

  def test22DeleteSubstruct(self) :
    query = Chem.MolFromSmarts('C(=O)O')
    mol = Chem.MolFromSmiles('CCC(=O)O')
    nmol = Chem.DeleteSubstructs(mol, query)
    
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CC')

    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    # now delete only fragments
    nmol = Chem.DeleteSubstructs(mol, query, 1)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCC(=O)O',Chem.MolToSmiles(nmol))
    
    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = Chem.DeleteSubstructs(mol, query, 0)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CC')
    
    mol = Chem.MolFromSmiles('CCCO')
    nmol = Chem.DeleteSubstructs(mol, query, 0)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCCO')

    # Issue 96 prevented this from working:
    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = Chem.DeleteSubstructs(mol, query, 1)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCC(=O)O')
    nmol = Chem.DeleteSubstructs(nmol, query, 1)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CCC(=O)O')
    nmol = Chem.DeleteSubstructs(nmol, query, 0)
    self.assertTrue(Chem.MolToSmiles(nmol) == 'CC')

    
  def test23MolFileParsing(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','triazine.mol')
    #fileN = "../FileParsers/test_data/triazine.mol"
    with open(fileN,'r') as inF:
      inD = inF.read()
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==9)

    m1 = Chem.MolFromMolFile(fileN)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==9)

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','triazine.mof')
    self.assertRaises(IOError,lambda :Chem.MolFromMolFile(fileN))

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','list-query.mol')
    query = Chem.MolFromMolFile(fileN)
    smi = Chem.MolToSmiles(query)
    self.assertEqual(smi,'c1ccccc1')
    smi = Chem.MolToSmarts(query)
    self.assertEqual(smi,'[#6]1:[#6]:[#6]:[#6]:[#6]:[#6,#7,#15]:1',smi)

    query = Chem.MolFromMolFile(fileN,sanitize=False)
    smi = Chem.MolToSmiles(query)
    self.assertEqual(smi,'C1=CC=CC=C1')
    query.UpdatePropertyCache()
    smi = Chem.MolToSmarts(query)
    self.assertEqual(smi,'[#6]1=[#6]-[#6]=[#6]-[#6]=[#6,#7,#15]-1')
    smi = "C1=CC=CC=C1"
    mol = Chem.MolFromSmiles(smi,0)
    self.assertTrue(mol.HasSubstructMatch(query))
    Chem.SanitizeMol(mol)
    self.assertTrue(not mol.HasSubstructMatch(query))

    mol = Chem.MolFromSmiles('N1=CC=CC=C1',0)
    self.assertTrue(mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('S1=CC=CC=C1',0)
    self.assertTrue(not mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('P1=CC=CC=C1',0)
    self.assertTrue(mol.HasSubstructMatch(query))


    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','issue123.mol')
    mol = Chem.MolFromMolFile(fileN)
    self.assertTrue(mol)
    self.assertEqual(mol.GetNumAtoms(),23)
    mol = Chem.MolFromMolFile(fileN,removeHs=False)
    self.assertTrue(mol)
    self.assertEqual(mol.GetNumAtoms(),39)

  # test23 was for Chem.DaylightFingerprint, which is deprecated
    
  def test24RDKFingerprint(self):
    from rdkit import DataStructs
    m1 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp1 = Chem.RDKFingerprint(m1)
    self.assertTrue(len(fp1)==2048)
    m2 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp2 = Chem.RDKFingerprint(m2)

    tmp = DataStructs.TanimotoSimilarity(fp1,fp2)
    self.assertTrue(tmp==1.0,tmp)

    m2 = Chem.MolFromSmiles('C1=CC=CC=N1')
    fp2 = Chem.RDKFingerprint(m2)
    self.assertTrue(len(fp2)==2048)
    tmp = DataStructs.TanimotoSimilarity(fp1,fp2)
    self.assertTrue(tmp<1.0,tmp)
    self.assertTrue(tmp>0.0,tmp)

    fp3 = Chem.RDKFingerprint(m1,tgtDensity=0.3)
    self.assertTrue(len(fp3)<2048)
    
    m1 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp1 = Chem.RDKFingerprint(m1)
    m2 = Chem.MolFromSmiles('C1=CC=CC=N1')
    fp2 = Chem.RDKFingerprint(m2)
    self.assertNotEqual(fp1,fp2)

    atomInvariants=[1]*6
    fp1 = Chem.RDKFingerprint(m1,atomInvariants=atomInvariants)
    fp2 = Chem.RDKFingerprint(m2,atomInvariants=atomInvariants)
    self.assertEqual(fp1,fp2)    

    m2 = Chem.MolFromSmiles('C1CCCCN1')
    fp1 = Chem.RDKFingerprint(m1,atomInvariants=atomInvariants,useBondOrder=False)
    fp2 = Chem.RDKFingerprint(m2,atomInvariants=atomInvariants,useBondOrder=False)
    self.assertEqual(fp1,fp2)    

    # rooted at atom
    m1 = Chem.MolFromSmiles('CCCCCO')
    fp1 = Chem.RDKFingerprint(m1,1,4,nBitsPerHash=1,fromAtoms=[0])
    self.assertEqual(fp1.GetNumOnBits(),4)    
    m1 = Chem.MolFromSmiles('CCCCCO')
    fp1 = Chem.RDKFingerprint(m1,1,4,nBitsPerHash=1,fromAtoms=[0,5])
    self.assertEqual(fp1.GetNumOnBits(),8)

    # test sf.net issue 270:
    fp1 = Chem.RDKFingerprint(m1,atomInvariants=[x.GetAtomicNum()+10 for x in m1.GetAtoms()])

    # atomBits
    m1 = Chem.MolFromSmiles('CCCO')
    l=[]
    fp1 = Chem.RDKFingerprint(m1,minPath=1,maxPath=2,nBitsPerHash=1,atomBits=l)
    self.assertEqual(fp1.GetNumOnBits(),4)
    self.assertEqual(len(l),m1.GetNumAtoms())
    self.assertEqual(len(l[0]),2)
    self.assertEqual(len(l[1]),3)
    self.assertEqual(len(l[2]),4)
    self.assertEqual(len(l[3]),2)
    
  def test25SDMolSupplier(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf')
    #fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
    sdSup = Chem.SDMolSupplier(fileN)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
            "192", "203", "210", "211", "213", "220", "229", "256"]

    chgs192 = {8:1, 11:1, 15:-1, 18:-1, 20:1, 21:1, 23:-1, 25:-1} 
    i = 0
    for mol in sdSup :
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
      if (mol.GetProp("_Name") == "192") :
        # test parsed charges on one of the molecules
        for id in chgs192.keys() :
          self.assertTrue(mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id])
    self.assertRaises(StopIteration,lambda:six.next(sdSup))
    sdSup.reset()
    
    ns = [mol.GetProp("_Name") for mol in sdSup]
    self.assertTrue(ns == molNames)

    sdSup = Chem.SDMolSupplier(fileN, 0)
    for mol in sdSup :
      self.assertTrue(not mol.HasProp("numArom"))

    sdSup = Chem.SDMolSupplier(fileN)
    self.assertTrue(len(sdSup) == 16)
    mol = sdSup[5]
    self.assertTrue(mol.GetProp("_Name") == "170")


    # test handling of H removal:
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','withHs.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    m = six.next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms()==23)
    m = six.next(sdSup)
    self.assertTrue(m)
    print(m.GetNumAtoms())
    self.assertTrue(m.GetNumAtoms()==28)

    sdSup = Chem.SDMolSupplier(fileN,removeHs=False)
    m = six.next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms()==39)
    m = six.next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms()==30)

    with open(fileN,'r') as dFile:
      d = dFile.read()
    sdSup.SetData(d)
    m = six.next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms()==23)
    m = six.next(sdSup)
    self.assertTrue(m)
    print(m.GetNumAtoms())
    self.assertTrue(m.GetNumAtoms()==28)

    sdSup.SetData(d,removeHs=False)
    m = six.next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms()==39)
    m = six.next(sdSup)
    self.assertTrue(m)
    self.assertTrue(m.GetNumAtoms()==30)

    
  def test26SmiMolSupplier(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','first_200.tpsa.csv')
    #fileN = "../FileParsers/test_data/first_200.tpsa.csv"
    smiSup = Chem.SmilesMolSupplier(fileN, ",", 0, -1)
    mol = smiSup[16];
    self.assertTrue(mol.GetProp("TPSA") == "46.25")

    mol = smiSup[8];
    self.assertTrue(mol.GetProp("TPSA") == "65.18")

    self.assertTrue(len(smiSup) == 200)

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','fewSmi.csv')
    #fileN = "../FileParsers/test_data/fewSmi.csv"
    smiSup = Chem.SmilesMolSupplier(fileN, delimiter=",",
                                      smilesColumn=1, nameColumn=0,
                                      titleLine=0)
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    i = 0
    for mol in smiSup:
      self.assertTrue(mol.GetProp("_Name") == names[i])
      i += 1
      
    mol = smiSup[3]
    
    self.assertTrue(mol.GetProp("_Name") == "4")
    self.assertTrue(mol.GetProp("Column_2") == "82.78")

    # and test doing a supplier from text:
    with open(fileN,'r') as inF:
      inD = inF.read()
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=1, nameColumn=0,
                   titleLine=0)
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    i = 0
    # iteration interface:
    for mol in smiSup:
      self.assertTrue(mol.GetProp("_Name") == names[i])
      i += 1
    self.assertTrue(i==10)
    # random access:
    mol = smiSup[3]
    self.assertTrue(len(smiSup)==10)
    self.assertTrue(mol.GetProp("_Name") == "4")
    self.assertTrue(mol.GetProp("Column_2") == "82.78")

    # issue 113:
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=1, nameColumn=0,
                   titleLine=0)
    self.assertTrue(len(smiSup)==10)

    # and test failure handling:
    inD = """mol-1,CCC
mol-2,CCCC
mol-3,fail
mol-4,CCOC
    """
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=1, nameColumn=0,
                   titleLine=0)
    # there are 4 entries in the supplier:
    self.assertTrue(len(smiSup)==4)
    # but the 3rd is a None:
    self.assertTrue(smiSup[2] is None)


    text="Id SMILES Column_2\n"+\
    "mol-1 C 1.0\n"+\
    "mol-2 CC 4.0\n"+\
    "mol-4 CCCC 16.0"
    smiSup.SetData(text, delimiter=" ",
                   smilesColumn=1, nameColumn=0,
                   titleLine=1)
    self.assertTrue(len(smiSup)==3)
    self.assertTrue(smiSup[0])
    self.assertTrue(smiSup[1])
    self.assertTrue(smiSup[2])
    m = [x for x in smiSup]
    self.assertTrue(smiSup[2])
    self.assertTrue(len(m)==3)
    self.assertTrue(m[0].GetProp("Column_2")=="1.0")
    
    # test simple parsing and Issue 114:
    smis = ['CC','CCC','CCOC','CCCOCC','CCCOCCC']
    inD = '\n'.join(smis)
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    self.assertTrue(len(smiSup)==5)
    m = [x for x in smiSup]
    self.assertTrue(smiSup[4])
    self.assertTrue(len(m)==5)
    
    # order dependance:
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    self.assertTrue(smiSup[4])
    self.assertTrue(len(smiSup)==5)

    # this was a nasty BC:
    # asking for a particular entry with a higher index than what we've
    # already seen resulted in a duplicate:
    smis = ['CC','CCC','CCOC','CCCCOC']
    inD = '\n'.join(smis)
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    m = six.next(smiSup)
    m = smiSup[3]
    self.assertTrue(len(smiSup)==4)

    try:
      smiSup[4]
    except:
      ok=1
    else:
      ok=0
    self.assertTrue(ok)

    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    try:
      smiSup[4]
    except:
      ok=1
    else:
      ok=0
    self.assertTrue(ok)
    sys.stderr.write('>>> This may result in an infinite loop.  It should finish almost instantly\n')
    self.assertTrue(len(smiSup)==4)
    sys.stderr.write('<<< OK, it finished.\n')


  def test27SmilesWriter(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','fewSmi.csv')
    #fileN = "../FileParsers/test_data/fewSmi.csv"

    smiSup = Chem.SmilesMolSupplier(fileN, delimiter=",",
                                      smilesColumn=1, nameColumn=0,
                                      titleLine=0)
    propNames = []
    propNames.append("Column_2")
    ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap','test_data','outSmiles.txt')
    writer = Chem.SmilesWriter(ofile)
    writer.SetProps(propNames)
    for mol in smiSup :
      writer.write(mol)
    writer.flush()

  def test28SmilesReverse(self):
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    props = ["34.14","25.78","106.51","82.78","60.16","87.74","37.38","77.28","65.18","0.00"]
    ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap','test_data','outSmiles.txt')
    #ofile = "test_data/outSmiles.csv"
    smiSup = Chem.SmilesMolSupplier(ofile)
    i = 0
    for mol in smiSup:
      #print([repr(x) for x in mol.GetPropNames()])
      self.assertTrue(mol.GetProp("_Name") == names[i])
      self.assertTrue(mol.GetProp("Column_2") == props[i])
      i += 1
      
  def writerSDFile(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf')
    #fileN = "../FileParsers/test_data/NCI_aids_few.sdf"
    ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap','test_data','outNCI_few.sdf');
    writer = Chem.SDWriter(ofile);
    sdSup = Chem.SDMolSupplier(fileN)
    for mol in sdSup :
      writer.write(mol)
    writer.flush()

  def test29SDWriterLoop(self) :
    self.writerSDFile()
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap',
                                            'test_data','outNCI_few.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
            "192", "203", "210", "211", "213", "220", "229", "256"]
    chgs192 = {8:1, 11:1, 15:-1, 18:-1, 20:1, 21:1, 23:-1, 25:-1}
    i = 0
    
    for mol in sdSup :
      #print('mol:',mol)
      #print('\t',molNames[i])
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
      if (mol.GetProp("_Name") == "192") :
        # test parsed charges on one of the molecules
        for id in chgs192.keys() :
          self.assertTrue(mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id])

  def test30Issues109and110(self) :
    """ issues 110 and 109 were both related to handling of explicit Hs in
       SMILES input.
       
    """ 
    m1 = Chem.MolFromSmiles('N12[CH](SC(C)(C)[CH]1C(O)=O)[CH](C2=O)NC(=O)[CH](N)c3ccccc3')
    self.assertTrue(m1.GetNumAtoms()==24)
    m2 = Chem.MolFromSmiles('C1C=C([CH](N)C(=O)N[C]2([H])[C]3([H])SC(C)(C)[CH](C(=O)O)N3C(=O)2)C=CC=1')
    self.assertTrue(m2.GetNumAtoms()==24)

    smi1 = Chem.MolToSmiles(m1)
    smi2 = Chem.MolToSmiles(m2)
    self.assertTrue(smi1==smi2)

    m1 = Chem.MolFromSmiles('[H]CCl')
    self.assertTrue(m1.GetNumAtoms()==2)
    self.assertTrue(m1.GetAtomWithIdx(0).GetNumExplicitHs()==0)
    m1 = Chem.MolFromSmiles('[H][CH2]Cl')
    self.assertTrue(m1.GetNumAtoms()==2)
    self.assertTrue(m1.GetAtomWithIdx(0).GetNumExplicitHs()==3)
    m2 = Chem.AddHs(m1)
    self.assertTrue(m2.GetNumAtoms()==5)
    m2 = Chem.RemoveHs(m2)
    self.assertTrue(m2.GetNumAtoms()==2)
    
  def test31ChiralitySmiles(self) :
    m1 = Chem.MolFromSmiles('F[C@](Br)(I)Cl')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==5)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[C@](Cl)(Br)I',Chem.MolToSmiles(m1,1))

    m1 = Chem.MolFromSmiles('CC1C[C@@]1(Cl)F')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==6)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='CC1C[C@]1(F)Cl',Chem.MolToSmiles(m1,1))

    m1 = Chem.MolFromSmiles('CC1C[C@]1(Cl)F')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==6)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='CC1C[C@@]1(F)Cl',Chem.MolToSmiles(m1,1))


  def test31aChiralitySubstructs(self) :
    m1 = Chem.MolFromSmiles('CC1C[C@@]1(Cl)F')
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==6)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='CC1C[C@]1(F)Cl',Chem.MolToSmiles(m1,1))

    m2 = Chem.MolFromSmiles('CC1C[C@]1(Cl)F')
    self.assertTrue(m2 is not None)
    self.assertTrue(m2.GetNumAtoms()==6)
    self.assertTrue(Chem.MolToSmiles(m2,1)=='CC1C[C@@]1(F)Cl',Chem.MolToSmiles(m2,1))

    self.assertTrue(m1.HasSubstructMatch(m1))
    self.assertTrue(m1.HasSubstructMatch(m2))
    self.assertTrue(m1.HasSubstructMatch(m1,useChirality=True))
    self.assertTrue(not m1.HasSubstructMatch(m2,useChirality=True))

    
    

  def _test32MolFilesWithChirality(self) :
    inD = """chiral1.mol
  ChemDraw10160313232D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  4  1  1      
  2  5  1  0      
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==5)
    self.assertTrue(smi=='F[C@](Cl)(Br)I',smi)

    inD = """chiral2.cdxml
  ChemDraw10160314052D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  4  1  6      
  2  5  1  0      
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==5)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[C@@](Cl)(Br)I')

    inD = """chiral1.mol
  ChemDraw10160313232D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  4  1  1      
  2  5  1  0      
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==5)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[C@](Cl)(Br)I')

    inD = """chiral1.mol
  ChemDraw10160313232D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0553   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7697   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
   -0.6592   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.7697   -0.6188    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0
    0.0553    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  1  3  1  1      
  1  4  1  0     
  1  5  1  0      
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==5)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[C@](Cl)(Br)I')

    inD = """chiral3.mol
  ChemDraw10160314362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.4125    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.2062    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  1      
  2  4  1  0      
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==4)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[C@H](Cl)Br')

    inD = """chiral4.mol
  ChemDraw10160314362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.4125    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.2062    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  1      
  2  4  1  0      
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==4)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='FN(Cl)Br')

    inD = """chiral5.mol
  ChemDraw10160314362D

  4  3  0  0  0  0  0  0  0  0999 V2000
    0.4125    0.6188    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.4125   -0.2062    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.3020   -0.6188    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
   -0.4125   -0.2062    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  1      
  2  4  1  0      
M  CHG  1   2   1
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==4)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[N@H+](Cl)Br')

    inD="""Case 10-14-3
  ChemDraw10140308512D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.8250   -0.4125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250   -0.4125    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.4125    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  2  4  1  1      
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==4)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[C@H](Cl)Br')

    inD="""Case 10-14-4
  ChemDraw10140308512D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -0.8250   -0.4125    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.4125    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8250   -0.4125    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.4125    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  1      
  2  4  1  0      
M  END

"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==4)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='F[C@H](Cl)Br')

    inD="""chiral4.mol
  ChemDraw10160315172D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.4422    0.1402    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4422   -0.6848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2723   -0.2723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8547    0.8547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6848    0.4422    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8547   -0.8547    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  1  1  0      
  1  4  1  0      
  3  5  1  1      
  3  6  1  0      
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==6)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='CC1C[C@@]1(F)Cl',Chem.MolToSmiles(m1,1))

    inD="""chiral4.mol
  ChemDraw10160315172D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.4422    0.1402    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4422   -0.6848    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.2723   -0.2723    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8547    0.8547    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6848    0.4422    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
    0.8547   -0.8547    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0      
  2  3  1  0      
  3  1  1  0      
  1  4  1  0      
  3  5  1  6      
  3  6  1  0      
M  END
"""
    m1 = Chem.MolFromMolBlock(inD)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==6)
    self.assertTrue(Chem.MolToSmiles(m1,1)=='CC1C[C@]1(F)Cl',Chem.MolToSmiles(m1,1))

  def test33Issue65(self) :
    """ issue 65 relates to handling of [H] in SMARTS
       
    """ 
    m1 = Chem.MolFromSmiles('OC(O)(O)O')
    m2 = Chem.MolFromSmiles('OC(O)O')
    m3 = Chem.MolFromSmiles('OCO')
    q1 = Chem.MolFromSmarts('OC[H]',1)
    q2 = Chem.MolFromSmarts('O[C;H1]',1)
    q3 = Chem.MolFromSmarts('O[C;H1][H]',1)

    self.assertTrue(not m1.HasSubstructMatch(q1))
    self.assertTrue(not m1.HasSubstructMatch(q2))
    self.assertTrue(not m1.HasSubstructMatch(q3))

    self.assertTrue(m2.HasSubstructMatch(q1))
    self.assertTrue(m2.HasSubstructMatch(q2))
    self.assertTrue(m2.HasSubstructMatch(q3))

    self.assertTrue(m3.HasSubstructMatch(q1))
    self.assertTrue(not m3.HasSubstructMatch(q2))
    self.assertTrue(not m3.HasSubstructMatch(q3))

    m1H = Chem.AddHs(m1)
    m2H = Chem.AddHs(m2)
    m3H = Chem.AddHs(m3)
    q1 = Chem.MolFromSmarts('OC[H]')
    q2 = Chem.MolFromSmarts('O[C;H1]')
    q3 = Chem.MolFromSmarts('O[C;H1][H]')

    self.assertTrue(not m1H.HasSubstructMatch(q1))
    self.assertTrue(not m1H.HasSubstructMatch(q2))
    self.assertTrue(not m1H.HasSubstructMatch(q3))

    #m2H.Debug()
    self.assertTrue(m2H.HasSubstructMatch(q1))
    self.assertTrue(m2H.HasSubstructMatch(q2))
    self.assertTrue(m2H.HasSubstructMatch(q3))

    self.assertTrue(m3H.HasSubstructMatch(q1))
    self.assertTrue(not m3H.HasSubstructMatch(q2))
    self.assertTrue(not m3H.HasSubstructMatch(q3))

  def test34Issue124(self) :
    """ issue 124 relates to calculation of the distance matrix
       
    """ 
    m = Chem.MolFromSmiles('CC=C')
    d = Chem.GetDistanceMatrix(m,0)
    self.assertTrue(feq(d[0,1],1.0))
    self.assertTrue(feq(d[0,2],2.0))
    # force an update:
    d = Chem.GetDistanceMatrix(m,1,0,1)
    self.assertTrue(feq(d[0,1],1.0))
    self.assertTrue(feq(d[0,2],1.5))

  def test35ChiralityPerception(self) :
    """ Test perception of chirality and CIP encoding       
    """ 
    m = Chem.MolFromSmiles('F[C@]([C@])(Cl)Br')
    Chem.AssignStereochemistry(m,1)
    self.assertTrue(m.GetAtomWithIdx(1).HasProp('_CIPCode'))
    self.assertFalse(m.GetAtomWithIdx(2).HasProp('_CIPCode'))
    Chem.RemoveStereochemistry(m)
    self.assertFalse(m.GetAtomWithIdx(1).HasProp('_CIPCode'))

    
    m = Chem.MolFromSmiles('F[C@H](C)C')
    Chem.AssignStereochemistry(m,1)
    self.assertTrue(m.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED)
    self.assertFalse(m.GetAtomWithIdx(1).HasProp('_CIPCode'))

    m = Chem.MolFromSmiles('F\\C=C/Cl')
    self.assertTrue(m.GetBondWithIdx(0).GetStereo()==Chem.BondStereo.STEREONONE)
    self.assertTrue(m.GetBondWithIdx(1).GetStereo()==Chem.BondStereo.STEREOZ)
    atoms = m.GetBondWithIdx(1).GetStereoAtoms()
    self.assertTrue(0 in atoms)
    self.assertTrue(3 in atoms)
    self.assertTrue(m.GetBondWithIdx(2).GetStereo()==Chem.BondStereo.STEREONONE)
    Chem.RemoveStereochemistry(m)
    self.assertTrue(m.GetBondWithIdx(1).GetStereo()==Chem.BondStereo.STEREONONE)

    
    m = Chem.MolFromSmiles('F\\C=CCl')
    self.assertTrue(m.GetBondWithIdx(1).GetStereo()==Chem.BondStereo.STEREONONE)

  def test36SubstructMatchStr(self):
    """ test the _SubstructMatchStr function """
    query = Chem.MolFromSmarts('[n,p]1ccccc1')
    self.assertTrue(query)
    mol = Chem.MolFromSmiles('N1=CC=CC=C1')
    self.assertTrue(mol.HasSubstructMatch(query))
    self.assertTrue(Chem._HasSubstructMatchStr(mol.ToBinary(),query))
    mol = Chem.MolFromSmiles('S1=CC=CC=C1')
    self.assertTrue(not Chem._HasSubstructMatchStr(mol.ToBinary(),query))
    self.assertTrue(not mol.HasSubstructMatch(query))
    mol = Chem.MolFromSmiles('P1=CC=CC=C1')
    self.assertTrue(mol.HasSubstructMatch(query))
    self.assertTrue(Chem._HasSubstructMatchStr(mol.ToBinary(),query))

    
  def test37SanitException(self):
    mol = Chem.MolFromSmiles('CC(C)(C)(C)C',0)
    self.assertTrue(mol)
    self.assertRaises(ValueError,lambda:Chem.SanitizeMol(mol))

  def test38TDTSuppliers(self):
    data="""$SMI<Cc1nnc(N)nc1C>
CAS<17584-12-2>
|
$SMI<Cc1n[nH]c(=O)nc1N>
CAS<~>
|
$SMI<Cc1n[nH]c(=O)[nH]c1=O>
CAS<932-53-6>
|
$SMI<Cc1nnc(NN)nc1O>
CAS<~>
|"""
    suppl = Chem.TDTMolSupplier()
    suppl.SetData(data,"CAS")
    i=0;
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetNumAtoms())
      self.assertTrue(mol.HasProp("CAS"))
      self.assertTrue(mol.HasProp("_Name"))
      self.assertTrue(mol.GetProp("CAS")==mol.GetProp("_Name"))
      self.assertTrue(mol.GetNumConformers()==0)
      i+=1
    self.assertTrue(i==4)
    self.assertTrue(len(suppl)==4)

  def test38Issue266(self):
    """ test issue 266: generation of kekulized smiles"""
    mol = Chem.MolFromSmiles('c1ccccc1')
    Chem.Kekulize(mol)
    smi = Chem.MolToSmiles(mol)
    self.assertTrue(smi=='c1ccccc1')
    smi = Chem.MolToSmiles(mol,kekuleSmiles=True)
    self.assertTrue(smi=='C1=CC=CC=C1')
    
  def test39Issue273(self):
    """ test issue 273: MolFileComments and MolFileInfo props ending up in SD files

    """
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap',
                                            'test_data','outNCI_few.sdf')
    suppl = Chem.SDMolSupplier(fileN)
    ms = [x for x in suppl]
    for m in ms:
      self.assertTrue(m.HasProp('_MolFileInfo'))
      self.assertTrue(m.HasProp('_MolFileComments'))
    fName = tempfile.mktemp('.sdf')
    w = Chem.SDWriter(fName)
    w.SetProps(ms[0].GetPropNames())
    for m in ms:
      w.write(m)
    w = None

    with open(fName, 'r') as txtFile:
      txt= txtFile.read()
    os.unlink(fName)
    self.assertTrue(txt.find('MolFileInfo')==-1)
    self.assertTrue(txt.find('MolFileComments')==-1)
    

  def test40SmilesRootedAtAtom(self):
    """ test the rootAtAtom functionality

    """
    smi = 'CN(C)C'
    m = Chem.MolFromSmiles(smi)

    self.assertTrue(Chem.MolToSmiles(m)=='CN(C)C')
    self.assertTrue(Chem.MolToSmiles(m,rootedAtAtom=1)=='N(C)(C)C')
    

  def test41SetStreamIndices(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
            "192", "203", "210", "211", "213", "220", "229", "256"]
    indices=[0, 2136, 6198, 8520, 11070, 12292, 14025, 15313, 17313, 20125, 22231,
	     23729, 26147, 28331, 32541, 33991]
    sdSup._SetStreamIndices(indices)
    self.assertTrue(len(sdSup) == 16)
    mol = sdSup[5]
    self.assertTrue(mol.GetProp("_Name") == "170")
    
    i = 0
    for mol in sdSup :
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
          
    ns = [mol.GetProp("_Name") for mol in sdSup]
    self.assertTrue(ns == molNames)

    # this can also be used to skip molecules in the file:
    indices=[0, 6198, 12292]
    sdSup._SetStreamIndices(indices)
    self.assertTrue(len(sdSup) == 3)
    mol = sdSup[2]
    self.assertTrue(mol.GetProp("_Name") == "170")

    # or to reorder them:
    indices=[0, 12292, 6198]
    sdSup._SetStreamIndices(indices)
    self.assertTrue(len(sdSup) == 3)
    mol = sdSup[1]
    self.assertTrue(mol.GetProp("_Name") == "170")

  def test42LifeTheUniverseAndEverything(self) :
    self.assertTrue(True)
    
  def test43TplFileParsing(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','cmpd2.tpl')
    m1 = Chem.MolFromTPLFile(fileN)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==12)
    self.assertTrue(m1.GetNumConformers()==2)
    
    m1 = Chem.MolFromTPLFile(fileN,skipFirstConf=True)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==12)
    self.assertTrue(m1.GetNumConformers()==1)

    with open(fileN, 'r') as blockFile:
      block = blockFile.read()
    m1 = Chem.MolFromTPLBlock(block)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==12)
    self.assertTrue(m1.GetNumConformers()==2)
    
  def test44TplFileWriting(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','cmpd2.tpl')
    m1 = Chem.MolFromTPLFile(fileN)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==12)
    self.assertTrue(m1.GetNumConformers()==2)

    block = Chem.MolToTPLBlock(m1)
    m1 = Chem.MolFromTPLBlock(block)
    self.assertTrue(m1 is not None)
    self.assertTrue(m1.GetNumAtoms()==12)
    self.assertTrue(m1.GetNumConformers()==2)

  def test45RingInfo(self):
    """ test the RingInfo class

    """
    smi = 'CNC'
    m = Chem.MolFromSmiles(smi)
    ri = m.GetRingInfo()
    self.assertTrue(ri)
    self.assertTrue(ri.NumRings()==0)
    self.assertFalse(ri.IsAtomInRingOfSize(0,3))
    self.assertFalse(ri.IsAtomInRingOfSize(1,3))
    self.assertFalse(ri.IsAtomInRingOfSize(2,3))
    self.assertFalse(ri.IsBondInRingOfSize(1,3))
    self.assertFalse(ri.IsBondInRingOfSize(2,3))

    smi = 'C1CC2C1C2'
    m = Chem.MolFromSmiles(smi)
    ri = m.GetRingInfo()
    self.assertTrue(ri)
    self.assertTrue(ri.NumRings()==2)
    self.assertFalse(ri.IsAtomInRingOfSize(0,3))
    self.assertTrue(ri.IsAtomInRingOfSize(0,4))
    self.assertFalse(ri.IsBondInRingOfSize(0,3))
    self.assertTrue(ri.IsBondInRingOfSize(0,4))
    self.assertTrue(ri.IsAtomInRingOfSize(2,4))
    self.assertTrue(ri.IsAtomInRingOfSize(2,3))
    self.assertTrue(ri.IsBondInRingOfSize(2,3))
    self.assertTrue(ri.IsBondInRingOfSize(2,4))

  def test46ReplaceCore(self):
    """ test the ReplaceCore functionality

    """

    core = Chem.MolFromSmiles('C=O')

    smi = 'CCC=O'
    m = Chem.MolFromSmiles(smi)
    r = Chem.ReplaceCore(m,core)
    self.assertTrue(r)
    self.assertEqual(Chem.MolToSmiles(r,True),'[1*]CC')

    smi = 'C1CC(=O)CC1'
    m = Chem.MolFromSmiles(smi)
    r = Chem.ReplaceCore(m,core)
    self.assertTrue(r)
    self.assertEqual(Chem.MolToSmiles(r,True),'[1*]CCCC[2*]')

    smi = 'C1CC(=N)CC1'
    m = Chem.MolFromSmiles(smi)
    r = Chem.ReplaceCore(m,core)
    self.assertFalse(r)

  def test47RWMols(self):
    """ test the RWMol class

    """
    mol = Chem.MolFromSmiles('C1CCC1')
    self.assertTrue(type(mol)==Chem.Mol)
    
    rwmol = Chem.EditableMol(mol)
    self.assertTrue(type(rwmol)==Chem.EditableMol)
    newAt = Chem.Atom(8)
    rwmol.ReplaceAtom(0,newAt)
    self.assertTrue(Chem.MolToSmiles(rwmol.GetMol())=='C1COC1')

    rwmol.RemoveBond(0,1)
    self.assertTrue(Chem.MolToSmiles(rwmol.GetMol())=='CCCO')
    a = Chem.Atom(7)
    idx=rwmol.AddAtom(a)
    self.assertEqual(rwmol.GetMol().GetNumAtoms(),5)
    self.assertEqual(idx,4)

    idx=rwmol.AddBond(0,4,order=Chem.BondType.SINGLE)
    self.assertEqual(idx,4)

    self.assertTrue(Chem.MolToSmiles(rwmol.GetMol())=='CCCON')
    rwmol.AddBond(4,1,order=Chem.BondType.SINGLE)
    self.assertTrue(Chem.MolToSmiles(rwmol.GetMol())=='C1CNOC1')
    
    rwmol.RemoveAtom(3)
    self.assertTrue(Chem.MolToSmiles(rwmol.GetMol())=='CCNO')
    
    # practice shooting ourselves in the foot:
    m = Chem.MolFromSmiles('c1ccccc1')
    em=Chem.EditableMol(m)
    em.RemoveAtom(0)
    m2 = em.GetMol()
    self.assertRaises(ValueError,lambda : Chem.SanitizeMol(m2))
    m = Chem.MolFromSmiles('c1ccccc1')
    em=Chem.EditableMol(m)
    em.RemoveBond(0,1)
    m2 = em.GetMol()
    self.assertRaises(ValueError,lambda : Chem.SanitizeMol(m2))
                    
    # boundary cases: 
    
    # removing non-existent bonds:
    m = Chem.MolFromSmiles('c1ccccc1')
    em=Chem.EditableMol(m)
    em.RemoveBond(0,2)
    m2 = em.GetMol()
    Chem.SanitizeMol(m2)
    self.assertTrue(Chem.MolToSmiles(m2)=='c1ccccc1')
    
    # removing non-existent atoms:
    m = Chem.MolFromSmiles('c1ccccc1')
    em=Chem.EditableMol(m)
    self.assertRaises(RuntimeError,lambda:em.RemoveAtom(12))

  def test47SmartsPieces(self):
    """ test the GetAtomSmarts and GetBondSmarts functions

    """
    m =Chem.MolFromSmarts("[C,N]C")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts()=='[C,N]')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts()=='C')
    self.assertTrue(m.GetBondBetweenAtoms(0,1).GetSmarts()=='-,:')
    
    m =Chem.MolFromSmarts("[$(C=O)]-O")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts()=='[$(C=O)]')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts()=='O')
    self.assertTrue(m.GetBondBetweenAtoms(0,1).GetSmarts()=='-')

    m =Chem.MolFromSmiles("CO")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts()=='C')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts()=='O')
    self.assertTrue(m.GetBondBetweenAtoms(0,1).GetSmarts()=='')
    self.assertTrue(m.GetBondBetweenAtoms(0,1).GetSmarts(allBondsExplicit=True)=='-')

    m =Chem.MolFromSmiles("C=O")
    self.assertTrue(m.GetAtomWithIdx(0).GetSmarts()=='C')
    self.assertTrue(m.GetAtomWithIdx(1).GetSmarts()=='O')
    self.assertTrue(m.GetBondBetweenAtoms(0,1).GetSmarts()=='=')

  def test48Issue1928819(self):
    """ test a crash involving looping directly over mol suppliers
    """
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf')
    ms = [x for x in Chem.SDMolSupplier(fileN)]
    self.assertEqual(len(ms),16)
    count=0
    for m in Chem.SDMolSupplier(fileN): count+=1
    self.assertEqual(count,16)

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','fewSmi.csv')
    count=0
    for m in Chem.SmilesMolSupplier(fileN,titleLine=False,smilesColumn=1,delimiter=','): count+=1
    self.assertEqual(count,10)
    
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','acd_few.tdt')
    count=0
    for m in Chem.TDTMolSupplier(fileN): count+=1
    self.assertEqual(count,10)

  def test49Issue1932365(self):
    """ test aromatic Se and Te from smiles/smarts
    """
    m = Chem.MolFromSmiles('c1ccc[se]1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(4).GetIsAromatic())
    m = Chem.MolFromSmiles('c1ccc[te]1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(4).GetIsAromatic())
    m = Chem.MolFromSmiles('C1=C[Se]C=C1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(2).GetIsAromatic())
    m = Chem.MolFromSmiles('C1=C[Te]C=C1')
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(0).GetIsAromatic())
    self.assertTrue(m.GetAtomWithIdx(2).GetIsAromatic())

    p = Chem.MolFromSmarts('[se]')
    self.assertTrue(Chem.MolFromSmiles('c1ccc[se]1').HasSubstructMatch(p))
    self.assertFalse(Chem.MolFromSmiles('C1=CCC[Se]1').HasSubstructMatch(p))
    
    p = Chem.MolFromSmarts('[te]')
    self.assertTrue(Chem.MolFromSmiles('c1ccc[te]1').HasSubstructMatch(p))
    self.assertFalse(Chem.MolFromSmiles('C1=CCC[Te]1').HasSubstructMatch(p))

  def test50Issue1968608(self):
    """ test sf.net issue 1968608
    """
    smarts = Chem.MolFromSmarts("[r5]")
    mol = Chem.MolFromSmiles("N12CCC36C1CC(C(C2)=CCOC4CC5=O)C4C3N5c7ccccc76")
    count = len(mol.GetSubstructMatches(smarts, uniquify=0))
    self.assertTrue(count==9)

  def test51RadicalHandling(self):
    """ test handling of atoms with radicals
    """
    mol = Chem.MolFromSmiles("[C]C")
    self.assertTrue(mol)
    atom=mol.GetAtomWithIdx(0)
    self.assertTrue(atom.GetNumRadicalElectrons()==3)
    self.assertTrue(atom.GetNoImplicit())
    atom.SetNoImplicit(False)
    atom.SetNumRadicalElectrons(1)
    mol.UpdatePropertyCache()
    self.assertTrue(atom.GetNumRadicalElectrons()==1)
    self.assertTrue(atom.GetNumImplicitHs()==2)
    
    mol = Chem.MolFromSmiles("[c]1ccccc1")
    self.assertTrue(mol)
    atom=mol.GetAtomWithIdx(0)
    self.assertTrue(atom.GetNumRadicalElectrons()==1)
    self.assertTrue(atom.GetNoImplicit())

    mol = Chem.MolFromSmiles("[n]1ccccc1")
    self.assertTrue(mol)
    atom=mol.GetAtomWithIdx(0)
    self.assertTrue(atom.GetNumRadicalElectrons()==0)
    self.assertTrue(atom.GetNoImplicit())


  def test52MolFrags(self):
    """ test GetMolFrags functionality
    """
    mol = Chem.MolFromSmiles("C.CC")
    self.assertTrue(mol)
    fs = Chem.GetMolFrags(mol)
    self.assertTrue(len(fs)==2)
    self.assertTrue(len(fs[0])==1)
    self.assertTrue(tuple(fs[0])==(0,))
    self.assertTrue(len(fs[1])==2)
    self.assertTrue(tuple(fs[1])==(1,2))

    fs = Chem.GetMolFrags(mol,True)
    self.assertTrue(len(fs)==2)
    self.assertTrue(fs[0].GetNumAtoms()==1)
    self.assertTrue(fs[1].GetNumAtoms()==2)
    
    mol = Chem.MolFromSmiles("CCC")
    self.assertTrue(mol)
    fs = Chem.GetMolFrags(mol)
    self.assertTrue(len(fs)==1)
    self.assertTrue(len(fs[0])==3)
    self.assertTrue(tuple(fs[0])==(0,1,2))
    fs = Chem.GetMolFrags(mol,True)
    self.assertTrue(len(fs)==1)
    self.assertTrue(fs[0].GetNumAtoms()==3)

    mol = Chem.MolFromSmiles("CO")
    em = Chem.EditableMol(mol)
    em.RemoveBond(0,1)
    nm = em.GetMol()
    fs = Chem.GetMolFrags(nm,asMols=True)
    self.assertEqual([x.GetNumAtoms(onlyExplicit=False) for x in fs],[5,3])
    fs = Chem.GetMolFrags(nm,asMols=True,sanitizeFrags=False)
    self.assertEqual([x.GetNumAtoms(onlyExplicit=False) for x in fs],[4,2])
    



  def test53Matrices(self) :
    """ test adjacency and distance matrices
       
    """ 
    m = Chem.MolFromSmiles('CC=C')
    d = Chem.GetDistanceMatrix(m,0)
    self.assertTrue(feq(d[0,1],1.0))
    self.assertTrue(feq(d[0,2],2.0))
    self.assertTrue(feq(d[1,0],1.0))
    self.assertTrue(feq(d[2,0],2.0))
    a = Chem.GetAdjacencyMatrix(m,0)
    self.assertTrue(a[0,1]==1)
    self.assertTrue(a[0,2]==0)
    self.assertTrue(a[1,2]==1)
    self.assertTrue(a[1,0]==1)
    self.assertTrue(a[2,0]==0)

    m = Chem.MolFromSmiles('C1CC1')
    d = Chem.GetDistanceMatrix(m,0)
    self.assertTrue(feq(d[0,1],1.0))
    self.assertTrue(feq(d[0,2],1.0))
    a = Chem.GetAdjacencyMatrix(m,0)
    self.assertTrue(a[0,1]==1)
    self.assertTrue(a[0,2]==1)
    self.assertTrue(a[1,2]==1)

    m = Chem.MolFromSmiles('CC.C')
    d = Chem.GetDistanceMatrix(m,0)
    self.assertTrue(feq(d[0,1],1.0))
    self.assertTrue(d[0,2]>1000)
    self.assertTrue(d[1,2]>1000)
    a = Chem.GetAdjacencyMatrix(m,0)
    self.assertTrue(a[0,1]==1)
    self.assertTrue(a[0,2]==0)
    self.assertTrue(a[1,2]==0)

  def test54Mol2Parser(self):
    """ test the mol2 parser
    """
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','pyrazole_pyridine.mol2')
    m = Chem.MolFromMol2File(fileN)
    self.assertTrue(m.GetNumAtoms()==5)
    self.assertTrue(Chem.MolToSmiles(m)=='c1cn[nH]c1',Chem.MolToSmiles(m))

  def test55LayeredFingerprint(self):
    m1 = Chem.MolFromSmiles('CC(C)C')
    fp1 = Chem.LayeredFingerprint(m1)
    self.assertEqual(len(fp1),2048)
    atomCounts=[0]*m1.GetNumAtoms()
    fp2 = Chem.LayeredFingerprint(m1,atomCounts=atomCounts)
    self.assertEqual(fp1,fp2)
    self.assertEqual(atomCounts,[4,7,4,4])

    fp2 = Chem.LayeredFingerprint(m1,atomCounts=atomCounts)
    self.assertEqual(fp1,fp2)
    self.assertEqual(atomCounts,[8,14,8,8])

    pbv=DataStructs.ExplicitBitVect(2048)
    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=pbv)
    self.assertEqual(fp3.GetNumOnBits(),0)

    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=fp2)
    self.assertEqual(fp3,fp2)

    m2=Chem.MolFromSmiles('CC')
    fp4 = Chem.LayeredFingerprint(m2)
    atomCounts=[0]*m1.GetNumAtoms()
    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=fp4,atomCounts=atomCounts)
    self.assertEqual(atomCounts,[1,3,1,1])

    m2=Chem.MolFromSmiles('CCC')
    fp4 = Chem.LayeredFingerprint(m2)
    atomCounts=[0]*m1.GetNumAtoms()
    fp3 = Chem.LayeredFingerprint(m1,setOnlyBits=fp4,atomCounts=atomCounts)
    self.assertEqual(atomCounts,[3,6,3,3])


  def test56LazySDMolSupplier(self) :
    if not hasattr(Chem,'CompressedSDMolSupplier'): return

    self.assertRaises(ValueError,lambda : Chem.CompressedSDMolSupplier('nosuchfile.sdf.gz'))

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf.gz')
    sdSup = Chem.CompressedSDMolSupplier(fileN)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
            "192", "203", "210", "211", "213", "220", "229", "256"]

    chgs192 = {8:1, 11:1, 15:-1, 18:-1, 20:1, 21:1, 23:-1, 25:-1} 
    i = 0
    for mol in sdSup :
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
      if (mol.GetProp("_Name") == "192") :
        # test parsed charges on one of the molecules
        for id in chgs192.keys() :
          self.assertTrue(mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id])
    self.assertEqual(i,16)
          
    sdSup = Chem.CompressedSDMolSupplier(fileN)
    ns = [mol.GetProp("_Name") for mol in sdSup]
    self.assertTrue(ns == molNames)

    sdSup = Chem.CompressedSDMolSupplier(fileN, 0)
    for mol in sdSup :
      self.assertTrue(not mol.HasProp("numArom"))

  def test57AddRecursiveQuery(self):
    q1 = Chem.MolFromSmiles('CC')
    q2 = Chem.MolFromSmiles('CO')
    Chem.AddRecursiveQuery(q1,q2,1)

    m1 = Chem.MolFromSmiles('OCC')
    self.assertTrue(m1.HasSubstructMatch(q2))
    self.assertTrue(m1.HasSubstructMatch(q1))
    self.assertTrue(m1.HasSubstructMatch(q1))
    self.assertTrue(m1.GetSubstructMatch(q1)==(2,1))

    q3 = Chem.MolFromSmiles('CS')
    Chem.AddRecursiveQuery(q1,q3,1)
    
    self.assertFalse(m1.HasSubstructMatch(q3))
    self.assertFalse(m1.HasSubstructMatch(q1))

    m2 = Chem.MolFromSmiles('OC(S)C')
    self.assertTrue(m2.HasSubstructMatch(q1))
    self.assertTrue(m2.GetSubstructMatch(q1)==(3,1))
    
    m3 = Chem.MolFromSmiles('SCC')
    self.assertTrue(m3.HasSubstructMatch(q3))
    self.assertFalse(m3.HasSubstructMatch(q1))

    q1 = Chem.MolFromSmiles('CC')
    Chem.AddRecursiveQuery(q1,q2,1)
    Chem.AddRecursiveQuery(q1,q3,1,False)
    self.assertTrue(m3.HasSubstructMatch(q1))
    self.assertTrue(m3.GetSubstructMatch(q1)==(2,1))

  def test58Issue2983794(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap',
                                            'test_data','issue2983794.sdf')
    m1 = Chem.MolFromMolFile(fileN)
    self.assertTrue(m1)
    em = Chem.EditableMol(m1)
    em.RemoveAtom(0)
    m2 = em.GetMol()
    Chem.Kekulize(m2)
    
  def test59Issue3007178(self) :
    m = Chem.MolFromSmiles('CCC')
    a = m.GetAtomWithIdx(0)
    m=None
    self.assertEqual(Chem.MolToSmiles(a.GetOwningMol()),'CCC')
    a=None
    m = Chem.MolFromSmiles('CCC')
    b = m.GetBondWithIdx(0)
    m=None
    self.assertEqual(Chem.MolToSmiles(b.GetOwningMol()),'CCC')

  def test60SmilesWriterClose(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','fewSmi.csv')
    smiSup = Chem.SmilesMolSupplier(fileN, delimiter=",",
                                      smilesColumn=1, nameColumn=0,
                                      titleLine=0)
    ms = [x for x in smiSup]
    
    ofile = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap','test_data','outSmiles.txt')
    writer = Chem.SmilesWriter(ofile)
    for mol in ms:
      writer.write(mol)
    writer.close()

    newsup=Chem.SmilesMolSupplier(ofile)
    newms = [x for x  in newsup]
    self.assertEqual(len(ms),len(newms))

  def test61PathToSubmol(self):
    m = Chem.MolFromSmiles('CCCCCC1C(O)CC(O)N1C=CCO')
    env = Chem.FindAtomEnvironmentOfRadiusN(m,2,11)
    self.assertEqual(len(env),8)
    amap={}
    submol = Chem.PathToSubmol(m,env,atomMap=amap)
    self.assertEqual(submol.GetNumAtoms(),len(amap.keys()))
    self.assertEqual(submol.GetNumAtoms(),9)
    smi=Chem.MolToSmiles(submol,rootedAtAtom=amap[11])
    self.assertEqual(smi[0],'N')
    refsmi = Chem.MolToSmiles(Chem.MolFromSmiles('N(C=C)(C(C)C)C(O)C'))
    csmi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    self.assertEqual(refsmi,csmi)
    
  def test62SmilesAndSmartsReplacements(self):
    mol = Chem.MolFromSmiles('C{branch}C',replacements={'{branch}':'C1(CC1)'})
    self.assertEqual(mol.GetNumAtoms(),5)
    mol = Chem.MolFromSmarts('C{branch}C',replacements={'{branch}':'C1(CC1)'})
    self.assertEqual(mol.GetNumAtoms(),5)
    mol = Chem.MolFromSmiles('C{branch}C{acid}',replacements={'{branch}':'C1(CC1)',
                                                              '{acid}':"C(=O)O"})
    self.assertEqual(mol.GetNumAtoms(),8)

  def test63Issue3313539(self):
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','rgroups1.mol')
    m = Chem.MolFromMolFile(fileN)
    self.assertTrue(m is not None)
    at = m.GetAtomWithIdx(3)
    self.assertTrue(at is not None)
    self.assertTrue(at.HasProp('_MolFileRLabel'))
    p = at.GetProp('_MolFileRLabel')
    self.assertEqual(p,'2')

    at = m.GetAtomWithIdx(4)
    self.assertTrue(at is not None)
    self.assertTrue(at.HasProp('_MolFileRLabel'))
    p = at.GetProp('_MolFileRLabel')
    self.assertEqual(p,'1')
    
  def test64MoleculeCleanup(self):
    m = Chem.MolFromSmiles('CN(=O)=O',False)
    self.assertTrue(m)
    self.assertTrue(m.GetAtomWithIdx(1).GetFormalCharge()==0 and \
                      m.GetAtomWithIdx(2).GetFormalCharge()==0 and \
                      m.GetAtomWithIdx(3).GetFormalCharge()==0)
    self.assertTrue(m.GetBondBetweenAtoms(1,3).GetBondType()==Chem.BondType.DOUBLE and \
                      m.GetBondBetweenAtoms(1,2).GetBondType()==Chem.BondType.DOUBLE )
    Chem.Cleanup(m)
    m.UpdatePropertyCache()
    self.assertTrue(m.GetAtomWithIdx(1).GetFormalCharge()==1 and \
                      (m.GetAtomWithIdx(2).GetFormalCharge()==-1 or \
                         m.GetAtomWithIdx(3).GetFormalCharge()==-1))
    self.assertTrue(m.GetBondBetweenAtoms(1,3).GetBondType()==Chem.BondType.SINGLE or \
                      m.GetBondBetweenAtoms(1,2).GetBondType()==Chem.BondType.SINGLE )

  def test65StreamSupplier(self):
    import gzip
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf.gz')
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
                "192", "203", "210", "211", "213", "220", "229", "256"]
    inf = gzip.open(fileN)
    if 0:
      sb = Chem.streambuf(inf)
      suppl = Chem.ForwardSDMolSupplier(sb)
    else:
      suppl = Chem.ForwardSDMolSupplier(inf)
      
    i = 0
    while not suppl.atEnd():
      mol = six.next(suppl)
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i,16)

    # make sure we have object ownership preserved
    inf = gzip.open(fileN)
    suppl = Chem.ForwardSDMolSupplier(inf)
    inf=None
    i = 0
    while not suppl.atEnd():
      mol = six.next(suppl)
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i,16)

  def test66StreamSupplierIter(self):
    import gzip
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf.gz')
    inf = gzip.open(fileN)
    if 0:
      sb = Chem.streambuf(inf)
      suppl = Chem.ForwardSDMolSupplier(sb)
    else:
      suppl = Chem.ForwardSDMolSupplier(inf)
      
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
                "192", "203", "210", "211", "213", "220", "229", "256"]
    i = 0
    for mol in suppl :
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i,16)

  def test67StreamSupplierStringIO(self):
    import gzip
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf.gz')
    if six.PY3:
      from io import BytesIO
      sio = BytesIO(gzip.open(fileN).read())
    else:
      import StringIO
      sio = StringIO.StringIO(gzip.open(fileN).read())
    suppl = Chem.ForwardSDMolSupplier(sio)
      
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
                "192", "203", "210", "211", "213", "220", "229", "256"]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i,16)

  def test68ForwardSupplierUsingFilename(self):
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf')
    suppl = Chem.ForwardSDMolSupplier(fileN)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
                "192", "203", "210", "211", "213", "220", "229", "256"]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i,16)

    self.assertRaises(IOError,lambda : Chem.ForwardSDMolSupplier('nosuchfile.sdf'))

  def test69StreamSupplierStreambuf(self):
    import gzip
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf.gz')
    sb = rdBase.streambuf(gzip.open(fileN))
    suppl = Chem.ForwardSDMolSupplier(sb)
      
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
                "192", "203", "210", "211", "213", "220", "229", "256"]
    i = 0
    for mol in suppl:
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i,16)
    
  def test70StreamSDWriter(self):
    import gzip
    if six.PY3:
      from io import BytesIO,StringIO
    else:
      from StringIO import StringIO


    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf.gz')
    inf = gzip.open(fileN)
    suppl = Chem.ForwardSDMolSupplier(inf)
    osio=StringIO()
    w = Chem.SDWriter(osio)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
                "192", "203", "210", "211", "213", "220", "229", "256"]
    i = 0
    for mol in suppl :
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      w.write(mol)
      i += 1
    self.assertEqual(i,16)
    w.flush()
    w=None
    if six.PY3:
      txt = osio.getvalue().encode()
      isio = BytesIO(txt)
    else:
      isio = StringIO(osio.getvalue())
    suppl = Chem.ForwardSDMolSupplier(isio)
    i = 0
    for mol in suppl :
      self.assertTrue(mol)
      self.assertTrue(mol.GetProp("_Name") == molNames[i])
      i += 1
    self.assertEqual(i,16)

  def test71StreamSmilesWriter(self):
    from rdkit.six.moves import StringIO
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','esters.sdf')
    suppl = Chem.ForwardSDMolSupplier(fileN)
    osio=StringIO()
    w = Chem.SmilesWriter(osio)
    ms = [x for x in suppl]
    w.SetProps(ms[0].GetPropNames())
    i=0
    for mol in ms:
      self.assertTrue(mol)
      w.write(mol)
      i+=1
    self.assertEqual(i,6)
    w.flush()
    w=None
    txt = osio.getvalue()
    self.assertEqual(txt.count('ID'),1)
    self.assertEqual(txt.count('\n'),7)

  def test72StreamTDTWriter(self):
    from rdkit.six.moves import StringIO
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','esters.sdf')
    suppl = Chem.ForwardSDMolSupplier(fileN)
    osio=StringIO()
    w = Chem.TDTWriter(osio)
    ms = [x for x in suppl]
    w.SetProps(ms[0].GetPropNames())
    i=0
    for mol in ms:
      self.assertTrue(mol)
      w.write(mol)
      i+=1
    self.assertEqual(i,6)
    w.flush()
    w=None
    txt = osio.getvalue()
    self.assertEqual(txt.count('ID'),6)
    self.assertEqual(txt.count('NAME'),6)

  def test73SanitizationOptions(self):
    m = Chem.MolFromSmiles('c1ccccc1',sanitize=False)
    res = Chem.SanitizeMol(m,catchErrors=True)
    self.assertEqual(res,0)

    m = Chem.MolFromSmiles('c1cccc1',sanitize=False)
    res = Chem.SanitizeMol(m,catchErrors=True)
    self.assertEqual(res,Chem.SanitizeFlags.SANITIZE_KEKULIZE)

    m = Chem.MolFromSmiles('CC(C)(C)(C)C',sanitize=False)
    res = Chem.SanitizeMol(m,catchErrors=True)
    self.assertEqual(res,Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    
    m = Chem.MolFromSmiles('c1cccc1',sanitize=False)
    res = Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE,
                           catchErrors=True)
    self.assertEqual(res,Chem.SanitizeFlags.SANITIZE_NONE)
    
  def test74Issue3510149(self):
    mol = Chem.MolFromSmiles("CCC1CNCC1CC")
    atoms = mol.GetAtoms()
    mol=None
    for atom in atoms:
      idx=atom.GetIdx()
      p= atom.GetOwningMol().GetNumAtoms()

    mol = Chem.MolFromSmiles("CCC1CNCC1CC")
    bonds = mol.GetBonds()
    mol=None
    for bond in bonds:
      idx=bond.GetIdx()
      p= atom.GetOwningMol().GetNumAtoms()

    mol = Chem.MolFromSmiles("CCC1CNCC1CC")
    bond = mol.GetBondBetweenAtoms(0,1)
    mol=None
    idx=bond.GetBeginAtomIdx()
    p= bond.GetOwningMol().GetNumAtoms()
      
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    mol = six.next(sdSup)
    nats = mol.GetNumAtoms()
    conf = mol.GetConformer()
    mol=None
    self.assertEqual(nats,conf.GetNumAtoms())
    conf.GetOwningMol().GetProp("_Name")

  def test75AllBondsExplicit(self):
    m = Chem.MolFromSmiles("CCC")
    smi = Chem.MolToSmiles(m)
    self.assertEqual(smi,"CCC")
    smi = Chem.MolToSmiles(m,allBondsExplicit=True)
    self.assertEqual(smi,"C-C-C")

    m = Chem.MolFromSmiles("c1ccccc1")
    smi = Chem.MolToSmiles(m)
    self.assertEqual(smi,"c1ccccc1")
    smi = Chem.MolToSmiles(m,allBondsExplicit=True)
    self.assertEqual(smi,"c1:c:c:c:c:c:1")
    
  def test76VeryLargeMolecule(self):
    # this is sf.net issue 3524984
    smi = '[C@H](F)(Cl)'+'c1cc[nH]c1'*500+'[C@H](F)(Cl)'
    m = Chem.MolFromSmiles(smi)
    self.assertTrue(m)
    self.assertEqual(m.GetNumAtoms(),2506)
    scs = Chem.FindMolChiralCenters(m)
    self.assertEqual(len(scs),2)

  def test77MolFragmentToSmiles(self):
    smi="OC1CC1CC"
    m = Chem.MolFromSmiles(smi)
    fsmi = Chem.MolFragmentToSmiles(m,[1,2,3])
    self.assertEqual(fsmi,"C1CC1")
    fsmi = Chem.MolFragmentToSmiles(m,[1,2,3],bondsToUse=[1,2,5])
    self.assertEqual(fsmi,"C1CC1")
    fsmi = Chem.MolFragmentToSmiles(m,[1,2,3],bondsToUse=[1,2])
    self.assertEqual(fsmi,"CCC")
    fsmi = Chem.MolFragmentToSmiles(m,[1,2,3],atomSymbols=["","[A]","[C]","[B]","",""])
    self.assertEqual(fsmi,"[C]1[B][A]1")
    fsmi = Chem.MolFragmentToSmiles(m,[1,2,3],bondSymbols=["","%","%","","","%"])
    self.assertEqual(fsmi,"C1%C%C%1")
    
    smi="c1ccccc1C"
    m = Chem.MolFromSmiles(smi)
    fsmi = Chem.MolFragmentToSmiles(m,range(6))
    self.assertEqual(fsmi,"c1ccccc1")
    Chem.Kekulize(m)
    fsmi = Chem.MolFragmentToSmiles(m,range(6),kekuleSmiles=True)
    self.assertEqual(fsmi,"C1=CC=CC=C1")    
    fsmi = Chem.MolFragmentToSmiles(m,range(6),atomSymbols=["[C]"]*7,kekuleSmiles=True)
    self.assertEqual(fsmi,"[C]1=[C][C]=[C][C]=[C]1")

    self.assertRaises(ValueError,lambda : Chem.MolFragmentToSmiles(m,[]))

  def test78AtomAndBondProps(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    at = m.GetAtomWithIdx(0)
    self.assertFalse(at.HasProp('foo'))
    at.SetProp('foo','bar')
    self.assertTrue(at.HasProp('foo'))
    self.assertEqual(at.GetProp('foo'),'bar')
    bond = m.GetBondWithIdx(0)
    self.assertFalse(bond.HasProp('foo'))    
    bond.SetProp('foo','bar')
    self.assertTrue(bond.HasProp('foo'))
    self.assertEqual(bond.GetProp('foo'),'bar')
    
  def test79AddRecursiveStructureQueries(self):
    qs = {'carbonyl':Chem.MolFromSmiles('CO'),
          'amine':Chem.MolFromSmiles('CN')}
    q = Chem.MolFromSmiles('CCC')
    q.GetAtomWithIdx(0).SetProp('query','carbonyl,amine')
    Chem.MolAddRecursiveQueries(q,qs,'query')
    m = Chem.MolFromSmiles('CCCO')
    self.assertTrue(m.HasSubstructMatch(q))
    m = Chem.MolFromSmiles('CCCN')
    self.assertTrue(m.HasSubstructMatch(q))
    m = Chem.MolFromSmiles('CCCC')
    self.assertFalse(m.HasSubstructMatch(q))

  def test80ParseMolQueryDefFile(self):
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','ChemTransforms',
                                            'testData','query_file1.txt')
    d = Chem.ParseMolQueryDefFile(fileN,standardize=False)
    self.assertTrue('CarboxylicAcid' in d)
    m = Chem.MolFromSmiles('CC(=O)O')
    self.assertTrue(m.HasSubstructMatch(d['CarboxylicAcid']))
    self.assertFalse(m.HasSubstructMatch(d['CarboxylicAcid.Aromatic']))

    d = Chem.ParseMolQueryDefFile(fileN)
    self.assertTrue('carboxylicacid' in d)
    self.assertFalse('CarboxylicAcid' in d)
                    
  def test81Issue275(self):
    smi = Chem.MolToSmiles(Chem.MurckoDecompose(Chem.MolFromSmiles('CCCCC[C@H]1CC[C@H](C(=O)O)CC1')))
    self.assertEqual(smi,'C1CCCCC1')

  def test82Issue288(self):
    m = Chem.MolFromSmiles('CC*')
    m.GetAtomWithIdx(2).SetProp('molAtomMapNumber','30')
    smi=Chem.MolToSmiles(m)
    self.assertEqual(smi,'[*:30]CC')

  def test83GitHubIssue19(self):
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','empty2.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    self.assertTrue(sdSup.atEnd())
    self.assertRaises(IndexError,lambda : sdSup[0])

    sdSup.SetData('')
    self.assertTrue(sdSup.atEnd())
    self.assertRaises(IndexError,lambda : sdSup[0])

    sdSup = Chem.SDMolSupplier(fileN)
    self.assertRaises(IndexError,lambda : sdSup[0])

    sdSup.SetData('')
    self.assertRaises(IndexError,lambda : sdSup[0])

    sdSup = Chem.SDMolSupplier(fileN)
    self.assertEqual(len(sdSup),0)

    sdSup.SetData('')
    self.assertEqual(len(sdSup),0)

  def test84PDBBasics(self):
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','1CRN.pdb')
    m = Chem.MolFromPDBFile(fileN)
    self.assertTrue(m is not None)
    self.assertEqual(m.GetNumAtoms(),327)
    self.assertEqual(m.GetNumBonds(),337)
    self.assertTrue(m.GetAtomWithIdx(0).GetPDBResidueInfo())
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetName()," N  ")
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName(),"THR")
    self.assertAlmostEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetTempFactor(),13.79,2)
    m = Chem.MolFromPDBBlock(Chem.MolToPDBBlock(m))
    self.assertEqual(m.GetNumAtoms(),327)
    self.assertEqual(m.GetNumBonds(),337)
    self.assertTrue(m.GetAtomWithIdx(0).GetPDBResidueInfo())
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetName()," N  ")
    self.assertEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetResidueName(),"THR")
    self.assertAlmostEqual(m.GetAtomWithIdx(0).GetPDBResidueInfo().GetTempFactor(),13.79,2)
    
  def test85MolCopying(self):
    m = Chem.MolFromSmiles('C1CC1[C@H](F)Cl')
    m.SetProp('foo','bar')
    m2 = Chem.Mol(m)
    self.assertEqual(Chem.MolToSmiles(m,True),Chem.MolToSmiles(m2,True))
    self.assertTrue(m2.HasProp('foo'))
    self.assertEqual(m2.GetProp('foo'),'bar')    
    ri = m2.GetRingInfo()
    self.assertTrue(ri)
    self.assertEqual(ri.NumRings(),1)
    
  def test86MolRenumbering(self):
    import random
    m = Chem.MolFromSmiles('C[C@H]1CC[C@H](C/C=C/[C@H](F)Cl)CC1')
    cSmi = Chem.MolToSmiles(m,True)
    for i in range(m.GetNumAtoms()):
      ans = list(range(m.GetNumAtoms()))
      random.shuffle(ans)
      m2 = Chem.RenumberAtoms(m,ans)
      nSmi = Chem.MolToSmiles(m2,True)
      self.assertEqual(cSmi,nSmi)

  def test87FragmentOnBonds(self):
    m = Chem.MolFromSmiles('CC1CC(O)C1CCC1CC1')
    bis = m.GetSubstructMatches(Chem.MolFromSmarts('[!R][R]'))
    bs = []
    labels=[]
    for bi in bis:
        b = m.GetBondBetweenAtoms(bi[0],bi[1])
        if b.GetBeginAtomIdx()==bi[0]:
            labels.append((10,1))
        else:
            labels.append((1,10))
        bs.append(b.GetIdx())
    nm = Chem.FragmentOnBonds(m,bs)
    frags = Chem.GetMolFrags(nm)
    self.assertEqual(len(frags),5)
    self.assertEqual(frags,((0, 12), (1, 2, 3, 5, 11, 14, 16), (4, 13), (6, 7, 15, 18), (8, 9, 10, 17)))
    smi = Chem.MolToSmiles(nm,True)
    self.assertEqual(smi,'[*]C1CC([4*])C1[6*].[1*]C.[3*]O.[5*]CC[8*].[7*]C1CC1')

    nm = Chem.FragmentOnBonds(m,bs,dummyLabels=labels)
    frags = Chem.GetMolFrags(nm)
    self.assertEqual(len(frags),5)
    self.assertEqual(frags,((0, 12), (1, 2, 3, 5, 11, 14, 16), (4, 13), (6, 7, 15, 18), (8, 9, 10, 17)))
    smi = Chem.MolToSmiles(nm,True)
    self.assertEqual(smi,'[1*]C.[1*]O.[1*]CC[1*].[10*]C1CC1.[10*]C1CC([10*])C1[10*]')

    m = Chem.MolFromSmiles('CCC(=O)CC(=O)C')
    bis = m.GetSubstructMatches(Chem.MolFromSmarts('C=O'))
    bs = []
    for bi in bis:
        b = m.GetBondBetweenAtoms(bi[0],bi[1])
        bs.append(b.GetIdx())
    bts = [Chem.BondType.DOUBLE]*len(bs)
    nm = Chem.FragmentOnBonds(m,bs,bondTypes=bts)
    frags = Chem.GetMolFrags(nm)
    self.assertEqual(len(frags),3)
    smi = Chem.MolToSmiles(nm,True)
    self.assertEqual(smi,'[2*]=O.[3*]=C(CC)CC(=[6*])C.[5*]=O')

  def test88QueryAtoms(self):
    from rdkit.Chem import rdqueries
    m = Chem.MolFromSmiles('c1nc(C)n(CC)c1')
    
    qa = rdqueries.ExplicitDegreeEqualsQueryAtom(3)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l,(2,4))

    qa.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(6,negate=True))
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l,(4,))
    
    qa = rdqueries.ExplicitDegreeEqualsQueryAtom(3)
    qa.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(6,negate=True),
                   how=Chem.CompositeQueryType.COMPOSITE_OR)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l,(1,2,4))
    
    qa = rdqueries.ExplicitDegreeEqualsQueryAtom(3)
    qa.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(6,negate=True),
                   how=Chem.CompositeQueryType.COMPOSITE_XOR)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l,(1,2))
    
    qa = rdqueries.ExplicitDegreeGreaterQueryAtom(2)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l,(2,4))

    qa = rdqueries.ExplicitDegreeLessQueryAtom(2)
    l = tuple([x.GetIdx() for x in m.GetAtomsMatchingQuery(qa)])
    self.assertEqual(l,(3,6))

    
  def test89UnicodeInput(self):
    m = Chem.MolFromSmiles(u'c1ccccc1')
    self.assertTrue(m is not None)
    self.assertEqual(m.GetNumAtoms(),6)
    m = Chem.MolFromSmarts(u'c1ccccc1')
    self.assertTrue(m is not None)
    self.assertEqual(m.GetNumAtoms(),6)
    
  def test90FragmentOnSomeBonds(self):
    m = Chem.MolFromSmiles('OCCCCN')
    pieces = Chem.FragmentOnSomeBonds(m,(0,2,4),2)
    self.assertEqual(len(pieces),3)
    
    frags = Chem.GetMolFrags(pieces[0])
    self.assertEqual(len(frags),3)
    self.assertEqual(len(frags[0]),2)
    self.assertEqual(len(frags[1]),4)
    self.assertEqual(len(frags[2]),4)

    frags = Chem.GetMolFrags(pieces[1])
    self.assertEqual(len(frags),3)
    self.assertEqual(len(frags[0]),2)
    self.assertEqual(len(frags[1]),6)
    self.assertEqual(len(frags[2]),2)

    frags = Chem.GetMolFrags(pieces[2])
    self.assertEqual(len(frags),3)
    self.assertEqual(len(frags[0]),4)
    self.assertEqual(len(frags[1]),4)
    self.assertEqual(len(frags[2]),2)

    
if __name__ == '__main__':
  unittest.main()

