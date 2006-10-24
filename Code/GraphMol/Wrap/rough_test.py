# $Id$
#
#  Copyright (C) 2003-2005  Rational Discovery LLC
#         All Rights Reserved
#
""" This is a rough coverage test of the python wrapper

it's intended to be shallow, but broad

"""
import RDConfig
import os,sys,tempfile
import unittest
import DataStructs
import Chem

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
    assert tbl

    assert feq(tbl.GetAtomicWeight(6),12.011)
    assert feq(tbl.GetAtomicWeight("C"),12.011)
    assert tbl.GetAtomicNumber('C')==6
    assert feq(tbl.GetRvdw(6),1.950)
    assert feq(tbl.GetRvdw("C"),1.950)
    assert feq(tbl.GetRcovalent(6),0.680)
    assert feq(tbl.GetRcovalent("C"),0.680)
    assert tbl.GetDefaultValence(6)==4
    assert tbl.GetDefaultValence("C")==4
    assert tuple(tbl.GetValenceList(6))==(4,)
    assert tuple(tbl.GetValenceList("C"))==(4,)
    assert tuple(tbl.GetValenceList(16))==(2,4,6)
    assert tuple(tbl.GetValenceList("S"))==(2,4,6)
    assert tbl.GetNOuterElecs(6)==4
    assert tbl.GetNOuterElecs("C")==4

  def test2Atom(self):
    atom = Chem.Atom(6)
    assert atom
    assert atom.GetAtomicNum()==6
    atom.SetAtomicNum(8)
    assert atom.GetAtomicNum()==8

    atom = Chem.Atom("C")
    assert atom
    assert atom.GetAtomicNum()==6


  def test3Bond(self):
    # No longer relevant, bonds are not constructible from Python
    pass

  def test4Mol(self):
    mol = Chem.Mol()
    assert mol

  def test5Smiles(self):
    mol = Chem.MolFromSmiles('n1ccccc1')
    assert mol
    assert mol.GetNumAtoms()==6
    assert mol.GetNumAtoms(1)==6
    assert mol.GetNumAtoms(0)==11
    at = mol.GetAtomWithIdx(2)
    assert at.GetAtomicNum()==6
    at = mol.GetAtomWithIdx(0)
    assert at.GetAtomicNum()==7
    
  def _test6Bookmarks(self):
    mol = Chem.MolFromSmiles('n1ccccc1')
    assert mol

    assert not mol.HasAtomBookmark(0)
    mol.SetAtomBookmark(mol.GetAtomWithIdx(0),0)
    mol.SetAtomBookmark(mol.GetAtomWithIdx(1),1)
    assert mol.HasAtomBookmark(0)
    assert mol.HasAtomBookmark(1)

    if 1:
      assert not mol.HasBondBookmark(0)
      assert not mol.HasBondBookmark(1)
      mol.SetBondBookmark(mol.GetBondWithIdx(0),0)
      mol.SetBondBookmark(mol.GetBondWithIdx(1),1)
      assert mol.HasBondBookmark(0)
      assert mol.HasBondBookmark(1)


    at = mol.GetAtomWithBookmark(0)
    assert at
    assert at.GetAtomicNum()==7
    mol.ClearAtomBookmark(0)
    assert not mol.HasAtomBookmark(0)
    assert mol.HasAtomBookmark(1)
    mol.ClearAllAtomBookmarks()
    assert not mol.HasAtomBookmark(0)
    assert not mol.HasAtomBookmark(1)
    
    mol.SetAtomBookmark(mol.GetAtomWithIdx(1),1)

    if 1:
      assert mol.HasBondBookmark(0)
      assert mol.HasBondBookmark(1)
      bond = mol.GetBondWithBookmark(0)
      assert bond
      mol.ClearBondBookmark(0)
      assert not mol.HasBondBookmark(0)
      assert mol.HasBondBookmark(1)
      mol.ClearAllBondBookmarks()
      assert not mol.HasBondBookmark(0)
      assert not mol.HasBondBookmark(1)

      assert mol.HasAtomBookmark(1)
    
  def test7Atom(self):
    mol = Chem.MolFromSmiles('n1ccccc1C[CH2-]')
    assert mol
    Chem.SanitizeMol(mol)
    a0 = mol.GetAtomWithIdx(0)
    a1 = mol.GetAtomWithIdx(1)
    a6 = mol.GetAtomWithIdx(6)
    a7 = mol.GetAtomWithIdx(7)
    
    assert a0.GetAtomicNum()==7
    assert a0.GetSymbol()=='N'
    assert a0.GetIdx()==0

    aList = [a0,a1,a6,a7]
    assert a0.GetDegree()==2
    assert a1.GetDegree()==2
    assert a6.GetDegree()==2
    assert a7.GetDegree()==1
    assert [x.GetDegree() for x in aList]==[2,2,2,1]

    assert [x.GetTotalNumHs() for x in aList]==[0,1,2,2]
    assert [x.GetNumImplicitHs() for x in aList]==[0,1,2,0]
    assert [x.GetExplicitValence() for x in aList]==[3,3,2,3]
    assert [x.GetImplicitValence() for x in aList]==[0,1,2,0]
    assert [x.GetFormalCharge() for x in aList]==[0,0,0,-1]
    assert [x.GetNoImplicit() for x in aList]==[0,0,0,1]
    assert [x.GetNumExplicitHs() for x in aList]==[0,0,0,2]
    assert [x.GetIsAromatic() for x in aList]==[1,1,0,0]
    assert [x.GetHybridization() for x in aList]==[Chem.HybridizationType.SP2,Chem.HybridizationType.SP2,
                                                   Chem.HybridizationType.SP3,Chem.HybridizationType.SP3],\
                                                   [x.GetHybridization() for x in aList]
    


  def test8Bond(self):
    mol = Chem.MolFromSmiles('n1ccccc1CC(=O)O')
    assert mol
    Chem.SanitizeMol(mol)
    # note bond numbering is funny because of ring closure
    b0 = mol.GetBondWithIdx(0)
    b6 = mol.GetBondWithIdx(6)
    b7 = mol.GetBondWithIdx(7)
    b8 = mol.GetBondWithIdx(8)

    bList = [b0,b6,b7,b8]
    assert [x.GetBondType() for x in bList] == \
           [Chem.BondType.AROMATIC,Chem.BondType.SINGLE,
            Chem.BondType.DOUBLE,Chem.BondType.SINGLE]
    assert [x.GetIsAromatic() for x in bList] == \
           [1,0,0,0]
    
    assert [x.GetIsConjugated()!=0 for x in bList] == \
           [1,0,1,1],[x.GetIsConjugated()!=0 for x in bList]
    assert [x.GetBeginAtomIdx() for x in bList] == \
           [0,6,7,7],[x.GetBeginAtomIdx() for x in bList]
    assert [x.GetBeginAtom().GetIdx() for x in bList] == \
           [0,6,7,7]
    assert [x.GetEndAtomIdx() for x in bList] == \
           [1,7,8,9]
    assert [x.GetEndAtom().GetIdx() for x in bList] == \
           [1,7,8,9]


  def test9Smarts(self):
    query1 = Chem.MolFromSmarts('C(=O)O')
    assert query1
    query2 = Chem.MolFromSmarts('C(=O)[O,N]')
    assert query2
    query3 = Chem.MolFromSmarts('[$(C(=O)O)]')
    assert query3


    mol = Chem.MolFromSmiles('CCC(=O)O')
    assert mol

    assert mol.HasSubstructMatch(query1)
    assert mol.HasSubstructMatch(query2)
    assert mol.HasSubstructMatch(query3)

    mol = Chem.MolFromSmiles('CCC(=O)N')
    assert mol

    assert not mol.HasSubstructMatch(query1)
    assert mol.HasSubstructMatch(query2)
    assert not mol.HasSubstructMatch(query3)

  def test10Iterators(self):
    mol = Chem.MolFromSmiles('CCOC')
    assert mol

    for atom in mol.GetAtoms():
      assert atom
    ats = mol.GetAtoms()
    try:
      ats[1]
    except:
      ok = 0
    else:
      ok = 1
    assert ok
    try:
      ats[12]
    except IndexError:
      ok = 1
    else:
      ok = 0
    assert ok

    if 0:
      for atom in mol.GetHeteros():
        assert atom
      ats = mol.GetHeteros()
      try:
        ats[0]
      except:
        ok = 0
      else:
        ok = 1
      assert ok
      assert ats[0].GetIdx()==2
      try:
        ats[12]
      except IndexError:
        ok = 1
      else:
        ok = 0
      assert ok


    for bond in mol.GetBonds():
      assert bond
    bonds = mol.GetBonds()
    try:
      bonds[1]
    except:
      ok = 0
    else:
      ok = 1
    assert ok
    try:
      bonds[12]
    except IndexError:
      ok = 1
    else:
      ok = 0
    assert ok
      
    if 0:
      mol = Chem.MolFromSmiles('c1ccccc1C')
      for atom in mol.GetAromaticAtoms():
        assert atom
      ats = mol.GetAromaticAtoms()
      try:
        ats[0]
      except:
        ok = 0
      else:
        ok = 1
      assert ok
      assert ats[0].GetIdx()==0
      assert ats[1].GetIdx()==1
      assert ats[2].GetIdx()==2
      try:
        ats[12]
      except IndexError:
        ok = 1
      else:
        ok = 0
      assert ok




  def test11MolOps(self) :
    mol = Chem.MolFromSmiles('C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=CC=CC=C3')
    assert mol
    smi = Chem.MolToSmiles(mol)
    Chem.SanitizeMol(mol)
    nr = Chem.GetSymmSSSR(mol)

    assert (len(nr) == 3)
    
  def test12Smarts(self):
    query1 = Chem.MolFromSmarts('C(=O)O')
    assert query1
    query2 = Chem.MolFromSmarts('C(=O)[O,N]')
    assert query2
    query3 = Chem.MolFromSmarts('[$(C(=O)O)]')
    assert query3


    mol = Chem.MolFromSmiles('CCC(=O)O')
    assert mol

    assert mol.HasSubstructMatch(query1)
    assert mol.GetSubstructMatch(query1)==(2,3,4)
    assert mol.HasSubstructMatch(query2)
    assert mol.GetSubstructMatch(query2)==(2,3,4)
    assert mol.HasSubstructMatch(query3)
    assert mol.GetSubstructMatch(query3)==(2,)

    mol = Chem.MolFromSmiles('CCC(=O)N')
    assert mol

    assert not mol.HasSubstructMatch(query1)
    assert not mol.GetSubstructMatch(query1)
    assert mol.HasSubstructMatch(query2)
    assert mol.GetSubstructMatch(query2)==(2,3,4)
    assert not mol.HasSubstructMatch(query3)

    mol = Chem.MolFromSmiles('OC(=O)CC(=O)O')
    assert mol
    assert mol.HasSubstructMatch(query1)
    assert mol.GetSubstructMatch(query1)==(1,2,0)
    assert mol.GetSubstructMatches(query1)==((1,2,0),(4,5,6))
    assert mol.HasSubstructMatch(query2)
    assert mol.GetSubstructMatch(query2)==(1,2,0)
    assert mol.GetSubstructMatches(query2)==((1,2,0),(4,5,6))
    assert mol.HasSubstructMatch(query3)
    assert mol.GetSubstructMatches(query3)==((1,),(4,))
    
  def test13Smarts(self):
    # previous smarts problems:
    query = Chem.MolFromSmarts('N(=,-C)')
    assert query
    mol = Chem.MolFromSmiles('N#C')
    assert not mol.HasSubstructMatch(query)
    mol = Chem.MolFromSmiles('N=C')
    assert mol.HasSubstructMatch(query)
    mol = Chem.MolFromSmiles('NC')
    assert mol.HasSubstructMatch(query)

    query = Chem.MolFromSmarts('[Cl,$(O)]')
    mol = Chem.MolFromSmiles('C(=O)O')
    assert len(mol.GetSubstructMatches(query))==2
    mol = Chem.MolFromSmiles('C(=N)N')
    assert len(mol.GetSubstructMatches(query))==0
    
    query = Chem.MolFromSmarts('[$([O,S]-[!$(*=O)])]')
    mol = Chem.MolFromSmiles('CC(S)C(=O)O')
    assert len(mol.GetSubstructMatches(query))==1
    mol = Chem.MolFromSmiles('C(=O)O')
    assert len(mol.GetSubstructMatches(query))==0
    
    
  def test14Hs(self):
    m = Chem.MolFromSmiles('CC(=O)[OH]')
    assert m.GetNumAtoms()==4

    m2 = Chem.AddHs(m,1)
    assert m2.GetNumAtoms()==5
    m2 = Chem.RemoveHs(m2,1)
    assert m2.GetNumAtoms()==5
    m2 = Chem.RemoveHs(m2,0)
    assert m2.GetNumAtoms()==4
    
    m2 = Chem.AddHs(m,0)
    assert m2.GetNumAtoms()==8
    m2 = Chem.RemoveHs(m2,1)
    assert m2.GetNumAtoms()==5
    m2 = Chem.RemoveHs(m2)
    assert m2.GetNumAtoms()==4
    
  def test15Neighbors(self):
    m = Chem.MolFromSmiles('CC(=O)[OH]')
    assert m.GetNumAtoms()==4
    
    a = m.GetAtomWithIdx(1)
    ns = a.GetNeighbors()
    assert len(ns)==3

    bs = a.GetBonds()
    assert len(bs)==3

    for b in bs:
      try:
        a2 = b.GetOtherAtom(a)
      except:
        a2=None
      assert a2
    assert len(bs)==3

  def test16Pickle(self):
    import cPickle
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    pkl = cPickle.dumps(m)
    m2 = cPickle.loads(pkl)
    smi1 = Chem.MolToSmiles(m)
    smi2 = Chem.MolToSmiles(m2)
    assert smi1==smi2
    
  def test16Props(self):
    m = Chem.MolFromSmiles('C1=CN=CC=C1')
    assert not m.HasProp('prop1')
    assert not m.HasProp('prop2')
    assert not m.HasProp('prop2')
    m.SetProp('prop1','foob')
    assert not m.HasProp('prop2')
    assert m.HasProp('prop1')
    assert m.GetProp('prop1')=='foob'
    assert not m.HasProp('propo')
    try:
      m.GetProp('prop2')
    except KeyError:
      ok=1
    else:
      ok=0
    assert ok

    # test computed properties
    m.SetProp('cprop1', 'foo', 1)
    m.SetProp('cprop2', 'foo2', 1)

    m.ClearComputedProps()
    assert not m.HasProp('cprop1')
    assert not m.HasProp('cprop2')

  def test17Kekulize(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    smi = Chem.MolToSmiles(m)
    assert smi=='c1ccccc1'

    Chem.Kekulize(m)
    smi = Chem.MolToSmiles(m)
    assert smi=='c1ccccc1'

    m = Chem.MolFromSmiles('c1ccccc1')
    smi = Chem.MolToSmiles(m)
    assert smi=='c1ccccc1'

    Chem.Kekulize(m,1)
    smi = Chem.MolToSmiles(m)
    assert smi=='C1=CC=CC=C1', smi

  def test18Paths(self):


    m = Chem.MolFromSmiles("C1CC2C1CC2")
    #assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7
    #print Chem.FindAllPathsOfLengthN(m,3,useBonds=0)
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==10,\
           Chem.FindAllPathsOfLengthN(m,2,useBonds=1)
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==14
    

    m = Chem.MolFromSmiles('C1CC1C')
    assert m
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==4
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3,Chem.FindAllPathsOfLengthN(m,3,useBonds=1)
    assert len(Chem.FindAllPathsOfLengthN(m,4,useBonds=1))==1,Chem.FindAllPathsOfLengthN(m,4,useBonds=1)
    assert len(Chem.FindAllPathsOfLengthN(m,5,useBonds=1))==0,Chem.FindAllPathsOfLengthN(m,5,useBonds=1)

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==4
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3
    
    m = Chem.MolFromSmiles("CCC(C)CC")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==4,Chem.FindAllPathsOfLengthN(m,3,useBonds=1)
    
    m = Chem.MolFromSmiles("CCCC(C)C")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3
    
    m = Chem.MolFromSmiles("CC(C)C(C)C")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==6
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==4
    
    m = Chem.MolFromSmiles("CC(C)(C)CC")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==7
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3,Chem.FindAllPathsOfLengthN(m,3,useBonds=1)
    
    m = Chem.MolFromSmiles("C1CCCCC1")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==6
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==6
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==6
    
    m = Chem.MolFromSmiles("C1CC2C1CC2")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==10,\
           Chem.FindAllPathsOfLengthN(m,2,useBonds=1)
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==14
    
    
    m = Chem.MolFromSmiles("CC2C1CCC12")
    assert len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7
    assert len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==11
    # FIX: this result disagrees with the paper (which says 13),
    #   but it seems right
    assert len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==15,\
           Chem.FindAllPathsOfLengthN(m,3,useBonds=1)
    
    
    
  def test19Subgraphs(self):
    m = Chem.MolFromSmiles('C1CC1C')
    assert m
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1,0))==4
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==4
    assert len(Chem.FindAllSubgraphsOfLengthN(m,4))==1
    assert len(Chem.FindAllSubgraphsOfLengthN(m,5))==0

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==4
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==3
    
    m = Chem.MolFromSmiles("CCC(C)CC")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==5
    
    m = Chem.MolFromSmiles("CCCC(C)C")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==4
    
    m = Chem.MolFromSmiles("CC(C)C(C)C")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==6
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==6
    
    m = Chem.MolFromSmiles("CC(C)(C)CC")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==5
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==7
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==7,Chem.FindAllSubgraphsOfLengthN(m,3)
    
    m = Chem.MolFromSmiles("C1CCCCC1")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==6
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==6
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==6
    #assert len(Chem.FindUniqueSubgraphsOfLengthN(m,1))==1
    assert len(Chem.FindUniqueSubgraphsOfLengthN(m,2))==1
    assert len(Chem.FindUniqueSubgraphsOfLengthN(m,3))==1
    
    m = Chem.MolFromSmiles("C1CC2C1CC2")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==7
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==10
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==16
    

    m = Chem.MolFromSmiles("CC2C1CCC12")
    assert len(Chem.FindAllSubgraphsOfLengthN(m,1))==7
    assert len(Chem.FindAllSubgraphsOfLengthN(m,2))==11
    assert len(Chem.FindAllSubgraphsOfLengthN(m,3))==18,\
           len(Chem.FindAllSubgraphsOfLengthN(m,3))
    
    
  def test20IsInRing(self):
    m = Chem.MolFromSmiles('C1CCC1C')
    assert m
    assert m.GetAtomWithIdx(0).IsInRingSize(4)
    assert m.GetAtomWithIdx(1).IsInRingSize(4)
    assert m.GetAtomWithIdx(2).IsInRingSize(4)
    assert m.GetAtomWithIdx(3).IsInRingSize(4)
    assert not m.GetAtomWithIdx(4).IsInRingSize(4)
    
    assert not m.GetAtomWithIdx(0).IsInRingSize(3)
    assert not m.GetAtomWithIdx(1).IsInRingSize(3)
    assert not m.GetAtomWithIdx(2).IsInRingSize(3)
    assert not m.GetAtomWithIdx(3).IsInRingSize(3)
    assert not m.GetAtomWithIdx(4).IsInRingSize(3)
    
    assert m.GetBondWithIdx(0).IsInRingSize(4)
    assert not m.GetBondWithIdx(3).IsInRingSize(4)
    assert not m.GetBondWithIdx(0).IsInRingSize(3)
    assert not m.GetBondWithIdx(3).IsInRingSize(3)
    
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
    assert ok  

  def test22DeleteSubstruct(self) :
    query = Chem.MolFromSmarts('C(=O)O')
    mol = Chem.MolFromSmiles('CCC(=O)O')
    nmol = Chem.DeleteSubstructs(mol, query)
    
    assert Chem.MolToSmiles(nmol) == 'CC'

    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    # now delete only fragments
    nmol = Chem.DeleteSubstructs(mol, query, 1)
    assert Chem.MolToSmiles(nmol) == 'CCC(O)=O',Chem.MolToSmiles(nmol)
    
    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = Chem.DeleteSubstructs(mol, query, 0)
    assert Chem.MolToSmiles(nmol) == 'CC'
    
    mol = Chem.MolFromSmiles('CCCO')
    nmol = Chem.DeleteSubstructs(mol, query, 0)
    assert Chem.MolToSmiles(nmol) == 'CCCO'

    # Issue 96 prevented this from working:
    mol = Chem.MolFromSmiles('CCC(=O)O.O=CO')
    nmol = Chem.DeleteSubstructs(mol, query, 1)
    assert Chem.MolToSmiles(nmol) == 'CCC(O)=O'
    nmol = Chem.DeleteSubstructs(nmol, query, 1)
    assert Chem.MolToSmiles(nmol) == 'CCC(O)=O'
    nmol = Chem.DeleteSubstructs(nmol, query, 0)
    assert Chem.MolToSmiles(nmol) == 'CC'

    
  def test23MolFileParsing(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','triazine.mol')
    #fileN = "../FileParsers/test_data/triazine.mol"
    inD = open(fileN,'r').read()
    m1 = Chem.MolFromMolBlock(inD)
    assert m1 is not None
    assert m1.GetNumAtoms()==9

    m1 = Chem.MolFromMolFile(fileN)
    assert m1 is not None
    assert m1.GetNumAtoms()==9

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','triazine.mof')
    m1 = Chem.MolFromMolFile(fileN)
    assert m1 is None

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','list-query.mol')
    #fileN = "../FileParsers/test_data/list-query.mol"
    query = Chem.MolFromMolFile(fileN)
    smi = Chem.MolToSmiles(query)
    assert smi=='[Du]1=CC=CC=C1'
    smi = "C1=CC=CC=C1"
    mol = Chem.MolFromSmiles(smi,0)
    assert mol.HasSubstructMatch(query)
    Chem.SanitizeMol(mol)
    assert not mol.HasSubstructMatch(query)

    mol = Chem.MolFromSmiles('N1=CC=CC=C1',0)
    assert mol.HasSubstructMatch(query)
    mol = Chem.MolFromSmiles('S1=CC=CC=C1',0)
    assert not mol.HasSubstructMatch(query)
    mol = Chem.MolFromSmiles('P1=CC=CC=C1',0)
    assert mol.HasSubstructMatch(query)

  def test24DaylightFingerprint(self):
    import DataStructs
    m1 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp1 = Chem.DaylightFingerprint(m1)
    m2 = Chem.MolFromSmiles('C1=CC=CC=C1')
    fp2 = Chem.DaylightFingerprint(m2)

    tmp = DataStructs.TanimotoSimilarity(fp1,fp2)
    assert tmp==1.0,tmp

    m2 = Chem.MolFromSmiles('C1=CC=CC=N1')
    fp2 = Chem.DaylightFingerprint(m2)
    tmp = DataStructs.TanimotoSimilarity(fp1,fp2)
    assert tmp<1.0,tmp
    assert tmp>0.0,tmp

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
      assert mol
      assert mol.GetProp("_Name") == molNames[i]
      i += 1
      if (mol.GetProp("_Name") == "192") :
        # test parsed charges on one of the molecules
        for id in chgs192.keys() :
          assert mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id]
          
    sdSup.reset()
    
    ns = [mol.GetProp("_Name") for mol in sdSup]
    assert ns == molNames

    sdSup = Chem.SDMolSupplier(fileN, 0)
    for mol in sdSup :
      assert not mol.HasProp("numArom")

    sdSup = Chem.SDMolSupplier(fileN)
    assert len(sdSup) == 16
    mol = sdSup[5]
    assert mol.GetProp("_Name") == "170"

  def test26SmiMolSupplier(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','first_200.tpsa.csv')
    #fileN = "../FileParsers/test_data/first_200.tpsa.csv"
    smiSup = Chem.SmilesMolSupplier(fileN, ",", 0, -1)
    mol = smiSup[16];
    assert mol.GetProp("TPSA") == "46.25"

    mol = smiSup[8];
    assert mol.GetProp("TPSA") == "65.18"

    assert len(smiSup) == 200

    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','fewSmi.csv')
    #fileN = "../FileParsers/test_data/fewSmi.csv"
    smiSup = Chem.SmilesMolSupplier(fileN, delimiter=",",
                                      smilesColumn=1, nameColumn=0,
                                      titleLine=0)
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    i = 0
    for mol in smiSup:
      assert mol.GetProp("_Name") == names[i]
      i += 1
      
    mol = smiSup[3]
    
    assert mol.GetProp("_Name") == "4"
    assert mol.GetProp("Column_2") == "82.78"

    # and test doing a supplier from text:
    inD = open(fileN,'r').read()
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=1, nameColumn=0,
                   titleLine=0)
    names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
    i = 0
    # iteration interface:
    for mol in smiSup:
      assert mol.GetProp("_Name") == names[i]
      i += 1
    assert i==10
    # random access:
    mol = smiSup[3]
    assert len(smiSup)==10
    assert mol.GetProp("_Name") == "4"
    assert mol.GetProp("Column_2") == "82.78"

    # issue 113:
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=1, nameColumn=0,
                   titleLine=0)
    assert len(smiSup)==10

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
    assert len(smiSup)==4
    # but the 3rd is a None:
    assert smiSup[2] is None


    text="Id SMILES Column_2\n"+"mol-1 C 1.0\n"+"mol-2 CC 4.0\n"+"mol-4 CCCC 16.0"
    smiSup.SetData(text, delimiter=" ",
                   smilesColumn=1, nameColumn=0,
                   titleLine=1)
    assert len(smiSup)==3
    m = [x for x in smiSup]
    assert smiSup[2]
    assert len(m)==3
    assert m[0].GetProp("Column_2")=="1.0"
    
    # test simple parsing and Issue 114:
    smis = ['CC','CCC','CCOC','CCCOCC','CCCOCCC']
    inD = '\n'.join(smis)
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    assert len(smiSup)==5
    m = [x for x in smiSup]
    assert smiSup[4]
    assert len(m)==5
    
    # order dependance:
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    assert smiSup[4]
    assert len(smiSup)==5

    # this was a nasty BC:
    # asking for a particular entry with a higher index than what we've
    # already seen resulted in a duplicate:
    smis = ['CC','CCC','CCOC','CCCCOC']
    inD = '\n'.join(smis)
    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    m = smiSup.next()
    m = smiSup[3]
    assert len(smiSup)==4

    try:
      smiSup[4]
    except:
      ok=1
    else:
      ok=0
    assert ok

    smiSup.SetData(inD, delimiter=",",
                   smilesColumn=0, nameColumn=-1,
                   titleLine=0)
    try:
      smiSup[4]
    except:
      ok=1
    else:
      ok=0
    assert ok  
    sys.stderr.write('>>> This may result in an infinite loop.  It should finish almost instantly\n')
    assert len(smiSup)==4
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
    ofile = "test_data/outSmiles.txt"
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
      #print [repr(x) for x in mol.GetPropNames()]
      assert mol.GetProp("_Name") == names[i]
      assert mol.GetProp("Column_2") == props[i]
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
      #print 'mol:',mol
      #print '\t',molNames[i]
      assert mol.GetProp("_Name") == molNames[i]
      i += 1
      if (mol.GetProp("_Name") == "192") :
        # test parsed charges on one of the molecules
        for id in chgs192.keys() :
          assert mol.GetAtomWithIdx(id).GetFormalCharge() == chgs192[id]

  def test30Issues109and110(self) :
    """ issues 110 and 109 were both related to handling of explicit Hs in
       SMILES input.
       
    """ 
    m1 = Chem.MolFromSmiles('N12[C@H](SC(C)(C)[C@@H]1C(O)=O)[C@@H](C2=O)NC(=O)[C@H](N)c3ccccc3')
    assert m1.GetNumAtoms()==24
    m2 = Chem.MolFromSmiles('C1C=C([C@H](N)C(=O)N[C@@]2([H])[C@]3([H])SC(C)(C)[C@@H](C(=O)O)N3C(=O)2)C=CC=1')
    assert m2.GetNumAtoms()==24

    smi1 = Chem.MolToSmiles(m1)
    smi2 = Chem.MolToSmiles(m2)
    assert smi1==smi2

    m1 = Chem.MolFromSmiles('[H]CCl')
    assert m1.GetNumAtoms()==2
    assert m1.GetAtomWithIdx(0).GetNumExplicitHs()==0
    m1 = Chem.MolFromSmiles('[H][CH2]Cl')
    assert m1.GetNumAtoms()==2
    assert m1.GetAtomWithIdx(0).GetNumExplicitHs()==3
    m2 = Chem.AddHs(m1)
    assert m2.GetNumAtoms()==5
    m2 = Chem.RemoveHs(m2)
    assert m2.GetNumAtoms()==2
    
  def test31ChiralitySmiles(self) :
    m1 = Chem.MolFromSmiles('F[C@](Br)(I)Cl')
    assert m1 is not None
    assert m1.GetNumAtoms()==5
    assert Chem.MolToSmiles(m1,1)=='F[C@](Cl)(Br)I',Chem.MolToSmiles(m1,1)

    m1 = Chem.MolFromSmiles('CC1C[C@@]1(Cl)F')
    assert m1 is not None
    assert m1.GetNumAtoms()==6
    assert Chem.MolToSmiles(m1,1)=='CC1C[C@]1(F)Cl',Chem.MolToSmiles(m1,1)

    m1 = Chem.MolFromSmiles('CC1C[C@]1(Cl)F')
    assert m1 is not None
    assert m1.GetNumAtoms()==6
    assert Chem.MolToSmiles(m1,1)=='CC1C[C@@]1(F)Cl',Chem.MolToSmiles(m1,1)


    

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
    assert m1 is not None
    assert m1.GetNumAtoms()==5
    assert smi=='F[C@](Cl)(Br)I',smi

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
    assert m1 is not None
    assert m1.GetNumAtoms()==5
    assert Chem.MolToSmiles(m1,1)=='F[C@@](Cl)(Br)I'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==5
    assert Chem.MolToSmiles(m1,1)=='F[C@](Cl)(Br)I'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==5
    assert Chem.MolToSmiles(m1,1)=='F[C@](Cl)(Br)I'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==4
    assert Chem.MolToSmiles(m1,1)=='F[C@H](Cl)Br'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==4
    assert Chem.MolToSmiles(m1,1)=='FN(Cl)Br'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==4
    assert Chem.MolToSmiles(m1,1)=='F[N@H+](Cl)Br'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==4
    assert Chem.MolToSmiles(m1,1)=='F[C@H](Cl)Br'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==4
    assert Chem.MolToSmiles(m1,1)=='F[C@H](Cl)Br'

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
    assert m1 is not None
    assert m1.GetNumAtoms()==6
    assert Chem.MolToSmiles(m1,1)=='CC1C[C@@]1(F)Cl',Chem.MolToSmiles(m1,1)

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
    assert m1 is not None
    assert m1.GetNumAtoms()==6
    assert Chem.MolToSmiles(m1,1)=='CC1C[C@]1(F)Cl',Chem.MolToSmiles(m1,1)

  def test33Issue65(self) :
    """ issue 65 relates to handling of [H] in SMARTS
       
    """ 
    m1 = Chem.MolFromSmiles('OC(O)(O)O')
    m2 = Chem.MolFromSmiles('OC(O)O')
    m3 = Chem.MolFromSmiles('OCO')
    q1 = Chem.MolFromSmarts('OC[H]',1)
    q2 = Chem.MolFromSmarts('O[C;H1]',1)
    q3 = Chem.MolFromSmarts('O[C;H1][H]',1)

    assert not m1.HasSubstructMatch(q1)
    assert not m1.HasSubstructMatch(q2)
    assert not m1.HasSubstructMatch(q3)

    assert m2.HasSubstructMatch(q1)
    assert m2.HasSubstructMatch(q2)
    assert m2.HasSubstructMatch(q3)

    assert m3.HasSubstructMatch(q1)
    assert not m3.HasSubstructMatch(q2)
    assert not m3.HasSubstructMatch(q3)

    m1H = Chem.AddHs(m1)
    m2H = Chem.AddHs(m2)
    m3H = Chem.AddHs(m3)
    q1 = Chem.MolFromSmarts('OC[H]')
    q2 = Chem.MolFromSmarts('O[C;H1]')
    q3 = Chem.MolFromSmarts('O[C;H1][H]')

    assert not m1H.HasSubstructMatch(q1)
    assert not m1H.HasSubstructMatch(q2)
    assert not m1H.HasSubstructMatch(q3)

    #m2H.Debug()
    assert m2H.HasSubstructMatch(q1)
    assert m2H.HasSubstructMatch(q2)
    assert m2H.HasSubstructMatch(q3)

    assert m3H.HasSubstructMatch(q1)
    assert not m3H.HasSubstructMatch(q2)
    assert not m3H.HasSubstructMatch(q3)



  def test34Issue124(self) :
    """ issue 124 relates to calculation of the distance matrix
       
    """ 
    m = Chem.MolFromSmiles('CC=C')
    d = Chem.GetDistanceMatrix(m,0)
    assert feq(d[0,1],1.0)
    assert feq(d[0,2],2.0)
    # force an update:
    d = Chem.GetDistanceMatrix(m,1,0,1)
    assert feq(d[0,1],1.0)
    assert feq(d[0,2],1.5)

  def test35ChiralityPerception(self) :
    """ Test perception of chirality and CIP encoding       
    """ 
    m = Chem.MolFromSmiles('F[C@]([C@])(Cl)Br')
    Chem.AssignAtomChiralCodes(m,1)
    self.failUnless(m.GetAtomWithIdx(1).HasProp('_CIPCode'))
    self.failIf(m.GetAtomWithIdx(2).HasProp('_CIPCode'))
    
    m = Chem.MolFromSmiles('F[C@H](C)C')
    Chem.AssignAtomChiralCodes(m,1)
    self.failUnless(m.GetAtomWithIdx(1).GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED)
    self.failIf(m.GetAtomWithIdx(1).HasProp('_CIPCode'))

    m = Chem.MolFromSmiles('F\\C=C/Cl')
    self.failUnless(m.GetBondWithIdx(0).GetStereo()==Chem.BondStereo.STEREONONE)
    self.failUnless(m.GetBondWithIdx(1).GetStereo()==Chem.BondStereo.STEREOZ)
    atoms = m.GetBondWithIdx(1).GetStereoAtoms()
    self.failUnless(0 in atoms)
    self.failUnless(3 in atoms)
    self.failUnless(m.GetBondWithIdx(2).GetStereo()==Chem.BondStereo.STEREONONE)
    
    m = Chem.MolFromSmiles('F\\C=CCl')
    self.failUnless(m.GetBondWithIdx(1).GetStereo()==Chem.BondStereo.STEREONONE)

  def test36SubstructMatchStr(self):
    """ test the _SubstructMatchStr function """
    query = Chem.MolFromSmarts('[n,p]1ccccc1')
    assert query
    mol = Chem.MolFromSmiles('N1=CC=CC=C1')
    assert mol.HasSubstructMatch(query)
    assert Chem._HasSubstructMatchStr(mol.ToBinary(),query)
    mol = Chem.MolFromSmiles('S1=CC=CC=C1')
    assert not Chem._HasSubstructMatchStr(mol.ToBinary(),query)
    assert not mol.HasSubstructMatch(query)
    mol = Chem.MolFromSmiles('P1=CC=CC=C1')
    assert mol.HasSubstructMatch(query)
    assert Chem._HasSubstructMatchStr(mol.ToBinary(),query)

    
  def test37SanitException(self):
    mol = Chem.MolFromSmiles('CC(C)(C)(C)C',0)
    assert mol
    self.failUnlessRaises(ValueError,lambda:Chem.SanitizeMol(mol))

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
      self.failUnless(mol)
      self.failUnless(mol.GetNumAtoms())
      self.failUnless(mol.HasProp("CAS"))
      self.failUnless(mol.HasProp("_Name"))
      self.failUnless(mol.GetProp("CAS")==mol.GetProp("_Name"))
      self.failUnless(mol.GetNumConformers()==0)
      i+=1
    self.failUnless(i==4)
    self.failUnless(len(suppl)==4)

  def test38Issue266(self):
    """ test issue 266: generation of kekulized smiles"""
    mol = Chem.MolFromSmiles('c1ccccc1')
    Chem.Kekulize(mol)
    smi = Chem.MolToSmiles(mol)
    self.failUnless(smi=='c1ccccc1')
    smi = Chem.MolToSmiles(mol,kekuleSmiles=True)
    self.failUnless(smi=='C1=CC=CC=C1')
    
  def test39Issue273(self):
    """ test issue 273: MolFileComments and MolFileInfo props ending up in SD files

    """
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','Wrap',
                                            'test_data','outNCI_few.sdf')
    suppl = Chem.SDMolSupplier(fileN)
    ms = [x for x in suppl]
    for m in ms:
      self.failUnless(m.HasProp('_MolFileInfo'))
      self.failUnless(m.HasProp('_MolFileComments'))
    fName = tempfile.mktemp('.sdf')
    w = Chem.SDWriter(fName)
    w.SetProps(ms[0].GetPropNames())
    for m in ms:
      w.write(m)
    w = None

    txt= file(fName,'r').read()
    os.unlink(fName)
    self.failUnless(txt.find('MolFileInfo')==-1)
    self.failUnless(txt.find('MolFileComments')==-1)
    

  def test40SmilesRootedAtAtom(self):
    """ test the rootAtAtom functionality

    """
    smi = 'CN(C)C'
    m = Chem.MolFromSmiles(smi)

    self.failUnless(Chem.MolToSmiles(m)=='CN(C)C')
    self.failUnless(Chem.MolToSmiles(m,rootedAtAtom=1)=='N(C)(C)C')
    

  def test41SetStreamIndices(self) :
    fileN = os.path.join(RDConfig.RDBaseDir,'Code','GraphMol','FileParsers',
                                            'test_data','NCI_aids_few.sdf')
    sdSup = Chem.SDMolSupplier(fileN)
    molNames = ["48", "78", "128", "163", "164", "170", "180", "186",
            "192", "203", "210", "211", "213", "220", "229", "256"]
    indices=[0, 2136, 6198, 8520, 11070, 12292, 14025, 15313, 17313, 20125, 22231,
	     23729, 26147, 28331, 32541, 33991]
    sdSup._SetStreamIndices(indices)
    assert len(sdSup) == 16
    mol = sdSup[5]
    assert mol.GetProp("_Name") == "170"
    
    i = 0
    for mol in sdSup :
      assert mol
      assert mol.GetProp("_Name") == molNames[i]
      i += 1
          
    ns = [mol.GetProp("_Name") for mol in sdSup]
    assert ns == molNames

    # this can also be used to skip molecules in the file:
    indices=[0, 6198, 12292]
    sdSup._SetStreamIndices(indices)
    assert len(sdSup) == 3
    mol = sdSup[2]
    assert mol.GetProp("_Name") == "170"

    # or to reorder them:
    indices=[0, 12292, 6198]
    sdSup._SetStreamIndices(indices)
    assert len(sdSup) == 3
    mol = sdSup[1]
    assert mol.GetProp("_Name") == "170"


if __name__ == '__main__':
  unittest.main()


