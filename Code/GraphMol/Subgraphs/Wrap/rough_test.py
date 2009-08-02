# $Id$
#
#  Copyright (C) 2009 Greg Landrum
#         All Rights Reserved
#
from rdkit import RDConfig
import os,sys,tempfile
import unittest
from rdkit import Chem

def feq(v1,v2,tol2=1e-4):
  return abs(v1-v2)<=tol2

class TestCase(unittest.TestCase):
  def setUp(self):
    pass
  def test19Subgraphs(self):
    m = Chem.MolFromSmiles('C1CC1C')
    self.failUnless(m)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1,0))==4)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==4)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,4))==1)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,5))==0)

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==4)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==3)
    
    m = Chem.MolFromSmiles("CCC(C)CC")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==5)
    
    m = Chem.MolFromSmiles("CCCC(C)C")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==4)
    
    m = Chem.MolFromSmiles("CC(C)C(C)C")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==6)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==6)
    
    m = Chem.MolFromSmiles("CC(C)(C)CC")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==5)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==7)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==7,Chem.FindAllSubgraphsOfLengthN(m,3))
    
    m = Chem.MolFromSmiles("C1CCCCC1")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==6)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==6)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==6)
    #self.failUnless(len(Chem.FindUniqueSubgraphsOfLengthN(m,1))==1)
    self.failUnless(len(Chem.FindUniqueSubgraphsOfLengthN(m,2))==1)
    self.failUnless(len(Chem.FindUniqueSubgraphsOfLengthN(m,3))==1)
    
    m = Chem.MolFromSmiles("C1CC2C1CC2")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==7)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==10)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==16)
    

    m = Chem.MolFromSmiles("CC2C1CCC12")
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,1))==7)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,2))==11)
    self.failUnless(len(Chem.FindAllSubgraphsOfLengthN(m,3))==18,
           len(Chem.FindAllSubgraphsOfLengthN(m,3)))
    
  def test18Paths(self):
    m = Chem.MolFromSmiles("C1CC2C1CC2")
    #self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7)
    #print Chem.FindAllPathsOfLengthN(m,3,useBonds=0)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==10,
           Chem.FindAllPathsOfLengthN(m,2,useBonds=1))
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==14)
    

    m = Chem.MolFromSmiles('C1CC1C')
    self.failUnless(m)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==4)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3,Chem.FindAllPathsOfLengthN(m,3,useBonds=1))
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,4,useBonds=1))==1,Chem.FindAllPathsOfLengthN(m,4,useBonds=1))
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,5,useBonds=1))==0,Chem.FindAllPathsOfLengthN(m,5,useBonds=1))

    #
    #  Hexane example from Hall-Kier Rev.Comp.Chem. paper
    #  Rev. Comp. Chem. vol 2, 367-422, (1991)
    #
    m = Chem.MolFromSmiles("CCCCCC")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==4)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3)
    
    m = Chem.MolFromSmiles("CCC(C)CC")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==4,Chem.FindAllPathsOfLengthN(m,3,useBonds=1))
    
    m = Chem.MolFromSmiles("CCCC(C)C")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3)
    
    m = Chem.MolFromSmiles("CC(C)C(C)C")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==6)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==4)
    
    m = Chem.MolFromSmiles("CC(C)(C)CC")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==5)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==7)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==3,Chem.FindAllPathsOfLengthN(m,3,useBonds=1))
    
    m = Chem.MolFromSmiles("C1CCCCC1")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==6)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==6)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==6)
    
    m = Chem.MolFromSmiles("C1CC2C1CC2")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==10,
           Chem.FindAllPathsOfLengthN(m,2,useBonds=1))
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==14)
    
    
    m = Chem.MolFromSmiles("CC2C1CCC12")
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,1,useBonds=1))==7)
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,2,useBonds=1))==11)
    # FIX: this result disagrees with the paper (which says 13),
    #   but it seems right
    self.failUnless(len(Chem.FindAllPathsOfLengthN(m,3,useBonds=1))==15,
           Chem.FindAllPathsOfLengthN(m,3,useBonds=1))

if __name__ == '__main__':
  unittest.main()
