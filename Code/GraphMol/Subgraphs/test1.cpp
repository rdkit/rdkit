// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <iostream>
using namespace std;
using namespace RDKit;


void testSubgraphs()
{
  std::cout << "-----------------------\n Subgraph retrieval" << std::endl;
  // build: CCC(C)CC
  RWMol mol;
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addBond(0,1,Bond::SINGLE);
  mol.addBond(1,2,Bond::SINGLE);
  mol.addBond(2,3,Bond::SINGLE);
  mol.addBond(2,4,Bond::SINGLE);
  mol.addBond(3,5,Bond::SINGLE);


  PATH_LIST tmp;
  PATH_LIST::iterator i;
  PATH_TYPE::iterator j;

  int totPs = 0;
  tmp = findAllSubgraphsOfLengthN(mol,1);
  CHECK_INVARIANT(tmp.size()==5,"");
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol,2);
  CHECK_INVARIANT(tmp.size()==5,"");
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol,3);
  CHECK_INVARIANT(tmp.size()==5,"");
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==3,"");
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol,5);
  CHECK_INVARIANT(tmp.size()==1,"");
  totPs += tmp.size();
  tmp = findAllSubgraphsOfLengthN(mol,6);
  CHECK_INVARIANT(tmp.size()==0,"");
  totPs += tmp.size();

  // now use the direct range function and check that we get the 
  // same anwswer
  INT_PATH_LIST_MAP tmpm;
  tmpm = findAllSubgraphsOfLengthsMtoN(mol, 1, 6);
  int newTot, idx;
  newTot = 0;
  for (idx = 1; idx <= 6; idx++) {
    newTot += tmpm[idx].size();
  }
  CHECK_INVARIANT(totPs==newTot, "");

  // add an H and make sure things don't change:
  mol.addAtom(new Atom(1));
  mol.addBond(5,6,Bond::SINGLE);
  
  tmp = findAllSubgraphsOfLengthN(mol,1);
  
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findAllSubgraphsOfLengthN(mol,2);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findAllSubgraphsOfLengthN(mol,3);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findAllSubgraphsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==3,"");
  tmp = findAllSubgraphsOfLengthN(mol,5);
  CHECK_INVARIANT(tmp.size()==1,"");
  tmp = findAllSubgraphsOfLengthN(mol,6);
  CHECK_INVARIANT(tmp.size()==0,"");


  
  std::cout << "Finished" << std::endl;
}

void testSubgraphs2a() {
  // these have been moved here from test2.cpp
  std::cout << "-----------------------\n testSubgraphs2a" << std::endl;
  RWMol *mol=SmilesToMol("C1CC2C1N2");
  CHECK_INVARIANT(mol,"");
  int nAll=0;
  int nUnique=0;
  for(unsigned int i=1;i<7;++i){
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol,i);
    nUnique += tmp.size();
  }
  
  CHECK_INVARIANT(nAll==49,"");
  CHECK_INVARIANT(nUnique==21,"");
  delete mol;

  std::cout << "Finished" << std::endl;
}

void testSubgraphs2 () {
  // these have been moved here from test2.cpp
  std::cout << "-----------------------\n testSubgraphs2" << std::endl;
  RWMol *mol=SmilesToMol("C12C3C4C1C1C2C3N41");
  CHECK_INVARIANT(mol,"");
  
  int nAll=0;
  int nUnique=0;
  int i;
  for(i=1;i<13;i++){
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol,i);
    nUnique += tmp.size();
  }
  
  CHECK_INVARIANT(nAll==2433,"");
  CHECK_INVARIANT(nUnique==300,"");
  delete mol;

  mol=SmilesToMol("CCC(O)C(c1ccccc1)CC(C)N(C)C");
  CHECK_INVARIANT(mol,"");

  nAll=0;
  nUnique=0;
  for(i=1;i<18;i++){
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,i);
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol,i);
    nUnique += tmp.size();
  }
  CHECK_INVARIANT(nAll==1990,"");
  
  CHECK_INVARIANT(nUnique==907,"");

  std::cout << "Finished" << std::endl;
}

void dumpVIV(PATH_LIST v){
  PATH_LIST::iterator i;
  PATH_TYPE::iterator j;
  for(i=v.begin();i!=v.end();i++){
    for(j=i->begin();j!=i->end();j++){
      std::cout << *j << " ";
    }
    std::cout << std::endl;
  }

}

void testPaths()
{
  std::cout << "-----------------------\n Path retrieval" << std::endl;
  // build: CCC(C)CC
  RWMol mol;
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addBond(0,1,Bond::SINGLE);
  mol.addBond(1,2,Bond::SINGLE);
  mol.addBond(2,3,Bond::SINGLE);
  mol.addBond(2,4,Bond::SINGLE);
  mol.addBond(3,5,Bond::SINGLE);


  PATH_LIST tmp;

  //
  //  Retrieve using bonds
  //
  tmp = findAllPathsOfLengthN(mol,1);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findAllPathsOfLengthN(mol,2);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findAllPathsOfLengthN(mol,3);
  CHECK_INVARIANT(tmp.size()==4,"");
  tmp = findAllPathsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==1,"");
  tmp = findAllPathsOfLengthN(mol,5);
  CHECK_INVARIANT(tmp.size()==0,"");
  tmp = findAllPathsOfLengthN(mol,6);
  CHECK_INVARIANT(tmp.size()==0,"");
  
  //
  //  Retrieve using atoms, which gives the results shifted by
  //  one (it takes two atoms to make one bond)
  //
  tmp = findAllPathsOfLengthN(mol,1,false);
  CHECK_INVARIANT(tmp.size()==6,"");
  tmp = findAllPathsOfLengthN(mol,2,false);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findAllPathsOfLengthN(mol,3,false);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findAllPathsOfLengthN(mol,4,false);
  CHECK_INVARIANT(tmp.size()==4,"");
  tmp = findAllPathsOfLengthN(mol,5,false);
  CHECK_INVARIANT(tmp.size()==1,"");
  tmp = findAllPathsOfLengthN(mol,6,false);
  CHECK_INVARIANT(tmp.size()==0,"");
  

  //
  //  try m->n
  //
  INT_PATH_LIST_MAP pths;
  pths = findAllPathsOfLengthsMtoN(mol,1,6);
  CHECK_INVARIANT(pths[1].size()==5,"");
  CHECK_INVARIANT(pths[2].size()==5,"");
  CHECK_INVARIANT(pths[3].size()==4,"");
  CHECK_INVARIANT(pths[4].size()==1,"");
  CHECK_INVARIANT(pths[5].size()==0,"");
  CHECK_INVARIANT(pths[6].size()==0,"");


  pths = findAllPathsOfLengthsMtoN(mol,1,6,false);
  CHECK_INVARIANT(pths[1].size()==6,"");
  CHECK_INVARIANT(pths[2].size()==5,"");
  CHECK_INVARIANT(pths[3].size()==5,"");
  CHECK_INVARIANT(pths[4].size()==4,"");
  CHECK_INVARIANT(pths[5].size()==1,"");
  CHECK_INVARIANT(pths[6].size()==0,"");


  //
  //  add an atom, close the ring and re-check a couple indices:
  //   (leaves us with CC1CCCCC1)
  //
  mol.addAtom(new Atom(6));
  mol.addBond(5,6,Bond::SINGLE);
  mol.addBond(0,6,Bond::SINGLE);
  tmp = findAllPathsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==8,"");
  tmp = findAllPathsOfLengthN(mol,5,false);
  CHECK_INVARIANT(tmp.size()==8,"");


  std::cout << "Finished" << std::endl;
}

void testPaths2()
{
  std::cout << "-----------------------\n Path retrieval2" << std::endl;
  // build: CCC(C)CC
  RWMol mol;
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addBond(0,1,Bond::SINGLE);
  mol.addBond(1,2,Bond::SINGLE);
  mol.addBond(2,3,Bond::SINGLE);
  mol.addBond(2,0,Bond::SINGLE);

  //
  //  Retrieve using bonds
  //
  PATH_LIST tmp = findAllPathsOfLengthN(mol,3);
  //std::cout << "\n3:" << std::endl;
  //dumpVIV(tmp);
  CHECK_INVARIANT(tmp.size()==3,"");
  

  std::cout << "Finished" << std::endl;
}

void testUniqueSubgraphs()
{
  std::cout << "-----------------------\n Unique Subgraph retrieval" << std::endl;
  // build: CCC(C)CC
  RWMol mol;
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addAtom(new Atom(6));
  mol.addBond(0,1,Bond::SINGLE);
  mol.addBond(1,2,Bond::SINGLE);
  mol.addBond(2,3,Bond::SINGLE);
  mol.addBond(2,4,Bond::SINGLE);
  mol.addBond(3,5,Bond::SINGLE);


  PATH_LIST tmp;
  PATH_LIST::iterator i;
  PATH_TYPE::iterator j;

  tmp = findAllSubgraphsOfLengthN(mol,1);
  CHECK_INVARIANT(tmp.size()==5,"");

  tmp = findAllSubgraphsOfLengthN(mol,2);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,2);
  CHECK_INVARIANT(tmp.size()==1,"");

  tmp = findAllSubgraphsOfLengthN(mol,3);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,3);
  CHECK_INVARIANT(tmp.size()==2,"");

  tmp = findAllSubgraphsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==3,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==2,"");

  tmp = findAllSubgraphsOfLengthN(mol,5);
  CHECK_INVARIANT(tmp.size()==1,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,5);
  CHECK_INVARIANT(tmp.size()==1,"");

  tmp = findAllSubgraphsOfLengthN(mol,6);
  CHECK_INVARIANT(tmp.size()==0,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,6);
  CHECK_INVARIANT(tmp.size()==0,"");

  // add an H and make sure things don't change:
  mol.addAtom(new Atom(1));
  mol.addBond(5,6,Bond::SINGLE);
  
  tmp = findAllSubgraphsOfLengthN(mol,2);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,2);
  CHECK_INVARIANT(tmp.size()==1,"");

  tmp = findAllSubgraphsOfLengthN(mol,3);
  CHECK_INVARIANT(tmp.size()==5,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,3);
  CHECK_INVARIANT(tmp.size()==2,"");

  tmp = findAllSubgraphsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==3,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,4);
  CHECK_INVARIANT(tmp.size()==2,"");

  tmp = findAllSubgraphsOfLengthN(mol,5);
  CHECK_INVARIANT(tmp.size()==1,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,5);
  CHECK_INVARIANT(tmp.size()==1,"");

  tmp = findAllSubgraphsOfLengthN(mol,6);
  CHECK_INVARIANT(tmp.size()==0,"");
  tmp = findUniqueSubgraphsOfLengthN(mol,6);
  CHECK_INVARIANT(tmp.size()==0,"");
  
  std::cout << "Finished" << std::endl;
}

void testUniqueSubgraphs2() {
  // moved here from test2.cpp
  std::cout << "-----------------------\n testUniqueSubgraphs2" << std::endl;
  RWMol *mol=SmilesToMol("O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))");
  CHECK_INVARIANT(mol,"");

  int nAll=0;
  int nUnique=0;
  for(int i=1;i<26;i++){
    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,i);
    //std::cout << i << "\t" << tmp.size();
    nAll += tmp.size();
    tmp = findUniqueSubgraphsOfLengthN(*mol,i);
    //std::cout << "\t" << tmp.size() << std::endl;;
    nUnique += tmp.size();
  }

  CHECK_INVARIANT(nAll==6435,"");
  CHECK_INVARIANT(nUnique==5618,"");
  delete mol;
  std::cout << "Finished" << std::endl;
}

void testLeak() {
  // moved here from test2.cpp
  std::cout << "-----------------------\n testing for a core leak (Issue 42) " << std::endl;
  RWMol *mol=SmilesToMol("O=C(O)CCCC=CC(C1C(O)CC(O)C1(C=CC(O)CCCCC))");
  CHECK_INVARIANT(mol,"");

  for(int rep=0;rep<100;rep++){
    int nAll=0;
    int nUnique=0;
    for(int i=1;i<26;i++){
      PATH_LIST tmp;
      tmp = findAllSubgraphsOfLengthN(*mol,i);
      //std::cout << i << "\t" << tmp.size();
      nAll += tmp.size();
      tmp = findUniqueSubgraphsOfLengthN(*mol,i);
      //std::cout << "\t" << tmp.size() << std::endl;;
      nUnique += tmp.size();
    }
    CHECK_INVARIANT(nAll==6435,"");
    CHECK_INVARIANT(nUnique==5618,"");
  }
  delete mol;
  std::cout << "Finished" << std::endl;
}

void testRootedSubgraphs () {
  std::cout << "-----------------------\n testRootedSubgraphs" << std::endl;
  {
    RWMol *mol=SmilesToMol("CC1CC1");
    TEST_ASSERT(mol);

    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,1,false,0);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllSubgraphsOfLengthN(*mol,2,false,0);
    TEST_ASSERT(tmp.size()==2);
    tmp = findAllSubgraphsOfLengthN(*mol,3,false,0);
    TEST_ASSERT(tmp.size()==3);
    tmp = findUniqueSubgraphsOfLengthN(*mol,2,false,false,0);
    TEST_ASSERT(tmp.size()==1);
    tmp = findUniqueSubgraphsOfLengthN(*mol,3,false,false,0);
    TEST_ASSERT(tmp.size()==2);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol,1,3,false,0);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==2);  
    TEST_ASSERT(tmpm[3].size()==3);  
  
    // edge case:
    tmp = findAllSubgraphsOfLengthN(*mol,1,false,10);
    TEST_ASSERT(tmp.size()==0);

    delete mol;
  }

  { // tests for sf.net issue 250
    RWMol *mol=SmilesToMol("C1CC1C");
    TEST_ASSERT(mol);

    PATH_LIST tmp;
    tmp = findAllSubgraphsOfLengthN(*mol,1,false,3);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllSubgraphsOfLengthN(*mol,2,false,3);
    TEST_ASSERT(tmp.size()==2);
    tmp = findAllSubgraphsOfLengthN(*mol,3,false,3);
    TEST_ASSERT(tmp.size()==3);
    tmp = findUniqueSubgraphsOfLengthN(*mol,2,false,false,3);
    TEST_ASSERT(tmp.size()==1);
    tmp = findUniqueSubgraphsOfLengthN(*mol,3,false,false,3);
    TEST_ASSERT(tmp.size()==2);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol,1,3,false,3);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==2);  
    TEST_ASSERT(tmpm[3].size()==3);  
  
    // edge case:
    tmp = findAllSubgraphsOfLengthN(*mol,1,false,10);
    TEST_ASSERT(tmp.size()==0);

    delete mol;
  }

  {
    RWMol *mol=SmilesToMol("CC1CC1");
    TEST_ASSERT(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol,1,2,false,0);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==2);
    TEST_ASSERT(tmpm[3].size()==0);    
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol,1,3,false,0);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==2);
    TEST_ASSERT(tmpm[3].size()==3);    
    delete mol;
  }
  { // tests for sf.net issue 250
    RWMol *mol=SmilesToMol("C1CC1C");
    TEST_ASSERT(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol,1,2,false,3);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==2);
    TEST_ASSERT(tmpm[3].size()==0);    
    tmpm = findAllSubgraphsOfLengthsMtoN(*mol,1,3,false,3);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==2);
    TEST_ASSERT(tmpm[3].size()==3);    
    delete mol;
  }

  std::cout << "Finished" << std::endl;
}

void testRootedPaths () {
  std::cout << "-----------------------\n testRootedPaths" << std::endl;
  {
    RWMol *mol=SmilesToMol("CC1CC1");
    TEST_ASSERT(mol);

    PATH_LIST tmp;

    // bond paths:
    tmp = findAllPathsOfLengthN(*mol,1,true,false,0);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllPathsOfLengthN(*mol,2,true,false,0);
    TEST_ASSERT(tmp.size()==2);
    tmp = findAllPathsOfLengthN(*mol,3,true,false,0);
    TEST_ASSERT(tmp.size()==2);

    // edge case:
    tmp = findAllPathsOfLengthN(*mol,1,true,false,10);
    TEST_ASSERT(tmp.size()==0);

    // atom paths:
    tmp = findAllPathsOfLengthN(*mol,1,false,false,0);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllPathsOfLengthN(*mol,2,false,false,0);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllPathsOfLengthN(*mol,3,false,false,0);
    TEST_ASSERT(tmp.size()==2);
    tmp = findAllPathsOfLengthN(*mol,4,false,false,0);
    TEST_ASSERT(tmp.size()==2);
  
    delete mol;
  }
  
  { // tests for sf.net issue 250
    RWMol *mol=SmilesToMol("C1CC1C");
    TEST_ASSERT(mol);

    PATH_LIST tmp;

    // bond paths:
    tmp = findAllPathsOfLengthN(*mol,1,true,false,3);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllPathsOfLengthN(*mol,2,true,false,3);
    TEST_ASSERT(tmp.size()==2);
    tmp = findAllPathsOfLengthN(*mol,3,true,false,3);
    TEST_ASSERT(tmp.size()==2);

    // edge case:
    tmp = findAllPathsOfLengthN(*mol,1,true,false,10);
    TEST_ASSERT(tmp.size()==0);

    // atom paths:
    tmp = findAllPathsOfLengthN(*mol,1,false,false,3);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllPathsOfLengthN(*mol,2,false,false,3);
    TEST_ASSERT(tmp.size()==1);
    tmp = findAllPathsOfLengthN(*mol,3,false,false,3);
    TEST_ASSERT(tmp.size()==2);
    tmp = findAllPathsOfLengthN(*mol,4,false,false,3);
    TEST_ASSERT(tmp.size()==2);
  
    delete mol;
  }
  
  {
    RWMol *mol=SmilesToMol("CC1CC1");
    TEST_ASSERT(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllPathsOfLengthsMtoN(*mol,1,2,false,false,0);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==1);
    TEST_ASSERT(tmpm[3].size()==0);    
    tmpm = findAllPathsOfLengthsMtoN(*mol,1,3,false,false,0);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==1);
    TEST_ASSERT(tmpm[3].size()==2);    
    delete mol;
  }
  { // tests for sf.net issue 250
    RWMol *mol=SmilesToMol("C1CC1C");
    TEST_ASSERT(mol);

    INT_PATH_LIST_MAP tmpm;
    tmpm = findAllPathsOfLengthsMtoN(*mol,1,2,false,false,3);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==1);
    TEST_ASSERT(tmpm[3].size()==0);    
    tmpm = findAllPathsOfLengthsMtoN(*mol,1,3,false,false,3);
    TEST_ASSERT(tmpm[1].size()==1);
    TEST_ASSERT(tmpm[2].size()==1);
    TEST_ASSERT(tmpm[3].size()==2);    
    delete mol;
  }

  std::cout << "Finished" << std::endl;
}



// -------------------------------------------------------------------
int main()
{
#if 1
  testSubgraphs();
  testSubgraphs2a();
  testSubgraphs2();
  testPaths();
  testPaths2();
  testUniqueSubgraphs();
  testUniqueSubgraphs2();
  testRootedSubgraphs();
  testRootedPaths();
#endif
  //testLeak();
  return 0;
}
