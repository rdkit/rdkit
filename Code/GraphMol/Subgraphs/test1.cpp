// $Id$
//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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



  PATH_LIST tmp;

  //
  //  Retrieve using bonds
  //
  tmp = findAllPathsOfLengthN(mol,3);
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


// -------------------------------------------------------------------
int main()
{
  testSubgraphs();
  testSubgraphs2();
  testPaths();
  testPaths2();
  testUniqueSubgraphs();
  testUniqueSubgraphs2();
  //testLeak();
  return 0;
}
