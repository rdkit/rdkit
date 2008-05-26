// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
//

// std bits
#include <iostream>

// RD bits
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include "SubstructMatch.h"
#include "SubstructUtils.h"

using namespace RDKit;

void test1(){
  std::cout << " ----------------- Test 1" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  unsigned int n;

  RWMol *m,*q1;
  m = new RWMol();
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(8));
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);

  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6));
  q1->addAtom(new QueryAtom(6));
  q1->addBond(0,1,Bond::SINGLE);
  n = SubstructMatch(*m,*q1,matches,false);
  CHECK_INVARIANT(n==2,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");
  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");

  CHECK_INVARIANT(SubstructMatch(*m,*q1,matchV),"");
  CHECK_INVARIANT(matchV.size()==2,"");
  std::cout << "Done\n" << std::endl;
}

void test2(){
  std::cout << " ----------------- Test 2" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  unsigned int n;

  RWMol *m,*q1;
  m = new RWMol();
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(8));
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);

  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6));
  q1->addAtom(new QueryAtom(8));
  q1->addBond(0,1,Bond::SINGLE);
  n = SubstructMatch(*m,*q1,matches,false);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");
  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");

  CHECK_INVARIANT(SubstructMatch(*m,*q1,matchV),"");
  CHECK_INVARIANT(matchV.size()==2,"");

  delete m;
  m = new RWMol();
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(8));
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::DOUBLE);

  matches.clear();
  n = SubstructMatch(*m,*q1,matches,false);
  CHECK_INVARIANT(n==0,"");
  CHECK_INVARIANT(matches.size()==n,"");
  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==0,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(!SubstructMatch(*m,*q1,matchV),"");

  std::cout << "Done\n" << std::endl;

}

void test3(){
  std::cout << " ----------------- Test 3" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  unsigned int n;

  RWMol *m,*q1;
  m = new RWMol();
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(8));
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);

  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6));
  q1->addAtom(new QueryAtom(8));
  q1->addBond(0,1,Bond::UNSPECIFIED);
  n = SubstructMatch(*m,*q1,matches,false);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");
  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");

  CHECK_INVARIANT(SubstructMatch(*m,*q1,matchV),"");
  CHECK_INVARIANT(matchV.size()==2,"");

  delete m;
  m = new RWMol();
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(8));
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::DOUBLE);

  matches.clear();
  n = SubstructMatch(*m,*q1,matches,false);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");
  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");

  CHECK_INVARIANT(SubstructMatch(*m,*q1,matchV),"");
  CHECK_INVARIANT(matchV.size()==2,"");

  std::cout << "Done\n" << std::endl;
}




void test4(){
  std::cout << " ----------------- Test 4" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  int n;

  RWMol *m,*q1,*q2;
  Atom *a6 = new Atom(6);
  Atom *a8 = new Atom(8);
  m = new RWMol();
  m->addAtom(a6);
  m->addAtom(a6);
  m->addAtom(a8);
  m->addAtom(a6);
  m->addAtom(a6);
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);
  m->addBond(1,3,Bond::SINGLE);
  m->addBond(2,4,Bond::SINGLE);

  // this will be the recursive query
  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(8),true);
  q1->addBond(0,1,Bond::UNSPECIFIED);

  // here's the main query
  q2 = new RWMol();
  QueryAtom *qA = new QueryAtom(6);
  RecursiveStructureQuery *rsq = new RecursiveStructureQuery(q1);
  qA->expandQuery(rsq,Queries::COMPOSITE_AND);
  //std::cout << "post expand: " << qA->getQuery() << std::endl;
  q2->addAtom(qA,true,true);
  //std::cout << "mol: " << q2->getAtomWithIdx(0)->getQuery() << std::endl;
  q2->addAtom(new QueryAtom(6),true,true);
  q2->addBond(0,1,Bond::UNSPECIFIED);

  bool found = SubstructMatch(*m,*q2,matchV);
  CHECK_INVARIANT(found,"");
  CHECK_INVARIANT(matchV.size()==2,"");
  n = SubstructMatch(*m,*q2,matches,true);
  CHECK_INVARIANT(n==2,"");
  CHECK_INVARIANT(matches[0].size()==2,"");


  std::cout << "Done\n" << std::endl;
}

void test5(){
  std::cout << " ----------------- Test 5" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  int n;

  RWMol *m,*q1,*q2;
  Atom *a6 = new Atom(6);
  Atom *a8 = new Atom(8);
  // CC(OC)C
  m = new RWMol();
  m->addAtom(a6);
  m->addAtom(a6);
  m->addAtom(a8);
  m->addAtom(a6);
  m->addAtom(a6);
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);
  m->addBond(1,4,Bond::SINGLE);
  m->addBond(2,3,Bond::SINGLE);

  // this will be the recursive query
  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(8),true);
  q1->addBond(0,1,Bond::UNSPECIFIED);

  // here's the main query
  q2 = new RWMol();
  QueryAtom *qA = new QueryAtom();
  RecursiveStructureQuery *rsq = new RecursiveStructureQuery(q1);
  qA->setQuery(rsq);
  //std::cout << "post expand: " << qA->getQuery() << std::endl;
  q2->addAtom(qA,true,true);
  //std::cout << "mol: " << q2->getAtomWithIdx(0)->getQuery() << std::endl;
  q2->addAtom(new QueryAtom(6),true,true);
  q2->addBond(0,1,Bond::UNSPECIFIED);

  bool found = SubstructMatch(*m,*q2,matchV);
  CHECK_INVARIANT(found,"");
  CHECK_INVARIANT(matchV.size()==2,"");
  n = SubstructMatch(*m,*q2,matches,true);
  CHECK_INVARIANT(n==2,"");
  CHECK_INVARIANT(matches[0].size()==2,"");


  std::cout << "Done\n" << std::endl;
}

void test6(){
  std::cout << " ----------------- Test 6 (Issue71 related)" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  int n;

  RWMol *m,*q1;
  Atom *a6 = new Atom(6);

  m = new RWMol();
  m->addAtom(a6);
  m->addAtom(a6);
  m->addAtom(a6);
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);
  m->addBond(0,2,Bond::SINGLE);

  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(6),true);
  q1->addBond(0,1,Bond::UNSPECIFIED);
  q1->addBond(1,2,Bond::UNSPECIFIED);

  bool found = SubstructMatch(*m,*q1,matchV);
  CHECK_INVARIANT(found,"");
  CHECK_INVARIANT(matchV.size()==3,"");
  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches[0].size()==3,"");

  // close the loop and try again (we should still match)
  q1->addBond(0,2,Bond::UNSPECIFIED);
  found = SubstructMatch(*m,*q1,matchV);
  CHECK_INVARIANT(found,"");
  CHECK_INVARIANT(matchV.size()==3,"");
  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches[0].size()==3,"");

  std::cout << "Done\n" << std::endl;
}

void test7(){
  std::cout << " ----------------- Test 7 (leak check)" << std::endl;
  MatchVectType matchV;
  int n;

  RWMol *m,*q1;
  Atom *a6 = new Atom(6);

  m = new RWMol();
  m->addAtom(a6);
  m->addAtom(a6);
  m->addAtom(a6);
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);
  m->addBond(0,2,Bond::SINGLE);

  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(6),true);
  q1->addBond(0,1,Bond::UNSPECIFIED);
  q1->addBond(1,2,Bond::UNSPECIFIED);


  bool found = SubstructMatch(*m,*q1,matchV);
  CHECK_INVARIANT(found,"");
  CHECK_INVARIANT(matchV.size()==3,"");
  std::vector< MatchVectType > matches;
  for(int i=0;i<300000;i++){
    n = SubstructMatch(*m,*q1,matches,true,true);
    CHECK_INVARIANT(n==1,"");
    CHECK_INVARIANT(matches[0].size()==3,"");
    if(! (i%500) ) std::cout << i << std::endl;
  }
  std::cout << "Done\n" << std::endl;
}

#ifdef CACHE_ARMOLGRAPHS
void test8(){
  std::cout << " ----------------- Test 8 (molgraph cache)" << std::endl;
  MatchVectType matchV;
  int n;

  RWMol *m,*q1;
  Atom *a6 = new Atom(6);

  m = new RWMol();
  m->addAtom(a6);
  m->addAtom(a6);
  m->addAtom(a6);
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);
  m->addBond(0,2,Bond::SINGLE);

  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(6),true);
  q1->addAtom(new QueryAtom(6),true);
  q1->addBond(0,1,Bond::UNSPECIFIED);
  q1->addBond(1,2,Bond::UNSPECIFIED);


  bool found = SubstructMatch(*m,*q1,matchV);
  CHECK_INVARIANT(found,"");
  CHECK_INVARIANT(matchV.size()==3,"");
  std::vector< MatchVectType > matches;
  for(int i=0;i<30000;i++){
    n = SubstructMatch(*m,*q1,matches,true,true);
    CHECK_INVARIANT(n==1,"");
    CHECK_INVARIANT(matches[0].size()==3,"");
    if(! (i%500) ) std::cout << i << std::endl;
  }
  std::cout << "Done\n" << std::endl;
}
#endif

void test9(){
  std::cout << " ----------------- Test 9 (chiral searches)" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  int n;

  RWMol *m,*q1;
  Atom *a6 = new Atom(6);

  m = new RWMol();
  m->addAtom(a6);
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(7));
  m->addAtom(new Atom(8));
  m->addAtom(new Atom(9));
  m->addBond(0,1,Bond::SINGLE);
  m->addBond(0,2,Bond::SINGLE);
  m->addBond(0,3,Bond::SINGLE);
  m->addBond(0,4,Bond::SINGLE);
  m->getAtomWithIdx(0)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
  
  q1 = new RWMol();
  q1->addAtom(a6);
  q1->addAtom(new Atom(6));
  q1->addAtom(new Atom(7));
  q1->addAtom(new Atom(8));
  q1->addAtom(new Atom(9));
  q1->addBond(0,1,Bond::SINGLE);
  q1->addBond(0,2,Bond::SINGLE);
  q1->addBond(0,3,Bond::SINGLE);
  q1->addBond(0,4,Bond::SINGLE);
  q1->getAtomWithIdx(0)->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);

  MolOps::sanitizeMol(*m);
  MolOps::assignAtomChiralCodes(*m);
  MolOps::sanitizeMol(*q1);
  MolOps::assignAtomChiralCodes(*q1);

  bool found;
  // test with default options (no chirality):
  found = SubstructMatch(*m,*q1,matchV);
  TEST_ASSERT(found);
  n = SubstructMatch(*m,*q1,matches,true);
  TEST_ASSERT(n==1);

  // test with chirality
  found = SubstructMatch(*m,*q1,matchV,true,true);
  TEST_ASSERT(!found);
  n = SubstructMatch(*m,*q1,matches,true,true,true);
  TEST_ASSERT(n==0);

  // self matches:
  found = SubstructMatch(*m,*m,matchV,true,true);
  TEST_ASSERT(found);
  n = SubstructMatch(*m,*m,matches,true,true,true);
  TEST_ASSERT(n==1);
  found = SubstructMatch(*q1,*q1,matchV,true,true);
  TEST_ASSERT(found);
  n = SubstructMatch(*q1,*q1,matches,true,true,true);
  TEST_ASSERT(n==1);



  
  std::cout << "Done\n" << std::endl;
}


int main(int argc,char *argv[])
{
  test1();  
  test2();  
  test3();  
  test4();  
  test5();  
  test6();  
  if(argc>1 && !strcmp(argv[1],"-l"))
    test7();  
#ifdef CACHE_ARMOLGRAPHS
  test8();
#endif
  test9();  
  return 0;
}

  
