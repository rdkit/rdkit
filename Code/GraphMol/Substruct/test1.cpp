// $Id$
//
//  Copyright (C) 2001-2010 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//

// std bits
#include <iostream>

// RD bits
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include "SubstructMatch.h"
#include "SubstructUtils.h"

#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/MolSupplier.h>

using namespace RDKit;

void test1(){
  std::cout << " ----------------- Test 1" << std::endl;
  MatchVectType matchV;
  std::vector< MatchVectType > matches;
  unsigned int n;

  RWMol *m,*q1;
  m = new RWMol();
  m->addAtom(new Atom(8));
  m->addAtom(new Atom(6));
  m->addAtom(new Atom(6));
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
  TEST_ASSERT(matches[0][0].first==0);
  TEST_ASSERT(matches[0][0].second==1||matches[0][0].second==2);
  TEST_ASSERT(matches[0][1].first==1);
  TEST_ASSERT(matches[0][1].second!=matches[0][0].second);

  TEST_ASSERT(matches[1][1].second==1||matches[0][1].second==2);
  TEST_ASSERT(matches[1][0].first==0);
  TEST_ASSERT(matches[1][0].second==1||matches[1][0].second==2);
  TEST_ASSERT(matches[1][0].second!=matches[0][0].second);
  TEST_ASSERT(matches[1][0].second==matches[0][1].second);
  TEST_ASSERT(matches[1][1].first==1);
  TEST_ASSERT(matches[1][1].second!=matches[1][0].second);
  TEST_ASSERT(matches[1][1].second==matches[0][0].second);

  n = SubstructMatch(*m,*q1,matches,true);
  CHECK_INVARIANT(n==1,"");
  CHECK_INVARIANT(matches.size()==n,"");
  CHECK_INVARIANT(matches[0].size()==2,"");
  TEST_ASSERT(matches[0][0].first==0);
  TEST_ASSERT(matches[0][0].second==1||matches[0][0].second==2);
  TEST_ASSERT(matches[0][1].first==1);
  TEST_ASSERT(matches[0][1].second!=matches[0][0].second);
  TEST_ASSERT(matches[0][1].second==1||matches[0][1].second==2);
  
  CHECK_INVARIANT(SubstructMatch(*m,*q1,matchV),"");
  CHECK_INVARIANT(matchV.size()==2,"");

  // make sure we reset the match vectors.
  // build a query we won't match:
  q1->addAtom(new QueryAtom(6));
  q1->addBond(1,2,Bond::SINGLE);
  q1->addAtom(new QueryAtom(6));
  q1->addBond(2,3,Bond::SINGLE);

  TEST_ASSERT(!SubstructMatch(*m,*q1,matchV));
  TEST_ASSERT(matchV.size()==0);

  n = SubstructMatch(*m,*q1,matches,false);
  TEST_ASSERT(n==0);
  TEST_ASSERT(matches.size()==0);
  
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

  n = SubstructMatch(*m,*q1,matchV);
  TEST_ASSERT(n);
  TEST_ASSERT(matchV.size()==2);
  TEST_ASSERT(matchV[0].first==0);
  TEST_ASSERT(matchV[0].second==1);
  TEST_ASSERT(matchV[1].first==1);
  TEST_ASSERT(matchV[1].second==2);

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


  delete q1;
  q1 = new RWMol();
  q1->addAtom(new QueryAtom(6));
  q1->addAtom(new QueryAtom(6));
  q1->addBond(0,1,Bond::UNSPECIFIED);
  n = SubstructMatch(*m,*q1,matches,false);
  TEST_ASSERT(n==2);
  TEST_ASSERT(matches.size()==n);
  TEST_ASSERT(matches[0].size()==2);
  TEST_ASSERT(matches[1].size()==2);
  TEST_ASSERT(matches[0][0].second!=matches[1][0].second);
  TEST_ASSERT(matches[0][1].second!=matches[1][1].second);
  n = SubstructMatch(*m,*q1,matches,true);
  TEST_ASSERT(n==1);
  TEST_ASSERT(matches.size()==n);

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
  m->addBond(1,0,Bond::SINGLE);
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
  TEST_ASSERT(matchV[0].first==0);
  TEST_ASSERT(matchV[0].second==1);
  TEST_ASSERT(matchV[1].first==1);
  TEST_ASSERT(matchV[1].second==0||matchV[1].second==3);
  n = SubstructMatch(*m,*q2,matches,true);
  TEST_ASSERT(n==2);
  TEST_ASSERT(matches.size()==n);
  TEST_ASSERT(matches[0].size()==2);
  TEST_ASSERT(matches[1].size()==2);
  TEST_ASSERT(matches[0][0].second==matches[1][0].second);
  TEST_ASSERT(matches[0][1].second!=matches[1][1].second);

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
  q2->addAtom(qA,true,true);
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
void test5QueryRoot(){
  std::cout << " ----------------- Test 5 QueryRoot" << std::endl;
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
  q1->addAtom(new QueryAtom(8),true);
  q1->addAtom(new QueryAtom(6),true);
  q1->addBond(0,1,Bond::UNSPECIFIED);
  q1->setProp("_queryRootAtom",1);
  
  // here's the main query
  q2 = new RWMol();
  QueryAtom *qA = new QueryAtom();
  RecursiveStructureQuery *rsq = new RecursiveStructureQuery(q1);
  qA->setQuery(rsq);
  q2->addAtom(qA,true,true);
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
  MolOps::assignStereochemistry(*m);
  MolOps::sanitizeMol(*q1);
  MolOps::assignStereochemistry(*q1);

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


void testRecursiveSerialNumbers(){
  std::cout << " ----------------- Testing serial numbers on recursive queries" << std::endl;
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
  m->addBond(1,0,Bond::SINGLE);
  m->addBond(1,2,Bond::SINGLE);
  m->addBond(1,3,Bond::SINGLE);
  m->addBond(2,4,Bond::SINGLE);

  {
    // this will be the recursive query
    q1 = new RWMol();
    q1->addAtom(new QueryAtom(6),true);
    q1->addAtom(new QueryAtom(8),true);
    q1->addBond(0,1,Bond::UNSPECIFIED);

    // here's the main query
    q2 = new RWMol();
    QueryAtom *qA = new QueryAtom(6);
    RecursiveStructureQuery *rsq = new RecursiveStructureQuery(new RWMol(*q1),1);
    qA->expandQuery(rsq,Queries::COMPOSITE_AND);
    //std::cout << "post expand: " << qA->getQuery() << std::endl;
    q2->addAtom(qA,true,true);
    //std::cout << "mol: " << q2->getAtomWithIdx(0)->getQuery() << std::endl;
    q2->addAtom(new QueryAtom(8),true,true);
    q2->addBond(0,1,Bond::UNSPECIFIED);

    qA = new QueryAtom(6);
    rsq = new RecursiveStructureQuery(new RWMol(*q1),1);
    qA->expandQuery(rsq,Queries::COMPOSITE_AND);
    q2->addAtom(qA,true,true);
    q2->addBond(1,2,Bond::UNSPECIFIED);

    bool found = SubstructMatch(*m,*q2,matchV);
    CHECK_INVARIANT(found,"");
    CHECK_INVARIANT(matchV.size()==3,"");
    n = SubstructMatch(*m,*q2,matches,true);
    TEST_ASSERT(n==1);
    TEST_ASSERT(matches.size()==1);
    TEST_ASSERT(matches[0].size()==3);

    delete q1;
    delete q2;
  }
  delete m;
  std::cout << "Done\n" << std::endl;
}

#ifdef RDK_TEST_MULTITHREADED
#include <boost/thread.hpp>  
#include <boost/dynamic_bitset.hpp>
namespace {
  void runblock(const std::vector<ROMol *> &mols,const ROMol *query,
                const boost::dynamic_bitset<> &hits,unsigned int count,unsigned int idx){
    for(unsigned int j=0;j<100;j++){
      for(unsigned int i=0;i<mols.size();++i){
        if(i%count != idx) continue;
        ROMol *mol = mols[i];

        MatchVectType matchV;
        bool found=SubstructMatch(*mols[i],*query,matchV);
        
        TEST_ASSERT(found==hits[i]);
      }
    }
  };
}
void testMultiThread(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test multithreading" << std::endl;

  std::string fName = getenv("RDBASE");
  fName += "/Data/NCI/first_200.props.sdf";
  SDMolSupplier suppl(fName);
  std::cerr<<"reading molecules"<<std::endl;
  std::vector<ROMol *> mols;
  while(!suppl.atEnd()&&mols.size()<100){
    ROMol *mol=0;
    try{
      mol=suppl.next();
    } catch(...){
      continue;
    }
    if(!mol) continue;
    mols.push_back(mol);
  }
  boost::thread_group tg;

  ROMol *query=SmartsToMol("[#6;$([#6]([#6])[!#6])]");
  boost::dynamic_bitset<> hits(mols.size());
  for(unsigned int i=0;i<mols.size();++i){
    MatchVectType matchV;
    hits[i]=SubstructMatch(*mols[i],*query,matchV);
  }
  unsigned int count=4;
#if 1 
  std::cerr<<" hits: "<<hits<<std::endl;
  std::cerr<<"processing"<<std::endl;
  for(unsigned int i=0;i<count;++i){
    std::cerr<<" launch :"<<i<<std::endl;std::cerr.flush();
    tg.add_thread(new boost::thread(runblock,mols,query,hits,count,i));
  }
  tg.join_all();
  std::cerr<<" done"<<std::endl;
  delete query;

  query=SmartsToMol("[#6]([#6])[!#6]");
  for(unsigned int i=0;i<mols.size();++i){
    MatchVectType matchV;
    hits[i]=SubstructMatch(*mols[i],*query,matchV);
  }
  std::cerr<<" hits2: "<<hits<<std::endl;
  std::cerr<<"processing2"<<std::endl;
  for(unsigned int i=0;i<count;++i){
    std::cerr<<" launch2 :"<<i<<std::endl;std::cerr.flush();
    tg.add_thread(new boost::thread(runblock,mols,query,hits,count,i));
  }
  tg.join_all();
  std::cerr<<" done"<<std::endl;
  delete query;
#endif
  
  std::cerr<<" preprocessing 3"<<std::endl;
  query=SmartsToMol("[$([O,S]-[!$(*=O)])]");
  for(unsigned int i=0;i<mols.size();++i){
    MatchVectType matchV;
    hits[i]=SubstructMatch(*mols[i],*query,matchV);
  }
  std::cerr<<" hits3: "<<hits<<std::endl;
  std::cerr<<"processing3"<<std::endl;
  for(unsigned int i=0;i<count;++i){
    std::cerr<<" launch3 :"<<i<<std::endl;std::cerr.flush();
    tg.add_thread(new boost::thread(runblock,mols,query,hits,count,i));
  }
  tg.join_all();
  std::cerr<<" done"<<std::endl;
  delete query;


  for(unsigned int i=0;i<mols.size();++i) delete mols[i];

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
#else
void testMultiThread(){
}
#endif

void testChiralMatch(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test chiral matching" << std::endl;

  {
    std::string qSmi="Cl[C@](C)(F)Br";
    std::string mSmi="Cl[C@](C)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="Cl[C@](C)(F)Br";
    std::string mSmi="Cl[C@@](C)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="Cl[C@](C)(F)Br";
    std::string mSmi="Cl[C@@](F)(C)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="Cl[C@](C)(F)Br";
    std::string mSmi="Cl[C@](F)(C)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="Cl[C@](C)(F)Br";
    std::string mSmi="Cl[C@@](Br)(C)F";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="Cl[C@](C)(F)Br";
    std::string mSmi="Cl[C@](Br)(C)F";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }

  {
    std::string qSmi="C[C@](O)(F)Br";
    std::string mSmi="O[C@](F)(Br)CC[C@](O)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }

  {
    std::string qSmi="C[C@](O)(F)Br";
    std::string mSmi="O[C@](F)(Br)CC[C@@](O)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }

  {
    std::string qSmi="C[C@](O)(F)Br";
    std::string mSmi="O[C@](F)(Br)CC(C[C@](O)(F)Br)C[C@](O)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    std::vector< MatchVectType > matches;
    int count=SubstructMatch(*mol,*query,matches,true,true,true);
    TEST_ASSERT(count==2);
  }

  {
    std::string qSmi="C[C@](O)(F)Br";
    std::string mSmi="O[C@@](F)(Br)CC(C[C@](O)(F)Br)C[C@](O)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    std::vector< MatchVectType > matches;
    int count=SubstructMatch(*mol,*query,matches,true,true,true);
    TEST_ASSERT(count==3);
  }
  {
    std::string qSmi="C[C@](O)(F)Br";
    std::string mSmi="O[C@](F)(Br)CC(C[C@@](O)(F)Br)C[C@](O)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    std::vector< MatchVectType > matches;
    int count=SubstructMatch(*mol,*query,matches,true,true,true);
    TEST_ASSERT(count==1);
  }
  {
    std::string qSmi="C[C@](O)(F)Br";
    std::string mSmi="O[C@](F)(Br)CC(C[C@@](O)(F)Br)C[C@@](O)(F)Br";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    std::vector< MatchVectType > matches;
    //std::cerr<<"\n\n------------------------------------------\n"<<qSmi<<" "<<mSmi<<"\n"<<std::endl;
    int count=SubstructMatch(*mol,*query,matches,true,true,true);
    //std::cerr<<"res: "<<count<<std::endl;
    TEST_ASSERT(count==0);
  }

  {
    std::string qSmi="Cl[C@](*)(F)Br";
    std::string mSmi="Cl[C@](C)(F)Br";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="Cl[C@@](*)(F)Br";
    std::string mSmi="Cl[C@](C)(F)Br";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="Cl[C@](*)(*)Br";
    std::string mSmi="Cl[C@](C)(F)Br";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="Cl[C@@](*)(*)Br";
    std::string mSmi="Cl[C@](C)(F)Br";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="[C@](C)(F)Br";
    std::string mSmi="Cl[C@](C)(F)Br";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="[C@@](C)(F)Br";
    std::string mSmi="Cl[C@](C)(F)Br";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }

  
  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testCisTransMatch(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test cis/trans matching" << std::endl;

  {
    std::string qSmi="CC=CC";
    std::string mSmi="CC=CC";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="CC=CC";
    std::string mSmi="C/C=C/C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="CC=CC";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="C/C=C\\C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="C/C=C/C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="C/C=C(/F)C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="C/C=C(\\F)C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="C/C(F)=C(\\F)C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="CC(/F)=C(\\F)C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(matched);
  }
  {
    std::string qSmi="C/C=C/C";
    std::string mSmi="CC(\\F)=C(\\F)C";
    ROMol *query=SmilesToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true);
    TEST_ASSERT(!matched);
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}

void testGitHubIssue15(){
  BOOST_LOG(rdErrorLog) << "-------------------------------------" << std::endl;
  BOOST_LOG(rdErrorLog) << "    Test GitHub issue 15" << std::endl;

  {
    std::string qSmi="[R2]~[R1]~[R2]";
    std::string mSmi="CCC";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmilesToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true,true);
    TEST_ASSERT(!matched);
  }
  {
    std::string qSmi="[R2]~[R1]~[R2]";
    std::string mSmi="CCC";
    ROMol *query=SmartsToMol(qSmi);
    ROMol *mol = SmartsToMol(mSmi);
    MatchVectType matchV;
    bool matched=SubstructMatch(*mol,*query,matchV,true,true,true);
    TEST_ASSERT(!matched);
  }

  BOOST_LOG(rdErrorLog) << "  done" << std::endl;
}
int main(int argc,char *argv[])
{
#if 1
  test1();  
  test2();  
  test3();  
  test4();  
  test5();  
  test5QueryRoot();  
  test6();  
  if(argc>1 && !strcmp(argv[1],"-l"))
    test7();  
  //test9();  
  testRecursiveSerialNumbers();  
  testMultiThread();
  testChiralMatch();
  testCisTransMatch();
#endif
  testGitHubIssue15();
  return 0;
}

  
