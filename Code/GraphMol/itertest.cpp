// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>

#include <iostream>
#include <stdlib.h>
#include <vector>

using namespace RDKit;
using namespace std;

typedef class ROMol Mol;

void test1(){
  string smi="CCOC";
  Mol *m = SmilesToMol(smi);
  Mol::AtomIterator atIt;

  unsigned int idx=0;
  for(atIt=m->beginAtoms();atIt!=m->endAtoms();atIt++){
    CHECK_INVARIANT((*atIt)->getIdx()==idx,"bad idx");
    idx++;
  }
  atIt=m->beginAtoms();
  CHECK_INVARIANT((*(atIt+2))->getIdx()==2,"bad idx");

  atIt=m->beginAtoms();
  Mol::AtomIterator atIt2=m->beginAtoms();
  CHECK_INVARIANT(atIt==atIt2,"iterators don't compare equal");

  atIt++;
  CHECK_INVARIANT((*atIt)->getIdx()==1,"bad idx");
  atIt+=2;
  CHECK_INVARIANT((*atIt)->getIdx()==3,"bad idx");
  atIt-=1;
  CHECK_INVARIANT((*atIt)->getIdx()==2,"bad idx");

  CHECK_INVARIANT(atIt!=atIt2,"iterators don't compare different");
  CHECK_INVARIANT(atIt2<atIt,"iterator inequality failed");
  CHECK_INVARIANT(atIt2<=atIt,"iterator inequality failed");
  CHECK_INVARIANT(atIt>atIt2,"iterator inequality failed");
  CHECK_INVARIANT(atIt>=atIt2,"iterator inequality failed");

  atIt--;
  --atIt;
  CHECK_INVARIANT((*atIt)->getIdx()==0,"bad idx");
  CHECK_INVARIANT(atIt==atIt2,"iterator inequality failed");
  
  atIt++;
  ++atIt;
  CHECK_INVARIANT((*atIt)->getIdx()==2,"bad idx");

  atIt = m->beginAtoms();
  atIt = atIt + 2;
  atIt = atIt - 2;
  CHECK_INVARIANT((*atIt)->getIdx()==0,"bad idx");

  atIt2 = m->beginAtoms();
  atIt2 += 2;
  CHECK_INVARIANT(atIt2-atIt==2,"subtraction failed");
  CHECK_INVARIANT(atIt-atIt2==-2,"subtraction failed");

  // past the end stuff
  atIt2 = m->endAtoms();

  atIt = m->beginAtoms()+10;
  CHECK_INVARIANT(atIt>=atIt2,"past-the-end failed");

  // this is whack
  atIt = m->beginAtoms();
  atIt -= 10;
  CHECK_INVARIANT(atIt>=atIt2,"past-the-end failed");
  
  BOOST_LOG(rdInfoLog)<< "test1 done" << endl;
};

void test2(){
  string smi="C1COC1";
  Mol *m = SmilesToMol(smi);
  Mol::BondIterator bondIt;

  unsigned int idx=0;
  for(bondIt=m->beginBonds();bondIt!=m->endBonds();bondIt++){
    CHECK_INVARIANT((*bondIt)->getIdx()==idx,"bad idx");
    idx++;
  }
  bondIt=m->beginBonds();
  Mol::BondIterator bondIt2=m->beginBonds();
  CHECK_INVARIANT(bondIt==bondIt2,"iterators don't compare equal");

  bondIt++;
  CHECK_INVARIANT((*bondIt)->getIdx()==1,"bad idx");
  bondIt++;
  bondIt++;
  bondIt--;
  CHECK_INVARIANT((*bondIt)->getIdx()==2,"bad idx");

  CHECK_INVARIANT(bondIt!=bondIt2,"iterators don't compare different");

  bondIt--;
  --bondIt;
  CHECK_INVARIANT((*bondIt)->getIdx()==0,"bad idx");
  CHECK_INVARIANT(bondIt==bondIt2,"iterator inequality failed");
  
  bondIt++;
  ++bondIt;
  CHECK_INVARIANT((*bondIt)->getIdx()==2,"bad idx");

  // past the end stuff
  bondIt2 = m->endBonds();
  bondIt = m->beginBonds();
  bondIt--;


  BOOST_LOG(rdInfoLog)<< "test2 done" << endl;
};

void test3(){
  string smi="C1COCCNCOCNSCC1";
  unsigned char heteros[]={2,5,7,9,10};
  
  Mol *m = SmilesToMol(smi);
  {
    unsigned int nSeen=0;
    for(Mol::HeteroatomIterator heteroIt=m->beginHeteros();heteroIt!=m->endHeteros();heteroIt++){
      CHECK_INVARIANT((*heteroIt)->getIdx()==heteros[nSeen],"bad hetero");
      nSeen++;
    }
  }
  {
    unsigned int nSeen=0;
    for(Mol::HeteroatomIterator heteroIt=m->beginHeteros();heteroIt!=m->endHeteros();++heteroIt){
      CHECK_INVARIANT((*heteroIt)->getIdx()==heteros[nSeen],"bad hetero");
      nSeen++;
    }
  }
  {
    Mol::HeteroatomIterator heteroIt = m->beginHeteros();
    heteroIt++;
    heteroIt++;
    heteroIt--;
    CHECK_INVARIANT((*heteroIt)->getIdx()==heteros[1],"bad hetero");
    CHECK_INVARIANT((*--heteroIt)->getIdx()==heteros[0],"bad hetero");
    CHECK_INVARIANT((*heteroIt)->getIdx()==heteros[0],"bad hetero");
  }
  BOOST_LOG(rdInfoLog)<< "test3 done" << endl;
};

void test4(){
  string smi="C1COCCNCOCNSCC1";
  unsigned int heteros1[]={2,7};
  
  Mol *m = SmilesToMol(smi);
  QueryAtom *q= new QueryAtom();
  q->setQuery(makeAtomNumQuery(8));
  {
    unsigned int nSeen=0;
    for(Mol::QueryAtomIterator queryIt=m->beginQueryAtoms(q);queryIt!=m->endQueryAtoms();queryIt++){
      CHECK_INVARIANT((*queryIt)->getIdx()==heteros1[nSeen],"bad query");
      nSeen++;
    }
  }
  {
    Mol::QueryAtomIterator queryIt = m->beginQueryAtoms(q);
    queryIt++;
    queryIt--;
    CHECK_INVARIANT((*queryIt)->getIdx()==heteros1[0],"bad query");
    CHECK_INVARIANT((*++queryIt)->getIdx()==heteros1[1],"bad query");
    CHECK_INVARIANT((*queryIt)->getIdx()==heteros1[1],"bad query");
  }
  {
    Mol::QueryAtomIterator queryIt = m->beginQueryAtoms(q);
    queryIt++;
    queryIt--;
    Mol::QueryAtomIterator queryIt2 = queryIt;
    CHECK_INVARIANT((*queryIt2)->getIdx()==heteros1[0],"bad query");
    CHECK_INVARIANT((*++queryIt2)->getIdx()==heteros1[1],"bad query");
    CHECK_INVARIANT((*queryIt2)->getIdx()==heteros1[1],"bad query");
  }
  smi = "CC(C)CC(C)CC(C)CC(C)C";
  unsigned int heteros2[]={1,4,7,10};
  m = SmilesToMol(smi);
  //m->debugMol(cout);
  q->setQuery(makeAtomImplicitValenceQuery(1));
  {
    unsigned int nSeen = 0;
    for(Mol::QueryAtomIterator queryIt=m->beginQueryAtoms(q);queryIt!=m->endQueryAtoms();++queryIt){
      CHECK_INVARIANT((*queryIt)->getIdx()==heteros2[nSeen],"bad query");
      nSeen++;
    }
  }
  BOOST_LOG(rdInfoLog)<< "test4 done" << endl;
};

void test5(){
  string smi="CCCC";
  Mol *m = SmilesToMol(smi);
  Mol::BondIterator bondIt;
  unsigned int idx=0;
  for(bondIt=m->beginBonds();bondIt!=m->endBonds();bondIt++){
    CHECK_INVARIANT((*bondIt)->getIdx()==idx,"bad idx");
    idx++;
  }
  CHECK_INVARIANT(idx==3,"bad idx");
  idx = 0;
  for(bondIt=m->beginBonds();bondIt!=m->endBonds();bondIt++){
    CHECK_INVARIANT((*bondIt)->getIdx()==idx,"bad idx");
    idx++;
  }
  CHECK_INVARIANT(idx==3,"bad idx");

  idx = 0;
  Mol::BondIterator beginP(m->beginBonds());
  Mol::BondIterator endP(m->endBonds());
  for(bondIt=beginP;bondIt!=endP;bondIt++){
    CHECK_INVARIANT((*bondIt)->getIdx()==idx,"bad idx");
    idx++;
  }
  CHECK_INVARIANT(idx==3,"bad idx");


  BOOST_LOG(rdInfoLog)<< "test5 done" << endl;
}

#if 1
void _test6Help(const ROMol *m){
  Mol::ConstAtomIterator atIt;
  unsigned int idx=0;
  for(atIt=m->beginAtoms();atIt!=m->endAtoms();atIt++){
    CHECK_INVARIANT((*atIt)->getIdx()==idx,"bad idx");
    idx++;
  }
  atIt=m->beginAtoms();
  CHECK_INVARIANT((*(atIt+2))->getIdx()==2,"bad idx");

  atIt=m->beginAtoms();
  Mol::ConstAtomIterator atIt2=m->beginAtoms();
  CHECK_INVARIANT(atIt==atIt2,"iterators don't compare equal");

  atIt++;
  CHECK_INVARIANT((*atIt)->getIdx()==1,"bad idx");
  atIt+=2;
  CHECK_INVARIANT((*atIt)->getIdx()==3,"bad idx");
  atIt-=1;
  CHECK_INVARIANT((*atIt)->getIdx()==2,"bad idx");

  CHECK_INVARIANT(atIt!=atIt2,"iterators don't compare different");
  CHECK_INVARIANT(atIt2<atIt,"iterator inequality failed");
  CHECK_INVARIANT(atIt2<=atIt,"iterator inequality failed");
  CHECK_INVARIANT(atIt>atIt2,"iterator inequality failed");
  CHECK_INVARIANT(atIt>=atIt2,"iterator inequality failed");

  atIt--;
  --atIt;
  CHECK_INVARIANT((*atIt)->getIdx()==0,"bad idx");
  CHECK_INVARIANT(atIt==atIt2,"iterator inequality failed");
  
  atIt++;
  ++atIt;
  CHECK_INVARIANT((*atIt)->getIdx()==2,"bad idx");

  atIt = m->beginAtoms();
  atIt = atIt + 2;
  atIt = atIt - 2;
  CHECK_INVARIANT((*atIt)->getIdx()==0,"bad idx");

  atIt2 = m->beginAtoms();
  atIt2 += 2;
  CHECK_INVARIANT(atIt2-atIt==2,"subtraction failed");
  CHECK_INVARIANT(atIt-atIt2==-2,"subtraction failed");

  // past the end stuff
  atIt2 = m->endAtoms();

  atIt = m->beginAtoms()+10;
  CHECK_INVARIANT(atIt>=atIt2,"past-the-end failed");

  // this is whack
  atIt = m->beginAtoms();
  atIt -= 10;
  CHECK_INVARIANT(atIt>=atIt2,"past-the-end failed");

}
void test6(){
  string smi="CCOC";
  Mol *m = SmilesToMol(smi);
  _test6Help(m);
  
  BOOST_LOG(rdInfoLog)<< "test6 done" << endl;
};
#endif

void test7(){
  string smi="c1ccccc1C";
#if 1
  Mol *m = SmilesToMol(smi);
  Mol::AromaticAtomIterator atomIt;
  Mol::AromaticAtomIterator beginP(m->beginAromaticAtoms());
  Mol::AromaticAtomIterator endP(m->endAromaticAtoms());
  unsigned int idx=0;
  for(atomIt=beginP;atomIt!=endP;atomIt++){
    TEST_ASSERT((*atomIt)->getIdx()==idx);
    idx++;
  }
  TEST_ASSERT(idx==6);

  atomIt = beginP;
  atomIt++;
  atomIt--;
  TEST_ASSERT((*atomIt)->getIdx()==0);

  delete m;
  smi = "Cc1ccccc1";
  m = SmilesToMol(smi);
  beginP = m->beginAromaticAtoms();
  endP = m->endAromaticAtoms();
  idx=0;
  for(atomIt=beginP;atomIt!=endP;atomIt++){
    TEST_ASSERT((*atomIt)->getIdx()==idx+1);
    idx++;
  }
  TEST_ASSERT(idx==6);
#endif  
  BOOST_LOG(rdInfoLog)<< "test7 done" << endl;
}

void testIssue263(){
  string smi="c1ccccc1C";
#if 1
  Mol *m = SmilesToMol(smi);
  Mol::AtomIterator atomIt;
  unsigned int idx=0;
  for(atomIt=m->beginAtoms();atomIt!=m->endAtoms();++atomIt){
    TEST_ASSERT((*atomIt)->getIdx()==idx);
    idx++;
  }
  TEST_ASSERT(idx==7);

  Mol::BondIterator bondIt;
  idx = 0;
  for(bondIt=m->beginBonds();bondIt!=m->endBonds();++bondIt){
    CHECK_INVARIANT((*bondIt)->getIdx()==idx,"bad idx");
    idx++;
  }
  CHECK_INVARIANT(idx==7,"bad idx");

#endif  
  BOOST_LOG(rdInfoLog)<< "testIssue263 done" << endl;
}


void test8(){
  {
    string smi="CC1CC2CC1C2";
    Mol *m = SmilesToMol(smi);
    QueryAtom *q= new QueryAtom();
    q->setQuery(makeAtomExplicitDegreeQuery(3));
    q->expandQuery(makeAtomRingBondCountQuery(2));
    unsigned int nSeen=0;
    for(Mol::QueryAtomIterator queryIt=m->beginQueryAtoms(q);queryIt!=m->endQueryAtoms();++queryIt){
      TEST_ASSERT((*queryIt)->getIdx()==1);
      nSeen++;
    }
    TEST_ASSERT(nSeen=1);
    delete m;
    delete q;
  }

  BOOST_LOG(rdInfoLog)<< "test8 done" << endl;
};


int main(){
  RDLog::InitLogs();
  test1();
  test2();
  test3();
  test4();
  test5();
  test6();
  test7();
  test8();
  testIssue263();

  
  return 0;
}

