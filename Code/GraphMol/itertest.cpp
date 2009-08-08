// $Id$
//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
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




int main(){
  RDLog::InitLogs();
  test1();
  test2();
  test5();
  test6();
  testIssue263();

  
  return 0;
}

