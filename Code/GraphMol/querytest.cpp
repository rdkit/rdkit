//  $Id: querytest.cpp 4961 2006-02-18 00:14:47Z glandrum $
// 
//   Copyright (C) 2002-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <iostream>

using namespace RDKit;
using namespace std;
typedef class RWMol Mol;


void test1(){
  BOOST_LOG(rdErrorLog) << "---------------------- Test1" << std::endl;
  Mol qM;
  Mol m;

  Atom *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0,1,Bond::SINGLE);
  a = new Atom(8);
  m.addAtom(a);
  delete a;
  m.addBond(1,2,Bond::DOUBLE);
  MolOps::sanitizeMol(m);

  QueryAtom *qA = new QueryAtom(6);
  CHECK_INVARIANT(qA->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(qA->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(2)),"");
  qA->expandQuery(makeAtomImplicitValenceQuery(3));
  CHECK_INVARIANT(qA->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(2)),"");


  qM.addAtom(qA);
  qA = new QueryAtom(6);
  qA->expandQuery(makeAtomNumEqualsQuery(8),Queries::COMPOSITE_OR);
  qM.addAtom(qA);
  qM.addAtom(new QueryAtom(8));
  //Atom::ATOM_SPTR qA(new QueryAtom(6));
  
  QueryBond *qB;
  qB = new QueryBond(Bond::UNSPECIFIED);
  qB->setOwningMol(qM);
  qB->setBeginAtomIdx(0);
  qB->setEndAtomIdx(1);

  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(0)),"");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(1)),"");
  qM.addBond(qB);
  qB = new QueryBond(Bond::DOUBLE);
  qB->setOwningMol(qM);
  qB->setBeginAtomIdx(1);
  qB->setEndAtomIdx(2);
  qM.addBond(qB);


  CHECK_INVARIANT(qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!qM.getAtomWithIdx(0)->Match(m.getAtomWithIdx(2)),"");
  CHECK_INVARIANT(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(qM.getAtomWithIdx(1)->Match(m.getAtomWithIdx(2)),"");
  CHECK_INVARIANT(!qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(qM.getAtomWithIdx(2)->Match(m.getAtomWithIdx(2)),"");


  CHECK_INVARIANT(qM.getBondWithIdx(0)->Match(m.getBondWithIdx(0)),"");
  CHECK_INVARIANT(qM.getBondWithIdx(0)->Match(m.getBondWithIdx(1)),"");
  CHECK_INVARIANT(!qM.getBondWithIdx(1)->Match(m.getBondWithIdx(0)),"");
  CHECK_INVARIANT(qM.getBondWithIdx(1)->Match(m.getBondWithIdx(1)),"");


  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void test2(){
  BOOST_LOG(rdErrorLog) << "---------------------- Test2" << std::endl;
  Mol qM;
  Mol m;

  Atom *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0,1,Bond::TRIPLE);
  m.addBond(1,2,Bond::SINGLE);
  a = new Atom(8);
  m.addAtom(a);
  delete a;
  m.addBond(2,3,Bond::DOUBLE);
  MolOps::sanitizeMol(m);

  QueryBond *qB;
  qB = new QueryBond(Bond::SINGLE);

  CHECK_INVARIANT(!qB->Match(m.getBondWithIdx(0)),"");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(1)),"");
  CHECK_INVARIANT(!qB->Match(m.getBondWithIdx(2)),"");

  BOND_EQUALS_QUERY *newQ = makeBondOrderEqualsQuery(Bond::DOUBLE);
  qB->expandQuery(newQ,Queries::COMPOSITE_OR);
  CHECK_INVARIANT(!qB->Match(m.getBondWithIdx(0)),"");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(1)),"");
  CHECK_INVARIANT(qB->Match(m.getBondWithIdx(2)),"");
  
  
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}



void test3(){
  BOOST_LOG(rdErrorLog) << "---------------------- Test3" << std::endl;
  Mol m;

  Atom *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0,1,Bond::SINGLE);
  m.addBond(1,2,Bond::DOUBLE);
  m.addBond(2,3,Bond::SINGLE);
  m.addBond(3,4,Bond::DOUBLE);
  m.addBond(4,5,Bond::SINGLE);
  m.addBond(5,0,Bond::DOUBLE);
  m.addBond(5,6,Bond::SINGLE);
  MolOps::sanitizeMol(m);
  
  ATOM_EQUALS_QUERY *aeq = makeAtomExplicitDegreeQuery(3);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;
  aeq = makeAtomExplicitDegreeQuery(2);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomHCountQuery(1);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomInNRingsQuery(1);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomInNRingsQuery(0);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomAromaticQuery();
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomAliphaticQuery();
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomInRingQuery();
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  aeq->setNegation(true);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomInRingOfSizeQuery(6);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  aeq->setNegation(true);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  aeq = makeAtomInRingOfSizeQuery(5);
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)),"");
  aeq->setNegation(true);
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;
 
  BOND_EQUALS_QUERY *beq = makeBondIsInRingQuery();
  CHECK_INVARIANT(beq->Match(m.getBondWithIdx(0)),"");
  CHECK_INVARIANT(!beq->Match(m.getBondWithIdx(6)),"");
  CHECK_INVARIANT(beq->Match(m.getBondBetweenAtoms(0,1)),"");
  CHECK_INVARIANT(beq->Match(m.getBondBetweenAtoms(1,0)),"");
  CHECK_INVARIANT(!beq->Match(m.getBondBetweenAtoms(5,6)),"");
  CHECK_INVARIANT(!beq->Match(m.getBondBetweenAtoms(6,5)),"");
  
  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}

void test4(){
  BOOST_LOG(rdErrorLog) << "---------------------- Test4" << std::endl;
  Mol m;

  Atom *a = new Atom(6);
  // we copy in addAtom, so this is safe
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  m.addAtom(a);
  delete a;
  m.addBond(0,1,Bond::SINGLE);
  m.addBond(1,2,Bond::DOUBLE);
  m.addBond(2,3,Bond::SINGLE);
  m.addBond(3,4,Bond::DOUBLE);
  m.addBond(4,5,Bond::SINGLE);
  m.addBond(5,0,Bond::DOUBLE);
  m.addBond(5,6,Bond::SINGLE);
  MolOps::sanitizeMol(m);
  
  ATOM_NULL_QUERY *aeq = makeAtomNullQuery();
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)),"");
  CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)),"");
  delete aeq;

  BOOST_LOG(rdErrorLog) << "Done!" << std::endl;
}


int main(){
  RDLog::InitLogs();
  test1();
  test2();
  test3();
  test4();
  return 0;
}
