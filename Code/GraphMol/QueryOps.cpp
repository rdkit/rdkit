// $Id$
//
// Copyright (C) 2003-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "QueryOps.h"
#include <algorithm>
#include <RDGeneral/types.h>

namespace RDKit{

  // common general queries


//! returns a Query for matching atoms with a particular number of ring bonds
ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what) {
  ATOM_EQUALS_QUERY *res=new AtomRingQuery(what);
  res->setDescription("AtomRingBondCount");
  res->setDataFunc(queryAtomRingBondCount);
  return res;
};

#ifndef OLDWIN32
template <int sz>  
int queryIsAtomInRingOfSize(Atom const * at) {
  return at->getOwningMol().getRingInfo()->isAtomInRingOfSize(at->getIdx(),sz);
};

template <int sz>  
int queryIsBondInRingOfSize(Bond const * at) {
  return at->getOwningMol().getRingInfo()->isBondInRingOfSize(at->getIdx(),sz);
};

ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,20);
  ATOM_EQUALS_QUERY *res = new ATOM_EQUALS_QUERY;
  res->setVal(true);
  switch(tgt){
  case 3:
    res->setDataFunc(queryIsAtomInRingOfSize<3>);break;
  case 4:
    res->setDataFunc(queryIsAtomInRingOfSize<4>);break;
  case 5:
    res->setDataFunc(queryIsAtomInRingOfSize<5>);break;
  case 6:
    res->setDataFunc(queryIsAtomInRingOfSize<6>);break;
  case 7:
    res->setDataFunc(queryIsAtomInRingOfSize<7>);break;
  case 8:
    res->setDataFunc(queryIsAtomInRingOfSize<8>);break;
  case 9:
    res->setDataFunc(queryIsAtomInRingOfSize<9>);break;
  case 10:
    res->setDataFunc(queryIsAtomInRingOfSize<10>);break;
  case 11:
    res->setDataFunc(queryIsAtomInRingOfSize<11>);break;
  case 12:
    res->setDataFunc(queryIsAtomInRingOfSize<12>);break;
  case 13:
    res->setDataFunc(queryIsAtomInRingOfSize<13>);break;
  case 14:
    res->setDataFunc(queryIsAtomInRingOfSize<14>);break;
  case 15:
    res->setDataFunc(queryIsAtomInRingOfSize<15>);break;
  case 16:
    res->setDataFunc(queryIsAtomInRingOfSize<16>);break;
  case 17:
    res->setDataFunc(queryIsAtomInRingOfSize<17>);break;
  case 18:
    res->setDataFunc(queryIsAtomInRingOfSize<18>);break;
  case 19:
    res->setDataFunc(queryIsAtomInRingOfSize<19>);break;
  case 20:
    res->setDataFunc(queryIsAtomInRingOfSize<20>);break;
  }
  res->setDescription("AtomRingSize");
  return res;
}
BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,15);
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  switch(tgt){
  case 3:
    res->setDataFunc(queryIsBondInRingOfSize<3>);break;
  case 4:
    res->setDataFunc(queryIsBondInRingOfSize<4>);break;
  case 5:
    res->setDataFunc(queryIsBondInRingOfSize<5>);break;
  case 6:
    res->setDataFunc(queryIsBondInRingOfSize<6>);break;
  case 7:
    res->setDataFunc(queryIsBondInRingOfSize<7>);break;
  case 8:
    res->setDataFunc(queryIsBondInRingOfSize<8>);break;
  case 9:
    res->setDataFunc(queryIsBondInRingOfSize<9>);break;
  case 10:
    res->setDataFunc(queryIsBondInRingOfSize<10>);break;
  case 11:
    res->setDataFunc(queryIsBondInRingOfSize<11>);break;
  case 12:
    res->setDataFunc(queryIsBondInRingOfSize<12>);break;
  case 13:
    res->setDataFunc(queryIsBondInRingOfSize<13>);break;
  case 14:
    res->setDataFunc(queryIsBondInRingOfSize<14>);break;
  case 15:
    res->setDataFunc(queryIsBondInRingOfSize<15>);break;
  case 16:
    res->setDataFunc(queryIsBondInRingOfSize<16>);break;
  case 17:
    res->setDataFunc(queryIsBondInRingOfSize<17>);break;
  case 18:
    res->setDataFunc(queryIsBondInRingOfSize<18>);break;
  case 19:
    res->setDataFunc(queryIsBondInRingOfSize<19>);break;
  case 20:
    res->setDataFunc(queryIsBondInRingOfSize<20>);break;
  }
  res->setDescription("BondRingSize");
  return res;
}
#else
  // YA MSVC++ bug
template <typename T>
int queryIsInRingOfSize3(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),3)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize4(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),4)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize5(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),5)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize6(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),6)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize7(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),7)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize8(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),8)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize9(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),9)!=tmp.end();
  }
};

template <typename T>
int queryIsInRingOfSize10(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),10)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize11(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),11)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize12(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),12)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize13(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),13)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize14(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),14)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize15(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),15)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize16(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),16)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize17(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),17)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize18(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),18)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize19(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),19)!=tmp.end();
  }
};
template <typename T>
int queryIsInRingOfSize20(T const * at) {
  if(!(at->hasProp("ringMembership")))
    return 0;
  else{
    INT_VECT tmp;
    at->getProp("ringMembership",tmp);
    return std::find(tmp.begin(),tmp.end(),20)!=tmp.end();
  }
};

ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,20);
  ATOM_EQUALS_QUERY *res = new ATOM_EQUALS_QUERY;
  res->setVal(true);
  switch(tgt){
  case 3:
    res->setDataFunc(queryIsInRingOfSize3<Atom>);break;
  case 4:
    res->setDataFunc(queryIsInRingOfSize4<Atom>);break;
  case 5:
    res->setDataFunc(queryIsInRingOfSize5<Atom>);break;
  case 6:
    res->setDataFunc(queryIsInRingOfSize6<Atom>);break;
  case 7:
    res->setDataFunc(queryIsInRingOfSize7<Atom>);break;
  case 8:
    res->setDataFunc(queryIsInRingOfSize8<Atom>);break;
  case 9:
    res->setDataFunc(queryIsInRingOfSize9<Atom>);break;
  case 10:
    res->setDataFunc(queryIsInRingOfSize10<Atom>);break;
  case 11:
    res->setDataFunc(queryIsInRingOfSize11<Atom>);break;
  case 12:
    res->setDataFunc(queryIsInRingOfSize12<Atom>);break;
  case 13:
    res->setDataFunc(queryIsInRingOfSize13<Atom>);break;
  case 14:
    res->setDataFunc(queryIsInRingOfSize14<Atom>);break;
  case 15:
    res->setDataFunc(queryIsInRingOfSize15<Atom>);break;
  case 16:
    res->setDataFunc(queryIsInRingOfSize16<Atom>);break;
  case 17:
    res->setDataFunc(queryIsInRingOfSize17<Atom>);break;
  case 18:
    res->setDataFunc(queryIsInRingOfSize18<Atom>);break;
  case 19:
    res->setDataFunc(queryIsInRingOfSize19<Atom>);break;
  case 20:
    res->setDataFunc(queryIsInRingOfSize20<Atom>);break;
  }
  res->setDescription("AtomRingSize");
  return res;
}
BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,20);
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  switch(tgt){
  case 3:
    res->setDataFunc(queryIsInRingOfSize3<Bond>);break;
  case 4:
    res->setDataFunc(queryIsInRingOfSize4<Bond>);break;
  case 5:
    res->setDataFunc(queryIsInRingOfSize5<Bond>);break;
  case 6:
    res->setDataFunc(queryIsInRingOfSize6<Bond>);break;
  case 7:
    res->setDataFunc(queryIsInRingOfSize7<Bond>);break;
  case 8:
    res->setDataFunc(queryIsInRingOfSize8<Bond>);break;
  case 9:
    res->setDataFunc(queryIsInRingOfSize9<Bond>);break;
  case 10:
    res->setDataFunc(queryIsInRingOfSize10<Bond>);break;
  case 11:
    res->setDataFunc(queryIsInRingOfSize11<Bond>);break;
  case 12:
    res->setDataFunc(queryIsInRingOfSize12<Bond>);break;
  case 13:
    res->setDataFunc(queryIsInRingOfSize13<Bond>);break;
  case 14:
    res->setDataFunc(queryIsInRingOfSize14<Bond>);break;
  case 15:
    res->setDataFunc(queryIsInRingOfSize15<Bond>);break;
  case 16:
    res->setDataFunc(queryIsInRingOfSize16<Bond>);break;
  case 17:
    res->setDataFunc(queryIsInRingOfSize17<Bond>);break;
  case 18:
    res->setDataFunc(queryIsInRingOfSize18<Bond>);break;
  case 19:
    res->setDataFunc(queryIsInRingOfSize19<Bond>);break;
  case 20:
    res->setDataFunc(queryIsInRingOfSize20<Bond>);break;
  }
  res->setDescription("BondRingSize");
  return res;
}

#endif

ATOM_EQUALS_QUERY *makeAtomSimpleQuery(int what,int func(Atom const *)){
  ATOM_EQUALS_QUERY *res = new ATOM_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(func);
  res->setDescription("AtomSimple");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what){
  ATOM_EQUALS_QUERY *res = makeAtomSimpleQuery(what,queryAtomImplicitValence);
  res->setDescription("AtomImplicitValence");
  return res;
}
  
ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomTotalValence);
  res->setDescription("AtomTotalValence");
  return res;
}
  
ATOM_EQUALS_QUERY *makeAtomNumEqualsQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomNum);
  res->setDescription("AtomAtomicNum");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomExplicitDegree);
  res->setDescription("AtomExplicitDegree");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomTotalDegree);
  res->setDescription("AtomTotalDegree");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomHCount);
  res->setDescription("AtomHCount");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomAromaticQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(true,queryAtomAromatic);
  res->setDescription("AtomIsAromatic");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomAliphaticQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(true,queryAtomAliphatic);
  res->setDescription("AtomIsAliphatic");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomMassQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomMass);
  res->setDescription("AtomMass");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomFormalCharge);
  res->setDescription("AtomFormalCharge");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(what,queryAtomHybridization);
  res->setDescription("AtomHybridization");
  return res;
}

  
ATOM_EQUALS_QUERY *makeAtomInRingQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery(true,queryIsAtomInRing);
  res->setDescription("AtomInRing");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what){
  ATOM_EQUALS_QUERY *res;
  res = makeAtomSimpleQuery(what,queryIsAtomInNRings);
  res->setDescription("AtomInNRings");
  return res;
}



BOND_EQUALS_QUERY *makeBondOrderEqualsQuery(Bond::BondType what){
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryBondOrder);
  res->setDescription("BondOrder");
  return res;
}

BOND_EQUALS_QUERY *makeBondDirEqualsQuery(Bond::BondDir what){
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryBondDir);
  res->setDescription("BondDir");
  return res;
}

BOND_EQUALS_QUERY *makeBondIsInRingQuery(){
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryIsBondInRing);
  res->setDescription("BondInRing");
  return res;
}

BOND_EQUALS_QUERY *makeBondInNRingsQuery(int what){
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryIsBondInNRings);
  res->setDescription("BondInNRings");
  return res;
}


template <typename T>
int nullDataFun(T arg) { return 1; } 
bool nullQueryFun(int arg) { return true; } 

BOND_NULL_QUERY *makeBondNullQuery(){
  BOND_NULL_QUERY *res = new BOND_NULL_QUERY;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("BondNull");
  return res;
}

ATOM_NULL_QUERY *makeAtomNullQuery(){
  ATOM_NULL_QUERY *res = new ATOM_NULL_QUERY;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("AtomNull");
  return res;
}
  
};
