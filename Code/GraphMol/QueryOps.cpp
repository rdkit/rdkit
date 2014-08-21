// $Id$
//
// Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
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

ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,20);
  ATOM_EQUALS_QUERY *res = new ATOM_EQUALS_QUERY;
  res->setVal(tgt);
  switch(tgt){
  case 3:
    res->setDataFunc(queryAtomIsInRingOfSize<3>);break;
  case 4:
    res->setDataFunc(queryAtomIsInRingOfSize<4>);break;
  case 5:
    res->setDataFunc(queryAtomIsInRingOfSize<5>);break;
  case 6:
    res->setDataFunc(queryAtomIsInRingOfSize<6>);break;
  case 7:
    res->setDataFunc(queryAtomIsInRingOfSize<7>);break;
  case 8:
    res->setDataFunc(queryAtomIsInRingOfSize<8>);break;
  case 9:
    res->setDataFunc(queryAtomIsInRingOfSize<9>);break;
  case 10:
    res->setDataFunc(queryAtomIsInRingOfSize<10>);break;
  case 11:
    res->setDataFunc(queryAtomIsInRingOfSize<11>);break;
  case 12:
    res->setDataFunc(queryAtomIsInRingOfSize<12>);break;
  case 13:
    res->setDataFunc(queryAtomIsInRingOfSize<13>);break;
  case 14:
    res->setDataFunc(queryAtomIsInRingOfSize<14>);break;
  case 15:
    res->setDataFunc(queryAtomIsInRingOfSize<15>);break;
  case 16:
    res->setDataFunc(queryAtomIsInRingOfSize<16>);break;
  case 17:
    res->setDataFunc(queryAtomIsInRingOfSize<17>);break;
  case 18:
    res->setDataFunc(queryAtomIsInRingOfSize<18>);break;
  case 19:
    res->setDataFunc(queryAtomIsInRingOfSize<19>);break;
  case 20:
    res->setDataFunc(queryAtomIsInRingOfSize<20>);break;
  }

  res->setDescription("AtomRingSize");
  return res;
}
BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,20);
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(tgt);
  switch(tgt){
  case 3:
    res->setDataFunc(queryBondIsInRingOfSize<3>);break;
  case 4:
    res->setDataFunc(queryBondIsInRingOfSize<4>);break;
  case 5:
    res->setDataFunc(queryBondIsInRingOfSize<5>);break;
  case 6:
    res->setDataFunc(queryBondIsInRingOfSize<6>);break;
  case 7:
    res->setDataFunc(queryBondIsInRingOfSize<7>);break;
  case 8:
    res->setDataFunc(queryBondIsInRingOfSize<8>);break;
  case 9:
    res->setDataFunc(queryBondIsInRingOfSize<9>);break;
  case 10:
    res->setDataFunc(queryBondIsInRingOfSize<10>);break;
  case 11:
    res->setDataFunc(queryBondIsInRingOfSize<11>);break;
  case 12:
    res->setDataFunc(queryBondIsInRingOfSize<12>);break;
  case 13:
    res->setDataFunc(queryBondIsInRingOfSize<13>);break;
  case 14:
    res->setDataFunc(queryBondIsInRingOfSize<14>);break;
  case 15:
    res->setDataFunc(queryBondIsInRingOfSize<15>);break;
  case 16:
    res->setDataFunc(queryBondIsInRingOfSize<16>);break;
  case 17:
    res->setDataFunc(queryBondIsInRingOfSize<17>);break;
  case 18:
    res->setDataFunc(queryBondIsInRingOfSize<18>);break;
  case 19:
    res->setDataFunc(queryBondIsInRingOfSize<19>);break;
  case 20:
    res->setDataFunc(queryBondIsInRingOfSize<20>);break;
  }
  res->setDescription("BondRingSize");
  return res;
}


ATOM_EQUALS_QUERY *makeAtomMinRingSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,20);
  ATOM_EQUALS_QUERY *res = new ATOM_EQUALS_QUERY;
  res->setVal(tgt);
  res->setDataFunc(queryAtomMinRingSize);
  res->setDescription("AtomMinRingSize");
  return res;
}
BOND_EQUALS_QUERY *makeBondMinRingSizeQuery(int tgt){
  RANGE_CHECK(3,tgt,20);
  BOND_EQUALS_QUERY *res = new BOND_EQUALS_QUERY;
  res->setVal(tgt);
  res->setDataFunc(queryBondMinRingSize);
  res->setDescription("BondMinRingSize");
  return res;
}

unsigned int queryAtomBondProduct(Atom const * at) {
  ROMol::OEDGE_ITER beg,end;
  boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
  unsigned int prod=1;
  while(beg!=end){
    prod *= static_cast<unsigned int>(firstThousandPrimes[at->getOwningMol()[*beg]->getBondType()]);
    ++beg;
  }
  return prod;
}
unsigned int queryAtomAllBondProduct(Atom const * at) {
  ROMol::OEDGE_ITER beg,end;

  boost::tie(beg,end) = at->getOwningMol().getAtomBonds(at);
  unsigned int prod=1;
  while(beg!=end){
    prod *= static_cast<unsigned int>(firstThousandPrimes[at->getOwningMol()[*beg]->getBondType()]);
    ++beg;
  }
  for(unsigned int i=0;i<at->getTotalNumHs();i++){
    prod *= static_cast<unsigned int>(firstThousandPrimes[Bond::SINGLE]);
  }
  return prod;
}



ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what){
  ATOM_EQUALS_QUERY *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomImplicitValence);
  res->setDescription("AtomImplicitValence");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomExplicitValenceQuery(int what){
  ATOM_EQUALS_QUERY *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomExplicitValence);
  res->setDescription("AtomExplicitValence");
  return res;
}
  
ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomTotalValence);
  res->setDescription("AtomTotalValence");
  return res;
}
  
ATOM_EQUALS_QUERY *makeAtomNumQuery(int what){
  return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomNum,"AtomAtomicNum");
}

ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomExplicitDegree);
  res->setDescription("AtomExplicitDegree");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomTotalDegree);
  res->setDescription("AtomTotalDegree");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomHCount);
  res->setDescription("AtomHCount");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomImplicitHCountQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomImplicitHCount);
  res->setDescription("AtomImplicitHCount");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHasImplicitHQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true,queryAtomHasImplicitH);
  res->setDescription("AtomHasImplicitH");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomAromaticQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true,queryAtomAromatic);
  res->setDescription("AtomIsAromatic");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomAliphaticQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true,queryAtomAliphatic);
  res->setDescription("AtomIsAliphatic");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true,queryAtomUnsaturated);
  res->setDescription("AtomUnsaturated");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomMassQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(massIntegerConversionFactor*what,
                                             queryAtomMass);
  res->setDescription("AtomMass");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,
                                             queryAtomIsotope);
  res->setDescription("AtomIsotope");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomFormalCharge);
  res->setDescription("AtomFormalCharge");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryAtomHybridization);
  res->setDescription("AtomHybridization");
  return res;
}
  
ATOM_EQUALS_QUERY *makeAtomInRingQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true,queryIsAtomInRing);
  res->setDescription("AtomInRing");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHasRingBondQuery(){
  ATOM_EQUALS_QUERY *res=makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true,queryAtomHasRingBond);
  res->setDescription("AtomHasRingBond");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what){
  ATOM_EQUALS_QUERY *res;
  res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what,queryIsAtomInNRings);
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
