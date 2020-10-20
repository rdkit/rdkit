//
// Copyright (C) 2003-2017 Greg Landrum and Rational Discovery LLC
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
#include <GraphMol/QueryAtom.h>

namespace RDKit {

// common general queries

//! returns a Query for matching atoms with a particular number of ring bonds
ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what) {
  ATOM_EQUALS_QUERY *res = new AtomRingQuery(what);
  res->setDescription("AtomRingBondCount");
  res->setDataFunc(queryAtomRingBondCount);
  return res;
};

ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt) {
  RANGE_CHECK(3, tgt, 20);
  auto *res = new ATOM_EQUALS_QUERY;
  res->setVal(tgt);
  switch (tgt) {
    case 3:
      res->setDataFunc(queryAtomIsInRingOfSize<3>);
      break;
    case 4:
      res->setDataFunc(queryAtomIsInRingOfSize<4>);
      break;
    case 5:
      res->setDataFunc(queryAtomIsInRingOfSize<5>);
      break;
    case 6:
      res->setDataFunc(queryAtomIsInRingOfSize<6>);
      break;
    case 7:
      res->setDataFunc(queryAtomIsInRingOfSize<7>);
      break;
    case 8:
      res->setDataFunc(queryAtomIsInRingOfSize<8>);
      break;
    case 9:
      res->setDataFunc(queryAtomIsInRingOfSize<9>);
      break;
    case 10:
      res->setDataFunc(queryAtomIsInRingOfSize<10>);
      break;
    case 11:
      res->setDataFunc(queryAtomIsInRingOfSize<11>);
      break;
    case 12:
      res->setDataFunc(queryAtomIsInRingOfSize<12>);
      break;
    case 13:
      res->setDataFunc(queryAtomIsInRingOfSize<13>);
      break;
    case 14:
      res->setDataFunc(queryAtomIsInRingOfSize<14>);
      break;
    case 15:
      res->setDataFunc(queryAtomIsInRingOfSize<15>);
      break;
    case 16:
      res->setDataFunc(queryAtomIsInRingOfSize<16>);
      break;
    case 17:
      res->setDataFunc(queryAtomIsInRingOfSize<17>);
      break;
    case 18:
      res->setDataFunc(queryAtomIsInRingOfSize<18>);
      break;
    case 19:
      res->setDataFunc(queryAtomIsInRingOfSize<19>);
      break;
    case 20:
      res->setDataFunc(queryAtomIsInRingOfSize<20>);
      break;
  }

  res->setDescription("AtomRingSize");
  return res;
}

BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int tgt) {
  RANGE_CHECK(3, tgt, 20);
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(tgt);
  switch (tgt) {
    case 3:
      res->setDataFunc(queryBondIsInRingOfSize<3>);
      break;
    case 4:
      res->setDataFunc(queryBondIsInRingOfSize<4>);
      break;
    case 5:
      res->setDataFunc(queryBondIsInRingOfSize<5>);
      break;
    case 6:
      res->setDataFunc(queryBondIsInRingOfSize<6>);
      break;
    case 7:
      res->setDataFunc(queryBondIsInRingOfSize<7>);
      break;
    case 8:
      res->setDataFunc(queryBondIsInRingOfSize<8>);
      break;
    case 9:
      res->setDataFunc(queryBondIsInRingOfSize<9>);
      break;
    case 10:
      res->setDataFunc(queryBondIsInRingOfSize<10>);
      break;
    case 11:
      res->setDataFunc(queryBondIsInRingOfSize<11>);
      break;
    case 12:
      res->setDataFunc(queryBondIsInRingOfSize<12>);
      break;
    case 13:
      res->setDataFunc(queryBondIsInRingOfSize<13>);
      break;
    case 14:
      res->setDataFunc(queryBondIsInRingOfSize<14>);
      break;
    case 15:
      res->setDataFunc(queryBondIsInRingOfSize<15>);
      break;
    case 16:
      res->setDataFunc(queryBondIsInRingOfSize<16>);
      break;
    case 17:
      res->setDataFunc(queryBondIsInRingOfSize<17>);
      break;
    case 18:
      res->setDataFunc(queryBondIsInRingOfSize<18>);
      break;
    case 19:
      res->setDataFunc(queryBondIsInRingOfSize<19>);
      break;
    case 20:
      res->setDataFunc(queryBondIsInRingOfSize<20>);
      break;
  }
  res->setDescription("BondRingSize");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomMinRingSizeQuery(int tgt) {
  auto *res = new ATOM_EQUALS_QUERY;
  res->setVal(tgt);
  res->setDataFunc(queryAtomMinRingSize);
  res->setDescription("AtomMinRingSize");
  return res;
}
BOND_EQUALS_QUERY *makeBondMinRingSizeQuery(int tgt) {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(tgt);
  res->setDataFunc(queryBondMinRingSize);
  res->setDescription("BondMinRingSize");
  return res;
}

unsigned int queryAtomBondProduct(Atom const *at) {
  unsigned int prod = 1;
  for(auto *bond : at->bonds()) {
    prod *= static_cast<unsigned int>(
        firstThousandPrimes[bond->getBondType()]);
  }
  return prod;
}
unsigned int queryAtomAllBondProduct(Atom const *at) {
  unsigned int prod = 1;
  for(auto *bond : at->bonds()) {
    prod *= static_cast<unsigned int>(
        firstThousandPrimes[bond->getBondType()]);
  }
  for (unsigned int i = 0; i < at->getTotalNumHs(); i++) {
    prod *= static_cast<unsigned int>(firstThousandPrimes[Bond::SINGLE]);
  }
  return prod;
}

ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomImplicitValence);
  res->setDescription("AtomImplicitValence");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomExplicitValenceQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomExplicitValence);
  res->setDescription("AtomExplicitValence");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomTotalValence);
  res->setDescription("AtomTotalValence");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomNumQuery(int what) {
  return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNum,
                                                "AtomAtomicNum");
}

ATOM_EQUALS_QUERY *makeAtomTypeQuery(int num, int aromatic) {
  return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(makeAtomType(num, aromatic),
                                                queryAtomType, "AtomType");
}
ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomExplicitDegree);
  res->setDescription("AtomExplicitDegree");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomTotalDegree);
  res->setDescription("AtomTotalDegree");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHeavyAtomDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHeavyAtomDegree);
  res->setDescription("AtomHeavyAtomDegree");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHCount);
  res->setDescription("AtomHCount");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomImplicitHCountQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomImplicitHCount);
  res->setDescription("AtomImplicitHCount");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHasImplicitHQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasImplicitH);
  res->setDescription("AtomHasImplicitH");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomAromaticQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomAromatic);
  res->setDescription("AtomIsAromatic");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomAliphaticQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomAliphatic);
  res->setDescription("AtomIsAliphatic");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomUnsaturated);
  res->setDescription("AtomUnsaturated");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomMassQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      massIntegerConversionFactor * what, queryAtomMass);
  res->setDescription("AtomMass");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomIsotope);
  res->setDescription("AtomIsotope");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomFormalCharge);
  res->setDescription("AtomFormalCharge");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomNegativeFormalChargeQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      what, queryAtomNegativeFormalCharge);
  res->setDescription("AtomNegativeFormalCharge");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHybridization);
  res->setDescription("AtomHybridization");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomNumRadicalElectronsQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      what, queryAtomNumRadicalElectrons);
  res->setDescription("AtomNumRadicalElectrons");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHasChiralTagQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasChiralTag);
  res->setDescription("AtomHasChiralTag");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomMissingChiralTagQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomMissingChiralTag);
  res->setDescription("AtomMissingChiralTag");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomInRingQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryIsAtomInRing);
  res->setDescription("AtomInRing");
  return res;
}

ATOM_OR_QUERY *makeQAtomQuery() {
  auto *res = new ATOM_OR_QUERY;
  res->setDescription("AtomOr");
  res->setTypeLabel("Q");
  res->setNegation(true);
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(6)));
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
  return res;
}
ATOM_EQUALS_QUERY *makeQHAtomQuery() {
  ATOM_EQUALS_QUERY *res = makeAtomNumQuery(6);
  res->setNegation(true);
  res->setTypeLabel("QH");
  return res;
}
ATOM_EQUALS_QUERY *makeAAtomQuery() {
  ATOM_EQUALS_QUERY *res = makeAtomNumQuery(1);
  res->setNegation(true);
  res->setTypeLabel("A");
  return res;
}
ATOM_NULL_QUERY *makeAHAtomQuery() {
  auto *res = makeAtomNullQuery();
  res->setTypeLabel("AH");
  return res;
}

ATOM_OR_QUERY *makeXAtomQuery() {
  auto *res = new ATOM_OR_QUERY;
  res->setDescription("AtomOr");
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(9)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(17)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(35)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(53)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(85)));
  res->setTypeLabel("X");

  return res;
}
ATOM_OR_QUERY *makeXHAtomQuery() {
  ATOM_OR_QUERY *res = makeXAtomQuery();
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
  res->setTypeLabel("XH");
  return res;
}

ATOM_OR_QUERY *makeMAtomQuery() {
  // using the definition from Marvin Sketch, which produces the following
  // SMARTS:
  // !#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
  // it's easier to define what isn't a metal than what is. :-)
  ATOM_OR_QUERY *res = makeMHAtomQuery();
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
  res->setTypeLabel("M");

  return res;
}
ATOM_OR_QUERY *makeMHAtomQuery() {
  // using the definition from Marvin Sketch, which produces the following
  // SMARTS:
  // !#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
  // it's easier to define what isn't a metal than what is. :-)
  auto *res = new ATOM_OR_QUERY;
  res->setDescription("AtomOr");
  res->setNegation(true);
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(2)));
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(5)));
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(6)));
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(7)));
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(8)));
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(9)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(10)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(14)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(15)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(16)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(17)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(18)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(33)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(34)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(35)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(36)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(52)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(53)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(54)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(85)));
  res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
      makeAtomNumQuery(86)));
  res->setTypeLabel("MH");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what) {
  ATOM_EQUALS_QUERY *res;
  res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryIsAtomInNRings);
  res->setDescription("AtomInNRings");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHasRingBondQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasRingBond);
  res->setDescription("AtomHasRingBond");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomNumHeteroatomNbrsQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNumHeteroatomNbrs);
  res->setDescription("AtomNumHeteroatomNeighbors");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHasHeteroatomNbrsQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasHeteroatomNbrs);
  res->setDescription("AtomHasHeteroatomNeighbors");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomNumAliphaticHeteroatomNbrsQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      what, queryAtomNumAliphaticHeteroatomNbrs);
  res->setDescription("AtomNumAliphaticHeteroatomNeighbors");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomHasAliphaticHeteroatomNbrsQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      true, queryAtomHasAliphaticHeteroatomNbrs);
  res->setDescription("AtomHasAliphaticHeteroatomNeighbors");
  return res;
}

ATOM_EQUALS_QUERY *makeAtomNonHydrogenDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNonHydrogenDegree);
  res->setDescription("AtomNonHydrogenDegree");
  return res;
}

BOND_EQUALS_QUERY *makeBondOrderEqualsQuery(Bond::BondType what) {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryBondOrder);
  res->setDescription("BondOrder");
  return res;
}

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrAromaticBondQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrAromatic);
  res->setDescription("SingleOrAromaticBond");
  return res;
};

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *
makeSingleOrDoubleOrAromaticBondQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrDoubleOrAromatic);
  res->setDescription("SingleOrDoubleOrAromaticBond");
  return res;
};

BOND_EQUALS_QUERY *makeBondDirEqualsQuery(Bond::BondDir what) {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryBondDir);
  res->setDescription("BondDir");
  return res;
}

BOND_EQUALS_QUERY *makeBondHasStereoQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondHasStereo);
  res->setDescription("BondStereo");
  return res;
}

BOND_EQUALS_QUERY *makeBondIsInRingQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryIsBondInRing);
  res->setDescription("BondInRing");
  return res;
}

BOND_EQUALS_QUERY *makeBondInNRingsQuery(int what) {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryIsBondInNRings);
  res->setDescription("BondInNRings");
  return res;
}

BOND_NULL_QUERY *makeBondNullQuery() {
  auto *res = new BOND_NULL_QUERY;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("BondNull");
  return res;
}

ATOM_NULL_QUERY *makeAtomNullQuery() {
  auto *res = new ATOM_NULL_QUERY;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("AtomNull");
  return res;
}

bool isComplexQuery(const Bond *b) {
  if (!b->hasQuery()) {
    return false;
  }
  // negated things are always complex:
  if (b->getQuery()->getNegation()) {
    return true;
  }
  std::string descr = b->getQuery()->getDescription();
  if (descr == "BondOrder" || descr == "SingleOrAromaticBond") {
    return false;
  }
  if (descr == "BondAnd" || descr == "BondXor") {
    return true;
  }
  if (descr == "BondOr") {
    // detect the types of queries that appear for unspecified bonds in SMARTS:
    if (b->getQuery()->endChildren() - b->getQuery()->beginChildren() == 2) {
      for (auto child = b->getQuery()->beginChildren();
           child != b->getQuery()->endChildren(); ++child) {
        if ((*child)->getDescription() != "BondOrder" ||
            (*child)->getNegation()) {
          return true;
        }
        if (static_cast<BOND_EQUALS_QUERY *>(child->get())->getVal() !=
                Bond::SINGLE &&
            static_cast<BOND_EQUALS_QUERY *>(child->get())->getVal() !=
                Bond::AROMATIC) {
          return true;
        }
      }
      return false;
    }
  }

  return true;
}

bool _complexQueryHelper(Atom::QUERYATOM_QUERY const *query, bool &hasAtNum) {
  if (!query) {
    return false;
  }
  if (query->getNegation()) {
    return true;
  }
  std::string descr = query->getDescription();
  // std::cerr<<" |"<<descr;
  if (descr == "AtomAtomicNum" || descr == "AtomType") {
    hasAtNum = true;
    return false;
  }
  if (descr == "AtomOr" || descr == "AtomXor") {
    return true;
  }
  if (descr == "AtomAnd") {
    auto childIt = query->beginChildren();
    while (childIt != query->endChildren()) {
      if (_complexQueryHelper(childIt->get(), hasAtNum)) {
        return true;
      }
      ++childIt;
    }
  }
  return false;
}
bool isComplexQuery(const Atom *a) {
  if (!a->hasQuery()) {
    return false;
  }
  // std::cerr<<"\n"<<a->getIdx();
  // negated things are always complex:
  if (a->getQuery()->getNegation()) {
    return true;
  }
  std::string descr = a->getQuery()->getDescription();
  // std::cerr<<" "<<descr;
  if (descr == "AtomNull" || descr == "AtomAtomicNum" || descr == "AtomType") {
    return false;
  }
  if (descr == "AtomOr" || descr == "AtomXor") {
    return true;
  }
  if (descr == "AtomAnd") {
    bool hasAtNum = false;
    if (_complexQueryHelper(a->getQuery(), hasAtNum)) {
      return true;
    }
    if (hasAtNum) {
      return false;
    } else {
      return true;
    }
  }

  return true;
}
bool isAtomAromatic(const Atom *a) {
  bool res = false;
  if (!a->hasQuery()) {
    res = a->getIsAromatic();
  } else {
    std::string descr = a->getQuery()->getDescription();
    if (descr == "AtomAtomicNum") {
      res = a->getIsAromatic();
    } else if (descr == "AtomIsAromatic") {
      res = true;
      if (a->getQuery()->getNegation()) {
        res = !res;
      }
    } else if (descr == "AtomIsAliphatic") {
      res = false;
      if (a->getQuery()->getNegation()) {
        res = !res;
      }
    } else if (descr == "AtomType") {
      res = getAtomTypeIsAromatic(
          static_cast<ATOM_EQUALS_QUERY *>(a->getQuery())->getVal());
      if (a->getQuery()->getNegation()) {
        res = !res;
      }
    } else if (descr == "AtomAnd") {
      auto childIt = a->getQuery()->beginChildren();
      if ((*childIt)->getDescription() == "AtomAtomicNum") {
        if (a->getQuery()->getNegation()) {
          res = false;
        } else if ((*(childIt + 1))->getDescription() == "AtomIsAliphatic") {
          res = false;
        } else if ((*(childIt + 1))->getDescription() == "AtomIsAromatic") {
          res = true;
        }
      }
    }
  }
  return res;
}

namespace QueryOps {
namespace {
void completeQueryAndChildren(ATOM_EQUALS_QUERY *query, Atom *tgt,
                              unsigned int magicVal) {
  PRECONDITION(query, "no query");
  PRECONDITION(tgt, "no atom");
  if (static_cast<unsigned int>(query->getVal()) == magicVal) {
    int tgtVal = query->getDataFunc()(tgt);
    query->setVal(tgtVal);
  }
  for (auto childIt = query->beginChildren(); childIt != query->endChildren();
       ++childIt) {
    completeQueryAndChildren((ATOM_EQUALS_QUERY *)(childIt->get()), tgt,
                             magicVal);
  }
}
}  // namespace
void completeMolQueries(RWMol *mol, unsigned int magicVal) {
  PRECONDITION(mol, "bad molecule");
  for (auto atom : mol->atoms()) {
    if (atom->hasQuery()) {
      auto *query = static_cast<ATOM_EQUALS_QUERY *>(atom->getQuery());
      completeQueryAndChildren(query, atom, magicVal);
    }
  }
}

Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atom, "bad atom");
  if (atom->hasQuery()) {
    return atom;
  }

  QueryAtom qa(*atom);
  unsigned int idx = atom->getIdx();

  if (atom->hasProp(common_properties::_hasMassQuery)) {
    qa.expandQuery(makeAtomMassQuery(static_cast<int>(atom->getMass())));
  }
  mol->replaceAtom(idx, &qa);
  return mol->getAtomWithIdx(idx);
}
}  // namespace QueryOps
};  // namespace RDKit
