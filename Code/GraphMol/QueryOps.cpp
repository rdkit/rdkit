//
// Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
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
#include <boost/range/iterator_range.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>

namespace RDKit {

// common general queries

int queryIsAtomBridgehead(Atom const *at) {
  // at least three ring bonds, all ring bonds in a ring which shares at
  // least two bonds with another ring involving this atom
  //
  // We can't just go with "at least three ring bonds shared between multiple
  // rings" because of structures like CC12CCN(CC1)C2 where there are only two
  // SSSRs
  PRECONDITION(at, "no atom");
  if (at->getDegree() < 3) {
    return 0;
  }
  const auto &mol = at->getOwningMol();
  const auto ri = mol.getRingInfo();
  if (!ri || !ri->isInitialized()) {
    return 0;
  }
  // track which ring bonds involve this atom
  boost::dynamic_bitset<> atomRingBonds(mol.getNumBonds());
  for (const auto bnd : mol.atomBonds(at)) {
    if (ri->numBondRings(bnd->getIdx())) {
      atomRingBonds.set(bnd->getIdx());
    }
  }
  if (atomRingBonds.count() < 3) {
    return 0;
  }

  boost::dynamic_bitset<> bondsInRingI(mol.getNumBonds());
  boost::dynamic_bitset<> ringsOverlap(ri->numRings());
  for (unsigned int i = 0; i < ri->bondRings().size(); ++i) {
    bondsInRingI.reset();
    bool atomInRingI = false;
    for (const auto bidx : ri->bondRings()[i]) {
      bondsInRingI.set(bidx);
      if (atomRingBonds[bidx]) {
        atomInRingI = true;
      }
    }
    if (!atomInRingI) {
      continue;
    }
    for (unsigned int j = i + 1; j < ri->bondRings().size(); ++j) {
      unsigned int overlap = 0;
      bool atomInRingJ = false;
      for (const auto bidx : ri->bondRings()[j]) {
        if (atomRingBonds[bidx]) {
          atomInRingJ = true;
        }
        if (bondsInRingI[bidx]) {
          ++overlap;
        }
        if (overlap >= 2 && atomInRingJ) {
          // we have two rings containing the atom which share at least two
          // bonds:
          ringsOverlap.set(i);
          ringsOverlap.set(j);
          break;
        }
      }
    }
    if (!ringsOverlap[i]) {
      return 0;
    }
  }
  return 1;
}

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
  ROMol::OEDGE_ITER beg, end;
  boost::tie(beg, end) = at->getOwningMol().getAtomBonds(at);
  unsigned int prod = 1;
  while (beg != end) {
    prod *= static_cast<unsigned int>(
        firstThousandPrimes[at->getOwningMol()[*beg]->getBondType()]);
    ++beg;
  }
  return prod;
}
unsigned int queryAtomAllBondProduct(Atom const *at) {
  ROMol::OEDGE_ITER beg, end;

  boost::tie(beg, end) = at->getOwningMol().getAtomBonds(at);
  unsigned int prod = 1;
  while (beg != end) {
    prod *= static_cast<unsigned int>(
        firstThousandPrimes[at->getOwningMol()[*beg]->getBondType()]);
    ++beg;
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

ATOM_EQUALS_QUERY *makeAtomIsBridgeheadQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryIsAtomBridgehead);
  res->setDescription("AtomIsBridgehead");
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
  // We expanded this with !#0 as part of #6106
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
  // We expanded this with !#0 as part of #6106
  // it's easier to define what isn't a metal than what is. :-)
  auto *res = new ATOM_OR_QUERY;
  res->setDescription("AtomOr");
  res->setNegation(true);
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(0)));
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
  res->setTypeLabel("BondOrder");
  return res;
}

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrAromaticBondQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrAromatic);
  res->setDescription("SingleOrAromaticBond");
  res->setTypeLabel("BondOrder");
  return res;
};

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeDoubleOrAromaticBondQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondIsDoubleOrAromatic);
  res->setDescription("DoubleOrAromaticBond");
  res->setTypeLabel("BondOrder");
  return res;
};

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrDoubleBondQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrDouble);
  res->setDescription("SingleOrDoubleBond");
  res->setTypeLabel("BondOrder");
  return res;
};

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *
makeSingleOrDoubleOrAromaticBondQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrDoubleOrAromatic);
  res->setDescription("SingleOrDoubleOrAromaticBond");
  res->setTypeLabel("BondOrder");
  return res;
};

namespace QueryOps {
// we don't use these anymore but we need to keep them around for backwards
// compatibility with pickled queries. There's no reason to update this list.
const std::vector<std::string> bondOrderQueryFunctions{
    std::string("BondOrder"), std::string("SingleOrAromaticBond"),
    std::string("DoubleOrAromaticBond"), std::string("SingleOrDoubleBond"),
    std::string("SingleOrDoubleOrAromaticBond")};
RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    const Queries::Query<int, Bond const *, true> &qry) {
  const auto df = qry.getDescription();
  const auto dt = qry.getTypeLabel();
  // is this a bond order query?
  if (dt == "BondOrder" ||
      (dt.empty() &&
       std::find(bondOrderQueryFunctions.begin(), bondOrderQueryFunctions.end(),
                 df) != bondOrderQueryFunctions.end())) {
    return true;
  }
  for (const auto &child :
       boost::make_iterator_range(qry.beginChildren(), qry.endChildren())) {
    if (hasBondTypeQuery(*child)) {
      return true;
    }
  }
  return false;
}

namespace {
bool hasComplexBondTypeQueryHelper(
    const Queries::Query<int, Bond const *, true> &qry, bool seenBondOrder) {
  const auto df = qry.getDescription();
  bool isBondOrder = (df == "BondOrder");
  // is this a bond order query?
  if (std::find(bondOrderQueryFunctions.begin(), bondOrderQueryFunctions.end(),
                df) != bondOrderQueryFunctions.end()) {
    if (seenBondOrder || !isBondOrder || qry.getNegation()) {
      return true;
    }
  }
  for (const auto &child :
       boost::make_iterator_range(qry.beginChildren(), qry.endChildren())) {
    if (hasComplexBondTypeQueryHelper(*child, seenBondOrder | isBondOrder)) {
      return true;
    }
    if (child->getDescription() == "BondOrder") {
      seenBondOrder = true;
    }
  }
  return false;
}
}  // namespace

RDKIT_GRAPHMOL_EXPORT bool hasComplexBondTypeQuery(
    const Queries::Query<int, Bond const *, true> &qry) {
  return hasComplexBondTypeQueryHelper(qry, false);
}
}  // namespace QueryOps

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

void convertComplexNameToQuery(Atom *query, std::string_view symb) {
  if (symb == "Q") {
    query->setQuery(makeQAtomQuery());
  } else if (symb == "QH") {
    query->setQuery(makeQHAtomQuery());
  } else if (symb == "A") {
    query->setQuery(makeAAtomQuery());
  } else if (symb == "AH") {
    query->setQuery(makeAHAtomQuery());
  } else if (symb == "X") {
    query->setQuery(makeXAtomQuery());
  } else if (symb == "XH") {
    query->setQuery(makeXHAtomQuery());
  } else if (symb == "M") {
    query->setQuery(makeMAtomQuery());
  } else if (symb == "MH") {
    query->setQuery(makeMHAtomQuery());
  } else {
    // we control what this function gets called with, so we should never land
    // here
    ASSERT_INVARIANT(0, "bad complex query symbol");
  }
}

bool isComplexQuery(const Bond *b) {
  PRECONDITION(b, "bad bond");
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
    // detect the types of queries that appear for unspecified bonds in
    // SMARTS:
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

namespace {
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

template <typename T>
bool _atomListQueryHelper(const T query, bool ignoreNegation) {
  PRECONDITION(query, "no query");
  if (!ignoreNegation && query->getNegation()) {
    return false;
  }
  if (query->getDescription() == "AtomAtomicNum" ||
      query->getDescription() == "AtomType") {
    return true;
  }
  if (query->getDescription() == "AtomOr") {
    for (const auto &child : boost::make_iterator_range(query->beginChildren(),
                                                        query->endChildren())) {
      if (!_atomListQueryHelper(child, ignoreNegation)) {
        return false;
      }
    }
    return true;
  }
  return false;
}
}  // namespace
bool isAtomListQuery(const Atom *a) {
  PRECONDITION(a, "bad atom");
  if (!a->hasQuery()) {
    return false;
  }
  if (a->getQuery()->getDescription() == "AtomOr") {
    for (const auto &child : boost::make_iterator_range(
             a->getQuery()->beginChildren(), a->getQuery()->endChildren())) {
      if (!_atomListQueryHelper(child, false)) {
        return false;
      }
    }
    return true;
  } else if (a->getQuery()->getNegation() &&
             _atomListQueryHelper(a->getQuery(), true)) {
    // this was github #5930: negated list queries containing a single atom were
    // being lost on output
    return true;
  }
  return false;
}

void getAtomListQueryVals(const Atom::QUERYATOM_QUERY *q,
                          std::vector<int> &vals) {
  // list queries are series of nested ors of AtomAtomicNum queries
  PRECONDITION(q, "bad query");
  auto descr = q->getDescription();
  if (descr == "AtomOr") {
    for (const auto &child :
         boost::make_iterator_range(q->beginChildren(), q->endChildren())) {
      auto descr = child->getDescription();
      if (child->getNegation() ||
          (descr != "AtomOr" && descr != "AtomAtomicNum" &&
           descr != "AtomType")) {
        throw ValueErrorException("bad query type1");
      }
      // we don't allow negation of any children of the query:
      if (descr == "AtomOr") {
        getAtomListQueryVals(child.get(), vals);
      } else if (descr == "AtomAtomicNum") {
        vals.push_back(static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal());
      } else if (descr == "AtomType") {
        auto v = static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal();
        // aromatic AtomType queries add 1000 to the atomic number;
        // correct for that:
        if (v >= 1000) {
          v -= 1000;
        }
        vals.push_back(v);
      }
    }
  } else if (descr == "AtomAtomicNum") {
    vals.push_back(static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal());
  } else if (descr == "AtomType") {
    auto v = static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal();
    // aromatic AtomType queries add 1000 to the atomic number;
    // correct for that:
    if (v >= 1000) {
      v -= 1000;
    }
    vals.push_back(v);
  } else {
    CHECK_INVARIANT(0, "bad query type");
  }
}

bool isComplexQuery(const Atom *a) {
  PRECONDITION(a, "bad atom");
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
    return !hasAtNum;
  }

  return true;
}
bool isAtomAromatic(const Atom *a) {
  PRECONDITION(a, "bad atom");
  bool res = false;
  if (!a->hasQuery()) {
    res = isAromaticAtom(*a);
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
void completeQueryAndChildren(Atom::QUERYATOM_QUERY *query, Atom *tgt,
                              unsigned int magicVal) {
  PRECONDITION(query, "no query");
  PRECONDITION(tgt, "no atom");
  auto eqQuery = dynamic_cast<ATOM_EQUALS_QUERY *>(query);
  if (eqQuery) {
    if (static_cast<unsigned int>(eqQuery->getVal()) == magicVal) {
      int tgtVal = eqQuery->getDataFunc()(tgt);
      eqQuery->setVal(tgtVal);
    }
  }
  for (auto childIt = query->beginChildren(); childIt != query->endChildren();
       ++childIt) {
    completeQueryAndChildren(childIt->get(), tgt, magicVal);
  }
}
}  // namespace
void completeMolQueries(RWMol *mol, unsigned int magicVal) {
  PRECONDITION(mol, "bad molecule");
  for (auto atom : mol->atoms()) {
    if (atom->hasQuery()) {
      completeQueryAndChildren(atom->getQuery(), atom, magicVal);
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

void finalizeQueryFromDescription(
    Queries::Query<int, Atom const *, true> *query, Atom const *) {
  std::string descr = query->getDescription();

  if (boost::starts_with(descr, "range_")) {
    descr = descr.substr(6);
  } else if (boost::starts_with(descr, "less_")) {
    descr = descr.substr(5);
  } else if (boost::starts_with(descr, "greater_")) {
    descr = descr.substr(8);
  }

  Queries::Query<int, Atom const *, true> *tmpQuery;
  if (descr == "AtomRingBondCount") {
    query->setDataFunc(queryAtomRingBondCount);
  } else if (descr == "AtomHasRingBond") {
    query->setDataFunc(queryAtomHasRingBond);
  } else if (descr == "AtomRingSize") {
    tmpQuery = makeAtomInRingOfSizeQuery(
        static_cast<ATOM_EQUALS_QUERY *>(query)->getVal());
    query->setDataFunc(tmpQuery->getDataFunc());
    delete tmpQuery;
  } else if (descr == "AtomMinRingSize") {
    query->setDataFunc(queryAtomMinRingSize);
  } else if (descr == "AtomImplicitValence") {
    query->setDataFunc(queryAtomImplicitValence);
  } else if (descr == "AtomTotalValence") {
    query->setDataFunc(queryAtomTotalValence);
  } else if (descr == "AtomAtomicNum") {
    query->setDataFunc(queryAtomNum);
  } else if (descr == "AtomExplicitDegree") {
    query->setDataFunc(queryAtomExplicitDegree);
  } else if (descr == "AtomTotalDegree") {
    query->setDataFunc(queryAtomTotalDegree);
  } else if (descr == "AtomHeavyAtomDegree") {
    query->setDataFunc(queryAtomHeavyAtomDegree);
  } else if (descr == "AtomHCount") {
    query->setDataFunc(queryAtomHCount);
  } else if (descr == "AtomImplicitHCount") {
    query->setDataFunc(queryAtomImplicitHCount);
  } else if (descr == "AtomHasImplicitH") {
    query->setDataFunc(queryAtomHasImplicitH);
  } else if (descr == "AtomIsAromatic") {
    query->setDataFunc(queryAtomAromatic);
  } else if (descr == "AtomIsAliphatic") {
    query->setDataFunc(queryAtomAliphatic);
  } else if (descr == "AtomUnsaturated") {
    query->setDataFunc(queryAtomUnsaturated);
  } else if (descr == "AtomMass") {
    query->setDataFunc(queryAtomMass);
  } else if (descr == "AtomIsotope") {
    query->setDataFunc(queryAtomIsotope);
  } else if (descr == "AtomFormalCharge") {
    query->setDataFunc(queryAtomFormalCharge);
  } else if (descr == "AtomNegativeFormalCharge") {
    query->setDataFunc(queryAtomNegativeFormalCharge);
  } else if (descr == "AtomHybridization") {
    query->setDataFunc(queryAtomHybridization);
  } else if (descr == "AtomInRing") {
    query->setDataFunc(queryIsAtomInRing);
  } else if (descr == "AtomInNRings") {
    query->setDataFunc(queryIsAtomInNRings);
  } else if (descr == "AtomHasHeteroatomNeighbors") {
    query->setDataFunc(queryAtomHasHeteroatomNbrs);
  } else if (descr == "AtomNumHeteroatomNeighbors") {
    query->setDataFunc(queryAtomNumHeteroatomNbrs);
  } else if (descr == "AtomNonHydrogenDegree") {
    query->setDataFunc(queryAtomNonHydrogenDegree);
  } else if (descr == "AtomHasAliphaticHeteroatomNeighbors") {
    query->setDataFunc(queryAtomHasAliphaticHeteroatomNbrs);
  } else if (descr == "AtomNumAliphaticHeteroatomNeighbors") {
    query->setDataFunc(queryAtomNumAliphaticHeteroatomNbrs);
  } else if (descr == "AtomNull") {
    query->setDataFunc(nullDataFun);
    query->setMatchFunc(nullQueryFun);
  } else if (descr == "AtomType") {
    query->setDataFunc(queryAtomType);
  } else if (descr == "AtomNumRadicalElectrons") {
    query->setDataFunc(queryAtomNumRadicalElectrons);
  } else if (descr == "AtomInNRings" || descr == "RecursiveStructure") {
    // don't need to do anything here because the classes
    // automatically have everything set
  } else if (descr == "AtomAnd" || descr == "AtomOr" || descr == "AtomXor" ||
             descr == "HasProp") {
    // don't need to do anything here because the classes
    // automatically have everything set
  } else {
    throw ValueErrorException("Do not know how to finalize query: '" + descr +
                              "'");
  }
}

void finalizeQueryFromDescription(
    Queries::Query<int, Bond const *, true> *query, Bond const *) {
  std::string descr = query->getDescription();
  Queries::Query<int, Bond const *, true> *tmpQuery;
  if (descr == "BondRingSize") {
    tmpQuery = makeBondInRingOfSizeQuery(
        static_cast<BOND_EQUALS_QUERY *>(query)->getVal());
    query->setDataFunc(tmpQuery->getDataFunc());
    delete tmpQuery;
  } else if (descr == "BondMinRingSize") {
    query->setDataFunc(queryBondMinRingSize);
  } else if (descr == "BondOrder") {
    query->setDataFunc(queryBondOrder);
  } else if (descr == "BondDir") {
    query->setDataFunc(queryBondDir);
  } else if (descr == "BondInRing") {
    query->setDataFunc(queryIsBondInRing);
  } else if (descr == "BondInNRings") {
    query->setDataFunc(queryIsBondInNRings);
  } else if (descr == "SingleOrAromaticBond") {
    query->setDataFunc(queryBondIsSingleOrAromatic);
  } else if (descr == "SingleOrDoubleBond") {
    query->setDataFunc(queryBondIsSingleOrDouble);
  } else if (descr == "DoubleOrAromaticBond") {
    query->setDataFunc(queryBondIsDoubleOrAromatic);
  } else if (descr == "SingleOrDoubleOrAromaticBond") {
    query->setDataFunc(queryBondIsSingleOrDoubleOrAromatic);
  } else if (descr == "BondNull") {
    query->setDataFunc(nullDataFun);
    query->setMatchFunc(nullQueryFun);
  } else if (descr == "BondAnd" || descr == "BondOr" || descr == "BondXor" ||
             descr == "HasProp") {
    // don't need to do anything here because the classes
    // automatically have everything set
  } else {
    throw ValueErrorException("Do not know how to finalize query: '" + descr +
                              "'");
  }
}

}  // namespace QueryOps
};  // namespace RDKit
