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
#include "RDMol.h"
#include <algorithm>
#include <RDGeneral/types.h>
#include <GraphMol/QueryAtom.h>
#include <boost/range/iterator_range.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/algorithm/string.hpp>

#include <Query/NullQueryAlgebra.h>

namespace RDKit {

// common general queries

int queryIsAtomBridgehead(Atom const *at) {
  return queryIsAtomBridgeheadInternal(at->getRDMol(), at->getIdx(),
                                       at->getRDMol().getRingInfo()) ? 1 : 0;
}
bool queryIsAtomBridgeheadInternal(const RDMol &mol, atomindex_t atomIndex,
                           const RingInfoCache &rings) {
  // at least three ring bonds, at least one ring bond in a ring which shares at
  // least two bonds with another ring involving this atom
  //
  // We can't just go with "at least three ring bonds shared between multiple
  // rings" because of structures like CC12CCN(CC1)C2 where there are only two
  // SSSRs
  PRECONDITION(atomIndex < mol.getNumAtoms(), "no atom");
  uint32_t atomDegree = mol.getAtomDegree(atomIndex);
  if (atomDegree < 3) {
    return false;
  }
  if (!rings.isInitialized()) {
    // FIXME: Should there be some sort of warning here?
    return false;
  }
  // track which bonds involve this atom
  auto [beginBonds, endBonds] = mol.getAtomBonds(atomIndex);
  size_t numAtomRingBonds = 0;
  for (; beginBonds != endBonds; ++beginBonds) {
    const uint32_t bondIndex = *beginBonds;
    if (rings.numBondRings(bondIndex)) {
      ++numAtomRingBonds;
    }
  }
  if (numAtomRingBonds < 3) {
    return false;
  }

  // Check all pairs of rings that contain this atom.
  const uint32_t beginAtomRings = rings.atomMembershipBegins[atomIndex];
  const uint32_t endAtomRings = rings.atomMembershipBegins[atomIndex+1];
  boost::dynamic_bitset<> bondsInRing(mol.getNumBonds());
  boost::dynamic_bitset<> ringsOverlap(rings.numRings());
  for (uint32_t atomRingIndexI = beginAtomRings;
       atomRingIndexI != endAtomRings; ++atomRingIndexI) {
    const uint32_t ringIndexI = rings.atomMemberships[atomRingIndexI];

    bondsInRing.reset();
    uint32_t beginBondRingI = rings.ringBegins[ringIndexI];
    const uint32_t endBondRingI = rings.ringBegins[ringIndexI+1];
    for (; beginBondRingI != endBondRingI; ++beginBondRingI) {
      const uint32_t bondIndex = rings.bondsInRings[beginBondRingI];
      bondsInRing.set(bondIndex);
    }

    for (uint32_t atomRingIndexJ = atomRingIndexI + 1;
         atomRingIndexJ != endAtomRings; ++atomRingIndexJ) {
      const uint32_t ringIndexJ = rings.atomMemberships[atomRingIndexJ];

      uint32_t overlap = 0;
      uint32_t beginBondRingJ = rings.ringBegins[ringIndexJ];
      const uint32_t endBondRingJ = rings.ringBegins[ringIndexJ + 1];
      for (; beginBondRingJ != endBondRingJ; ++beginBondRingJ) {
        const uint32_t bondIndex = rings.bondsInRings[beginBondRingJ];
        if (bondsInRing[bondIndex]) {
          ++overlap;
          if (overlap >= 2) {
            // we have two rings containing the atom which share at least two
            // bonds:
            ringsOverlap.set(ringIndexI);
            ringsOverlap.set(ringIndexJ);
            break;
          }
        }
      }
    }
    if (!ringsOverlap[ringIndexI]) {
      return false;
    }
  }
  return true;
}

//! returns a Query for matching atoms with a particular number of ring bonds
ATOM_EQUALS_QUERY2 *makeAtomRingBondCountQuery2(int what) {
  ATOM_EQUALS_QUERY2 *res = new AtomRingQuery2(what);
  res->setDescription("AtomRingBondCount");
  res->setDataFunc(queryAtomRingBondCount2);
  return res;
};
//! returns a Query for matching atoms with a particular number of ring bonds
ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what) {
  ATOM_EQUALS_QUERY *res = new AtomRingQuery(what);
  res->setDescription("AtomRingBondCount");
  res->setDataFunc(queryAtomRingBondCount);
  return res;
};

template <bool isBond, bool newVersion>
using ParamType = std::conditional_t<isBond,
      std::conditional_t<newVersion, ConstRDMolBond, Bond const *>,
      std::conditional_t<newVersion, ConstRDMolAtom, Atom const *>>;
template <bool isBond, bool newVersion>
using EqualsType = std::conditional_t<isBond,
    std::conditional_t<newVersion, BOND_EQUALS_QUERY2, BOND_EQUALS_QUERY>,
    std::conditional_t<newVersion, ATOM_EQUALS_QUERY2, ATOM_EQUALS_QUERY>>;

template <bool isBond, bool newVersion>
using DataFuncPtrType = int (*)(ParamType<isBond, newVersion>);

template <bool isBond, bool newVersion, int tgt>
int queryIsInRingOfSizeTmpl(ParamType<isBond, newVersion> at) {
  if constexpr (isBond) {
    if constexpr (newVersion) {
      return queryBondIsInRingOfSize2<tgt>(at);
    } else {
      return queryBondIsInRingOfSize<tgt>(at);
    }
  } else {
    if constexpr (newVersion) {
      return queryAtomIsInRingOfSize2<tgt>(at);
    } else {
      return queryAtomIsInRingOfSize<tgt>(at);
    }
  }
}

template <bool isBond, bool newVersion>
DataFuncPtrType<isBond, newVersion> getQueryIsInRingOfSizeTmpl(int tgt) {
  switch (tgt) {
    case 3:  return queryIsInRingOfSizeTmpl<isBond, newVersion, 3>;
    case 4:  return queryIsInRingOfSizeTmpl<isBond, newVersion, 4>;
    case 5:  return queryIsInRingOfSizeTmpl<isBond, newVersion, 5>;
    case 6:  return queryIsInRingOfSizeTmpl<isBond, newVersion, 6>;
    case 7:  return queryIsInRingOfSizeTmpl<isBond, newVersion, 7>;
    case 8:  return queryIsInRingOfSizeTmpl<isBond, newVersion, 8>;
    case 9:  return queryIsInRingOfSizeTmpl<isBond, newVersion, 9>;
    case 10: return queryIsInRingOfSizeTmpl<isBond, newVersion, 10>;
    case 11: return queryIsInRingOfSizeTmpl<isBond, newVersion, 11>;
    case 12: return queryIsInRingOfSizeTmpl<isBond, newVersion, 12>;
    case 13: return queryIsInRingOfSizeTmpl<isBond, newVersion, 13>;
    case 14: return queryIsInRingOfSizeTmpl<isBond, newVersion, 14>;
    case 15: return queryIsInRingOfSizeTmpl<isBond, newVersion, 15>;
    case 16: return queryIsInRingOfSizeTmpl<isBond, newVersion, 16>;
    case 17: return queryIsInRingOfSizeTmpl<isBond, newVersion, 17>;
    case 18: return queryIsInRingOfSizeTmpl<isBond, newVersion, 18>;
    case 19: return queryIsInRingOfSizeTmpl<isBond, newVersion, 19>;
    case 20: return queryIsInRingOfSizeTmpl<isBond, newVersion, 20>;
  }
  return nullptr;
}

template <bool isBond, bool newVersion>
EqualsType<isBond, newVersion> *makeInRingOfSizeQueryTmpl(int tgt) {
  RANGE_CHECK(3, tgt, 20);
  auto *res = new EqualsType<isBond, newVersion>;
  res->setVal(tgt);

  auto *dataFunc = getQueryIsInRingOfSizeTmpl<isBond, newVersion>(tgt);
  if (dataFunc != nullptr) {
    res->setDataFunc(dataFunc);
  }

  res->setDescription(isBond ? "BondRingSize" : "AtomRingSize");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomInRingOfSizeQuery2(int tgt) {
  return makeInRingOfSizeQueryTmpl<false, true>(tgt);
}
ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt) {
  return makeInRingOfSizeQueryTmpl<false, false>(tgt);
}

BOND_EQUALS_QUERY2 *makeBondInRingOfSizeQuery2(int tgt) {
  return makeInRingOfSizeQueryTmpl<true, true>(tgt);
}
BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int tgt) {
  return makeInRingOfSizeQueryTmpl<true, false>(tgt);
}

ATOM_EQUALS_QUERY2 *makeAtomMinRingSizeQuery2(int tgt) {
  auto *res = new ATOM_EQUALS_QUERY2;
  res->setVal(tgt);
  res->setDataFunc(queryAtomMinRingSize2);
  res->setDescription("AtomMinRingSize");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomMinRingSizeQuery(int tgt) {
  auto *res = new ATOM_EQUALS_QUERY;
  res->setVal(tgt);
  res->setDataFunc(queryAtomMinRingSize);
  res->setDescription("AtomMinRingSize");
  return res;
}
BOND_EQUALS_QUERY2 *makeBondMinRingSizeQuery2(int tgt) {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(tgt);
  res->setDataFunc(queryBondMinRingSize2);
  res->setDescription("BondMinRingSize");
  return res;
}
BOND_EQUALS_QUERY *makeBondMinRingSizeQuery(int tgt) {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(tgt);
  res->setDataFunc(queryBondMinRingSize);
  res->setDescription("BondMinRingSize");
  return res;
}

unsigned int queryAtomBondProduct2(ConstRDMolAtom at) {
  auto [beg, end] = at.mol().getAtomBonds(at.index());
  unsigned int prod = 1;
  while (beg != end) {
    prod *= static_cast<unsigned int>(
        firstThousandPrimes[at.mol().getBond(*beg).getBondType()]);
    ++beg;
  }
  return prod;
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
unsigned int queryAtomAllBondProduct2(ConstRDMolAtom at) {
  auto [beg, end] = at.mol().getAtomBonds(at.index());
  unsigned int prod = 1;
  while (beg != end) {
    prod *= static_cast<unsigned int>(
        firstThousandPrimes[at.mol().getBond(*beg).getBondType()]);
    ++beg;
  }
  for (unsigned int i = 0, n = at.mol().getTotalNumHs(at.index()); i < n; i++) {
    prod *= static_cast<unsigned int>(firstThousandPrimes[Bond::SINGLE]);
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

ATOM_EQUALS_QUERY2 *makeAtomImplicitValenceQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomImplicitValence2);
  res->setDescription("AtomImplicitValence");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomImplicitValence);
  res->setDescription("AtomImplicitValence");
  return res;
}
ATOM_EQUALS_QUERY2 *makeAtomExplicitValenceQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomExplicitValence2);
  res->setDescription("AtomExplicitValence");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomExplicitValenceQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomExplicitValence);
  res->setDescription("AtomExplicitValence");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomTotalValenceQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomTotalValence2);
  res->setDescription("AtomTotalValence");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomTotalValence);
  res->setDescription("AtomTotalValence");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomNumQuery2(int what) {
  return makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomNum2,
                                                  "AtomAtomicNum");
}
ATOM_EQUALS_QUERY *makeAtomNumQuery(int what) {
  return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNum,
                                                "AtomAtomicNum");
}

ATOM_EQUALS_QUERY2 *makeAtomTypeQuery2(int num, int aromatic) {
  return makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(makeAtomType(num, aromatic),
                                                  queryAtomType2, "AtomType");
}
ATOM_EQUALS_QUERY *makeAtomTypeQuery(int num, int aromatic) {
  return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(makeAtomType(num, aromatic),
                                                queryAtomType, "AtomType");
}
ATOM_EQUALS_QUERY2 *makeAtomExplicitDegreeQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomExplicitDegree2);
  res->setDescription("AtomExplicitDegree");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomExplicitDegree);
  res->setDescription("AtomExplicitDegree");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomTotalDegreeQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomTotalDegree2);
  res->setDescription("AtomTotalDegree");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomTotalDegree);
  res->setDescription("AtomTotalDegree");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomHeavyAtomDegreeQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomHeavyAtomDegree2);
  res->setDescription("AtomHeavyAtomDegree");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHeavyAtomDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHeavyAtomDegree);
  res->setDescription("AtomHeavyAtomDegree");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomHCountQuery2(int what) {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomHCount2);
  res->setDescription("AtomHCount");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHCount);
  res->setDescription("AtomHCount");
  return res;
}
ATOM_EQUALS_QUERY2 *makeAtomImplicitHCountQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomImplicitHCount2);
  res->setDescription("AtomImplicitHCount");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomImplicitHCountQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomImplicitHCount);
  res->setDescription("AtomImplicitHCount");
  return res;
}
ATOM_EQUALS_QUERY2 *makeAtomHasImplicitHQuery2() {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomHasImplicitH2);
  res->setDescription("AtomHasImplicitH");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHasImplicitHQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasImplicitH);
  res->setDescription("AtomHasImplicitH");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomAromaticQuery2() {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomAromatic2);
  res->setDescription("AtomIsAromatic");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomAromaticQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomAromatic);
  res->setDescription("AtomIsAromatic");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomAliphaticQuery2() {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomAliphatic2);
  res->setDescription("AtomIsAliphatic");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomAliphaticQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomAliphatic);
  res->setDescription("AtomIsAliphatic");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomUnsaturatedQuery2() {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomUnsaturated2);
  res->setDescription("AtomUnsaturated");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomUnsaturated);
  res->setDescription("AtomUnsaturated");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomMassQuery2(int what) {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(
      massIntegerConversionFactor * what, queryAtomMass2);
  res->setDescription("AtomMass");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomMassQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      massIntegerConversionFactor * what, queryAtomMass);
  res->setDescription("AtomMass");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomIsotopeQuery2(int what) {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomIsotope2);
  res->setDescription("AtomIsotope");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomIsotope);
  res->setDescription("AtomIsotope");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomFormalChargeQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomFormalCharge2);
  res->setDescription("AtomFormalCharge");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomFormalCharge);
  res->setDescription("AtomFormalCharge");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomNegativeFormalChargeQuery2(int what) {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(
      what, queryAtomNegativeFormalCharge2);
  res->setDescription("AtomNegativeFormalCharge");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomNegativeFormalChargeQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      what, queryAtomNegativeFormalCharge);
  res->setDescription("AtomNegativeFormalCharge");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomHybridizationQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomHybridization2);
  res->setDescription("AtomHybridization");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHybridization);
  res->setDescription("AtomHybridization");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomNumRadicalElectronsQuery2(int what) {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(
      what, queryAtomNumRadicalElectrons2);
  res->setDescription("AtomNumRadicalElectrons");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomNumRadicalElectronsQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      what, queryAtomNumRadicalElectrons);
  res->setDescription("AtomNumRadicalElectrons");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomHasChiralTagQuery2() {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomHasChiralTag2);
  res->setDescription("AtomHasChiralTag");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHasChiralTagQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasChiralTag);
  res->setDescription("AtomHasChiralTag");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomMissingChiralTagQuery2() {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomMissingChiralTag2);
  res->setDescription("AtomMissingChiralTag");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomMissingChiralTagQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomMissingChiralTag);
  res->setDescription("AtomMissingChiralTag");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomInRingQuery2() {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryIsAtomInRing2);
  res->setDescription("AtomInRing");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomInRingQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryIsAtomInRing);
  res->setDescription("AtomInRing");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomIsBridgeheadQuery2() {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryIsAtomBridgehead2);
  res->setDescription("AtomIsBridgehead");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomIsBridgeheadQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryIsAtomBridgehead);
  res->setDescription("AtomIsBridgehead");
  return res;
}

ATOM_OR_QUERY2 *makeQAtomQuery2() {
  auto *res = new ATOM_OR_QUERY2;
  res->setDescription("AtomOr");
  res->setTypeLabel("Q");
  res->setNegation(true);
  res->addChild(Queries::Query<int, ConstRDMolAtom, true>::CHILD_TYPE(
      makeAtomNumQuery2(6)));
  res->addChild(Queries::Query<int, ConstRDMolAtom, true>::CHILD_TYPE(
      makeAtomNumQuery2(1)));
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
ATOM_EQUALS_QUERY2 *makeQHAtomQuery2() {
  ATOM_EQUALS_QUERY2 *res = makeAtomNumQuery2(6);
  res->setNegation(true);
  res->setTypeLabel("QH");
  return res;
}
ATOM_EQUALS_QUERY *makeQHAtomQuery() {
  ATOM_EQUALS_QUERY *res = makeAtomNumQuery(6);
  res->setNegation(true);
  res->setTypeLabel("QH");
  return res;
}
ATOM_EQUALS_QUERY2 *makeAAtomQuery2() {
  ATOM_EQUALS_QUERY2 *res = makeAtomNumQuery2(1);
  res->setNegation(true);
  res->setTypeLabel("A");
  return res;
}
ATOM_EQUALS_QUERY *makeAAtomQuery() {
  ATOM_EQUALS_QUERY *res = makeAtomNumQuery(1);
  res->setNegation(true);
  res->setTypeLabel("A");
  return res;
}
ATOM_NULL_QUERY2 *makeAHAtomQuery2() {
  auto *res = makeAtomNullQuery2();
  res->setTypeLabel("AH");
  return res;
}
ATOM_NULL_QUERY *makeAHAtomQuery() {
  auto *res = makeAtomNullQuery();
  res->setTypeLabel("AH");
  return res;
}

template <bool newVersion, size_t count>
void addAtomNumQueries(
    std::conditional_t<newVersion, ATOM_OR_QUERY2, ATOM_OR_QUERY> *res,
    const std::array<int, count> &nums) {
  for (size_t i = 0; i < count; ++i) {
    if constexpr (newVersion) {
      res->addChild(
          typename Queries::Query<int, ParamType<false, newVersion>, true>::CHILD_TYPE(
              makeAtomNumQuery2(nums[i])));
    } else {
      res->addChild(
          typename Queries::Query<int, ParamType<false, newVersion>, true>::CHILD_TYPE(
              makeAtomNumQuery(nums[i])));
    }
  }
}

ATOM_OR_QUERY2 *makeXAtomQuery2() {
  auto *res = new ATOM_OR_QUERY2;
  res->setDescription("AtomOr");
  const std::array nums = {9, 17, 35, 53, 85};
  addAtomNumQueries<true>(res, nums);
  res->setTypeLabel("X");
  return res;
}
ATOM_OR_QUERY *makeXAtomQuery() {
  auto *res = new ATOM_OR_QUERY;
  res->setDescription("AtomOr");
  const std::array nums = {9, 17, 35, 53, 85};
  addAtomNumQueries<false>(res, nums);
  res->setTypeLabel("X");
  return res;
}
ATOM_OR_QUERY2 *makeXHAtomQuery2() {
  ATOM_OR_QUERY2 *res = makeXAtomQuery2();
  res->addChild(
      Queries::Query<int, ConstRDMolAtom, true>::CHILD_TYPE(makeAtomNumQuery2(1)));
  res->setTypeLabel("XH");
  return res;
}
ATOM_OR_QUERY *makeXHAtomQuery() {
  ATOM_OR_QUERY *res = makeXAtomQuery();
  res->addChild(
      Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
  res->setTypeLabel("XH");
  return res;
}

ATOM_OR_QUERY2 *makeMAtomQuery2() {
  // using the definition from Marvin Sketch, which produces the following
  // SMARTS:
  // !#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
  // We expanded this with !#0 as part of #6106
  // it's easier to define what isn't a metal than what is. :-)
  ATOM_OR_QUERY2 *res = makeMHAtomQuery2();
  res->addChild(
      Queries::Query<int, ConstRDMolAtom, true>::CHILD_TYPE(makeAtomNumQuery2(1)));
  res->setTypeLabel("M");
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
ATOM_OR_QUERY2 *makeMHAtomQuery2() {
  // using the definition from Marvin Sketch, which produces the following
  // SMARTS:
  // !#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
  // We expanded this with !#0 as part of #6106
  // it's easier to define what isn't a metal than what is. :-)
  auto *res = new ATOM_OR_QUERY2;
  res->setDescription("AtomOr");
  res->setNegation(true);
  const std::array nums = {
      0, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86
  };
  addAtomNumQueries<true>(res, nums);
  res->setTypeLabel("MH");
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
  const std::array nums = {
      0, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86
  };
  addAtomNumQueries<false>(res, nums);
  res->setTypeLabel("MH");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomInNRingsQuery2(int what) {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryIsAtomInNRings2);
  res->setDescription("AtomInNRings");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryIsAtomInNRings);
  res->setDescription("AtomInNRings");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomHasRingBondQuery2() {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomHasRingBond2);
  res->setDescription("AtomHasRingBond");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHasRingBondQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasRingBond);
  res->setDescription("AtomHasRingBond");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomNumHeteroatomNbrsQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomNumHeteroatomNbrs2);
  res->setDescription("AtomNumHeteroatomNeighbors");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomNumHeteroatomNbrsQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNumHeteroatomNbrs);
  res->setDescription("AtomNumHeteroatomNeighbors");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomHasHeteroatomNbrsQuery2() {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(true, queryAtomHasHeteroatomNbrs2);
  res->setDescription("AtomHasHeteroatomNeighbors");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHasHeteroatomNbrsQuery() {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasHeteroatomNbrs);
  res->setDescription("AtomHasHeteroatomNeighbors");
  return res;
}
ATOM_EQUALS_QUERY2 *makeAtomNumAliphaticHeteroatomNbrsQuery2(int what) {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(
      what, queryAtomNumAliphaticHeteroatomNbrs2);
  res->setDescription("AtomNumAliphaticHeteroatomNeighbors");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomNumAliphaticHeteroatomNbrsQuery(int what) {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      what, queryAtomNumAliphaticHeteroatomNbrs);
  res->setDescription("AtomNumAliphaticHeteroatomNeighbors");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomHasAliphaticHeteroatomNbrsQuery2() {
  auto *res = makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(
      true, queryAtomHasAliphaticHeteroatomNbrs2);
  res->setDescription("AtomHasAliphaticHeteroatomNeighbors");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomHasAliphaticHeteroatomNbrsQuery() {
  auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
      true, queryAtomHasAliphaticHeteroatomNbrs);
  res->setDescription("AtomHasAliphaticHeteroatomNeighbors");
  return res;
}

ATOM_EQUALS_QUERY2 *makeAtomNonHydrogenDegreeQuery2(int what) {
  auto *res =
      makeAtomSimpleQuery2<ATOM_EQUALS_QUERY2>(what, queryAtomNonHydrogenDegree2);
  res->setDescription("AtomNonHydrogenDegree");
  return res;
}
ATOM_EQUALS_QUERY *makeAtomNonHydrogenDegreeQuery(int what) {
  auto *res =
      makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNonHydrogenDegree);
  res->setDescription("AtomNonHydrogenDegree");
  return res;
}

BOND_EQUALS_QUERY2 *makeBondOrderEqualsQuery2(Bond::BondType what) {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(what);
  res->setDataFunc(queryBondOrder2);
  res->setDescription("BondOrder");
  res->setTypeLabel("BondOrder");
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

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeSingleOrAromaticBondQuery2() {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrAromatic2);
  res->setDescription("SingleOrAromaticBond");
  res->setTypeLabel("BondOrder");
  return res;
};
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrAromaticBondQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrAromatic);
  res->setDescription("SingleOrAromaticBond");
  res->setTypeLabel("BondOrder");
  return res;
};

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeDoubleOrAromaticBondQuery2() {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(true);
  res->setDataFunc(queryBondIsDoubleOrAromatic2);
  res->setDescription("DoubleOrAromaticBond");
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

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeSingleOrDoubleBondQuery2() {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrDouble2);
  res->setDescription("SingleOrDoubleBond");
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

RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *
makeSingleOrDoubleOrAromaticBondQuery2() {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(true);
  res->setDataFunc(queryBondIsSingleOrDoubleOrAromatic2);
  res->setDescription("SingleOrDoubleOrAromaticBond");
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
namespace {
// we don't use these anymore but we need to keep them around for backwards
// compatibility with pickled queries. There's no reason to update this list.
const std::vector<std::string> bondOrderQueryFunctions{
    std::string("BondOrder"), std::string("SingleOrAromaticBond"),
    std::string("DoubleOrAromaticBond"), std::string("SingleOrDoubleBond"),
    std::string("SingleOrDoubleOrAromaticBond")};
template <typename BondType>
bool hasBondTypeQueryHelper(const Queries::Query<int, BondType, true> &qry) {
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
}  // namespace

RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    const Queries::Query<int, Bond const *, true> &qry) {
  return hasBondTypeQueryHelper(qry);
}
RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    const Queries::Query<int, ConstRDMolBond, true> &qry) {
  return hasBondTypeQueryHelper(qry);
}

namespace {
template <typename BondType>
bool hasComplexBondTypeQueryHelper(
    const Queries::Query<int, BondType, true> &qry, bool seenBondOrder) {
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
RDKIT_GRAPHMOL_EXPORT bool hasComplexBondTypeQuery(
    const Queries::Query<int, ConstRDMolBond, true> &qry) {
  return hasComplexBondTypeQueryHelper(qry, false);
}
}  // namespace QueryOps

BOND_EQUALS_QUERY2 *makeBondDirEqualsQuery2(Bond::BondDir what) {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(what);
  res->setDataFunc(queryBondDir2);
  res->setDescription("BondDir");
  return res;
}
BOND_EQUALS_QUERY *makeBondDirEqualsQuery(Bond::BondDir what) {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryBondDir);
  res->setDescription("BondDir");
  return res;
}

BOND_EQUALS_QUERY2 *makeBondHasStereoQuery2() {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(true);
  res->setDataFunc(queryBondHasStereo2);
  res->setDescription("BondStereo");
  return res;
}
BOND_EQUALS_QUERY *makeBondHasStereoQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryBondHasStereo);
  res->setDescription("BondStereo");
  return res;
}

BOND_EQUALS_QUERY2 *makeBondIsInRingQuery2() {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(true);
  res->setDataFunc(queryIsBondInRing2);
  res->setDescription("BondInRing");
  return res;
}
BOND_EQUALS_QUERY *makeBondIsInRingQuery() {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(true);
  res->setDataFunc(queryIsBondInRing);
  res->setDescription("BondInRing");
  return res;
}

BOND_EQUALS_QUERY2 *makeBondInNRingsQuery2(int what) {
  auto *res = new BOND_EQUALS_QUERY2;
  res->setVal(what);
  res->setDataFunc(queryIsBondInNRings2);
  res->setDescription("BondInNRings");
  return res;
}
BOND_EQUALS_QUERY *makeBondInNRingsQuery(int what) {
  auto *res = new BOND_EQUALS_QUERY;
  res->setVal(what);
  res->setDataFunc(queryIsBondInNRings);
  res->setDescription("BondInNRings");
  return res;
}

BOND_NULL_QUERY2 *makeBondNullQuery2() {
  auto *res = new BOND_NULL_QUERY2;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("BondNull");
  return res;
}
BOND_NULL_QUERY *makeBondNullQuery() {
  auto *res = new BOND_NULL_QUERY;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("BondNull");
  return res;
}

ATOM_NULL_QUERY2 *makeAtomNullQuery2() {
  auto *res = new ATOM_NULL_QUERY2;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("AtomNull");
  return res;
}
ATOM_NULL_QUERY *makeAtomNullQuery() {
  auto *res = new ATOM_NULL_QUERY;
  res->setDataFunc(nullDataFun);
  res->setMatchFunc(nullQueryFun);
  res->setDescription("AtomNull");
  return res;
}

void convertComplexNameToQuery2(RDMolAtom atom, std::string_view symb) {
  if (symb == "Q") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeQAtomQuery2()));
  } else if (symb == "QH") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeQHAtomQuery2()));
  } else if (symb == "A") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeAAtomQuery2()));
  } else if (symb == "AH") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeAHAtomQuery2()));
  } else if (symb == "X") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeXAtomQuery2()));
  } else if (symb == "XH") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeXHAtomQuery2()));
  } else if (symb == "M") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeMAtomQuery2()));
  } else if (symb == "MH") {
    atom.mol().setAtomQuery(atom.index(), std::unique_ptr<RDMol::QUERYATOM_QUERY>(makeMHAtomQuery2()));
  } else {
    // we control what this function gets called with, so we should never land
    // here
    ASSERT_INVARIANT(0, "bad complex query symbol");
  }
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

namespace {
template<typename QUERY_T>
bool _complexBondQueryHelper(const QUERY_T* query) {
  if (query == nullptr) {
    return false;
  }
  // negated things are always complex:
  if (query->getNegation()) {
    return true;
  }
  std::string descr = query->getDescription();
  if (descr == "BondOrder" || descr == "SingleOrAromaticBond") {
    return false;
  }
  if (descr == "BondAnd" || descr == "BondXor") {
    return true;
  }
  if (descr == "BondOr") {
    // detect the types of queries that appear for unspecified bonds in
    // SMARTS:
    if (query->endChildren() - query->beginChildren() == 2) {
      for (auto child = query->beginChildren(); child != query->endChildren();
           ++child) {
        if ((*child)->getDescription() != "BondOrder" ||
            (*child)->getNegation()) {
          return true;
        }
        using BOND_EQUALS_QUERY_TYPE = Queries::EqualityQuery<int, typename QUERY_T::DATA_FUNC_ARG_TYPE, true>;
        if (static_cast<BOND_EQUALS_QUERY_TYPE *>(child->get())->getVal() !=
                Bond::SINGLE &&
            static_cast<BOND_EQUALS_QUERY_TYPE *>(child->get())->getVal() !=
                Bond::AROMATIC) {
          return true;
        }
      }
      return false;
    }
  }

  return true;
}
}  // namespace

bool isComplexQuery(ConstRDMolBond b) {
  PRECONDITION(b.index() < b.mol().getNumBonds(), "bad bond");
  const auto* query = b.mol().getBondQuery(b.index());
  return _complexBondQueryHelper(query);
}
bool isComplexQuery(const Bond *b) {
  PRECONDITION(b, "bad bond");
  if (!b->hasQuery()) {
    return false;
  }
  return _complexBondQueryHelper(b->getQuery());
}

namespace {
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

template <typename QUERY_T>
bool _isAtomListQuery(const QUERY_T *query, int atomicNum) {
  if (query == nullptr) {
    return false;
  }
  using EQUALS_QUERY_T = std::conditional_t<
      std::is_same_v<typename QUERY_T::DATA_FUNC_ARG_TYPE, ConstRDMolAtom>,
      ATOM_EQUALS_QUERY2,
      ATOM_EQUALS_QUERY>;
  if (query->getDescription() == "AtomOr") {
    for (const auto &child : boost::make_iterator_range(query->beginChildren(),
                                                        query->endChildren())) {
      if (!_atomListQueryHelper(child, false)) {
        return false;
      }
    }
    return true;
  } else if (query->getNegation() && _atomListQueryHelper(query, true)) {
    // this was github #5930: negated list queries containing a single atom were
    // being lost on output
    return true;
  } else if (query->getDescription() == "AtomAtomicNum" &&
             static_cast<const EQUALS_QUERY_T *>(query)->getVal() != atomicNum) {
    // when reading single-member atom lists from CTABs we end up with simple
    // AtomAtomicNum queries where the atomic number of the atom itself is zero.
    // Recognize this case.
    return true;
  }
  return false;
}
}  // namespace
bool isAtomListQuery(ConstRDMolAtom a) {
  PRECONDITION(a.index() < a.mol().getNumAtoms(), "bad atom");
  const auto *query = a.mol().getAtomQuery(a.index());
  return _isAtomListQuery(query, a.data().getAtomicNum());
}
bool isAtomListQuery(const Atom *a) {
  PRECONDITION(a, "bad atom");
  if (!a->hasQuery()) {
    return false;
  }
  return _isAtomListQuery(a->getQuery(), a->getAtomicNum());
}

namespace {
template <typename QUERY_T>
void _getAtomListQueryValsHelper(const QUERY_T *q, std::vector<int> &vals) {
  // list queries are series of nested ors of AtomAtomicNum queries
  PRECONDITION(q, "bad query");
  using ATOM_EQUALS_QUERY_TYPE = Queries::EqualityQuery<int, typename QUERY_T::DATA_FUNC_ARG_TYPE, true>;
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
        _getAtomListQueryValsHelper(child.get(), vals);
      } else if (descr == "AtomAtomicNum") {
        vals.push_back(static_cast<ATOM_EQUALS_QUERY_TYPE *>(child.get())->getVal());
      } else if (descr == "AtomType") {
        auto v = static_cast<ATOM_EQUALS_QUERY_TYPE *>(child.get())->getVal();
        // aromatic AtomType queries add 1000 to the atomic number;
        // correct for that:
        if (v >= 1000) {
          v -= 1000;
        }
        vals.push_back(v);
      }
    }
  } else if (descr == "AtomAtomicNum") {
    vals.push_back(static_cast<const ATOM_EQUALS_QUERY_TYPE *>(q)->getVal());
  } else if (descr == "AtomType") {
    auto v = static_cast<const ATOM_EQUALS_QUERY_TYPE *>(q)->getVal();
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
}  // namespace

void getAtomListQueryVals(const RDMol::QUERYATOM_QUERY *q,
                          std::vector<int> &vals) {
  _getAtomListQueryValsHelper(q, vals);
}
void getAtomListQueryVals(const Atom::QUERYATOM_QUERY *q,
                          std::vector<int> &vals) {
  _getAtomListQueryValsHelper(q, vals);
}

namespace {
template <typename QUERY_T>
bool _complexQueryHelper(const QUERY_T *query, bool &hasAtNum) {
  if (!query) {
    return false;
  }
  if (query->getNegation()) {
    return true;
  }
  const std::string &descr = query->getDescription();
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

template <typename QUERY_T>
bool _isComplexQuery(const QUERY_T *query) {
  if (query == nullptr) {
    return false;
  }
  // std::cerr<<"\n"<<a->getIdx();
  // negated things are always complex:
  if (query->getNegation()) {
    return true;
  }
  const std::string &descr = query->getDescription();
  // std::cerr<<" "<<descr;
  if (descr == "AtomNull" || descr == "AtomAtomicNum" || descr == "AtomType") {
    return false;
  }
  if (descr == "AtomOr" || descr == "AtomXor") {
    return true;
  }
  if (descr == "AtomAnd") {
    bool hasAtNum = false;
    if (_complexQueryHelper(query, hasAtNum)) {
      return true;
    }
    return !hasAtNum;
  }

  return true;
}
}  // namespace
bool isComplexQuery(ConstRDMolAtom a) {
  PRECONDITION(a.index() < a.mol().getNumAtoms(), "bad atom");
  const auto *query = a.mol().getAtomQuery(a.index());
  return _isComplexQuery(query);
}
bool isComplexQuery(const Atom *a) {
  PRECONDITION(a, "bad atom");
  if (!a->hasQuery()) {
    return false;
  }
  return _isComplexQuery(a->getQuery());
}

namespace {
template <typename QUERY_T, typename ATOM_T>
bool _isAtomAromatic(const QUERY_T *query, const ATOM_T &atom, const bool isAromaticDataMember) {
  if (query == nullptr) {
    return isAromaticAtom(atom);
  }
  const std::string &descr = query->getDescription();
  if (descr == "AtomAtomicNum") {
    return isAromaticDataMember;
  }
  if (descr == "AtomIsAromatic") {
    return !query->getNegation();
  }
  if (descr == "AtomIsAliphatic") {
    return query->getNegation();
  }
  if (descr == "AtomType") {
    using ATOM_EQUALS_QUERY_TYPE =
        Queries::EqualityQuery<int, typename QUERY_T::DATA_FUNC_ARG_TYPE, true>;
    bool res = getAtomTypeIsAromatic(
        static_cast<const ATOM_EQUALS_QUERY_TYPE *>(query)->getVal());
    if (query->getNegation()) {
      res = !res;
    }
    return res;
  }
  if (descr == "AtomAnd") {
    auto childIt = query->beginChildren();
    if ((*childIt)->getDescription() == "AtomAtomicNum") {
      if (query->getNegation()) {
        return false;
      } else if ((*(childIt + 1))->getDescription() == "AtomIsAliphatic") {
        return false;
      } else if ((*(childIt + 1))->getDescription() == "AtomIsAromatic") {
        return true;
      }
    }
  }

  return false;
}
}  // namespace
bool isAtomAromatic(ConstRDMolAtom a) {
  PRECONDITION(a.index() < a.mol().getNumAtoms(), "bad atom");
  const auto *query = a.mol().getAtomQuery(a.index());
  return _isAtomAromatic(query, a, a.data().getIsAromatic());
}
bool isAtomAromatic(const Atom *a) {
  PRECONDITION(a, "bad atom");
  if (!a->hasQuery()) {
    return isAromaticAtom(*a);
  }
  return _isAtomAromatic(a->getQuery(), *a, a->getIsAromatic());
}

namespace QueryOps {
namespace {
template <typename QUERY_T, typename ATOM_T>
void completeQueryAndChildren(QUERY_T *query, ATOM_T tgt,
                              unsigned int magicVal) {
  PRECONDITION(query, "no query");
  using ATOM_EQUALS_QUERY_TYPE = Queries::EqualityQuery<int, typename QUERY_T::DATA_FUNC_ARG_TYPE, true>;
  auto eqQuery = dynamic_cast<ATOM_EQUALS_QUERY_TYPE *>(query);
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
void completeMolQueries(RDMol *mol, unsigned int magicVal) {
  PRECONDITION(mol, "bad molecule");
  for (size_t atomIndex = 0, numAtoms = mol->getNumAtoms();
       atomIndex < numAtoms; ++atomIndex) {
    if (mol->hasAtomQuery(atomIndex)) {
      completeQueryAndChildren(mol->getAtomQuery(atomIndex),
                               ConstRDMolAtom(mol, atomIndex), magicVal);
    }
  }
}
void completeMolQueries(RWMol *mol, unsigned int magicVal) {
  PRECONDITION(mol, "bad molecule");
  for (auto atom : mol->atoms()) {
    if (atom->hasQuery()) {
      completeQueryAndChildren(atom->getQuery(), atom, magicVal);
    }
  }
}

void replaceAtomWithQueryAtom(RDMol *mol, atomindex_t atomIndex) {
  PRECONDITION(mol, "bad molecule");
  PRECONDITION(atomIndex < mol->getNumAtoms(), "bad atom");
  if (mol->hasAtomQuery(atomIndex)) {
    return;
  }

  // This is based on the QueryAtom constructor that accepts (const Atom &)
  auto &atom = mol->getAtom(atomIndex);
  std::unique_ptr<RDMol::QUERYATOM_QUERY> query(makeAtomNumQuery2(atom.getAtomicNum()));
  if (atom.getIsotope()) {
    std::unique_ptr<RDMol::QUERYATOM_QUERY> newQuery(
        makeAtomIsotopeQuery2(atom.getIsotope()));
    expandQuery(query, std::move(newQuery),
                      Queries::CompositeQueryType::COMPOSITE_AND);
  }
  if (atom.getFormalCharge()) {
    std::unique_ptr<RDMol::QUERYATOM_QUERY> newQuery(
        makeAtomFormalChargeQuery2(atom.getFormalCharge()));
    expandQuery(query, std::move(newQuery),
                      Queries::CompositeQueryType::COMPOSITE_AND);
  }
  if (atom.getNumRadicalElectrons()) {
    std::unique_ptr<RDMol::QUERYATOM_QUERY> newQuery(
        makeAtomNumRadicalElectronsQuery2(atom.getNumRadicalElectrons()));
    expandQuery(query, std::move(newQuery),
        Queries::CompositeQueryType::COMPOSITE_AND);
  }

  // This is based on the original replaceAtomWithQueryAtom
  bool massQueryPropValue = false;
  bool hasMassQueryProp = mol->getAtomPropIfPresent(
      common_properties::_hasMassQueryToken, atomIndex, massQueryPropValue);
  if (hasMassQueryProp) {
    std::unique_ptr<RDMol::QUERYATOM_QUERY> newQuery(
        makeAtomMassQuery2(static_cast<int>(atom.getMass())));
    expandQuery(query, std::move(newQuery),
                Queries::CompositeQueryType::COMPOSITE_AND);
  }
  mol->setAtomQuery(atomIndex, std::move(query));
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

static const std::unordered_map<std::string,
    std::pair<DataFuncPtrType<false, false>, DataFuncPtrType<false, true>>>
    atom_descr_to_func = {
        {"AtomRingBondCount", {queryAtomRingBondCount, queryAtomRingBondCount2}},
        {"AtomHasRingBond", {queryAtomHasRingBond, queryAtomHasRingBond2}},
        {"AtomRingSize", {nullptr, nullptr}},
        {"AtomMinRingSize", {queryAtomMinRingSize, queryAtomMinRingSize2}},
        {"AtomImplicitValence", {queryAtomImplicitValence, queryAtomImplicitValence2}},
        {"AtomTotalValence", {queryAtomTotalValence, queryAtomTotalValence2}},
        {"AtomAtomicNum", {queryAtomNum, queryAtomNum2}},
        {"AtomExplicitDegree", {queryAtomExplicitDegree, queryAtomExplicitDegree2}},
        {"AtomTotalDegree", {queryAtomTotalDegree, queryAtomTotalDegree2}},
        {"AtomHeavyAtomDegree", {queryAtomHeavyAtomDegree, queryAtomHeavyAtomDegree2}},
        {"AtomHCount", {queryAtomHCount, queryAtomHCount2}},
        {"AtomImplicitHCount", {queryAtomImplicitHCount, queryAtomImplicitHCount2}},
        {"AtomHasImplicitH", {queryAtomHasImplicitH, queryAtomHasImplicitH2}},
        {"AtomIsAromatic", {queryAtomAromatic, queryAtomAromatic2}},
        {"AtomIsAliphatic", {queryAtomAliphatic, queryAtomAliphatic2}},
        {"AtomUnsaturated", {queryAtomUnsaturated, queryAtomUnsaturated2}},
        {"AtomMass", {queryAtomMass, queryAtomMass2}},
        {"AtomIsotope", {queryAtomIsotope, queryAtomIsotope2}},
        {"AtomFormalCharge", {queryAtomFormalCharge, queryAtomFormalCharge2}},
        {"AtomNegativeFormalCharge", {queryAtomNegativeFormalCharge, queryAtomNegativeFormalCharge2}},
        {"AtomHybridization", {queryAtomHybridization, queryAtomHybridization2}},
        {"AtomInRing", {queryIsAtomInRing, queryIsAtomInRing2}},
        {"AtomInNRings", {queryIsAtomInNRings, queryIsAtomInNRings2}},
        {"AtomHasHeteroatomNeighbors", {queryAtomHasHeteroatomNbrs, queryAtomHasHeteroatomNbrs2}},
        {"AtomNumHeteroatomNeighbors", {queryAtomNumHeteroatomNbrs, queryAtomNumHeteroatomNbrs2}},
        {"AtomNonHydrogenDegree", {queryAtomNonHydrogenDegree, queryAtomNonHydrogenDegree2}},
        {"AtomHasAliphaticHeteroatomNeighbors", {queryAtomHasAliphaticHeteroatomNbrs, queryAtomHasAliphaticHeteroatomNbrs2}},
        {"AtomNumAliphaticHeteroatomNeighbors", {queryAtomNumAliphaticHeteroatomNbrs, queryAtomNumAliphaticHeteroatomNbrs2}},
        {"AtomNull", {nullDataFun, nullDataFun}},
        {"AtomType", {queryAtomType, queryAtomType2}},
        {"AtomNumRadicalElectrons", {queryAtomNumRadicalElectrons, queryAtomNumRadicalElectrons2}},
        {"RecursiveStructure", {nullptr, nullptr}},
        {"AtomAnd", {nullptr, nullptr}},
        {"AtomOr", {nullptr, nullptr}},
        {"AtomXor", {nullptr, nullptr}},
        {"HasProp", {nullptr, nullptr}},
        {"HasPropWithValue", {nullptr, nullptr}},
};

void finalizeQueryFromDescription(
    Queries::Query<int, ConstRDMolAtom, true> *query, ConstRDMolAtom) {
  std::string descr = query->getDescription();

  if (boost::starts_with(descr, "range_")) {
    descr = descr.substr(6);
  } else if (boost::starts_with(descr, "less_")) {
    descr = descr.substr(5);
  } else if (boost::starts_with(descr, "greater_")) {
    descr = descr.substr(8);
  }

  auto it = atom_descr_to_func.find(descr);
  if (it != atom_descr_to_func.end()) {
    if (it->second.second != nullptr) {
      query->setDataFunc(it->second.second);
      if (descr == "AtomNull") {
        query->setMatchFunc(nullQueryFun);
      }
    } else if (descr == "AtomRingSize") {
      auto *dataFunc = getQueryIsInRingOfSizeTmpl<false, true>(
          static_cast<ATOM_EQUALS_QUERY2 *>(query)->getVal());
      query->setDataFunc(dataFunc);
    } else {
      // don't need to do anything here because the classes
      // automatically have everything set
    }
  } else {
    throw ValueErrorException("Do not know how to finalize query: '" + descr +
                              "'");
  }
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

  auto it = atom_descr_to_func.find(descr);
  if (it != atom_descr_to_func.end()) {
    if (it->second.first != nullptr) {
      query->setDataFunc(it->second.first);
      if (descr == "AtomNull") {
        query->setMatchFunc(nullQueryFun);
      }
    } else if (descr == "AtomRingSize") {
      auto *dataFunc = getQueryIsInRingOfSizeTmpl<false, false>(
          static_cast<ATOM_EQUALS_QUERY *>(query)->getVal());
      query->setDataFunc(dataFunc);
    } else {
      // don't need to do anything here because the classes
      // automatically have everything set
    }
  } else {
    throw ValueErrorException("Do not know how to finalize query: '" + descr +
                              "'");
  }
}

static const std::unordered_map<std::string,
    std::pair<DataFuncPtrType<true, false>, DataFuncPtrType<true, true>>>
    bond_descr_to_func = {
        {"BondRingSize", {nullptr, nullptr}},
        {"BondMinRingSize", {queryBondMinRingSize, queryBondMinRingSize2}},
        {"BondOrder", {queryBondOrder, queryBondOrder2}},
        {"BondDir", {queryBondDir, queryBondDir2}},
        {"BondInRing", {queryIsBondInRing, queryIsBondInRing2}},
        {"BondInNRings", {queryIsBondInNRings, queryIsBondInNRings2}},
        {"SingleOrAromaticBond", {queryBondIsSingleOrAromatic, queryBondIsSingleOrAromatic2}},
        {"SingleOrDoubleBond", {queryBondIsSingleOrDouble, queryBondIsSingleOrDouble2}},
        {"DoubleOrAromaticBond", {queryBondIsDoubleOrAromatic, queryBondIsDoubleOrAromatic2}},
        {"SingleOrDoubleOrAromaticBond", {queryBondIsSingleOrDoubleOrAromatic, queryBondIsSingleOrDoubleOrAromatic2}},
        {"BondNull", {nullDataFun, nullDataFun}},
        {"BondAnd", {nullptr, nullptr}},
        {"BondOr", {nullptr, nullptr}},
        {"BondXor", {nullptr, nullptr}},
        {"HasProp", {nullptr, nullptr}},
        {"HasPropWithValue", {nullptr, nullptr}},
};

void finalizeQueryFromDescription(
    Queries::Query<int, ConstRDMolBond, true> *query, ConstRDMolBond) {
  const std::string &descr = query->getDescription();
  auto it = bond_descr_to_func.find(descr);
  if (it != bond_descr_to_func.end()) {
    if (it->second.second != nullptr) {
      query->setDataFunc(it->second.second);
      if (descr == "BondNull") {
        query->setMatchFunc(nullQueryFun);
      }
    } else if (descr == "BondRingSize") {
      auto *dataFunc = getQueryIsInRingOfSizeTmpl<true, true>(
          static_cast<BOND_EQUALS_QUERY2 *>(query)->getVal());
      query->setDataFunc(dataFunc);
    } else {
      // don't need to do anything here because the classes
      // automatically have everything set
    }
  } else {
    throw ValueErrorException("Do not know how to finalize query: '" + descr +
                              "'");
  }
}
void finalizeQueryFromDescription(
    Queries::Query<int, Bond const *, true> *query, Bond const *) {
  const std::string &descr = query->getDescription();
  auto it = bond_descr_to_func.find(descr);
  if (it != bond_descr_to_func.end()) {
    if (it->second.first != nullptr) {
      query->setDataFunc(it->second.first);
      if (descr == "BondNull") {
        query->setMatchFunc(nullQueryFun);
      }
    } else if (descr == "BondRingSize") {
      auto *dataFunc = getQueryIsInRingOfSizeTmpl<true, false>(
          static_cast<BOND_EQUALS_QUERY *>(query)->getVal());
      query->setDataFunc(dataFunc);
    } else {
      // don't need to do anything here because the classes
      // automatically have everything set
    }
  } else {
    throw ValueErrorException("Do not know how to finalize query: '" + descr +
                              "'");
  }
}

bool isMetal(const Atom &atom) {
  static const std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
  return q->Match(&atom);
}

}  // namespace QueryOps
};  // namespace RDKit
