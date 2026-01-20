//
//  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

//! \file QueryOps.h
/*!
    \brief Includes a bunch of functionality for handling Atom and Bond queries.
*/
#include <RDGeneral/export.h>
#include <RDGeneral/Dict.h>
#ifndef RD_QUERY_OPS_H
#define RD_QUERY_OPS_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDMol.h>
#include <Query/QueryObjects.h>
#include <Query/Query.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#include <utility>
#endif

namespace RDKit {
class RDMol;
class RingInfoCache;

typedef Queries::Query<bool, ConstRDMolAtom, true> ATOM_BOOL_QUERY2;
typedef Queries::Query<bool, ConstRDMolBond, true> BOND_BOOL_QUERY2;
typedef Queries::Query<bool, Atom const *, true> ATOM_BOOL_QUERY;
typedef Queries::Query<bool, Bond const *, true> BOND_BOOL_QUERY;

typedef Queries::AndQuery<int, ConstRDMolAtom, true> ATOM_AND_QUERY2;
typedef Queries::AndQuery<int, ConstRDMolBond, true> BOND_AND_QUERY2;
typedef Queries::AndQuery<int, Atom const *, true> ATOM_AND_QUERY;
typedef Queries::AndQuery<int, Bond const *, true> BOND_AND_QUERY;

typedef Queries::OrQuery<int, ConstRDMolAtom, true> ATOM_OR_QUERY2;
typedef Queries::OrQuery<int, ConstRDMolBond, true> BOND_OR_QUERY2;
typedef Queries::OrQuery<int, Atom const *, true> ATOM_OR_QUERY;
typedef Queries::OrQuery<int, Bond const *, true> BOND_OR_QUERY;

typedef Queries::XOrQuery<int, ConstRDMolAtom, true> ATOM_XOR_QUERY2;
typedef Queries::XOrQuery<int, ConstRDMolBond, true> BOND_XOR_QUERY2;
typedef Queries::XOrQuery<int, Atom const *, true> ATOM_XOR_QUERY;
typedef Queries::XOrQuery<int, Bond const *, true> BOND_XOR_QUERY;

typedef Queries::EqualityQuery<int, ConstRDMolAtom, true> ATOM_EQUALS_QUERY2;
typedef Queries::EqualityQuery<int, ConstRDMolBond, true> BOND_EQUALS_QUERY2;
typedef Queries::EqualityQuery<int, Atom const *, true> ATOM_EQUALS_QUERY;
typedef Queries::EqualityQuery<int, Bond const *, true> BOND_EQUALS_QUERY;

typedef Queries::GreaterQuery<int, ConstRDMolAtom, true> ATOM_GREATER_QUERY2;
typedef Queries::GreaterQuery<int, ConstRDMolBond, true> BOND_GREATER_QUERY2;
typedef Queries::GreaterQuery<int, Atom const *, true> ATOM_GREATER_QUERY;
typedef Queries::GreaterQuery<int, Bond const *, true> BOND_GREATER_QUERY;

typedef Queries::GreaterEqualQuery<int, ConstRDMolAtom, true>
    ATOM_GREATEREQUAL_QUERY2;
typedef Queries::GreaterEqualQuery<int, ConstRDMolBond, true>
    BOND_GREATEREQUAL_QUERY2;
typedef Queries::GreaterEqualQuery<int, Atom const *, true>
    ATOM_GREATEREQUAL_QUERY;
typedef Queries::GreaterEqualQuery<int, Bond const *, true>
    BOND_GREATEREQUAL_QUERY;

typedef Queries::LessQuery<int, ConstRDMolAtom, true> ATOM_LESS_QUERY2;
typedef Queries::LessQuery<int, ConstRDMolBond, true> BOND_LESS_QUERY2;
typedef Queries::LessQuery<int, Atom const *, true> ATOM_LESS_QUERY;
typedef Queries::LessQuery<int, Bond const *, true> BOND_LESS_QUERY;

typedef Queries::LessEqualQuery<int, ConstRDMolAtom, true>
    ATOM_LESSEQUAL_QUERY2;
typedef Queries::LessEqualQuery<int, ConstRDMolBond, true>
    BOND_LESSEQUAL_QUERY2;
typedef Queries::LessEqualQuery<int, Atom const *, true> ATOM_LESSEQUAL_QUERY;
typedef Queries::LessEqualQuery<int, Bond const *, true> BOND_LESSEQUAL_QUERY;

typedef Queries::RangeQuery<int, ConstRDMolAtom, true> ATOM_RANGE_QUERY2;
typedef Queries::RangeQuery<int, ConstRDMolBond, true> BOND_RANGE_QUERY2;
typedef Queries::RangeQuery<int, Atom const *, true> ATOM_RANGE_QUERY;
typedef Queries::RangeQuery<int, Bond const *, true> BOND_RANGE_QUERY;

typedef Queries::SetQuery<int, ConstRDMolAtom, true> ATOM_SET_QUERY2;
typedef Queries::SetQuery<int, ConstRDMolBond, true> BOND_SET_QUERY2;
typedef Queries::SetQuery<int, Atom const *, true> ATOM_SET_QUERY;
typedef Queries::SetQuery<int, Bond const *, true> BOND_SET_QUERY;

typedef Queries::Query<int, ConstRDMolAtom, true> ATOM_NULL_QUERY2;
typedef Queries::Query<int, ConstRDMolBond, true> BOND_NULL_QUERY2;
typedef Queries::Query<int, Atom const *, true> ATOM_NULL_QUERY;
typedef Queries::Query<int, Bond const *, true> BOND_NULL_QUERY;

// -------------------------------------------------
// common atom queries

static inline int queryAtomAromatic2(ConstRDMolAtom at) {
  return at.data().getIsAromatic();
};
static inline int queryAtomAromatic(Atom const *at) {
  return at->getIsAromatic();
};
static inline int queryAtomAliphatic2(ConstRDMolAtom at) {
  return !(at.data().getIsAromatic());
};
static inline int queryAtomAliphatic(Atom const *at) {
  return !(at->getIsAromatic());
};
static inline int queryAtomExplicitDegree2(ConstRDMolAtom at) {
  return at.mol().getAtomDegree(at.index());
};
static inline int queryAtomExplicitDegree(Atom const *at) {
  return at->getDegree();
};
static inline int queryAtomTotalDegree2(ConstRDMolAtom at) {
  return at.mol().getAtomTotalDegree(at.index());
};
static inline int queryAtomTotalDegree(Atom const *at) {
  return at->getTotalDegree();
};
//! D and T are treated as "non-hydrogen" here
static inline int queryAtomNonHydrogenDegree2(ConstRDMolAtom at) {
  int res = 0;
  for (const auto nbri :
       boost::make_iterator_range(at.mol().getAtomNeighbors(at.index()))) {
    const auto nbr = at.mol().getAtom(nbri);
    if (nbr.getAtomicNum() != 1 || nbr.getIsotope() > 1) {
      res++;
    }
  }
  return res;
};
//! D and T are treated as "non-hydrogen" here
static inline int queryAtomNonHydrogenDegree(Atom const *at) {
  int res = 0;
  for (const auto nbri :
       boost::make_iterator_range(at->getOwningMol().getAtomNeighbors(at))) {
    const auto nbr = at->getOwningMol()[nbri];
    if (nbr->getAtomicNum() != 1 || nbr->getIsotope() > 1) {
      res++;
    }
  }

  return res;
};
//! D and T are not treated as heavy atoms here
static inline int queryAtomHeavyAtomDegree2(ConstRDMolAtom at) {
  int heavyDegree = 0;
  for (const auto nbri :
       boost::make_iterator_range(at.mol().getAtomNeighbors(at.index()))) {
    const auto nbr = at.mol().getAtom(nbri);
    if (nbr.getAtomicNum() > 1) {
      heavyDegree++;
    }
  }

  return heavyDegree;
};
//! D and T are not treated as heavy atoms here
static inline int queryAtomHeavyAtomDegree(Atom const *at) {
  int heavyDegree = 0;
  for (const auto nbri :
       boost::make_iterator_range(at->getOwningMol().getAtomNeighbors(at))) {
    const auto nbr = at->getOwningMol()[nbri];
    if (nbr->getAtomicNum() > 1) {
      heavyDegree++;
    }
  }

  return heavyDegree;
};
static inline int queryAtomHCount2(ConstRDMolAtom at) {
  return at.mol().getTotalNumHs(at.index(), true);
};
static inline int queryAtomHCount(Atom const *at) {
  return at->getTotalNumHs(true);
};
static inline int queryAtomImplicitHCount2(ConstRDMolAtom at) {
  return at.mol().getTotalNumHs(at.index(), false);
};
static inline int queryAtomImplicitHCount(Atom const *at) {
  return at->getTotalNumHs(false);
};
static inline int queryAtomHasImplicitH2(ConstRDMolAtom at) {
  return int(at.mol().getTotalNumHs(at.index(), false) > 0);
};
static inline int queryAtomHasImplicitH(Atom const *at) {
  return int(at->getTotalNumHs(false) > 0);
};
static inline int queryAtomImplicitValence2(ConstRDMolAtom at) {
  return at.data().getImplicitValence();
};
static inline int queryAtomImplicitValence(Atom const *at) {
  return at->getValence(Atom::ValenceType::IMPLICIT);
};
static inline int queryAtomExplicitValence2(ConstRDMolAtom at) {
  const AtomData &data = at.data();
  return data.getExplicitValence() - data.getNumExplicitHs();
};
static inline int queryAtomExplicitValence(Atom const *at) {
  return at->getValence(Atom::ValenceType::EXPLICIT) - at->getNumExplicitHs();
};
static inline int queryAtomTotalValence2(ConstRDMolAtom at) {
  return at.data().getTotalValence();
};
static inline int queryAtomTotalValence(Atom const *at) {
  return at->getTotalValence();
};
static inline int queryAtomUnsaturated2(ConstRDMolAtom at) {
  const AtomData &data = at.data();
  return at.mol().getAtomTotalDegree(at.index()) < data.getTotalValence();
};
static inline int queryAtomUnsaturated(Atom const *at) {
  return at->getTotalDegree() < at->getTotalValence();
};
static inline int queryAtomNum2(ConstRDMolAtom at) {
  return at.data().getAtomicNum();
}
static inline int queryAtomNum(Atom const *at) { return at->getAtomicNum(); }
static inline int makeAtomType(int atomic_num, bool aromatic) {
  return atomic_num + 1000 * static_cast<int>(aromatic);
}
static inline void parseAtomType(int val, int &atomic_num, bool &aromatic) {
  if (val > 1000) {
    aromatic = true;
    atomic_num = val - 1000;
  } else {
    aromatic = false;
    atomic_num = val;
  }
}
static inline bool getAtomTypeIsAromatic(int val) { return val > 1000; }
static inline int getAtomTypeAtomicNum(int val) {
  if (val > 1000) {
    return val - 1000;
  }
  return val;
}

static inline int queryAtomType2(ConstRDMolAtom at) {
  const AtomData &data = at.data();
  return makeAtomType(data.getAtomicNum(), data.getIsAromatic());
};
static inline int queryAtomType(Atom const *at) {
  return makeAtomType(at->getAtomicNum(), at->getIsAromatic());
};
const int massIntegerConversionFactor = 1000;
static inline int queryAtomMass2(ConstRDMolAtom at) {
  return static_cast<int>(
      std::round(massIntegerConversionFactor * at.data().getMass()));
};
static inline int queryAtomMass(Atom const *at) {
  return static_cast<int>(
      std::round(massIntegerConversionFactor * at->getMass()));
};
static inline int queryAtomIsotope2(ConstRDMolAtom at) {
  return static_cast<int>(at.data().getIsotope());
};
static inline int queryAtomIsotope(Atom const *at) {
  return static_cast<int>(at->getIsotope());
};
static inline int queryAtomFormalCharge2(ConstRDMolAtom at) {
  return static_cast<int>(at.data().getFormalCharge());
};
static inline int queryAtomFormalCharge(Atom const *at) {
  return static_cast<int>(at->getFormalCharge());
};
static inline int queryAtomNegativeFormalCharge2(ConstRDMolAtom at) {
  return static_cast<int>(-1 * at.data().getFormalCharge());
};
static inline int queryAtomNegativeFormalCharge(Atom const *at) {
  return static_cast<int>(-1 * at->getFormalCharge());
};
static inline int queryAtomHybridization2(ConstRDMolAtom at) {
  return at.data().getHybridization();
};
static inline int queryAtomHybridization(Atom const *at) {
  return at->getHybridization();
};
static inline int queryAtomNumRadicalElectrons2(ConstRDMolAtom at) {
  return at.data().getNumRadicalElectrons();
};
static inline int queryAtomNumRadicalElectrons(Atom const *at) {
  return at->getNumRadicalElectrons();
};
static inline int queryAtomHasChiralTag2(ConstRDMolAtom at) {
  return at.data().getChiralTag() != AtomEnums::ChiralType::CHI_UNSPECIFIED;
};
static inline int queryAtomHasChiralTag(Atom const *at) {
  return at->getChiralTag() != Atom::CHI_UNSPECIFIED;
};
static inline int queryAtomMissingChiralTag2(ConstRDMolAtom at) {
  if (at.data().getChiralTag() != AtomEnums::ChiralType::CHI_UNSPECIFIED) {
    return false;
  }
  bool value = false;
  bool present = at.mol().getAtomPropIfPresent(
      common_properties::_ChiralityPossibleToken, at.index(), value);
  return present && value;
};
static inline int queryAtomMissingChiralTag(Atom const *at) {
  return at->getChiralTag() == Atom::CHI_UNSPECIFIED &&
         at->hasProp(common_properties::_ChiralityPossible);
};

static inline int queryAtomHasHeteroatomNbrs2(ConstRDMolAtom at) {
  auto [nbrIdx, endNbrs] = at.mol().getAtomNeighbors(at.index());
  while (nbrIdx != endNbrs) {
    const AtomData &nbr = at.mol().getAtom(*nbrIdx);
    if (nbr.getAtomicNum() != 6 && nbr.getAtomicNum() != 1) {
      return 1;
    }
    ++nbrIdx;
  }
  return 0;
};
static inline int queryAtomHasHeteroatomNbrs(Atom const *at) {
  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
  while (nbrIdx != endNbrs) {
    const Atom *nbr = at->getOwningMol()[*nbrIdx];
    if (nbr->getAtomicNum() != 6 && nbr->getAtomicNum() != 1) {
      return 1;
    }
    ++nbrIdx;
  }
  return 0;
};

static inline int queryAtomNumHeteroatomNbrs2(ConstRDMolAtom at) {
  int res = 0;
  auto [nbrIdx, endNbrs] = at.mol().getAtomNeighbors(at.index());
  while (nbrIdx != endNbrs) {
    const AtomData &nbr = at.mol().getAtom(*nbrIdx);
    if (nbr.getAtomicNum() != 6 && nbr.getAtomicNum() != 1) {
      ++res;
    }
    ++nbrIdx;
  }
  return res;
};
static inline int queryAtomNumHeteroatomNbrs(Atom const *at) {
  int res = 0;
  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
  while (nbrIdx != endNbrs) {
    const Atom *nbr = at->getOwningMol()[*nbrIdx];
    if (nbr->getAtomicNum() != 6 && nbr->getAtomicNum() != 1) {
      ++res;
    }
    ++nbrIdx;
  }
  return res;
};

static inline int queryAtomHasAliphaticHeteroatomNbrs2(ConstRDMolAtom at) {
  auto [nbrIdx, endNbrs] = at.mol().getAtomNeighbors(at.index());
  while (nbrIdx != endNbrs) {
    const AtomData &nbr = at.mol().getAtom(*nbrIdx);
    if ((!nbr.getIsAromatic()) && nbr.getAtomicNum() != 6 &&
        nbr.getAtomicNum() != 1) {
      return 1;
    }
    ++nbrIdx;
  }
  return 0;
};
static inline int queryAtomHasAliphaticHeteroatomNbrs(Atom const *at) {
  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
  while (nbrIdx != endNbrs) {
    const Atom *nbr = at->getOwningMol()[*nbrIdx];
    if ((!nbr->getIsAromatic()) && nbr->getAtomicNum() != 6 &&
        nbr->getAtomicNum() != 1) {
      return 1;
    }
    ++nbrIdx;
  }
  return 0;
};

static inline int queryAtomNumAliphaticHeteroatomNbrs2(ConstRDMolAtom at) {
  int res = 0;
  auto [nbrIdx, endNbrs] = at.mol().getAtomNeighbors(at.index());
  while (nbrIdx != endNbrs) {
    const AtomData &nbr = at.mol().getAtom(*nbrIdx);
    if ((!nbr.getIsAromatic()) && nbr.getAtomicNum() != 6 &&
        nbr.getAtomicNum() != 1) {
      ++res;
    }
    ++nbrIdx;
  }
  return res;
};
static inline int queryAtomNumAliphaticHeteroatomNbrs(Atom const *at) {
  int res = 0;
  ROMol::ADJ_ITER nbrIdx, endNbrs;
  boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
  while (nbrIdx != endNbrs) {
    const Atom *nbr = at->getOwningMol()[*nbrIdx];
    if ((!nbr->getIsAromatic()) && nbr->getAtomicNum() != 6 &&
        nbr->getAtomicNum() != 1) {
      ++res;
    }
    ++nbrIdx;
  }
  return res;
};

RDKIT_GRAPHMOL_EXPORT unsigned int queryAtomBondProduct2(ConstRDMolAtom at);
RDKIT_GRAPHMOL_EXPORT unsigned int queryAtomBondProduct(Atom const *at);
RDKIT_GRAPHMOL_EXPORT unsigned int queryAtomAllBondProduct2(ConstRDMolAtom at);
RDKIT_GRAPHMOL_EXPORT unsigned int queryAtomAllBondProduct(Atom const *at);

// -------------------------------------------------
// common bond queries

static inline int queryBondOrder2(ConstRDMolBond bond) {
  return static_cast<int>(bond.data().getBondType());
};
static inline int queryBondOrder(Bond const *bond) {
  return static_cast<int>(bond->getBondType());
};
static inline int queryBondIsSingleOrAromatic2(ConstRDMolBond bond) {
  auto type = bond.data().getBondType();
  return static_cast<int>(type == BondEnums::BondType::SINGLE ||
                          type == BondEnums::BondType::AROMATIC);
};
static inline int queryBondIsSingleOrAromatic(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
                          bond->getBondType() == Bond::AROMATIC);
};
static inline int queryBondIsDoubleOrAromatic2(ConstRDMolBond bond) {
  auto type = bond.data().getBondType();
  return static_cast<int>(type == BondEnums::BondType::DOUBLE ||
                          type == BondEnums::BondType::AROMATIC);
};
static inline int queryBondIsDoubleOrAromatic(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::DOUBLE ||
                          bond->getBondType() == Bond::AROMATIC);
};
static inline int queryBondIsSingleOrDouble2(ConstRDMolBond bond) {
  auto type = bond.data().getBondType();
  return static_cast<int>(type == BondEnums::BondType::SINGLE ||
                          type == BondEnums::BondType::DOUBLE);
};
static inline int queryBondIsSingleOrDouble(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
                          bond->getBondType() == Bond::DOUBLE);
};
static inline int queryBondIsSingleOrDoubleOrAromatic2(ConstRDMolBond bond) {
  auto type = bond.data().getBondType();
  return static_cast<int>(type == BondEnums::BondType::SINGLE ||
                          type == BondEnums::BondType::DOUBLE ||
                          type == BondEnums::BondType::AROMATIC);
};
static inline int queryBondIsSingleOrDoubleOrAromatic(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
                          bond->getBondType() == Bond::DOUBLE ||
                          bond->getBondType() == Bond::AROMATIC);
};
static inline int queryBondDir2(ConstRDMolBond bond) {
  return static_cast<int>(bond.data().getBondDir());
};
static inline int queryBondDir(Bond const *bond) {
  return static_cast<int>(bond->getBondDir());
};
static inline int queryIsBondInNRings2(ConstRDMolBond at) {
  return at.mol().getRingInfo().numBondRings(at.index());
};
static inline int queryIsBondInNRings(Bond const *at) {
  return at->getOwningMol().getRingInfo()->numBondRings(at->getIdx());
};
static inline int queryBondHasStereo2(ConstRDMolBond bnd) {
  return bnd.data().getStereo() > BondEnums::BondStereo::STEREONONE;
};
static inline int queryBondHasStereo(Bond const *bnd) {
  return bnd->getStereo() > Bond::STEREONONE;
};

// -------------------------------------------------
// ring queries

static inline int queryIsAtomInNRings2(ConstRDMolAtom at) {
  return at.mol().getRingInfo().numAtomRings(at.index());
};
static inline int queryIsAtomInNRings(Atom const *at) {
  return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx());
};
static inline int queryIsAtomInRing2(ConstRDMolAtom at) {
  return at.mol().getRingInfo().numAtomRings(at.index()) != 0;
};
static inline int queryIsAtomInRing(Atom const *at) {
  return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx()) != 0;
};
static inline int queryAtomHasRingBond2(ConstRDMolAtom at) {
  auto [begin, end] = at.mol().getAtomBonds(at.index());
  while (begin != end) {
    unsigned int bondIdx = *begin;
    if (at.mol().getRingInfo().numBondRings(bondIdx)) {
      return 1;
    }
    ++begin;
  }
  return 0;
};
static inline int queryAtomHasRingBond(Atom const *at) {
  ROMol::OBOND_ITER_PAIR atomBonds = at->getOwningMol().getAtomBonds(at);
  while (atomBonds.first != atomBonds.second) {
    unsigned int bondIdx =
        at->getOwningMol().getTopology()[*atomBonds.first]->getIdx();
    if (at->getOwningMol().getRingInfo()->numBondRings(bondIdx)) {
      return 1;
    }
    ++atomBonds.first;
  }
  return 0;
};
RDKIT_GRAPHMOL_EXPORT bool queryIsAtomBridgeheadInternal(
    const RDMol &mol, atomindex_t atomIndex, const RingInfoCache &rings);
static inline int queryIsAtomBridgehead2(ConstRDMolAtom at) {
  return queryIsAtomBridgeheadInternal(at.mol(), at.index(),
                                       at.mol().getRingInfo())
             ? 1
             : 0;
}
RDKIT_GRAPHMOL_EXPORT int queryIsAtomBridgehead(Atom const *at);

static inline int queryIsBondInRing2(ConstRDMolBond bond) {
  return bond.mol().getRingInfo().numBondRings(bond.index()) != 0;
};
static inline int queryIsBondInRing(Bond const *bond) {
  return bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx()) != 0;
};
static inline int queryAtomMinRingSize2(ConstRDMolAtom at) {
  return at.mol().getRingInfo().minAtomRingSize(at.index());
};
static inline int queryAtomMinRingSize(Atom const *at) {
  return at->getOwningMol().getRingInfo()->minAtomRingSize(at->getIdx());
};
static inline int queryBondMinRingSize2(ConstRDMolBond bond) {
  return bond.mol().getRingInfo().minBondRingSize(bond.index());
};
static inline int queryBondMinRingSize(Bond const *bond) {
  return bond->getOwningMol().getRingInfo()->minBondRingSize(bond->getIdx());
};

static inline int queryAtomRingBondCount2(ConstRDMolAtom at) {
  // EFF: cache this result
  int res = 0;
  auto [begin, end] = at.mol().getAtomBonds(at.index());
  while (begin != end) {
    unsigned int bondIdx = *begin;
    if (at.mol().getRingInfo().numBondRings(bondIdx)) {
      res++;
    }
    ++begin;
  }
  return res;
}
static inline int queryAtomRingBondCount(Atom const *at) {
  // EFF: cache this result
  int res = 0;
  ROMol::OBOND_ITER_PAIR atomBonds = at->getOwningMol().getAtomBonds(at);
  while (atomBonds.first != atomBonds.second) {
    unsigned int bondIdx =
        at->getOwningMol().getTopology()[*atomBonds.first]->getIdx();
    if (at->getOwningMol().getRingInfo()->numBondRings(bondIdx)) {
      res++;
    }
    ++atomBonds.first;
  }
  return res;
}

template <int tgt>
int queryAtomIsInRingOfSize2(ConstRDMolAtom at) {
  if (at.mol().getRingInfo().isAtomInRingOfSize(at.index(), tgt)) {
    return tgt;
  } else {
    return 0;
  }
};
template <int tgt>
int queryAtomIsInRingOfSize(Atom const *at) {
  if (at->getOwningMol().getRingInfo()->isAtomInRingOfSize(at->getIdx(), tgt)) {
    return tgt;
  } else {
    return 0;
  }
};
template <int tgt>
int queryBondIsInRingOfSize2(ConstRDMolBond bond) {
  if (bond.mol().getRingInfo().isBondInRingOfSize(bond.index(), tgt)) {
    return tgt;
  } else {
    return 0;
  }
};
template <int tgt>
int queryBondIsInRingOfSize(Bond const *bond) {
  if (bond->getOwningMol().getRingInfo()->isBondInRingOfSize(bond->getIdx(),
                                                             tgt)) {
    return tgt;
  } else {
    return 0;
  }
};

template <class T>
T *makeAtomSimpleQuery2(int what, int func(ConstRDMolAtom),
                        const std::string &description = "Atom Simple") {
  T *res = new T;
  res->setVal(what);
  res->setDataFunc(func);
  res->setDescription(description);
  return res;
}
template <class T>
T *makeAtomSimpleQuery(int what, int func(Atom const *),
                       const std::string &description = "Atom Simple") {
  T *res = new T;
  res->setVal(what);
  res->setDataFunc(func);
  res->setDescription(description);
  return res;
}

static inline ATOM_RANGE_QUERY2 *makeAtomRangeQuery2(
    int lower, int upper, bool lowerOpen, bool upperOpen,
    int func(ConstRDMolAtom), const std::string &description = "Atom Range") {
  ATOM_RANGE_QUERY2 *res = new ATOM_RANGE_QUERY2(lower, upper);
  res->setDataFunc(func);
  res->setDescription(description);
  res->setEndsOpen(lowerOpen, upperOpen);
  return res;
}
static inline ATOM_RANGE_QUERY *makeAtomRangeQuery(
    int lower, int upper, bool lowerOpen, bool upperOpen,
    int func(Atom const *), const std::string &description = "Atom Range") {
  ATOM_RANGE_QUERY *res = new ATOM_RANGE_QUERY(lower, upper);
  res->setDataFunc(func);
  res->setDescription(description);
  res->setEndsOpen(lowerOpen, upperOpen);
  return res;
}

//! returns a Query for matching atomic number
template <class T>
T *makeAtomNumQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomNum2, descr);
}
//! returns a Query for matching atomic number
template <class T>
T *makeAtomNumQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNum, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomNumQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNumQuery(int what);

//! returns a Query for matching atomic number and aromaticity
template <class T>
T *makeAtomTypeQuery2(int num, int aromatic, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(makeAtomType(num, aromatic), queryAtomType2,
                                 descr);
}
//! returns a Query for matching atomic number and aromaticity
template <class T>
T *makeAtomTypeQuery(int num, int aromatic, const std::string &descr) {
  return makeAtomSimpleQuery<T>(makeAtomType(num, aromatic), queryAtomType,
                                descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomTypeQuery2(int num,
                                                             int aromatic);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomTypeQuery(int num,
                                                           int aromatic);

//! returns a Query for matching implicit valence
template <class T>
T *makeAtomImplicitValenceQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomImplicitValence2, descr);
}
//! returns a Query for matching implicit valence
template <class T>
T *makeAtomImplicitValenceQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomImplicitValence, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomImplicitValenceQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what);

//! returns a Query for matching explicit valence
template <class T>
T *makeAtomExplicitValenceQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomExplicitValence2, descr);
}
//! returns a Query for matching explicit valence
template <class T>
T *makeAtomExplicitValenceQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomExplicitValence, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomExplicitValenceQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomExplicitValenceQuery(int what);

//! returns a Query for matching total valence
template <class T>
T *makeAtomTotalValenceQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomTotalValence2, descr);
}
//! returns a Query for matching total valence
template <class T>
T *makeAtomTotalValenceQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomTotalValence, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomTotalValenceQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what);

//! returns a Query for matching explicit degree
template <class T>
T *makeAtomExplicitDegreeQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomExplicitDegree2, descr);
}
//! returns a Query for matching explicit degree
template <class T>
T *makeAtomExplicitDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomExplicitDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomExplicitDegreeQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what);

//! returns a Query for matching atomic degree
template <class T>
T *makeAtomTotalDegreeQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomTotalDegree2, descr);
}
//! returns a Query for matching atomic degree
template <class T>
T *makeAtomTotalDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomTotalDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomTotalDegreeQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what);

//! returns a Query for matching heavy atom degree
template <class T>
T *makeAtomHeavyAtomDegreeQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomHeavyAtomDegree2, descr);
}
//! returns a Query for matching heavy atom degree
template <class T>
T *makeAtomHeavyAtomDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomHeavyAtomDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomHeavyAtomDegreeQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHeavyAtomDegreeQuery(int what);

//! returns a Query for matching hydrogen count
template <class T>
T *makeAtomHCountQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomHCount2, descr);
}
//! returns a Query for matching hydrogen count
template <class T>
T *makeAtomHCountQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomHCount, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomHCountQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what);

//! returns a Query for matching ring atoms
template <class T>
T *makeAtomHasImplicitHQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryAtomHasImplicitH2, descr);
}
//! returns a Query for matching ring atoms
template <class T>
T *makeAtomHasImplicitHQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomHasImplicitH, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomHasImplicitHQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasImplicitHQuery();

//! returns a Query for matching implicit hydrogen count
template <class T>
T *makeAtomImplicitHCountQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomImplicitHCount2, descr);
}
//! returns a Query for matching implicit hydrogen count
template <class T>
T *makeAtomImplicitHCountQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomImplicitHCount, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomImplicitHCountQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomImplicitHCountQuery(int what);

//! returns a Query for matching the \c isAromatic flag
template <class T>
T *makeAtomAromaticQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryAtomAromatic2, descr);
}
//! returns a Query for matching the \c isAromatic flag
template <class T>
T *makeAtomAromaticQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomAromatic, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomAromaticQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomAromaticQuery();

//! returns a Query for matching aliphatic atoms
template <class T>
T *makeAtomAliphaticQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryAtomAliphatic2, descr);
}
//! returns a Query for matching aliphatic atoms
template <class T>
T *makeAtomAliphaticQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomAliphatic, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomAliphaticQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomAliphaticQuery();

//! returns a Query for matching atoms with a particular mass
template <class T>
T *makeAtomMassQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(massIntegerConversionFactor * what,
                                 queryAtomMass2, descr);
}
//! returns a Query for matching atoms with a particular mass
template <class T>
T *makeAtomMassQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(massIntegerConversionFactor * what,
                                queryAtomMass, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomMassQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomMassQuery(int what);

//! returns a Query for matching atoms with a particular isotope
template <class T>
T *makeAtomIsotopeQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomIsotope2, descr);
}
//! returns a Query for matching atoms with a particular isotope
template <class T>
T *makeAtomIsotopeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomIsotope, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomIsotopeQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what);

//! returns a Query for matching formal charge
template <class T>
T *makeAtomFormalChargeQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomFormalCharge2, descr);
}
//! returns a Query for matching formal charge
template <class T>
T *makeAtomFormalChargeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomFormalCharge, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomFormalChargeQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what);

//! returns a Query for matching negative formal charges (i.e. a query val of 1
//! matches a formal charge of -1)
template <class T>
T *makeAtomNegativeFormalChargeQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomNegativeFormalCharge2, descr);
}
//! returns a Query for matching negative formal charges (i.e. a query val of 1
//! matches a formal charge of -1)
template <class T>
T *makeAtomNegativeFormalChargeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNegativeFormalCharge, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomNegativeFormalChargeQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNegativeFormalChargeQuery(
    int what);

//! returns a Query for matching hybridization
template <class T>
T *makeAtomHybridizationQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomHybridization2, descr);
}
//! returns a Query for matching hybridization
template <class T>
T *makeAtomHybridizationQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomHybridization, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomHybridizationQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what);

//! returns a Query for matching the number of radical electrons
template <class T>
T *makeAtomNumRadicalElectronsQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomNumRadicalElectrons2, descr);
}
//! returns a Query for matching the number of radical electrons
template <class T>
T *makeAtomNumRadicalElectronsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNumRadicalElectrons, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomNumRadicalElectronsQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNumRadicalElectronsQuery(
    int what);

//! returns a Query for matching whether or not chirality has been set on the
//! atom
template <class T>
T *makeAtomHasChiralTagQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryAtomHasChiralTag2, descr);
}
//! returns a Query for matching whether or not chirality has been set on the
//! atom
template <class T>
T *makeAtomHasChiralTagQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomHasChiralTag, descr);
}
//! \overloadquery
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomHasChiralTagQuery2();
//! \overloadquery
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasChiralTagQuery();

//! returns a Query for matching whether or not a potentially chiral atom is
//! missing a chiral tag
template <class T>
T *makeAtomMissingChiralTagQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryAtomMissingChiralTag2, descr);
}
//! returns a Query for matching whether or not a potentially chiral atom is
//! missing a chiral tag
template <class T>
T *makeAtomMissingChiralTagQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomMissingChiralTag, descr);
}
//! \overloadquery
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomMissingChiralTagQuery2();
//! \overloadquery
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomMissingChiralTagQuery();

//! returns a Query for matching atoms with unsaturation:
template <class T>
T *makeAtomUnsaturatedQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryAtomUnsaturated2, descr);
}
//! returns a Query for matching atoms with unsaturation:
template <class T>
T *makeAtomUnsaturatedQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomUnsaturated, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomUnsaturatedQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery();

//! returns a Query for matching ring atoms
template <class T>
T *makeAtomInRingQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryIsAtomInRing2, descr);
}
//! returns a Query for matching ring atoms
template <class T>
T *makeAtomInRingQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryIsAtomInRing, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomInRingQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomInRingQuery();

//! returns a Query for matching atoms in a particular number of rings
template <class T>
T *makeAtomInNRingsQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryIsAtomInNRings2, descr);
}
//! returns a Query for matching atoms in a particular number of rings
template <class T>
T *makeAtomInNRingsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryIsAtomInNRings, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomInNRingsQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what);

//! returns a Query for matching atoms in rings of a particular size
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomInRingOfSizeQuery2(int tgt);
//! returns a Query for matching atoms in rings of a particular size
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt);

//! returns a Query for matching an atom's minimum ring size
template <class T>
T *makeAtomMinRingSizeQuery2(int tgt, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(tgt, queryAtomMinRingSize2, descr);
}
//! returns a Query for matching an atom's minimum ring size
template <class T>
T *makeAtomMinRingSizeQuery(int tgt, const std::string &descr) {
  return makeAtomSimpleQuery<T>(tgt, queryAtomMinRingSize, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomMinRingSizeQuery2(int tgt);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomMinRingSizeQuery(int tgt);

//! returns a Query for matching atoms with a particular number of ring bonds
template <class T>
T *makeAtomRingBondCountQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomRingBondCount2, descr);
}
//! returns a Query for matching atoms with a particular number of ring bonds
template <class T>
T *makeAtomRingBondCountQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomRingBondCount, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomRingBondCountQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what);

//! returns a Query for matching generic A atoms (heavy atoms)
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAAtomQuery2();
//! returns a Query for matching generic A atoms (heavy atoms)
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAAtomQuery();
//! returns a Query for matching generic AH atoms (any atom)
RDKIT_GRAPHMOL_EXPORT ATOM_NULL_QUERY2 *makeAHAtomQuery2();
//! returns a Query for matching generic AH atoms (any atom)
RDKIT_GRAPHMOL_EXPORT ATOM_NULL_QUERY *makeAHAtomQuery();
//! returns a Query for matching generic Q atoms (heteroatoms)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY2 *makeQAtomQuery2();
//! returns a Query for matching generic Q atoms (heteroatoms)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeQAtomQuery();
//! returns a Query for matching generic QH atoms (heteroatom or H)
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeQHAtomQuery2();
//! returns a Query for matching generic QH atoms (heteroatom or H)
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeQHAtomQuery();
//! returns a Query for matching generic X atoms (halogens)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY2 *makeXAtomQuery2();
//! returns a Query for matching generic X atoms (halogens)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeXAtomQuery();
//! returns a Query for matching generic XH atoms (halogen or H)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY2 *makeXHAtomQuery2();
//! returns a Query for matching generic XH atoms (halogen or H)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeXHAtomQuery();
//! returns a Query for matching generic M atoms (metals)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY2 *makeMAtomQuery2();
//! returns a Query for matching generic M atoms (metals)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeMAtomQuery();
//! returns a Query for matching generic MH atoms (metals or H)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY2 *makeMHAtomQuery2();
//! returns a Query for matching generic MH atoms (metals or H)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeMHAtomQuery();

// We support the same special atom queries that we can read from
// CXSMILES
const std::vector<std::string> complexQueries = {"A", "AH", "Q", "QH",
                                                 "X", "XH", "M", "MH"};
RDKIT_GRAPHMOL_EXPORT void convertComplexNameToQuery2(RDMolAtom query,
                                                      std::string_view symb);
RDKIT_GRAPHMOL_EXPORT void convertComplexNameToQuery(Atom *query,
                                                     std::string_view symb);

//! returns a Query for matching atoms that have ring bonds
template <class T>
T *makeAtomHasRingBondQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(1, queryAtomHasRingBond2, descr);
}
//! returns a Query for matching atoms that have ring bonds
template <class T>
T *makeAtomHasRingBondQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(1, queryAtomHasRingBond, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomHasRingBondQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasRingBondQuery();

//! returns a Query for matching the number of heteroatom neighbors
template <class T>
T *makeAtomNumHeteroatomNbrsQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomNumHeteroatomNbrs2, descr);
}
//! returns a Query for matching the number of heteroatom neighbors
template <class T>
T *makeAtomNumHeteroatomNbrsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNumHeteroatomNbrs, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomNumHeteroatomNbrsQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNumHeteroatomNbrsQuery(
    int what);

//! returns a Query for matching atoms that have heteroatom neighbors
template <class T>
T *makeAtomHasHeteroatomNbrsQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(1, queryAtomHasHeteroatomNbrs2, descr);
}
//! returns a Query for matching atoms that have heteroatom neighbors
template <class T>
T *makeAtomHasHeteroatomNbrsQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(1, queryAtomHasHeteroatomNbrs, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomHasHeteroatomNbrsQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasHeteroatomNbrsQuery();

//! returns a Query for matching the number of aliphatic heteroatom neighbors
template <class T>
T *makeAtomNumAliphaticHeteroatomNbrsQuery2(int what,
                                            const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomNumAliphaticHeteroatomNbrs2,
                                 descr);
}
//! returns a Query for matching the number of aliphatic heteroatom neighbors
template <class T>
T *makeAtomNumAliphaticHeteroatomNbrsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNumAliphaticHeteroatomNbrs,
                                descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *
makeAtomNumAliphaticHeteroatomNbrsQuery2(int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *
makeAtomNumAliphaticHeteroatomNbrsQuery(int what);

//! returns a Query for matching atoms that have heteroatom neighbors
template <class T>
T *makeAtomHasAliphaticHeteroatomNbrsQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(1, queryAtomHasAliphaticHeteroatomNbrs2,
                                 descr);
}
//! returns a Query for matching atoms that have heteroatom neighbors
template <class T>
T *makeAtomHasAliphaticHeteroatomNbrsQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(1, queryAtomHasAliphaticHeteroatomNbrs, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *
makeAtomHasAliphaticHeteroatomNbrsQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *
makeAtomHasAliphaticHeteroatomNbrsQuery();

//! returns a Query for matching the number of non-hydrogen neighbors
template <class T>
T *makeAtomNonHydrogenDegreeQuery2(int what, const std::string &descr) {
  return makeAtomSimpleQuery2<T>(what, queryAtomNonHydrogenDegree2, descr);
}
//! returns a Query for matching the number of non-hydrogen neighbors
template <class T>
T *makeAtomNonHydrogenDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNonHydrogenDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomNonHydrogenDegreeQuery2(
    int what);
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNonHydrogenDegreeQuery(
    int what);

//! returns a Query for matching bridgehead atoms
template <class T>
T *makeAtomIsBridgeheadQuery2(const std::string &descr) {
  return makeAtomSimpleQuery2<T>(true, queryIsAtomBridgehead2, descr);
}
//! returns a Query for matching bridgehead atoms
template <class T>
T *makeAtomIsBridgeheadQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryIsAtomBridgehead, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY2 *makeAtomIsBridgeheadQuery2();
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomIsBridgeheadQuery();

//! returns a Query for matching bond orders
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeBondOrderEqualsQuery2(
    BondEnums::BondType what);
//! returns a Query for matching bond orders
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondOrderEqualsQuery(
    Bond::BondType what);
//! returns a Query for unspecified SMARTS bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeSingleOrAromaticBondQuery2();
//! returns a Query for unspecified SMARTS bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrAromaticBondQuery();
//! returns a Query for double|aromatic bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeDoubleOrAromaticBondQuery2();
//! returns a Query for single|double bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeDoubleOrAromaticBondQuery();
//! returns a Query for single|double bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeSingleOrDoubleBondQuery2();
//! returns a Query for tautomeric bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrDoubleBondQuery();
//! returns a Query for tautomeric bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *
makeSingleOrDoubleOrAromaticBondQuery2();
//! returns a Query for tautomeric bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *
makeSingleOrDoubleOrAromaticBondQuery();

//! returns a Query for matching bond directions
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeBondDirEqualsQuery2(
    BondEnums::BondDir what);
//! returns a Query for matching bond directions
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondDirEqualsQuery(
    Bond::BondDir what);
//! returns a Query for matching bonds with stereo set
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeBondHasStereoQuery2();
//! returns a Query for matching bonds with stereo set
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondHasStereoQuery();
//! returns a Query for matching ring bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeBondIsInRingQuery2();
//! returns a Query for matching ring bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondIsInRingQuery();
//! returns a Query for matching bonds in rings of a particular size
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeBondInRingOfSizeQuery2(int what);
//! returns a Query for matching bonds in rings of a particular size
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int what);
//! returns a Query for matching a bond's minimum ring size
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeBondMinRingSizeQuery2(int what);
//! returns a Query for matching a bond's minimum ring size
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondMinRingSizeQuery(int what);
//! returns a Query for matching bonds in a particular number of rings
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY2 *makeBondInNRingsQuery2(int tgt);
//! returns a Query for matching bonds in a particular number of rings
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondInNRingsQuery(int tgt);

//! returns a Query for matching any bond
RDKIT_GRAPHMOL_EXPORT BOND_NULL_QUERY2 *makeBondNullQuery2();
//! returns a Query for matching any bond
RDKIT_GRAPHMOL_EXPORT BOND_NULL_QUERY *makeBondNullQuery();
//! returns a Query for matching any atom
RDKIT_GRAPHMOL_EXPORT ATOM_NULL_QUERY2 *makeAtomNullQuery2();
//! returns a Query for matching any atom
RDKIT_GRAPHMOL_EXPORT ATOM_NULL_QUERY *makeAtomNullQuery();

static inline int queryAtomRingMembership2(ConstRDMolAtom at) {
  return static_cast<int>(at.mol().getRingInfo().numAtomRings(at.index()));
}
static inline int queryAtomRingMembership(Atom const *at) {
  return static_cast<int>(
      at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx()));
}
// I'm pretty sure that this typedef shouldn't be necessary,
// but VC++ generates a warning about const Atom const * in
// the definition of Match, then complains about an override
// that differs only by const/volatile (c4301), then generates
// incorrect code if we don't do this... so let's do it.
typedef Atom const *ConstAtomPtr;
class AtomRingQuery;
class AtomRingQuery2;

template <bool newVersion>
class AtomRingQueryBase
    : public Queries::EqualityQuery<
          int, std::conditional_t<newVersion, ConstRDMolAtom, ConstAtomPtr>,
          true> {
  using DataFuncArgType =
      std::conditional_t<newVersion, ConstRDMolAtom, ConstAtomPtr>;

 public:
  AtomRingQueryBase() : Queries::EqualityQuery<int, DataFuncArgType, true>(-1) {
    // default is to just do a number of rings query:
    this->setDescription("AtomInNRings");
    if constexpr (newVersion) {
      this->setDataFunc(queryAtomRingMembership2);
    } else {
      this->setDataFunc(queryAtomRingMembership);
    }
  }
  explicit AtomRingQueryBase(int v)
      : Queries::EqualityQuery<int, DataFuncArgType, true>(v) {
    // default is to just do a number of rings query:
    this->setDescription("AtomInNRings");
    if constexpr (newVersion) {
      this->setDataFunc(queryAtomRingMembership2);
    } else {
      this->setDataFunc(queryAtomRingMembership);
    }
  }

  bool Match(const DataFuncArgType what) const override {
    int v = this->TypeConvert(what, Queries::Int2Type<true>());
    bool res;
    if (this->d_val < 0) {
      res = v != 0;
    } else {
      res = !Queries::queryCmp(v, this->d_val, this->d_tol);
    }
    if (this->getNegation()) {
      res = !res;
    }
    return res;
  }

  //! returns a copy of this query
  Queries::Query<int, DataFuncArgType, true> *copy() const override {
    using RingQueryType =
        std::conditional_t<newVersion, AtomRingQuery2, AtomRingQuery>;
    RingQueryType *res = new RingQueryType(this->d_val);
    res->setNegation(this->getNegation());
    res->setTol(this->getTol());
    res->d_description = this->d_description;
    res->d_dataFunc = this->d_dataFunc;
    return res;
  }
};
class RDKIT_GRAPHMOL_EXPORT AtomRingQuery2 : public AtomRingQueryBase<true> {
 public:
  using AtomRingQueryBase::AtomRingQueryBase;
};
class RDKIT_GRAPHMOL_EXPORT AtomRingQuery : public AtomRingQueryBase<false> {
 public:
  using AtomRingQueryBase::AtomRingQueryBase;
};

//! allows use of recursive structure queries (e.g. recursive SMARTS)
class RDKIT_GRAPHMOL_EXPORT RecursiveStructureQuery2
    : public Queries::SetQuery<int, ConstRDMolAtom, true> {
 public:
  RecursiveStructureQuery2() : Queries::SetQuery<int, ConstRDMolAtom, true>() {
    setDataFunc(getAtIdx);
    setDescription("RecursiveStructure");
  }
  //! initialize from an RDMol pointer
  /*!
    <b>Notes</b>
      - this takes over ownership of the pointer
  */
  RecursiveStructureQuery2(RDMol *query, unsigned int serialNumber = 0)
      : Queries::SetQuery<int, ConstRDMolAtom, true>(),
        d_serialNumber(serialNumber) {
    setQueryMol(query);
    setDataFunc(getAtIdx);
    setDescription("RecursiveStructure");
  }
  //! returns the index of an atom
  static inline int getAtIdx(ConstRDMolAtom at) {
    PRECONDITION(at.index() < at.mol().getNumAtoms(), "bad atom argument");
    return at.index();
  }

  //! sets the molecule we'll use recursively
  /*!
    <b>Notes</b>
      - this takes over ownership of the pointer
  */
  void setQueryMol(RDMol *query) { dp_queryMol.reset(query); }
  //! returns a pointer to our query molecule
  RDMol const *getQueryMol() const { return dp_queryMol.get(); }

  //! returns a copy of this query
  Queries::Query<int, ConstRDMolAtom, true> *copy() const override {
    RecursiveStructureQuery2 *res = new RecursiveStructureQuery2();
    res->dp_queryMol.reset(new RDMol(*dp_queryMol, true));

    std::set<int>::const_iterator i;
    for (i = d_set.begin(); i != d_set.end(); i++) {
      res->insert(*i);
    }
    res->setNegation(getNegation());
    res->d_description = d_description;
    res->d_serialNumber = d_serialNumber;
    return res;
  }
  unsigned int getSerialNumber() const { return d_serialNumber; }

#ifdef RDK_BUILD_THREADSAFE_SSS
  std::mutex d_mutex;
#endif
 private:
  std::unique_ptr<const RDMol> dp_queryMol;
  unsigned int d_serialNumber{0};
};

//! allows use of recursive structure queries (e.g. recursive SMARTS)
class RDKIT_GRAPHMOL_EXPORT RecursiveStructureQuery
    : public Queries::SetQuery<int, Atom const *, true> {
 public:
  RecursiveStructureQuery() : Queries::SetQuery<int, Atom const *, true>() {
    setDataFunc(getAtIdx);
    setDescription("RecursiveStructure");
  }
  //! initialize from an ROMol pointer
  /*!
    <b>Notes</b>
      - this takes over ownership of the pointer
  */
  RecursiveStructureQuery(ROMol const *query, unsigned int serialNumber = 0)
      : Queries::SetQuery<int, Atom const *, true>(),
        d_serialNumber(serialNumber) {
    setQueryMol(query);
    setDataFunc(getAtIdx);
    setDescription("RecursiveStructure");
  }
  //! returns the index of an atom
  static inline int getAtIdx(Atom const *at) {
    PRECONDITION(at, "bad atom argument");
    return at->getIdx();
  }

  //! sets the molecule we'll use recursively
  /*!
    <b>Notes</b>
      - this takes over ownership of the pointer
  */
  void setQueryMol(ROMol const *query) { dp_queryMol.reset(query); }
  //! returns a pointer to our query molecule
  ROMol const *getQueryMol() const { return dp_queryMol.get(); }

  //! returns a copy of this query
  Queries::Query<int, Atom const *, true> *copy() const override {
    RecursiveStructureQuery *res = new RecursiveStructureQuery();
    res->dp_queryMol.reset(new ROMol(*dp_queryMol, true));

    std::set<int>::const_iterator i;
    for (i = d_set.begin(); i != d_set.end(); i++) {
      res->insert(*i);
    }
    res->setNegation(getNegation());
    res->d_description = d_description;
    res->d_serialNumber = d_serialNumber;
    return res;
  }
  unsigned int getSerialNumber() const { return d_serialNumber; }

#ifdef RDK_BUILD_THREADSAFE_SSS
  std::mutex d_mutex;
#endif
 private:
  boost::shared_ptr<const ROMol> dp_queryMol;
  unsigned int d_serialNumber{0};
};

template <typename T>
int nullDataFun(T) {
  return 1;
}
template <typename T>
bool nullQueryFun(T) {
  return true;
}

typedef Bond const *ConstBondPtr;

// ! Query whether an atom has a property
template <class TargetPtr>
class HasPropQuery : public Queries::EqualityQuery<int, TargetPtr, true> {
  PropToken propname;

 public:
  HasPropQuery() : Queries::EqualityQuery<int, TargetPtr, true>(), propname() {
    this->setDescription("HasProp");
    this->setDataFunc(nullptr);
  }
  explicit HasPropQuery(std::string v)
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(std::move(v)) {
    this->setDescription("HasProp");
    this->setDataFunc(nullptr);
  }
  explicit HasPropQuery(PropToken v)
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(std::move(v)) {
    this->setDescription("HasProp");
    this->setDataFunc(nullptr);
  }

  bool Match(const TargetPtr what) const override {
    bool res;
    if constexpr (std::is_same_v<TargetPtr, ConstRDMolAtom>) {
      res = what.mol().hasAtomProp(propname, what.index());
    } else if constexpr (std::is_same_v<TargetPtr, ConstRDMolBond>) {
      res = what.mol().hasBondProp(propname, what.index());
    } else {
      res = what->hasProp(propname.getString());
    }
    if (this->getNegation()) {
      res = !res;
    }
    return res;
  }

  //! returns a copy of this query
  Queries::Query<int, TargetPtr, true> *copy() const override {
    HasPropQuery *res = new HasPropQuery(this->propname);
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  }

  const std::string &getPropName() const { return propname.getString(); }
};

typedef Queries::EqualityQuery<int, ConstRDMolAtom, true> ATOM_PROP_QUERY2;
typedef Queries::EqualityQuery<int, ConstRDMolBond, true> BOND_PROP_QUERY2;
typedef Queries::EqualityQuery<int, Atom const *, true> ATOM_PROP_QUERY;
typedef Queries::EqualityQuery<int, Bond const *, true> BOND_PROP_QUERY;

//! returns a Query for matching atoms or bonds that have a particular property
template <class Target>
Queries::EqualityQuery<int, Target, true> *makeHasPropQuery2(
    const std::string &property) {
  return new HasPropQuery<Target>(property);
}
//! returns a Query for matching atoms or bonds that have a particular property
template <class Target>
Queries::EqualityQuery<int, const Target *, true> *makeHasPropQuery(
    const std::string &property) {
  return new HasPropQuery<const Target *>(property);
}

// ! Query whether an atom has a property with a value
class HasPropWithValueQueryBase {
 public:
  HasPropWithValueQueryBase() = default;
  virtual ~HasPropWithValueQueryBase() = default;
  virtual PairHolder getPair() const = 0;
  virtual double getTolerance() const = 0;
};

template <class TargetPtr, class T>
class HasPropWithValueQuery
    : public HasPropWithValueQueryBase,
      public Queries::EqualityQuery<int, TargetPtr, true> {
  PropToken propname;
  T val;
  double tolerance{0.0};

 public:
  HasPropWithValueQuery()
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(), val() {
    // default is to just do a number of rings query:
    this->setDescription("HasPropWithValue");
    this->setDataFunc(0);
  }

  PairHolder getPair() const override {
    return PairHolder(Dict::Pair(propname.getString(), val));
  }

  double getTolerance() const override { return tolerance; }

  explicit HasPropWithValueQuery(std::string prop, const T &v,
                                 const T &tol = 0.0)
      : Queries::EqualityQuery<int, TargetPtr, true>(),
        propname(std::move(prop)),
        val(v),
        tolerance(tol) {
    // default is to just do a number of rings query:
    this->setDescription("HasPropWithValue");
    this->setDataFunc(nullptr);
  }

  bool Match(const TargetPtr what) const override {
    bool res;
    if constexpr (std::is_same_v<TargetPtr, ConstRDMolAtom>) {
      res = what.mol().hasAtomProp(propname, what.index());
    } else if constexpr (std::is_same_v<TargetPtr, ConstRDMolBond>) {
      res = what.mol().hasBondProp(propname, what.index());
    } else {
      res = what->hasProp(propname.getString());
    }
    if (res) {
      try {
        T atom_val;
        if constexpr (std::is_same_v<TargetPtr, ConstRDMolAtom>) {
          atom_val =
              what.mol().template getAtomPropValue<T>(propname, what.index());
        } else if constexpr (std::is_same_v<TargetPtr, ConstRDMolBond>) {
          atom_val =
              what.mol().template getBondPropValue<T>(propname, what.index());
        } else {
          atom_val = what->template getProp<T>(propname.getString());
        }
        res = Queries::queryCmp(atom_val, this->val,
                                static_cast<T>(this->tolerance)) == 0;
      } catch (KeyErrorException &) {
        res = false;
      } catch (std::bad_any_cast &) {
        res = false;
      }
#ifdef __GNUC__
#if (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 2))
      catch (...) {
        // catch all -- this is currently necessary to
        //  trap some bugs in boost+gcc configurations
        //  Normally, this is not the correct thing to
        //  do, but the only exception above is due
        //  to the boost any_cast which is trapped
        //  by the Boost python wrapper when it shouldn't
        //  be.
        res = false;
      }
#endif
#endif
    }
    if (this->getNegation()) {
      res = !res;
    }
    return res;
  }

  //! returns a copy of this query
  Queries::Query<int, TargetPtr, true> *copy() const override {
    HasPropWithValueQuery *res = new HasPropWithValueQuery(
        this->propname.getString(), this->val, this->tolerance);
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  }
};

template <class TargetPtr>
class HasPropWithValueQuery<TargetPtr, std::string>
    : public HasPropWithValueQueryBase,
      public Queries::EqualityQuery<int, TargetPtr, true> {
  PropToken propname;
  PropToken val;

 public:
  HasPropWithValueQuery()
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(), val() {
    // default is to just do a number of rings query:
    this->setDescription("HasPropWithValue");
    this->setDataFunc(0);
  }
  explicit HasPropWithValueQuery(std::string prop, std::string v,
                                 const double /*tol*/ = 0.0)
      : Queries::EqualityQuery<int, TargetPtr, true>(),
        propname(std::move(prop)),
        val(std::move(v)) {
    // default is to just do a number of rings query:
    this->setDescription("HasPropWithValue");
    this->setDataFunc(nullptr);
  }

  PairHolder getPair() const override {
    return PairHolder(Dict::Pair(propname.getString(), val.getString()));
  }

  double getTolerance() const override { return 0.0; }

  bool Match(const TargetPtr what) const override {
    bool res;
    if constexpr (std::is_same_v<TargetPtr, ConstRDMolAtom>) {
      res = what.mol().hasAtomProp(propname, what.index());
    } else if constexpr (std::is_same_v<TargetPtr, ConstRDMolBond>) {
      res = what.mol().hasBondProp(propname, what.index());
    } else {
      res = what->hasProp(propname.getString());
    }
    if (res) {
      try {
        if constexpr (std::is_same_v<TargetPtr, ConstRDMolAtom>) {
          PropToken token = what.mol().template getAtomProp<PropToken>(
              propname, what.index());
          res = token == this->val;
        } else if constexpr (std::is_same_v<TargetPtr, ConstRDMolBond>) {
          PropToken token = what.mol().template getBondProp<PropToken>(
              propname, what.index());
          res = token == this->val;
        } else {
          std::string atom_val =
              what->template getProp<std::string>(propname.getString());
          res = atom_val == this->val.getString();
        }
      } catch (KeyErrorException &) {
        res = false;
      } catch (std::bad_any_cast &) {
        res = false;
      }
#ifdef __GNUC__
#if (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 2))
      catch (...) {
        // catch all -- this is currently necessary to
        //  trap some bugs in boost+gcc configurations
        //  Normally, this is not the correct thing to
        //  do, but the only exception above is due
        //  to the boost any_cast which is trapped
        //  by the Boost python wrapper when it shouldn't
        //  be.
        res = false;
      }
#endif
#endif
    }
    if (this->getNegation()) {
      res = !res;
    }
    return res;
  }

  //! returns a copy of this query
  Queries::Query<int, TargetPtr, true> *copy() const override {
    HasPropWithValueQuery<TargetPtr, std::string> *res =
        new HasPropWithValueQuery<TargetPtr, std::string>(
            this->propname.getString(), this->val.getString());
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  }
};

template <class TargetPtr>
class HasPropWithValueQuery<TargetPtr, ExplicitBitVect>
    : public HasPropWithValueQueryBase,
      public Queries::EqualityQuery<int, TargetPtr, true> {
  PropToken propname;
  ExplicitBitVect val;
  double tol{0.0};

 public:
  HasPropWithValueQuery()
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(), val() {
    this->setDescription("HasPropWithValue");
    this->setDataFunc(0);
  }

  explicit HasPropWithValueQuery(std::string prop, const ExplicitBitVect &v,
                                 double tol = 0.0)
      : Queries::EqualityQuery<int, TargetPtr, true>(),
        propname(std::move(prop)),
        val(v),
        tol(tol) {
    this->setDescription("HasPropWithValue");
    this->setDataFunc(nullptr);
  }

  PairHolder getPair() const override {
    return PairHolder(Dict::Pair(propname.getString(), val));
  }

  double getTolerance() const override { return tol; }

  bool Match(const TargetPtr what) const override {
    bool res;
    if constexpr (std::is_same_v<TargetPtr, ConstRDMolAtom>) {
      res = what.mol().hasAtomProp(propname, what.index());
    } else if constexpr (std::is_same_v<TargetPtr, ConstRDMolBond>) {
      res = what.mol().hasBondProp(propname, what.index());
    } else {
      res = what->hasProp(propname.getString());
    }
    if (res) {
      try {
        const ExplicitBitVect &bv =
            what->template getProp<const ExplicitBitVect &>(
                propname.getString());
        const double tani = TanimotoSimilarity(val, bv);
        res = (1.0 - tani) <= tol;
      } catch (KeyErrorException &) {
        res = false;
      } catch (std::bad_any_cast &) {
        res = false;
      }
#ifdef __GNUC__
#if (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 2))
      catch (...) {
        // catch all -- this is currently necessary to
        //  trap some bugs in boost+gcc configurations
        //  Normally, this is not the correct thing to
        //  do, but the only exception above is due
        //  to the boost any_cast which is trapped
        //  by the Boost python wrapper when it shouldn't
        //  be.
        res = false;
      }
#endif
#endif
    }
    if (this->getNegation()) {
      res = !res;
    }
    return res;
  }

  //! returns a copy of this query
  Queries::Query<int, TargetPtr, true> *copy() const override {
    HasPropWithValueQuery<TargetPtr, ExplicitBitVect> *res =
        new HasPropWithValueQuery<TargetPtr, ExplicitBitVect>(
            this->propname.getString(), this->val, this->tol);
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  }
};

template <class Target, class T>
Queries::EqualityQuery<int, const Target *, true> *makePropQuery(
    const std::string &propname, const T &val, double tolerance = 0.0) {
  return new HasPropWithValueQuery<const Target *, T>(propname, val, tolerance);
}

template <class Target>
Queries::EqualityQuery<int, const Target *, true> *makePropQuery(
    const std::string &propname, const ExplicitBitVect &val,
    double tolerance = 0.0) {
  return new HasPropWithValueQuery<const Target *, ExplicitBitVect>(
      propname, val, tolerance);
}

RDKIT_GRAPHMOL_EXPORT bool isComplexQuery(ConstRDMolBond b);
RDKIT_GRAPHMOL_EXPORT bool isComplexQuery(const Bond *b);
RDKIT_GRAPHMOL_EXPORT bool isComplexQuery(ConstRDMolAtom a);
RDKIT_GRAPHMOL_EXPORT bool isComplexQuery(const Atom *a);
RDKIT_GRAPHMOL_EXPORT bool isAtomAromatic(ConstRDMolAtom a);
RDKIT_GRAPHMOL_EXPORT bool isAtomAromatic(const Atom *a);
RDKIT_GRAPHMOL_EXPORT bool isAtomListQuery(ConstRDMolAtom a);
RDKIT_GRAPHMOL_EXPORT bool isAtomListQuery(const Atom *a);
RDKIT_GRAPHMOL_EXPORT void getAtomListQueryVals(const RDMol::QUERYATOM_QUERY *q,
                                                std::vector<int> &vals);
RDKIT_GRAPHMOL_EXPORT void getAtomListQueryVals(const Atom::QUERYATOM_QUERY *q,
                                                std::vector<int> &vals);

//! Checks if an atom is dummy or not.
//! 1. A dummy non-query atom (e.g., "*" in SMILES) is defined by its zero
//! atomic
//!    number. This rule breaks for query atoms because a COMPOSITE_OR query
//!    atom also has a zero atomic number (#6349).
//! 2. A dummy query atom (e.g., "*" in SMARTS) is defined by its explicit
//!    description: "AtomNull".
inline bool isAtomDummy(ConstRDMolAtom a) {
  return (!a.mol().hasAtomQuery(a.index()) && a.data().getAtomicNum() == 0) ||
         (a.mol().hasAtomQuery(a.index()) &&
          !a.mol().getAtomQuery(a.index())->getNegation() &&
          a.mol().getAtomQuery(a.index())->getDescription() == "AtomNull");
}
//! \overload
inline bool isAtomDummy(const Atom *a) {
  return (!a->hasQuery() && a->getAtomicNum() == 0) ||
         (a->hasQuery() && !a->getQuery()->getNegation() &&
          a->getQuery()->getDescription() == "AtomNull");
}

namespace QueryOps {
RDKIT_GRAPHMOL_EXPORT void completeMolQueries(
    RDMol *mol, unsigned int magicVal = 0xDEADBEEF);
RDKIT_GRAPHMOL_EXPORT void completeMolQueries(
    RWMol *mol, unsigned int magicVal = 0xDEADBEEF);

template <typename DataFuncArgType>
void expandQuery(
    std::unique_ptr<Queries::Query<int, DataFuncArgType, true>> &firstQuery,
    std::unique_ptr<Queries::Query<int, DataFuncArgType, true>> newQuery,
    Queries::CompositeQueryType how = Queries::COMPOSITE_AND,
    bool maintainOrder = true) {
  constexpr bool isAtom = std::is_same_v<DataFuncArgType, const Atom *> ||
                          std::is_same_v<DataFuncArgType, ConstRDMolAtom>;
  constexpr const char *nullString = isAtom ? "AtomNull" : "BondNull";
  bool thisIsNullQuery = firstQuery->getDescription() == nullString;
  bool otherIsNullQuery = newQuery->getDescription() == nullString;

  if (thisIsNullQuery || otherIsNullQuery) {
    auto *query1 = firstQuery.release();
    auto *query2 = newQuery.release();
    mergeNullQueries(query1, thisIsNullQuery, query2, otherIsNullQuery, how);
    delete query2;
    firstQuery.reset(query1);
    return;
  }

  auto *origQ = firstQuery.release();
  std::string descrip;
  switch (how) {
    case Queries::COMPOSITE_AND:
      firstQuery.reset(new Queries::AndQuery<int, DataFuncArgType, true>);
      descrip = isAtom ? "AtomAnd" : "BondAnd";
      break;
    case Queries::COMPOSITE_OR:
      firstQuery.reset(new Queries::OrQuery<int, DataFuncArgType, true>);
      descrip = isAtom ? "AtomOr" : "BondOr";
      break;
    case Queries::COMPOSITE_XOR:
      firstQuery.reset(new Queries::XOrQuery<int, DataFuncArgType, true>);
      descrip = isAtom ? "AtomXor" : "BondXor";
      break;
    default:
      UNDER_CONSTRUCTION("unrecognized combination query");
  }
  firstQuery->setDescription(descrip);
  using CHILD_TYPE =
      typename Queries::Query<int, DataFuncArgType, true>::CHILD_TYPE;
  if (maintainOrder) {
    firstQuery->addChild(CHILD_TYPE(origQ));
    firstQuery->addChild(CHILD_TYPE(newQuery.release()));
  } else {
    firstQuery->addChild(CHILD_TYPE(newQuery.release()));
    firstQuery->addChild(CHILD_TYPE(origQ));
  }
}

//! Adds the default Query to the given atom in mol.
//! If the atom already has a query, nothing will be changed.
RDKIT_GRAPHMOL_EXPORT void replaceAtomWithQueryAtom(RDMol *mol,
                                                    atomindex_t atomIndex);
//! Replaces the given atom in the molecule with a QueryAtom that is otherwise
//! a copy of the given atom.  Returns a pointer to that atom.
//! if the atom already has a query, nothing will be changed
RDKIT_GRAPHMOL_EXPORT Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom);

RDKIT_GRAPHMOL_EXPORT void finalizeQueryFromDescription(
    Queries::Query<int, ConstRDMolAtom, true> *query, ConstRDMolAtom owner);
RDKIT_GRAPHMOL_EXPORT void finalizeQueryFromDescription(
    Queries::Query<int, Atom const *, true> *query, Atom const *owner);
RDKIT_GRAPHMOL_EXPORT void finalizeQueryFromDescription(
    Queries::Query<int, ConstRDMolBond, true> *query, ConstRDMolBond owner);
RDKIT_GRAPHMOL_EXPORT void finalizeQueryFromDescription(
    Queries::Query<int, Bond const *, true> *query, Bond const *owner);

RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    const Queries::Query<int, Bond const *, true> &qry);
RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    const Queries::Query<int, ConstRDMolBond, true> &qry);
inline bool hasBondTypeQuery(const Bond &bond) {
  if (!bond.hasQuery()) {
    return false;
  }
  return hasBondTypeQuery(*bond.getQuery());
}
inline bool hasBondTypeQuery(ConstRDMolBond bond) {
  auto *query = bond.mol().getBondQuery(bond.index());
  if (query == nullptr) {
    return false;
  }
  return hasBondTypeQuery(*query);
}
RDKIT_GRAPHMOL_EXPORT bool hasComplexBondTypeQuery(
    const Queries::Query<int, Bond const *, true> &qry);
RDKIT_GRAPHMOL_EXPORT bool hasComplexBondTypeQuery(
    const Queries::Query<int, ConstRDMolBond, true> &qry);
inline bool hasComplexBondTypeQuery(const Bond &bond) {
  if (!bond.hasQuery()) {
    return false;
  }
  return hasComplexBondTypeQuery(*bond.getQuery());
}
inline bool hasComplexBondTypeQuery(ConstRDMolBond bond) {
  auto *query = bond.mol().getBondQuery(bond.index());
  if (query == nullptr) {
    return false;
  }
  return hasComplexBondTypeQuery(*query);
}

RDKIT_GRAPHMOL_EXPORT bool isMetal(const Atom &atom);
}  // namespace QueryOps
}  // namespace RDKit
#endif
