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
#ifndef RD_QUERY_OPS_H
#define RD_QUERY_OPS_H

#include <GraphMol/RDKitBase.h>
#include <Query/QueryObjects.h>
#include <Query/Query.h>
#include <DataStructs/BitVects.h>
#include <DataStructs/BitOps.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#include <utility>
#endif

namespace RDKit {
typedef Queries::Query<bool, Atom const *, true> ATOM_BOOL_QUERY;
typedef Queries::Query<bool, Bond const *, true> BOND_BOOL_QUERY;

typedef Queries::AndQuery<int, Atom const *, true> ATOM_AND_QUERY;
typedef Queries::AndQuery<int, Bond const *, true> BOND_AND_QUERY;

typedef Queries::OrQuery<int, Atom const *, true> ATOM_OR_QUERY;
typedef Queries::OrQuery<int, Bond const *, true> BOND_OR_QUERY;

typedef Queries::XOrQuery<int, Atom const *, true> ATOM_XOR_QUERY;
typedef Queries::XOrQuery<int, Bond const *, true> BOND_XOR_QUERY;

typedef Queries::EqualityQuery<int, Atom const *, true> ATOM_EQUALS_QUERY;
typedef Queries::EqualityQuery<int, Bond const *, true> BOND_EQUALS_QUERY;

typedef Queries::GreaterQuery<int, Atom const *, true> ATOM_GREATER_QUERY;
typedef Queries::GreaterQuery<int, Bond const *, true> BOND_GREATER_QUERY;

typedef Queries::GreaterEqualQuery<int, Atom const *, true>
    ATOM_GREATEREQUAL_QUERY;
typedef Queries::GreaterEqualQuery<int, Bond const *, true>
    BOND_GREATEREQUAL_QUERY;

typedef Queries::LessQuery<int, Atom const *, true> ATOM_LESS_QUERY;
typedef Queries::LessQuery<int, Bond const *, true> BOND_LESS_QUERY;

typedef Queries::LessEqualQuery<int, Atom const *, true> ATOM_LESSEQUAL_QUERY;
typedef Queries::LessEqualQuery<int, Bond const *, true> BOND_LESSEQUAL_QUERY;

typedef Queries::RangeQuery<int, Atom const *, true> ATOM_RANGE_QUERY;
typedef Queries::RangeQuery<int, Bond const *, true> BOND_RANGE_QUERY;

typedef Queries::SetQuery<int, Atom const *, true> ATOM_SET_QUERY;
typedef Queries::SetQuery<int, Bond const *, true> BOND_SET_QUERY;

typedef Queries::Query<int, Bond const *, true> BOND_NULL_QUERY;
typedef Queries::Query<int, Atom const *, true> ATOM_NULL_QUERY;

// -------------------------------------------------
// common atom queries

static inline int queryAtomAromatic(Atom const *at) {
  return at->getIsAromatic();
};
static inline int queryAtomAliphatic(Atom const *at) {
  return !(at->getIsAromatic());
};
static inline int queryAtomExplicitDegree(Atom const *at) {
  return at->getDegree();
};
static inline int queryAtomTotalDegree(Atom const *at) {
  return at->getTotalDegree();
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
static inline int queryAtomHCount(Atom const *at) {
  return at->getTotalNumHs(true);
};
static inline int queryAtomImplicitHCount(Atom const *at) {
  return at->getTotalNumHs(false);
};
static inline int queryAtomHasImplicitH(Atom const *at) {
  return int(at->getTotalNumHs(false) > 0);
};
static inline int queryAtomImplicitValence(Atom const *at) {
  return at->getImplicitValence();
};
static inline int queryAtomExplicitValence(Atom const *at) {
  return at->getExplicitValence() - at->getNumExplicitHs();
};
static inline int queryAtomTotalValence(Atom const *at) {
  return at->getTotalValence();
};
static inline int queryAtomUnsaturated(Atom const *at) {
  return at->getTotalDegree() < at->getTotalValence();
};
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

static inline int queryAtomType(Atom const *at) {
  return makeAtomType(at->getAtomicNum(), at->getIsAromatic());
};
const int massIntegerConversionFactor = 1000;
static inline int queryAtomMass(Atom const *at) {
  return static_cast<int>(
      std::round(massIntegerConversionFactor * at->getMass()));
};
static inline int queryAtomIsotope(Atom const *at) {
  return static_cast<int>(at->getIsotope());
};
static inline int queryAtomFormalCharge(Atom const *at) {
  return static_cast<int>(at->getFormalCharge());
};
static inline int queryAtomNegativeFormalCharge(Atom const *at) {
  return static_cast<int>(-1 * at->getFormalCharge());
};
static inline int queryAtomHybridization(Atom const *at) {
  return at->getHybridization();
};
static inline int queryAtomNumRadicalElectrons(Atom const *at) {
  return at->getNumRadicalElectrons();
};
static inline int queryAtomHasChiralTag(Atom const *at) {
  return at->getChiralTag() != Atom::CHI_UNSPECIFIED;
};
static inline int queryAtomMissingChiralTag(Atom const *at) {
  return at->getChiralTag() == Atom::CHI_UNSPECIFIED &&
         at->hasProp(common_properties::_ChiralityPossible);
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

RDKIT_GRAPHMOL_EXPORT unsigned int queryAtomBondProduct(Atom const *at);
RDKIT_GRAPHMOL_EXPORT unsigned int queryAtomAllBondProduct(Atom const *at);

// -------------------------------------------------
// common bond queries

static inline int queryBondOrder(Bond const *bond) {
  return static_cast<int>(bond->getBondType());
};
static inline int queryBondIsSingleOrAromatic(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
                          bond->getBondType() == Bond::AROMATIC);
};
static inline int queryBondIsDoubleOrAromatic(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::DOUBLE ||
                          bond->getBondType() == Bond::AROMATIC);
};
static inline int queryBondIsSingleOrDouble(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
                          bond->getBondType() == Bond::DOUBLE);
};
static inline int queryBondIsSingleOrDoubleOrAromatic(Bond const *bond) {
  return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
                          bond->getBondType() == Bond::DOUBLE ||
                          bond->getBondType() == Bond::AROMATIC);
};
static inline int queryBondDir(Bond const *bond) {
  return static_cast<int>(bond->getBondDir());
};
static inline int queryIsBondInNRings(Bond const *at) {
  return at->getOwningMol().getRingInfo()->numBondRings(at->getIdx());
};
static inline int queryBondHasStereo(Bond const *bnd) {
  return bnd->getStereo() > Bond::STEREONONE;
};

// -------------------------------------------------
// ring queries

static inline int queryIsAtomInNRings(Atom const *at) {
  return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx());
};
static inline int queryIsAtomInRing(Atom const *at) {
  return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx()) != 0;
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
RDKIT_GRAPHMOL_EXPORT int queryIsAtomBridgehead(Atom const *at);

static inline int queryIsBondInRing(Bond const *bond) {
  return bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx()) != 0;
};
static inline int queryAtomMinRingSize(Atom const *at) {
  return at->getOwningMol().getRingInfo()->minAtomRingSize(at->getIdx());
};
static inline int queryBondMinRingSize(Bond const *bond) {
  return bond->getOwningMol().getRingInfo()->minBondRingSize(bond->getIdx());
};

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
int queryAtomIsInRingOfSize(Atom const *at) {
  if (at->getOwningMol().getRingInfo()->isAtomInRingOfSize(at->getIdx(), tgt)) {
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
T *makeAtomSimpleQuery(int what, int func(Atom const *),
                       const std::string &description = "Atom Simple") {
  T *res = new T;
  res->setVal(what);
  res->setDataFunc(func);
  res->setDescription(description);
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
T *makeAtomNumQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNum, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNumQuery(int what);

//! returns a Query for matching atomic number and aromaticity
template <class T>
T *makeAtomTypeQuery(int num, int aromatic, const std::string &descr) {
  return makeAtomSimpleQuery<T>(makeAtomType(num, aromatic), queryAtomType,
                                descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomTypeQuery(int num,
                                                           int aromatic);

//! returns a Query for matching implicit valence
template <class T>
T *makeAtomImplicitValenceQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomImplicitValence, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what);

//! returns a Query for matching explicit valence
template <class T>
T *makeAtomExplicitValenceQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomExplicitValence, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomExplicitValenceQuery(int what);

//! returns a Query for matching total valence
template <class T>
T *makeAtomTotalValenceQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomTotalValence, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what);

//! returns a Query for matching explicit degree
template <class T>
T *makeAtomExplicitDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomExplicitDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what);

//! returns a Query for matching atomic degree
template <class T>
T *makeAtomTotalDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomTotalDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what);

//! returns a Query for matching heavy atom degree
template <class T>
T *makeAtomHeavyAtomDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomHeavyAtomDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHeavyAtomDegreeQuery(int what);

//! returns a Query for matching hydrogen count
template <class T>
T *makeAtomHCountQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomHCount, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what);

//! returns a Query for matching ring atoms
template <class T>
T *makeAtomHasImplicitHQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomHasImplicitH, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasImplicitHQuery();

//! returns a Query for matching implicit hydrogen count
template <class T>
T *makeAtomImplicitHCountQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomImplicitHCount, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomImplicitHCountQuery(int what);

//! returns a Query for matching the \c isAromatic flag
template <class T>
T *makeAtomAromaticQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomAromatic, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomAromaticQuery();

//! returns a Query for matching aliphatic atoms
template <class T>
T *makeAtomAliphaticQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomAliphatic, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomAliphaticQuery();

//! returns a Query for matching atoms with a particular mass
template <class T>
T *makeAtomMassQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(massIntegerConversionFactor * what,
                                queryAtomMass, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomMassQuery(int what);

//! returns a Query for matching atoms with a particular isotope
template <class T>
T *makeAtomIsotopeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomIsotope, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what);

//! returns a Query for matching formal charge
template <class T>
T *makeAtomFormalChargeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomFormalCharge, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what);

//! returns a Query for matching negative formal charges (i.e. a query val of 1
//! matches a formal charge of -1)
template <class T>
T *makeAtomNegativeFormalChargeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNegativeFormalCharge, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNegativeFormalChargeQuery(
    int what);

//! returns a Query for matching hybridization
template <class T>
T *makeAtomHybridizationQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomHybridization, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what);

//! returns a Query for matching the number of radical electrons
template <class T>
T *makeAtomNumRadicalElectronsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNumRadicalElectrons, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNumRadicalElectronsQuery(
    int what);

//! returns a Query for matching whether or not chirality has been set on the
//! atom
template <class T>
T *makeAtomHasChiralTagQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomHasChiralTag, descr);
}
//! \overloadquery
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasChiralTagQuery();

//! returns a Query for matching whether or not a potentially chiral atom is
//! missing a chiral tag
template <class T>
T *makeAtomMissingChiralTagQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomMissingChiralTag, descr);
}
//! \overloadquery
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomMissingChiralTagQuery();

//! returns a Query for matching atoms with unsaturation:
template <class T>
T *makeAtomUnsaturatedQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryAtomUnsaturated, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery();

//! returns a Query for matching ring atoms
template <class T>
T *makeAtomInRingQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryIsAtomInRing, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomInRingQuery();

//! returns a Query for matching atoms in a particular number of rings
template <class T>
T *makeAtomInNRingsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryIsAtomInNRings, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what);

//! returns a Query for matching atoms in rings of a particular size
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt);

//! returns a Query for matching an atom's minimum ring size
template <class T>
T *makeAtomMinRingSizeQuery(int tgt, const std::string &descr) {
  return makeAtomSimpleQuery<T>(tgt, queryAtomMinRingSize, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomMinRingSizeQuery(int tgt);

//! returns a Query for matching atoms with a particular number of ring bonds
template <class T>
T *makeAtomRingBondCountQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomRingBondCount, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what);

//! returns a Query for matching generic A atoms (heavy atoms)
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAAtomQuery();
//! returns a Query for matching generic AH atoms (any atom)
RDKIT_GRAPHMOL_EXPORT ATOM_NULL_QUERY *makeAHAtomQuery();
//! returns a Query for matching generic Q atoms (heteroatoms)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeQAtomQuery();
//! returns a Query for matching generic QH atoms (heteroatom or H)
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeQHAtomQuery();
//! returns a Query for matching generic X atoms (halogens)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeXAtomQuery();
//! returns a Query for matching generic XH atoms (halogen or H)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeXHAtomQuery();
//! returns a Query for matching generic M atoms (metals)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeMAtomQuery();
//! returns a Query for matching generic MH atoms (metals or H)
RDKIT_GRAPHMOL_EXPORT ATOM_OR_QUERY *makeMHAtomQuery();

// We support the same special atom queries that we can read from
// CXSMILES
const std::vector<std::string> complexQueries = {"A", "AH", "Q", "QH",
                                                 "X", "XH", "M", "MH"};
RDKIT_GRAPHMOL_EXPORT void convertComplexNameToQuery(Atom *query,
                                                     std::string_view symb);

//! returns a Query for matching atoms that have ring bonds
template <class T>
T *makeAtomHasRingBondQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(1, queryAtomHasRingBond, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasRingBondQuery();

//! returns a Query for matching the number of heteroatom neighbors
template <class T>
T *makeAtomNumHeteroatomNbrsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNumHeteroatomNbrs, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNumHeteroatomNbrsQuery(
    int what);

//! returns a Query for matching atoms that have heteroatom neighbors
template <class T>
T *makeAtomHasHeteroatomNbrsQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(1, queryAtomHasHeteroatomNbrs, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomHasHeteroatomNbrsQuery();

//! returns a Query for matching the number of aliphatic heteroatom neighbors
template <class T>
T *makeAtomNumAliphaticHeteroatomNbrsQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNumAliphaticHeteroatomNbrs,
                                descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *
makeAtomNumAliphaticHeteroatomNbrsQuery(int what);

//! returns a Query for matching atoms that have heteroatom neighbors
template <class T>
T *makeAtomHasAliphaticHeteroatomNbrsQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(1, queryAtomHasAliphaticHeteroatomNbrs, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *
makeAtomHasAliphaticHeteroatomNbrsQuery();

//! returns a Query for matching the number of non-hydrogen neighbors
template <class T>
T *makeAtomNonHydrogenDegreeQuery(int what, const std::string &descr) {
  return makeAtomSimpleQuery<T>(what, queryAtomNonHydrogenDegree, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomNonHydrogenDegreeQuery(
    int what);

//! returns a Query for matching bridgehead atoms
template <class T>
T *makeAtomIsBridgeheadQuery(const std::string &descr) {
  return makeAtomSimpleQuery<T>(true, queryIsAtomBridgehead, descr);
}
//! \overload
RDKIT_GRAPHMOL_EXPORT ATOM_EQUALS_QUERY *makeAtomIsBridgeheadQuery();

//! returns a Query for matching bond orders
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondOrderEqualsQuery(
    Bond::BondType what);
//! returns a Query for unspecified SMARTS bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrAromaticBondQuery();
//! returns a Query for double|aromatic bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeDoubleOrAromaticBondQuery();
//! returns a Query for single|double bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrDoubleBondQuery();
//! returns a Query for tautomeric bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *
makeSingleOrDoubleOrAromaticBondQuery();

//! returns a Query for matching bond directions
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondDirEqualsQuery(
    Bond::BondDir what);
//! returns a Query for matching bonds with stereo set
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondHasStereoQuery();
//! returns a Query for matching ring bonds
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondIsInRingQuery();
//! returns a Query for matching bonds in rings of a particular size
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int what);
//! returns a Query for matching a bond's minimum ring size
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondMinRingSizeQuery(int what);
//! returns a Query for matching bonds in a particular number of rings
RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeBondInNRingsQuery(int tgt);

//! returns a Query for matching any bond
RDKIT_GRAPHMOL_EXPORT BOND_NULL_QUERY *makeBondNullQuery();
//! returns a Query for matching any atom
RDKIT_GRAPHMOL_EXPORT ATOM_NULL_QUERY *makeAtomNullQuery();

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

class RDKIT_GRAPHMOL_EXPORT AtomRingQuery
    : public Queries::EqualityQuery<int, ConstAtomPtr, true> {
 public:
  AtomRingQuery() : Queries::EqualityQuery<int, ConstAtomPtr, true>(-1) {
    // default is to just do a number of rings query:
    this->setDescription("AtomInNRings");
    this->setDataFunc(queryAtomRingMembership);
  }
  explicit AtomRingQuery(int v)
      : Queries::EqualityQuery<int, ConstAtomPtr, true>(v) {
    // default is to just do a number of rings query:
    this->setDescription("AtomInNRings");
    this->setDataFunc(queryAtomRingMembership);
  }

  bool Match(const ConstAtomPtr what) const override {
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
  Queries::Query<int, ConstAtomPtr, true> *copy() const override {
    AtomRingQuery *res = new AtomRingQuery(this->d_val);
    res->setNegation(getNegation());
    res->setTol(this->getTol());
    res->d_description = this->d_description;
    res->d_dataFunc = this->d_dataFunc;
    return res;
  }
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
  std::string propname;

 public:
  HasPropQuery() : Queries::EqualityQuery<int, TargetPtr, true>(), propname() {
    // default is to just do a number of rings query:
    this->setDescription("AtomHasProp");
    this->setDataFunc(0);
  }
  explicit HasPropQuery(std::string v)
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(std::move(v)) {
    // default is to just do a number of rings query:
    this->setDescription("AtomHasProp");
    this->setDataFunc(nullptr);
  }

  bool Match(const TargetPtr what) const override {
    bool res = what->hasProp(propname);
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
};

typedef Queries::EqualityQuery<int, Atom const *, true> ATOM_PROP_QUERY;
typedef Queries::EqualityQuery<int, Bond const *, true> BOND_PROP_QUERY;

//! returns a Query for matching atoms that have a particular property
template <class Target>
Queries::EqualityQuery<int, const Target *, true> *makeHasPropQuery(
    const std::string &property) {
  return new HasPropQuery<const Target *>(property);
}

// ! Query whether an atom has a property with a value
template <class TargetPtr, class T>
class HasPropWithValueQuery
    : public Queries::EqualityQuery<int, TargetPtr, true> {
  std::string propname;
  T val;
  T tolerance;

 public:
  HasPropWithValueQuery()
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(), val() {
    // default is to just do a number of rings query:
    this->setDescription("HasPropWithValue");
    this->setDataFunc(0);
  }
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
    bool res = what->hasProp(propname);
    if (res) {
      try {
        T atom_val = what->template getProp<T>(propname);
        res = Queries::queryCmp(atom_val, this->val, this->tolerance) == 0;
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
    HasPropWithValueQuery *res =
        new HasPropWithValueQuery(this->propname, this->val, this->tolerance);
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  }
};

template <class TargetPtr>
class HasPropWithValueQuery<TargetPtr, std::string>
    : public Queries::EqualityQuery<int, TargetPtr, true> {
  std::string propname;
  std::string val;

 public:
  HasPropWithValueQuery()
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(), val() {
    // default is to just do a number of rings query:
    this->setDescription("HasPropWithValue");
    this->setDataFunc(0);
  }
  explicit HasPropWithValueQuery(std::string prop, std::string v,
                                 const std::string &tol = "")
      : Queries::EqualityQuery<int, TargetPtr, true>(),
        propname(std::move(prop)),
        val(std::move(v)) {
    RDUNUSED_PARAM(tol);
    // default is to just do a number of rings query:
    this->setDescription("HasPropWithValue");
    this->setDataFunc(nullptr);
  }

  bool Match(const TargetPtr what) const override {
    bool res = what->hasProp(propname);
    if (res) {
      try {
        std::string atom_val = what->template getProp<std::string>(propname);
        res = atom_val == this->val;
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
        new HasPropWithValueQuery<TargetPtr, std::string>(this->propname,
                                                          this->val);
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  }
};

template <class TargetPtr>
class HasPropWithValueQuery<TargetPtr, ExplicitBitVect>
    : public Queries::EqualityQuery<int, TargetPtr, true> {
  std::string propname;
  ExplicitBitVect val;
  float tol{0.0};

 public:
  HasPropWithValueQuery()
      : Queries::EqualityQuery<int, TargetPtr, true>(), propname(), val() {
    this->setDescription("HasPropWithValue");
    this->setDataFunc(0);
  }

  explicit HasPropWithValueQuery(std::string prop, const ExplicitBitVect &v,
                                 float tol = 0.0)
      : Queries::EqualityQuery<int, TargetPtr, true>(),
        propname(std::move(prop)),
        val(v),
        tol(tol) {
    this->setDescription("HasPropWithValue");
    this->setDataFunc(nullptr);
  }

  bool Match(const TargetPtr what) const override {
    bool res = what->hasProp(propname);
    if (res) {
      try {
        const ExplicitBitVect &bv =
            what->template getProp<const ExplicitBitVect &>(propname);
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
            this->propname, this->val, this->tol);
    res->setNegation(this->getNegation());
    res->d_description = this->d_description;
    return res;
  }
};

template <class Target, class T>
Queries::EqualityQuery<int, const Target *, true> *makePropQuery(
    const std::string &propname, const T &val, const T &tolerance = T()) {
  return new HasPropWithValueQuery<const Target *, T>(propname, val, tolerance);
}

template <class Target>
Queries::EqualityQuery<int, const Target *, true> *makePropQuery(
    const std::string &propname, const ExplicitBitVect &val,
    float tolerance = 0.0) {
  return new HasPropWithValueQuery<const Target *, ExplicitBitVect>(
      propname, val, tolerance);
}

RDKIT_GRAPHMOL_EXPORT bool isComplexQuery(const Bond *b);
RDKIT_GRAPHMOL_EXPORT bool isComplexQuery(const Atom *a);
RDKIT_GRAPHMOL_EXPORT bool isAtomAromatic(const Atom *a);
RDKIT_GRAPHMOL_EXPORT bool isAtomListQuery(const Atom *a);
RDKIT_GRAPHMOL_EXPORT void getAtomListQueryVals(const Atom::QUERYATOM_QUERY *q,
                                                std::vector<int> &vals);

// Checks if an atom is dummy or not.
// 1. A dummy non-query atom (e.g., "*" in SMILES) is defined by its zero atomic
//    number. This rule breaks for query atoms because a COMPOSITE_OR query atom
//    also has a zero atomic number (#6349).
// 2. A dummy query atom (e.g., "*" in SMARTS) is defined by its explicit
//    description: "AtomNull".
inline bool isAtomDummy(const Atom *a) {
  return (!a->hasQuery() && a->getAtomicNum() == 0) ||
         (a->hasQuery() && !a->getQuery()->getNegation() &&
          a->getQuery()->getDescription() == "AtomNull");
}

namespace QueryOps {
RDKIT_GRAPHMOL_EXPORT void completeMolQueries(
    RWMol *mol, unsigned int magicVal = 0xDEADBEEF);
// Replaces the given atom in the molecule with a QueryAtom that is otherwise
// a copy of the given atom.  Returns a pointer to that atom.
// if the atom already has a query, nothing will be changed
RDKIT_GRAPHMOL_EXPORT Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom);

RDKIT_GRAPHMOL_EXPORT void finalizeQueryFromDescription(
    Queries::Query<int, Atom const *, true> *query, Atom const *owner);
RDKIT_GRAPHMOL_EXPORT void finalizeQueryFromDescription(
    Queries::Query<int, Bond const *, true> *query, Bond const *owner);

RDKIT_GRAPHMOL_EXPORT bool hasBondTypeQuery(
    const Queries::Query<int, Bond const *, true> &qry);
inline bool hasBondTypeQuery(const Bond &bond) {
  if (!bond.hasQuery()) {
    return false;
  }
  return hasBondTypeQuery(*bond.getQuery());
}
RDKIT_GRAPHMOL_EXPORT bool hasComplexBondTypeQuery(
    const Queries::Query<int, Bond const *, true> &qry);
inline bool hasComplexBondTypeQuery(const Bond &bond) {
  if (!bond.hasQuery()) {
    return false;
  }
  return hasComplexBondTypeQuery(*bond.getQuery());
}

}  // namespace QueryOps
}  // namespace RDKit
#endif
