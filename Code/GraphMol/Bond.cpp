//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Bond.h"
#include "Atom.h"
#include "ROMol.h"
#include <RDGeneral/Invariant.h>

namespace RDKit {

Bond::Bond() : RDProps() { initBond(); };

Bond::Bond(BondType bT) : RDProps() {
  initBond();
  d_bondType = bT;
};

Bond::Bond(const Bond &other) : RDProps(other) {
  // NOTE: we do *not* copy ownership!
  dp_mol = nullptr;
  d_bondType = other.d_bondType;
  d_beginAtomIdx = other.d_beginAtomIdx;
  d_endAtomIdx = other.d_endAtomIdx;
  d_dirTag = other.d_dirTag;
  d_stereo = other.d_stereo;
  if (other.dp_stereoAtoms) {
    dp_stereoAtoms = new INT_VECT(*other.dp_stereoAtoms);
  } else {
    dp_stereoAtoms = nullptr;
  }
  df_isAromatic = other.df_isAromatic;
  df_isConjugated = other.df_isConjugated;
  d_index = other.d_index;
}

Bond::~Bond() { delete dp_stereoAtoms; }

Bond &Bond::operator=(const Bond &other) {
  if (this == &other) {
    return *this;
  }
  dp_mol = other.dp_mol;
  d_bondType = other.d_bondType;
  d_beginAtomIdx = other.d_beginAtomIdx;
  d_endAtomIdx = other.d_endAtomIdx;
  d_dirTag = other.d_dirTag;
  delete dp_stereoAtoms;
  if (other.dp_stereoAtoms) {
    dp_stereoAtoms = new INT_VECT(*other.dp_stereoAtoms);
  } else {
    dp_stereoAtoms = nullptr;
  }
  df_isAromatic = other.df_isAromatic;
  df_isConjugated = other.df_isConjugated;
  d_index = other.d_index;
  d_props = other.d_props;

  return *this;
}

Bond *Bond::copy() const {
  auto *res = new Bond(*this);
  return res;
}

void Bond::setOwningMol(ROMol *other) {
  // FIX: doesn't update topology
  dp_mol = other;
}

unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
  PRECONDITION(d_beginAtomIdx == thisIdx || d_endAtomIdx == thisIdx,
               "bad index");
  if (d_beginAtomIdx == thisIdx) {
    return d_endAtomIdx;
  } else if (d_endAtomIdx == thisIdx) {
    return d_beginAtomIdx;
  }
  // we cannot actually get down here
  return 0;
}

void Bond::setBeginAtomIdx(unsigned int what) {
  if (dp_mol) {
    URANGE_CHECK(what, getOwningMol().getNumAtoms());
  }
  d_beginAtomIdx = what;
};

void Bond::setEndAtomIdx(unsigned int what) {
  if (dp_mol) {
    URANGE_CHECK(what, getOwningMol().getNumAtoms());
  }
  d_endAtomIdx = what;
};

void Bond::setBeginAtom(Atom *at) {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  setBeginAtomIdx(at->getIdx());
}
void Bond::setEndAtom(Atom *at) {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  setEndAtomIdx(at->getIdx());
}

Atom *Bond::getBeginAtom() const {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  return dp_mol->getAtomWithIdx(d_beginAtomIdx);
};
Atom *Bond::getEndAtom() const {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");
  return dp_mol->getAtomWithIdx(d_endAtomIdx);
};
Atom *Bond::getOtherAtom(Atom const *what) const {
  PRECONDITION(dp_mol != nullptr, "no owning molecule for bond");

  return dp_mol->getAtomWithIdx(getOtherAtomIdx(what->getIdx()));
};

double Bond::getBondTypeAsDouble() const {
  double res;
  switch (getBondType()) {
    case UNSPECIFIED:
    case IONIC:
    case ZERO:
      res = 0;
      break;
    case SINGLE:
      res = 1;
      break;
    case DOUBLE:
      res = 2;
      break;
    case TRIPLE:
      res = 3;
      break;
    case QUADRUPLE:
      res = 4;
      break;
    case QUINTUPLE:
      res = 5;
      break;
    case HEXTUPLE:
      res = 6;
      break;
    case ONEANDAHALF:
      res = 1.5;
      break;
    case TWOANDAHALF:
      res = 2.5;
      break;
    case THREEANDAHALF:
      res = 3.5;
      break;
    case FOURANDAHALF:
      res = 4.5;
      break;
    case FIVEANDAHALF:
      res = 5.5;
      break;
    case AROMATIC:
      res = 1.5;
      break;
    case DATIVEONE:
      res = 1.0;
      break;  // FIX: this should probably be different
    case DATIVE:
      res = 1.0;
      break;  // FIX: again probably wrong
    case HYDROGEN:
      res = 0.0;
      break;
    default:
      UNDER_CONSTRUCTION("Bad bond type");
  }
  return res;
}

double Bond::getValenceContrib(const Atom *atom) const {
  if (atom != getBeginAtom() && atom != getEndAtom()) {
    return 0.0;
  }
  double res;
  if ((getBondType() == DATIVE || getBondType() == DATIVEONE) &&
      atom->getIdx() != getEndAtomIdx()) {
    res = 0.0;
  } else {
    res = getBondTypeAsDouble();
  }

  return res;
}

void Bond::setQuery(QUERYBOND_QUERY *) {
  //  Bonds don't have queries at the moment because I have not
  //  yet figured out what a good base query should be.
  //  It would be nice to be able to do substructure searches
  //  using molecules alone, so it'd be nice if we got this
  //  issue resolved ASAP.
  PRECONDITION(0, "plain bonds have no Query");
}

Bond::QUERYBOND_QUERY *Bond::getQuery() const {
  PRECONDITION(0, "plain bonds have no Query");
  return nullptr;
};

bool Bond::Match(Bond const *what) const {
  bool res;
  if (getBondType() == Bond::UNSPECIFIED ||
      what->getBondType() == Bond::UNSPECIFIED) {
    res = true;
  } else {
    res = getBondType() == what->getBondType();
  }
  return res;
};

void Bond::expandQuery(Bond::QUERYBOND_QUERY *, Queries::CompositeQueryType,
                       bool) {
  PRECONDITION(0, "plain bonds have no query");
};

void Bond::initBond() {
  d_bondType = UNSPECIFIED;
  d_dirTag = NONE;
  d_stereo = STEREONONE;
  dp_mol = nullptr;
  d_beginAtomIdx = 0;
  d_endAtomIdx = 0;
  df_isAromatic = 0;
  d_index = 0;
  df_isConjugated = 0;
  dp_stereoAtoms = nullptr;
};

void Bond::setStereoAtoms(unsigned int bgnIdx, unsigned int endIdx) {
  PRECONDITION(
      getOwningMol().getBondBetweenAtoms(getBeginAtomIdx(), bgnIdx) != nullptr,
      "bgnIdx not connected to begin atom of bond");
  PRECONDITION(
      getOwningMol().getBondBetweenAtoms(getEndAtomIdx(), endIdx) != nullptr,
      "endIdx not connected to end atom of bond");

  INT_VECT &atoms = getStereoAtoms();
  atoms.clear();
  atoms.push_back(bgnIdx);
  atoms.push_back(endIdx);
};

uint8_t getTwiceBondType(const Bond &b) {
  switch (b.getBondType()) {
    case Bond::UNSPECIFIED:
    case Bond::IONIC:
    case Bond::ZERO:
      return 0;
      break;
    case Bond::SINGLE:
      return 2;
      break;
    case Bond::DOUBLE:
      return 4;
      break;
    case Bond::TRIPLE:
      return 6;
      break;
    case Bond::QUADRUPLE:
      return 8;
      break;
    case Bond::QUINTUPLE:
      return 10;
      break;
    case Bond::HEXTUPLE:
      return 12;
      break;
    case Bond::ONEANDAHALF:
      return 3;
      break;
    case Bond::TWOANDAHALF:
      return 5;
      break;
    case Bond::THREEANDAHALF:
      return 7;
      break;
    case Bond::FOURANDAHALF:
      return 9;
      break;
    case Bond::FIVEANDAHALF:
      return 11;
      break;
    case Bond::AROMATIC:
      return 3;
      break;
    case Bond::DATIVEONE:
      return 2;
      break;  // FIX: this should probably be different
    case Bond::DATIVE:
      return 2;
      break;  // FIX: again probably wrong
    case Bond::HYDROGEN:
      return 0;
      break;
    default:
      UNDER_CONSTRUCTION("Bad bond type");
  }
}
};  // namespace RDKit

std::ostream &operator<<(std::ostream &target, const RDKit::Bond &bond) {
  target << bond.getIdx() << " ";
  target << bond.getBeginAtomIdx() << "->" << bond.getEndAtomIdx();
  target << " order: " << bond.getBondType();
  if (bond.getBondDir()) {
    target << " dir: " << bond.getBondDir();
  }
  if (bond.getStereo()) {
    target << " stereo: " << bond.getStereo();
    if (bond.getStereoAtoms().size() == 2) {
      const auto &ats = bond.getStereoAtoms();
      target << " stereoAts: (" << ats[0] << " " << ats[1] << ")";
    }
  }
  target << " conj?: " << bond.getIsConjugated();
  target << " aromatic?: " << bond.getIsAromatic();

  return target;
}
