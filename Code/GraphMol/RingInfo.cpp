//
//  Copyright (C) 2004-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RingInfo.h"
#include <RDGeneral/Invariant.h>
#include <algorithm>

namespace RDKit {
RingInfo::INT_VECT RingInfo::atomRingSizes(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_atomMembers.size()) {
    return d_atomMembers[idx];
  } else {
    return INT_VECT{0};
  }
}
bool RingInfo::isAtomInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_atomMembers.size()) {
    return std::find(d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
                     static_cast<int>(size)) != d_atomMembers[idx].end();
  } else {
    return false;
  }
}
unsigned int RingInfo::minAtomRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_atomMembers.size() && d_atomMembers[idx].size()) {
    return *std::min_element(d_atomMembers[idx].begin(),
                             d_atomMembers[idx].end());
  } else {
    return 0;
  }
}
unsigned int RingInfo::numAtomRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_atomMembers.size()) {
    return rdcast<unsigned int>(d_atomMembers[idx].size());
  } else {
    return 0;
  }
}
RingInfo::INT_VECT RingInfo::bondRingSizes(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size()) {
    return d_bondMembers[idx];
  } else {
    return INT_VECT{0};
  }
}
bool RingInfo::isBondInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size()) {
    return std::find(d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
                     static_cast<int>(size)) != d_bondMembers[idx].end();
  } else {
    return false;
  }
}
unsigned int RingInfo::minBondRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size() && d_bondMembers[idx].size()) {
    return *std::min_element(d_bondMembers[idx].begin(),
                             d_bondMembers[idx].end());
  } else {
    return 0;
  }
}
unsigned int RingInfo::numBondRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size()) {
    return rdcast<unsigned int>(d_bondMembers[idx].size());
  } else {
    return 0;
  }
}

unsigned int RingInfo::numRings() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(d_atomRings.size() == d_bondRings.size(), "length mismatch");
  return rdcast<unsigned int>(d_atomRings.size());
}

unsigned int RingInfo::addRing(const INT_VECT &atomIndices,
                               const INT_VECT &bondIndices) {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(atomIndices.size() == bondIndices.size(), "length mismatch");
  int sz = rdcast<int>(atomIndices.size());
  for (auto i : atomIndices) {
    if (i >= static_cast<int>(d_atomMembers.size()))
      d_atomMembers.resize((i) + 1);
    d_atomMembers[i].push_back(sz);
  }
  for (auto i : bondIndices) {
    if (i >= static_cast<int>(d_bondMembers.size()))
      d_bondMembers.resize((i) + 1);
    d_bondMembers[i].push_back(sz);
  }
  d_atomRings.push_back(atomIndices);
  d_bondRings.push_back(bondIndices);
  POSTCONDITION(d_atomRings.size() == d_bondRings.size(), "length mismatch");
  return rdcast<unsigned int>(d_atomRings.size());
}

#ifdef RDK_USE_URF
unsigned int RingInfo::numRingFamilies() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return d_atomRingFamilies.size();
};

unsigned int RingInfo::numRelevantCycles() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return rdcast<unsigned int>(RDL_getNofRC(dp_urfData.get()));
};

unsigned int RingInfo::addRingFamily(const INT_VECT &atomIndices,
                                     const INT_VECT &bondIndices) {
  PRECONDITION(df_init, "RingInfo not initialized");
  d_atomRingFamilies.push_back(atomIndices);
  d_bondRingFamilies.push_back(bondIndices);
  POSTCONDITION(d_atomRingFamilies.size() == d_bondRingFamilies.size(),
                "length mismatch");

  return rdcast<unsigned int>(d_atomRingFamilies.size());
}
#endif

void RingInfo::initialize() {
  PRECONDITION(!df_init, "already initialized");
  df_init = true;
};
void RingInfo::reset() {
  if (!df_init) return;
  df_init = false;
  d_atomMembers.clear();
  d_bondMembers.clear();
  d_atomRings.clear();
  d_bondRings.clear();
}
void RingInfo::preallocate(unsigned int numAtoms, unsigned int numBonds) {
  d_atomMembers.resize(numAtoms);
  d_bondMembers.resize(numBonds);
}
}  // namespace RDKit
