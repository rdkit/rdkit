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
    INT_VECT res(d_atomMembers[idx].size());
    std::transform(d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
                   res.begin(),
                   [this](int ri) { return d_atomRings.at(ri).size(); });
    return res;
  }
  return INT_VECT();
}
bool RingInfo::isAtomInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_atomMembers.size()) {
    return std::find_if(d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
                        [this, size](int ri) {
                          return d_atomRings.at(ri).size() == size;
                        }) != d_atomMembers[idx].end();
  }
  return false;
}
unsigned int RingInfo::minAtomRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_atomMembers.size() && !d_atomMembers[idx].empty()) {
    auto ri = *std::min_element(
        d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
        [this](int ri1, int ri2) {
          return d_atomRings.at(ri1).size() < d_atomRings.at(ri2).size();
        });
    return d_atomRings.at(ri).size();
  }
  return 0;
}
unsigned int RingInfo::numAtomRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_atomMembers.size()) {
    return rdcast<unsigned int>(d_atomMembers[idx].size());
  }
  return 0;
}
RingInfo::INT_VECT RingInfo::atomMembers(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  static const INT_VECT emptyVect;
  if (idx < d_atomMembers.size()) {
    return d_atomMembers[idx];
  }
  return INT_VECT();
}
bool RingInfo::areAtomsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx1 >= d_atomMembers.size() || idx2 >= d_atomMembers.size()) {
    return false;
  }
  auto it1 = d_atomMembers[idx1].begin();
  auto it2 = d_atomMembers[idx2].begin();
  while (it1 != d_atomMembers[idx1].end() && it2 != d_atomMembers[idx2].end()) {
    if (*it1 < *it2) {
      ++it1;
    } else if (*it1 > *it2) {
      ++it2;
    } else if (!size || d_atomRings.at(*it1).size() == size) {
      return true;
    } else {
      ++it1;
      ++it2;
    }
  }
  return false;
}
RingInfo::INT_VECT RingInfo::bondRingSizes(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size()) {
    INT_VECT res(d_bondMembers[idx].size());
    std::transform(d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
                   res.begin(),
                   [this](int ri) { return d_bondRings.at(ri).size(); });
    return res;
  }
  return INT_VECT();
}
bool RingInfo::isBondInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size()) {
    return std::find_if(d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
                        [this, size](int ri) {
                          return d_bondRings.at(ri).size() == size;
                        }) != d_bondMembers[idx].end();
  }
  return false;
}
unsigned int RingInfo::minBondRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size() && d_bondMembers[idx].size()) {
    return d_bondRings
        .at(*std::min_element(
            d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
            [this](int ri1, int ri2) {
              return d_bondRings.at(ri1).size() < d_bondRings.at(ri2).size();
            }))
        .size();
  }
  return 0;
}
unsigned int RingInfo::numBondRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size()) {
    return rdcast<unsigned int>(d_bondMembers[idx].size());
  }
  return 0;
}
RingInfo::INT_VECT RingInfo::bondMembers(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx < d_bondMembers.size()) {
    return d_bondMembers[idx];
  }
  return INT_VECT();
}
bool RingInfo::areBondsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (idx1 >= d_bondMembers.size() || idx2 >= d_bondMembers.size()) {
    return false;
  }
  auto it1 = d_bondMembers[idx1].begin();
  auto it2 = d_bondMembers[idx2].begin();
  while (it1 != d_bondMembers[idx1].end() && it2 != d_bondMembers[idx2].end()) {
    if (*it1 < *it2) {
      ++it1;
    } else if (*it1 > *it2) {
      ++it2;
    } else if (!size || d_bondRings.at(*it1).size() == size) {
      return true;
    } else {
      ++it1;
      ++it2;
    }
  }
  return false;
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
  for (auto i : atomIndices) {
    if (i >= static_cast<int>(d_atomMembers.size())) {
      d_atomMembers.resize((i) + 1);
    }
    d_atomMembers[i].push_back(d_atomRings.size());
  }
  for (auto i : bondIndices) {
    if (i >= static_cast<int>(d_bondMembers.size())) {
      d_bondMembers.resize((i) + 1);
    }
    d_bondMembers[i].push_back(d_bondRings.size());
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
  if (!df_init) {
    return;
  }
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
