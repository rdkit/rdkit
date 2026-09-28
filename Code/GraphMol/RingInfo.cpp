//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Rings.h"
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
const RingInfo::INT_VECT &RingInfo::atomMembers(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  static const INT_VECT emptyVect;
  if (idx < d_atomMembers.size()) {
    return d_atomMembers[idx];
  }
  return emptyVect;
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
const RingInfo::INT_VECT &RingInfo::bondMembers(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  static const INT_VECT emptyVect;
  if (idx < d_bondMembers.size()) {
    return d_bondMembers[idx];
  }
  return emptyVect;
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
  for (const auto &i : atomIndices) {
    if (i >= static_cast<int>(d_atomMembers.size())) {
      d_atomMembers.resize(i + 1);
    }
    d_atomMembers[i].push_back(d_atomRings.size());
  }
  for (const auto &i : bondIndices) {
    if (i >= static_cast<int>(d_bondMembers.size())) {
      d_bondMembers.resize(i + 1);
    }
    d_bondMembers[i].push_back(d_bondRings.size());
  }
  d_atomRings.push_back(atomIndices);
  d_bondRings.push_back(bondIndices);
  POSTCONDITION(d_atomRings.size() == d_bondRings.size(), "length mismatch");
  return rdcast<unsigned int>(d_atomRings.size());
}

bool RingInfo::isRingFused(unsigned int ringIdx) {
  initFusedRings();
  PRECONDITION(ringIdx < d_fusedRings.size(), "ringIdx out of bounds");
  return d_fusedRings[ringIdx].any();
}

bool RingInfo::areRingsFused(unsigned int ring1Idx, unsigned int ring2Idx) {
  initFusedRings();
  PRECONDITION(ring1Idx < d_fusedRings.size(), "ring1Idx out of bounds");
  PRECONDITION(ring2Idx < d_fusedRings.size(), "ring2Idx out of bounds");
  return d_fusedRings[ring1Idx].test(ring2Idx);
}

unsigned int RingInfo::numFusedBonds(unsigned int ringIdx) {
  PRECONDITION(ringIdx < d_bondRings.size(), "ringIdx out of bounds");
  if (d_numFusedBonds.size() != d_bondRings.size()) {
    d_numFusedBonds.clear();
    d_numFusedBonds.resize(d_bondRings.size(), 0);
    for (unsigned int ri = 0; ri < d_bondRings.size(); ++ri) {
      d_numFusedBonds[ri] += std::count_if(
          d_bondRings[ri].begin(), d_bondRings[ri].end(),
          [this](unsigned int bi) { return numBondRings(bi) > 1; });
    }
  }
  return d_numFusedBonds[ringIdx];
}

unsigned int RingInfo::numFusedRingNeighbors(unsigned int ringIdx) {
  initFusedRings();
  PRECONDITION(ringIdx < d_fusedRings.size(), "ringIdx out of bounds");
  return d_fusedRings[ringIdx].count();
}

std::vector<unsigned int> RingInfo::fusedRingNeighbors(unsigned int ringIdx) {
  initFusedRings();
  PRECONDITION(ringIdx < d_fusedRings.size(), "ringIdx out of bounds");
  std::vector<unsigned int> res;
  res.reserve(d_fusedRings[ringIdx].count());
  for (unsigned int i = 0; i < d_fusedRings[ringIdx].size(); ++i) {
    if (d_fusedRings[ringIdx].test(i)) {
      res.push_back(i);
    }
  }
  return res;
}

void RingInfo::initFusedRings() {
  if (d_fusedRings.size() == d_bondRings.size()) {
    return;
  }
  d_fusedRings.clear();
  if (d_bondRings.empty()) {
    return;
  }
  d_fusedRings.resize(d_bondRings.size());
  for (auto &fusedRing : d_fusedRings) {
    fusedRing.resize(d_bondRings.size());
  }
  for (const auto &ringIndices : d_bondMembers) {
    if (ringIndices.size() <= 1) {
      continue;
    }
    for (unsigned int i = 0; i < ringIndices.size() - 1; ++i) {
      unsigned int ringIdx1 = ringIndices[i];
      for (unsigned int j = i + 1; j < ringIndices.size(); ++j) {
        unsigned int ringIdx2 = ringIndices[j];
        d_fusedRings[ringIdx1].set(ringIdx2);
        d_fusedRings[ringIdx2].set(ringIdx1);
      }
    }
  }
}

unsigned int RingInfo::numRingFamilies() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return d_atomRingFamilies.size();
};

unsigned int RingInfo::numRelevantCycles() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(dp_urfData, "Ring families not initialized");
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

void RingInfo::resetRingFamilies() {
  d_atomRingFamilies.clear();
  d_bondRingFamilies.clear();
  dp_urfData.reset();
}

void RingInfo::initialize(RDKit::FIND_RING_TYPE ringType) {
  df_init = true;
  df_find_type_type = ringType;
};
void RingInfo::reset(bool doRingFamilies) {
  if (!df_init) {
    return;
  }
  df_init = false;
  df_find_type_type = RDKit::FIND_RING_TYPE_OTHER_OR_UNKNOWN;
  d_atomMembers.clear();
  d_bondMembers.clear();
  d_atomRings.clear();
  d_bondRings.clear();
  if (doRingFamilies) {
    resetRingFamilies();
  }
}
void RingInfo::preallocate(unsigned int numAtoms, unsigned int numBonds) {
  d_atomMembers.resize(numAtoms);
  d_bondMembers.resize(numBonds);
}

std::vector<std::vector<int>> RingInfo::atomRelevantCycles() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(dp_urfData, "Ring families not initialized");
  std::vector<std::vector<int>> res;
  if (dp_urfData) {
    res.reserve(RDL_getNofRC(dp_urfData.get()));

    RDL_cycleIterator *it = RDL_getRCyclesIterator(dp_urfData.get());
    while (!RDL_cycleIteratorAtEnd(it)) {
      auto cycle = RDL_cycleIteratorGetCycle(it);
      res.push_back(RingUtils::rdlCycleToAtomRing(cycle));
      RDL_deleteCycle(cycle);
      RDL_cycleIteratorNext(it);
    }
    RDL_deleteCycleIterator(it);
  }
  return res;
}

}  // namespace RDKit
namespace RingUtils {
// normalizes a ring by rotating/reversing it so that the first atom
// is the one with the smallest index, and the second atom is the neighbor
// to the first one that again has the smallest index.
// This change should have a small performance footprint while it helps
// keeping test results consistent when making changes to ring detection.
void normalizeRing(std::vector<int> &ring) {
  auto newStart = std::ranges::min_element(ring);
  std::ranges::rotate(ring, newStart);

  if (ring.back() < ring[1]) {
    // we don't need to move the central element!
    auto numPairsToMove = (ring.size() - 1) / 2;
    auto front = ring.begin() + 1;
    std::swap_ranges(front, front + numPairsToMove, ring.rbegin());
  }
}

std::vector<int> rdlCycleToAtomRing(RDL_cycle *cycle) {
  PRECONDITION(cycle, "cycle is null");
  std::vector<int> ring;
  ring.reserve(cycle->weight);

  // Edges in a cycle are not returned in iteration order.
  // so we need to take care of that while we convert them
  // into an atom ring.
  boost::dynamic_bitset<> unseen_edges(cycle->weight);
  unseen_edges.set();

  ring.push_back(cycle->edges[0][0]);
  ring.push_back(cycle->edges[0][1]);
  unseen_edges.set(0, false);

  while (ring.size() < cycle->weight) {
    // Note we don't want to close the cycle: that would
    // add the initial atom at the end too.
    for (auto edgeIdx = unseen_edges.find_first();
         edgeIdx != boost::dynamic_bitset<>::npos;
         edgeIdx = unseen_edges.find_next(edgeIdx)) {
      auto edge = cycle->edges[edgeIdx];
      for (auto j = 0; j < 2; ++j) {
        if (static_cast<unsigned int>(ring.back()) == edge[j]) {
          ring.push_back(edge[1 - j]);
          unseen_edges.set(edgeIdx, false);
          break;
        }
      }
      if (ring.size() == cycle->weight) {
        break;
      }
    }
  }

  // For consistency, normalize the ring
  normalizeRing(ring);

  return ring;
}
}  // namespace RingUtils
