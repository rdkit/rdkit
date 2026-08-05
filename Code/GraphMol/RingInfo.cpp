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
#include <numeric>

namespace RDKit {
namespace {
RingInfo::VECT_INT_VECT materializeRings(RingInfo::RingsView rings) {
  RingInfo::VECT_INT_VECT result;
  result.reserve(rings.size());
  for (const auto ring : rings) {
    result.emplace_back(ring.begin(), ring.end());
  }
  return result;
}
}  // namespace

RingInfo::VECT_INT_VECT RingInfo::atomRingsAsVectors() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return materializeRings(atomRings());
}

RingInfo::VECT_INT_VECT RingInfo::bondRingsAsVectors() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return materializeRings(bondRings());
}

RingInfo::INT_VECT RingInfo::atomRingSizes(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  const auto members = atomMembers(idx);
  if (!members.empty()) {
    INT_VECT res(members.size());
    std::transform(members.begin(), members.end(), res.begin(),
                   [this](int ri) { return atomRings().at(ri).size(); });
    return res;
  }
  return INT_VECT();
}
bool RingInfo::isAtomInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  const auto members = atomMembers(idx);
  if (!members.empty()) {
    return std::find_if(members.begin(), members.end(), [this, size](int ri) {
             return atomRings().at(ri).size() == size;
           }) != members.end();
  }
  return false;
}
unsigned int RingInfo::minAtomRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  const auto members = atomMembers(idx);
  if (!members.empty()) {
    auto ri = *std::min_element(
        members.begin(), members.end(), [this](int ri1, int ri2) {
          return atomRings().at(ri1).size() < atomRings().at(ri2).size();
        });
    return atomRings().at(ri).size();
  }
  return 0;
}
unsigned int RingInfo::numAtomRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  return rdcast<unsigned int>(atomMembers(idx).size());
}
RingInfo::IntView RingInfo::atomMembers(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (!d_atomMembershipBegins.empty() &&
      idx < d_atomMembershipBegins.size() - 1) {
    const auto begin = d_atomMembershipBegins[idx];
    const auto size = d_atomMembershipBegins[idx + 1] - begin;
    return {size ? d_atomMemberships.data() + begin : nullptr, size};
  }
  return {};
}
bool RingInfo::areAtomsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  const auto members1 = atomMembers(idx1);
  const auto members2 = atomMembers(idx2);
  if (members1.empty() || members2.empty()) {
    return false;
  }
  auto it1 = members1.begin();
  auto it2 = members2.begin();
  while (it1 != members1.end() && it2 != members2.end()) {
    if (*it1 < *it2) {
      ++it1;
    } else if (*it1 > *it2) {
      ++it2;
    } else if (!size || atomRings().at(*it1).size() == size) {
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

  const auto members = bondMembers(idx);
  if (!members.empty()) {
    INT_VECT res(members.size());
    std::transform(members.begin(), members.end(), res.begin(),
                   [this](int ri) { return bondRings().at(ri).size(); });
    return res;
  }
  return INT_VECT();
}
bool RingInfo::isBondInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  const auto members = bondMembers(idx);
  if (!members.empty()) {
    return std::find_if(members.begin(), members.end(), [this, size](int ri) {
             return bondRings().at(ri).size() == size;
           }) != members.end();
  }
  return false;
}
unsigned int RingInfo::minBondRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  const auto members = bondMembers(idx);
  if (!members.empty()) {
    return bondRings()
        .at(*std::min_element(members.begin(), members.end(),
                              [this](int ri1, int ri2) {
                                return bondRings().at(ri1).size() <
                                       bondRings().at(ri2).size();
                              }))
        .size();
  }
  return 0;
}
unsigned int RingInfo::numBondRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  return rdcast<unsigned int>(bondMembers(idx).size());
}
RingInfo::IntView RingInfo::bondMembers(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  if (!d_bondMembershipBegins.empty() &&
      idx < d_bondMembershipBegins.size() - 1) {
    const auto begin = d_bondMembershipBegins[idx];
    const auto size = d_bondMembershipBegins[idx + 1] - begin;
    return {size ? d_bondMemberships.data() + begin : nullptr, size};
  }
  return {};
}
bool RingInfo::areBondsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");

  const auto members1 = bondMembers(idx1);
  const auto members2 = bondMembers(idx2);
  if (members1.empty() || members2.empty()) {
    return false;
  }
  auto it1 = members1.begin();
  auto it2 = members2.begin();
  while (it1 != members1.end() && it2 != members2.end()) {
    if (*it1 < *it2) {
      ++it1;
    } else if (*it1 > *it2) {
      ++it2;
    } else if (!size || bondRings().at(*it1).size() == size) {
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
  PRECONDITION(atomRings().size() == bondRings().size(), "length mismatch");
  return rdcast<unsigned int>(atomRings().size());
}

unsigned int RingInfo::addRing(const INT_VECT &atomIndices,
                               const INT_VECT &bondIndices) {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(atomIndices.size() == bondIndices.size(), "length mismatch");
  PRECONDITION(
      std::ranges::none_of(atomIndices, [](int idx) { return idx < 0; }),
      "atom index must be non-negative");
  PRECONDITION(
      std::ranges::none_of(bondIndices, [](int idx) { return idx < 0; }),
      "bond index must be non-negative");
  if (d_atomRingBegins.empty()) {
    d_atomRingBegins.push_back(0);
  }
  if (d_bondRingBegins.empty()) {
    d_bondRingBegins.push_back(0);
  }
  d_atomsInRings.insert(d_atomsInRings.end(), atomIndices.begin(),
                        atomIndices.end());
  d_bondsInRings.insert(d_bondsInRings.end(), bondIndices.begin(),
                        bondIndices.end());
  d_atomRingBegins.push_back(rdcast<uint32_t>(d_atomsInRings.size()));
  d_bondRingBegins.push_back(rdcast<uint32_t>(d_bondsInRings.size()));
  rebuildMemberships();
  invalidateFusedRings();
  POSTCONDITION(atomRings().size() == bondRings().size(), "length mismatch");
  return rdcast<unsigned int>(atomRings().size());
}

unsigned int RingInfo::addRing(IntView atomIndices, IntView bondIndices) {
  return addRing(INT_VECT(atomIndices.begin(), atomIndices.end()),
                 INT_VECT(bondIndices.begin(), bondIndices.end()));
}

unsigned int RingInfo::addRings(const VECT_INT_VECT &atomRings,
                                const VECT_INT_VECT &bondRings) {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(atomRings.size() == bondRings.size(), "length mismatch");
  for (size_t i = 0; i < atomRings.size(); ++i) {
    PRECONDITION(atomRings[i].size() == bondRings[i].size(), "length mismatch");
    PRECONDITION(
        std::ranges::none_of(atomRings[i], [](int idx) { return idx < 0; }),
        "atom index must be non-negative");
    PRECONDITION(
        std::ranges::none_of(bondRings[i], [](int idx) { return idx < 0; }),
        "bond index must be non-negative");
  }
  if (d_atomRingBegins.empty()) {
    d_atomRingBegins.push_back(0);
  }
  if (d_bondRingBegins.empty()) {
    d_bondRingBegins.push_back(0);
  }
  size_t atomCount = 0, bondCount = 0;
  for (const auto &ring : atomRings) atomCount += ring.size();
  for (const auto &ring : bondRings) bondCount += ring.size();
  d_atomsInRings.reserve(d_atomsInRings.size() + atomCount);
  d_bondsInRings.reserve(d_bondsInRings.size() + bondCount);
  d_atomRingBegins.reserve(d_atomRingBegins.size() + atomRings.size());
  d_bondRingBegins.reserve(d_bondRingBegins.size() + bondRings.size());
  for (size_t i = 0; i < atomRings.size(); ++i) {
    d_atomsInRings.insert(d_atomsInRings.end(), atomRings[i].begin(),
                          atomRings[i].end());
    d_bondsInRings.insert(d_bondsInRings.end(), bondRings[i].begin(),
                          bondRings[i].end());
    d_atomRingBegins.push_back(rdcast<uint32_t>(d_atomsInRings.size()));
    d_bondRingBegins.push_back(rdcast<uint32_t>(d_bondsInRings.size()));
  }
  rebuildMemberships();
  invalidateFusedRings();
  return numRings();
}

void RingInfo::rebuildMemberships() {
  size_t numAtoms =
      d_atomMembershipBegins.empty() ? 0 : d_atomMembershipBegins.size() - 1;
  size_t numBonds =
      d_bondMembershipBegins.empty() ? 0 : d_bondMembershipBegins.size() - 1;
  for (const auto idx : d_atomsInRings)
    numAtoms = std::max(numAtoms, size_t(idx + 1));
  for (const auto idx : d_bondsInRings)
    numBonds = std::max(numBonds, size_t(idx + 1));
  d_atomMembershipBegins.assign(numAtoms + 1, 0);
  d_bondMembershipBegins.assign(numBonds + 1, 0);
  for (const auto idx : d_atomsInRings) ++d_atomMembershipBegins[idx + 1];
  for (const auto idx : d_bondsInRings) ++d_bondMembershipBegins[idx + 1];
  std::partial_sum(d_atomMembershipBegins.begin(), d_atomMembershipBegins.end(),
                   d_atomMembershipBegins.begin());
  std::partial_sum(d_bondMembershipBegins.begin(), d_bondMembershipBegins.end(),
                   d_bondMembershipBegins.begin());
  d_atomMemberships.resize(d_atomsInRings.size());
  d_bondMemberships.resize(d_bondsInRings.size());
  auto atomNext = d_atomMembershipBegins;
  auto bondNext = d_bondMembershipBegins;
  for (size_t ringIdx = 0; ringIdx < atomRings().size(); ++ringIdx) {
    for (const auto idx : atomRings()[ringIdx])
      d_atomMemberships[atomNext[idx]++] = rdcast<int>(ringIdx);
    for (const auto idx : bondRings()[ringIdx])
      d_bondMemberships[bondNext[idx]++] = rdcast<int>(ringIdx);
  }
}

void RingInfo::invalidateFusedRings() {
  d_fusedRings.clear();
  d_numFusedBonds.clear();
  d_fusedRingsInitialized = false;
}

bool RingInfo::isRingFused(unsigned int ringIdx) {
  initFusedRings();
  PRECONDITION(ringIdx < numRings(), "ringIdx out of bounds");
  const auto n = numRings();
  return std::find(d_fusedRings.begin() + ringIdx * n,
                   d_fusedRings.begin() + (ringIdx + 1) * n,
                   true) != d_fusedRings.begin() + (ringIdx + 1) * n;
}

bool RingInfo::areRingsFused(unsigned int ring1Idx, unsigned int ring2Idx) {
  initFusedRings();
  const auto n = numRings();
  PRECONDITION(ring1Idx < n, "ring1Idx out of bounds");
  PRECONDITION(ring2Idx < n, "ring2Idx out of bounds");
  return d_fusedRings[ring1Idx * n + ring2Idx];
}

unsigned int RingInfo::numFusedBonds(unsigned int ringIdx) {
  PRECONDITION(ringIdx < numRings(), "ringIdx out of bounds");
  initFusedRings();
  return d_numFusedBonds[ringIdx];
}

unsigned int RingInfo::numFusedRingNeighbors(unsigned int ringIdx) {
  initFusedRings();
  const auto n = numRings();
  PRECONDITION(ringIdx < n, "ringIdx out of bounds");
  return std::count(d_fusedRings.begin() + ringIdx * n,
                    d_fusedRings.begin() + (ringIdx + 1) * n, true);
}

std::vector<unsigned int> RingInfo::fusedRingNeighbors(unsigned int ringIdx) {
  initFusedRings();
  const auto n = numRings();
  PRECONDITION(ringIdx < n, "ringIdx out of bounds");
  std::vector<unsigned int> res;
  res.reserve(numFusedRingNeighbors(ringIdx));
  for (unsigned int i = 0; i < n; ++i) {
    if (d_fusedRings[ringIdx * n + i]) {
      res.push_back(i);
    }
  }
  return res;
}

void RingInfo::initFusedRings() {
  if (d_fusedRingsInitialized) {
    return;
  }
  const auto n = numRings();
  d_fusedRings.assign(n * n, false);
  d_numFusedBonds.assign(n, 0);
  for (unsigned int bondIdx = 0; bondIdx + 1 < d_bondMembershipBegins.size();
       ++bondIdx) {
    const auto ringIndices = bondMembers(bondIdx);
    if (ringIndices.size() <= 1) {
      continue;
    }
    for (unsigned int i = 0; i < ringIndices.size() - 1; ++i) {
      unsigned int ringIdx1 = ringIndices[i];
      for (unsigned int j = i + 1; j < ringIndices.size(); ++j) {
        unsigned int ringIdx2 = ringIndices[j];
        d_fusedRings[ringIdx1 * n + ringIdx2] = true;
        d_fusedRings[ringIdx2 * n + ringIdx1] = true;
      }
    }
    for (const auto ringIdx : ringIndices) ++d_numFusedBonds[ringIdx];
  }
  d_fusedRingsInitialized = true;
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
  d_atomRingBegins.assign(1, 0);
  d_bondRingBegins.assign(1, 0);
  d_atomsInRings.clear();
  d_bondsInRings.clear();
  d_atomMembershipBegins.assign(1, 0);
  d_bondMembershipBegins.assign(1, 0);
  d_atomMemberships.clear();
  d_bondMemberships.clear();
  invalidateFusedRings();
  if (doRingFamilies) {
    resetRingFamilies();
  }
}
void RingInfo::preallocate(unsigned int numAtoms, unsigned int numBonds) {
  if (d_atomMembershipBegins.empty()) {
    d_atomMembershipBegins.push_back(0);
  }
  if (d_bondMembershipBegins.empty()) {
    d_bondMembershipBegins.push_back(0);
  }
  if (d_atomMembershipBegins.size() < numAtoms + 1) {
    d_atomMembershipBegins.resize(numAtoms + 1, d_atomMembershipBegins.back());
  }
  if (d_bondMembershipBegins.size() < numBonds + 1) {
    d_bondMembershipBegins.resize(numBonds + 1, d_bondMembershipBegins.back());
  }
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
