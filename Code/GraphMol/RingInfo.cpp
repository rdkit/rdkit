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
#include <limits>
#include <numeric>
#include <utility>

namespace RDKit {
namespace {
void validateRing(std::span<const int> atomIndices,
                  std::span<const int> bondIndices) {
  PRECONDITION(atomIndices.size() == bondIndices.size(), "length mismatch");
  PRECONDITION(
      std::ranges::none_of(atomIndices, [](int idx) { return idx < 0; }),
      "atom index must be non-negative");
  PRECONDITION(
      std::ranges::none_of(bondIndices, [](int idx) { return idx < 0; }),
      "bond index must be non-negative");
}

unsigned int ringSizeUnchecked(const std::vector<uint32_t> &ringBegins,
                               unsigned int ringIdx) {
  return ringBegins[ringIdx + 1] - ringBegins[ringIdx];
}

std::pair<uint32_t, uint32_t> membershipRange(
    const std::vector<uint32_t> &membershipBegins, unsigned int idx) {
  if (membershipBegins.empty() || idx >= membershipBegins.size() - 1) {
    return {0, 0};
  }
  return {membershipBegins[idx], membershipBegins[idx + 1]};
}

std::vector<int> memberRingSizes(
    const std::vector<uint32_t> &membershipBegins,
    const std::vector<int> &memberships,
    const std::vector<uint32_t> &ringBegins, unsigned int idx) {
  const auto [begin, end] = membershipRange(membershipBegins, idx);
  std::vector<int> result(end - begin);
  for (auto pos = begin; pos < end; ++pos) {
    result[pos - begin] = ringSizeUnchecked(ringBegins, memberships[pos]);
  }
  return result;
}

bool isMemberInRingOfSize(const std::vector<uint32_t> &membershipBegins,
                          const std::vector<int> &memberships,
                          const std::vector<uint32_t> &ringBegins,
                          unsigned int idx, unsigned int size) {
  const auto [begin, end] = membershipRange(membershipBegins, idx);
  for (auto pos = begin; pos < end; ++pos) {
    if (ringSizeUnchecked(ringBegins, memberships[pos]) == size) {
      return true;
    }
  }
  return false;
}

unsigned int minMemberRingSize(
    const std::vector<uint32_t> &membershipBegins,
    const std::vector<int> &memberships,
    const std::vector<uint32_t> &ringBegins, unsigned int idx) {
  const auto [begin, end] = membershipRange(membershipBegins, idx);
  if (begin == end) {
    return 0;
  }
  auto result = std::numeric_limits<unsigned int>::max();
  for (auto pos = begin; pos < end; ++pos) {
    result =
        std::min(result, ringSizeUnchecked(ringBegins, memberships[pos]));
  }
  return result;
}

unsigned int numMemberRings(const std::vector<uint32_t> &membershipBegins,
                            unsigned int idx) {
  const auto [begin, end] = membershipRange(membershipBegins, idx);
  return end - begin;
}

size_t fusedRingPairIndex(unsigned int ring1Idx, unsigned int ring2Idx) {
  const auto high = std::max(ring1Idx, ring2Idx);
  const auto low = std::min(ring1Idx, ring2Idx);
  return size_t(high) * (high - 1) / 2 + low;
}

bool areRingsFusedUnchecked(const std::vector<bool> &fusedRings,
                            unsigned int ring1Idx, unsigned int ring2Idx) {
  return ring1Idx != ring2Idx &&
         fusedRings[fusedRingPairIndex(ring1Idx, ring2Idx)];
}

template <bool checkRingSize>
bool areMembersInSameRing(const std::vector<uint32_t> &membershipBegins,
                          const std::vector<int> &memberships,
                          const std::vector<uint32_t> &ringBegins,
                          unsigned int idx1, unsigned int idx2,
                          unsigned int ringSize = 0) {
  if (membershipBegins.empty() || idx1 >= membershipBegins.size() - 1 ||
      idx2 >= membershipBegins.size() - 1) {
    return false;
  }
  auto pos1 = membershipBegins[idx1];
  auto pos2 = membershipBegins[idx2];
  const auto end1 = membershipBegins[idx1 + 1];
  const auto end2 = membershipBegins[idx2 + 1];

  if (end1 - pos1 == 1 && end2 - pos2 == 1) {
    const auto ring1 = memberships[pos1];
    if (ring1 != memberships[pos2]) {
      return false;
    }
    if constexpr (checkRingSize) {
      return ringSizeUnchecked(ringBegins, ring1) == ringSize;
    }
    return true;
  }

  while (pos1 < end1 && pos2 < end2) {
    const auto ring1 = memberships[pos1];
    const auto ring2 = memberships[pos2];
    if (ring1 < ring2) {
      ++pos1;
    } else if (ring1 > ring2) {
      ++pos2;
    } else if constexpr (!checkRingSize) {
      return true;
    } else if (ringSizeUnchecked(ringBegins, ring1) == ringSize) {
      return true;
    } else {
      ++pos1;
      ++pos2;
    }
  }
  return false;
}

void rebuildRingTypeMemberships(
    std::vector<uint32_t> &membershipBegins,
    std::vector<int> &memberships, const std::vector<int> &membersInRings,
    const std::vector<uint32_t> &ringBegins) {
  size_t numMembers =
      membershipBegins.empty() ? 0 : membershipBegins.size() - 1;
  for (const auto idx : membersInRings) {
    numMembers = std::max(numMembers, size_t(idx) + 1);
  }
  membershipBegins.assign(numMembers + 1, 0);
  for (const auto idx : membersInRings) {
    ++membershipBegins[idx + 1];
  }
  std::partial_sum(membershipBegins.begin(), membershipBegins.end(),
                   membershipBegins.begin());

  memberships.resize(membersInRings.size());
  auto next = membershipBegins;
  for (size_t ringIdx = 0; ringIdx + 1 < ringBegins.size(); ++ringIdx) {
    for (auto pos = ringBegins[ringIdx]; pos < ringBegins[ringIdx + 1]; ++pos) {
      memberships[next[membersInRings[pos]]++] = rdcast<int>(ringIdx);
    }
  }
}
}  // namespace

void RingInfo::checkInitialized() const {
  PRECONDITION(df_init, "RingInfo not initialized");
}

RingInfo::INT_VECT RingInfo::atomRingSizes(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return memberRingSizes(d_atomMembershipBegins, d_atomMemberships,
                         d_atomRingBegins, idx);
}
bool RingInfo::isAtomInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return isMemberInRingOfSize(d_atomMembershipBegins, d_atomMemberships,
                              d_atomRingBegins, idx, size);
}
unsigned int RingInfo::minAtomRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return minMemberRingSize(d_atomMembershipBegins, d_atomMemberships,
                           d_atomRingBegins, idx);
}
unsigned int RingInfo::numAtomRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return numMemberRings(d_atomMembershipBegins, idx);
}
bool RingInfo::areAtomsInSameRing(unsigned int idx1,
                                  unsigned int idx2) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return areMembersInSameRing<false>(d_atomMembershipBegins,
                                     d_atomMemberships, d_atomRingBegins, idx1,
                                     idx2);
}
bool RingInfo::areAtomsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  if (!size) {
    return areMembersInSameRing<false>(d_atomMembershipBegins,
                                       d_atomMemberships, d_atomRingBegins,
                                       idx1, idx2);
  }
  return areMembersInSameRing<true>(d_atomMembershipBegins, d_atomMemberships,
                                    d_atomRingBegins, idx1, idx2, size);
}
RingInfo::INT_VECT RingInfo::bondRingSizes(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return memberRingSizes(d_bondMembershipBegins, d_bondMemberships,
                         d_bondRingBegins, idx);
}
bool RingInfo::isBondInRingOfSize(unsigned int idx, unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return isMemberInRingOfSize(d_bondMembershipBegins, d_bondMemberships,
                              d_bondRingBegins, idx, size);
}
unsigned int RingInfo::minBondRingSize(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return minMemberRingSize(d_bondMembershipBegins, d_bondMemberships,
                           d_bondRingBegins, idx);
}
unsigned int RingInfo::numBondRings(unsigned int idx) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return numMemberRings(d_bondMembershipBegins, idx);
}
bool RingInfo::areBondsInSameRing(unsigned int idx1,
                                  unsigned int idx2) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  return areMembersInSameRing<false>(d_bondMembershipBegins,
                                     d_bondMemberships, d_bondRingBegins, idx1,
                                     idx2);
}
bool RingInfo::areBondsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  PRECONDITION(df_init, "RingInfo not initialized");
  if (!size) {
    return areMembersInSameRing<false>(d_bondMembershipBegins,
                                       d_bondMemberships, d_bondRingBegins,
                                       idx1, idx2);
  }
  return areMembersInSameRing<true>(d_bondMembershipBegins, d_bondMemberships,
                                    d_bondRingBegins, idx1, idx2, size);
}

unsigned int RingInfo::numRings() const {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(d_atomRingBegins.size() == d_bondRingBegins.size(),
               "length mismatch");
  return d_atomRingBegins.empty()
             ? 0
             : rdcast<unsigned int>(d_atomRingBegins.size() - 1);
}

void RingInfo::appendRingUnchecked(std::span<const int> atomIndices,
                                   std::span<const int> bondIndices) {
  // A default-moved RingInfo may have lost its leading zero offsets.
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
}

unsigned int RingInfo::finishRingUpdates() {
  rebuildMemberships();
  invalidateFusedRings();
  POSTCONDITION(atomRings().size() == bondRings().size(), "length mismatch");
  return rdcast<unsigned int>(atomRings().size());
}

unsigned int RingInfo::addRing(std::span<const int> atomIndices,
                               std::span<const int> bondIndices) {
  PRECONDITION(df_init, "RingInfo not initialized");
  validateRing(atomIndices, bondIndices);
  appendRingUnchecked(atomIndices, bondIndices);
  return finishRingUpdates();
}

unsigned int RingInfo::addRings(const VECT_INT_VECT &atomRings,
                                const VECT_INT_VECT &bondRings) {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(atomRings.size() == bondRings.size(), "length mismatch");
  size_t atomCount = 0;
  size_t bondCount = 0;
  for (size_t i = 0; i < atomRings.size(); ++i) {
    validateRing(atomRings[i], bondRings[i]);
    atomCount += atomRings[i].size();
    bondCount += bondRings[i].size();
  }
  d_atomsInRings.reserve(d_atomsInRings.size() + atomCount);
  d_bondsInRings.reserve(d_bondsInRings.size() + bondCount);
  d_atomRingBegins.reserve(d_atomRingBegins.size() + atomRings.size());
  d_bondRingBegins.reserve(d_bondRingBegins.size() + bondRings.size());
  for (size_t i = 0; i < atomRings.size(); ++i) {
    appendRingUnchecked(atomRings[i], bondRings[i]);
  }
  return finishRingUpdates();
}

unsigned int RingInfo::addRings(const RingInfo &source,
                                const std::vector<unsigned int> &ringIndices) {
  PRECONDITION(df_init, "RingInfo not initialized");
  PRECONDITION(source.df_init, "source RingInfo not initialized");
  PRECONDITION(&source != this, "source and destination RingInfo must differ");
  size_t atomCount = 0;
  size_t bondCount = 0;
  for (const auto ringIdx : ringIndices) {
    PRECONDITION(ringIdx < source.numRings(), "ring index out of bounds");
    atomCount +=
        source.d_atomRingBegins[ringIdx + 1] - source.d_atomRingBegins[ringIdx];
    bondCount +=
        source.d_bondRingBegins[ringIdx + 1] - source.d_bondRingBegins[ringIdx];
  }
  d_atomsInRings.reserve(d_atomsInRings.size() + atomCount);
  d_bondsInRings.reserve(d_bondsInRings.size() + bondCount);
  d_atomRingBegins.reserve(d_atomRingBegins.size() + ringIndices.size());
  d_bondRingBegins.reserve(d_bondRingBegins.size() + ringIndices.size());
  for (const auto ringIdx : ringIndices) {
    const auto atomRing = source.atomRings()[ringIdx];
    const auto bondRing = source.bondRings()[ringIdx];
    appendRingUnchecked(atomRing, bondRing);
  }
  return finishRingUpdates();
}

void RingInfo::rebuildMemberships() {
  rebuildRingTypeMemberships(d_atomMembershipBegins, d_atomMemberships,
                             d_atomsInRings, d_atomRingBegins);
  rebuildRingTypeMemberships(d_bondMembershipBegins, d_bondMemberships,
                             d_bondsInRings, d_bondRingBegins);
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
  for (unsigned int otherRingIdx = 0; otherRingIdx < n; ++otherRingIdx) {
    if (areRingsFusedUnchecked(d_fusedRings, ringIdx, otherRingIdx)) {
      return true;
    }
  }
  return false;
}

bool RingInfo::areRingsFused(unsigned int ring1Idx, unsigned int ring2Idx) {
  initFusedRings();
  const auto n = numRings();
  PRECONDITION(ring1Idx < n, "ring1Idx out of bounds");
  PRECONDITION(ring2Idx < n, "ring2Idx out of bounds");
  return areRingsFusedUnchecked(d_fusedRings, ring1Idx, ring2Idx);
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
  unsigned int result = 0;
  for (unsigned int otherRingIdx = 0; otherRingIdx < n; ++otherRingIdx) {
    result += areRingsFusedUnchecked(d_fusedRings, ringIdx, otherRingIdx);
  }
  return result;
}

std::vector<unsigned int> RingInfo::fusedRingNeighbors(unsigned int ringIdx) {
  initFusedRings();
  const auto n = numRings();
  PRECONDITION(ringIdx < n, "ringIdx out of bounds");
  std::vector<unsigned int> res;
  res.reserve(numFusedRingNeighbors(ringIdx));
  for (unsigned int i = 0; i < n; ++i) {
    if (areRingsFusedUnchecked(d_fusedRings, ringIdx, i)) {
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
  const auto numRingPairs = n < 2 ? 0 : size_t(n) * (n - 1) / 2;
  d_fusedRings.assign(numRingPairs, false);
  d_numFusedBonds.assign(n, 0);
  for (unsigned int bondIdx = 0; bondIdx + 1 < d_bondMembershipBegins.size();
       ++bondIdx) {
    const auto ringIndices = bondMembersUnchecked(bondIdx);
    // A bond must belong to at least two rings to fuse them.
    if (ringIndices.size() <= 1) {
      continue;
    }
    for (unsigned int i = 0; i < ringIndices.size() - 1; ++i) {
      unsigned int ringIdx1 = ringIndices[i];
      for (unsigned int j = i + 1; j < ringIndices.size(); ++j) {
        unsigned int ringIdx2 = ringIndices[j];
        if (ringIdx1 != ringIdx2) {
          d_fusedRings[fusedRingPairIndex(ringIdx1, ringIdx2)] = true;
        }
      }
    }
    for (const auto ringIdx : ringIndices) {
      ++d_numFusedBonds[ringIdx];
    }
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
