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
#include <algorithm>

namespace RDKit {

RingInfo::INT_VECT RingInfo::atomRingSizes(unsigned int idx) const {
  checkInitialized();

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
  checkInitialized();

  if (idx < d_atomMembers.size()) {
    return std::find_if(d_atomMembers[idx].begin(), d_atomMembers[idx].end(),
                        [this, size](int ri) {
                          return d_atomRings.at(ri).size() == size;
                        }) != d_atomMembers[idx].end();
  }
  return false;
}
unsigned int RingInfo::minAtomRingSize(unsigned int idx) const {
  checkInitialized();

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
  checkInitialized();

  if (idx < d_atomMembers.size()) {
    return rdcast<unsigned int>(d_atomMembers[idx].size());
  }
  return 0;
}
RingInfo::INT_VECT RingInfo::atomMembers(unsigned int idx) const {
  checkInitialized();

  static const INT_VECT emptyVect;
  if (idx < d_atomMembers.size()) {
    return d_atomMembers[idx];
  }
  return INT_VECT();
}
bool RingInfo::areAtomsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  checkInitialized();

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
  checkInitialized();

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
  checkInitialized();

  if (idx < d_bondMembers.size()) {
    return std::find_if(d_bondMembers[idx].begin(), d_bondMembers[idx].end(),
                        [this, size](int ri) {
                          return d_bondRings.at(ri).size() == size;
                        }) != d_bondMembers[idx].end();
  }
  return false;
}
unsigned int RingInfo::minBondRingSize(unsigned int idx) const {
  checkInitialized();

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
  checkInitialized();

  if (idx < d_bondMembers.size()) {
    return rdcast<unsigned int>(d_bondMembers[idx].size());
  }
  return 0;
}
RingInfo::INT_VECT RingInfo::bondMembers(unsigned int idx) const {
  checkInitialized();

  if (idx < d_bondMembers.size()) {
    return d_bondMembers[idx];
  }
  return INT_VECT();
}
bool RingInfo::areBondsInSameRingOfSize(unsigned int idx1, unsigned int idx2,
                                        unsigned int size) const {
  checkInitialized();

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
  checkInitialized();
  PRECONDITION(d_atomRings.size() == d_bondRings.size(), "length mismatch");
  return rdcast<unsigned int>(d_atomRings.size());
}

unsigned int RingInfo::addRing(const INT_VECT &atomIndices,
                               const INT_VECT &bondIndices) {
  checkInitialized();
  PRECONDITION(atomIndices.size() == bondIndices.size(), "length mismatch");
  for (auto i : atomIndices) {
    if (i >= static_cast<int>(d_atomMembers.size())) {
      d_atomMembers.resize(i + 1);
    }
    d_atomMembers[i].push_back(d_atomRings.size());
  }
  for (auto i : bondIndices) {
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
  checkInitialized();
  PRECONDITION(ringIdx < d_bondRings.size(), "ringIdx out of bounds");
  initFusedRingInfo();
  return d_fusedRingInfo.isRingFused(ringIdx);
}

bool RingInfo::areRingsFused(unsigned int ring1Idx, unsigned int ring2Idx) {
  checkInitialized();
  PRECONDITION(ring1Idx < d_bondRings.size(), "ring1Idx out of bounds");
  PRECONDITION(ring2Idx < d_bondRings.size(), "ring2Idx out of bounds");
  initFusedRingInfo();
  return d_fusedRingInfo.areRingsFused(ring1Idx, ring2Idx);
}

RingInfo::INT_VECT RingInfo::fusedRings(unsigned int ringIdx) {
  checkInitialized();
  PRECONDITION(ringIdx < d_bondRings.size(), "ringIdx out of bounds");
  initFusedRingInfo();
  return d_fusedRingInfo.fusedRings(ringIdx);
}

unsigned int RingInfo::numFusedBonds(unsigned int ringIdx) {
  checkInitialized();
  PRECONDITION(ringIdx < d_bondRings.size(), "ringIdx out of bounds");
  initFusedRingInfo();
  return d_fusedRingInfo.numFusedBonds(ringIdx);
}

bool RingInfo::hasRingFusionInfoForBond(unsigned int idx) {
  checkInitialized();
  if (idx < d_bondMembers.size()) {
    initFusedRingInfo();
    return d_fusedRingInfo.hasRingFusionInfoForBond(idx);
  }
  return true;
}

bool RingInfo::isBondInFusedRingOfSize(unsigned int idx, unsigned int size) {
  checkInitialized();
  if (idx < d_bondMembers.size()) {
    initFusedRingInfo();
    const auto &ringSizes = d_fusedRingInfo.ringSizesForBond(idx);
    return ringSizes.test(size);
  }
  return false;
}

RingInfo::INT_VECT RingInfo::bondFusedRingSizes(unsigned int idx) {
  checkInitialized();
  INT_VECT res;
  if (idx < d_bondMembers.size()) {
    initFusedRingInfo();
    const auto &ringSizes = d_fusedRingInfo.ringSizesForBond(idx);
    for (unsigned int i = 0; i < ringSizes.size(); ++i) {
      if (ringSizes.test(i)) {
        res.push_back(i);
      }
    }
  }
  return res;
}

boost::dynamic_bitset<> RingInfo::bondFusedRingSizesAsBitset(unsigned int idx) {
  checkInitialized();
  boost::dynamic_bitset<> res;
  if (idx < d_bondMembers.size()) {
    initFusedRingInfo();
    res = d_fusedRingInfo.ringSizesForBond(idx);
  }
  return res;
}

RingInfo::RingNbrPermutations::RingNbrPermutations(
    const FusedRingInfo *fusedRingInfo, unsigned int ringIdx)
    : d_fusedRingInfo(fusedRingInfo) {
  constexpr unsigned int MAX_FUSED_NBRS =
      std::numeric_limits<perm_type>::digits;
  const auto &fusedRingSystem = fusedRingInfo->d_fusedRingSystems.at(ringIdx);
  unsigned int numFusedRings = fusedRingSystem.count();
  if (numFusedRings > MAX_FUSED_NBRS) {
    return;
  }
  d_toLocalIndex.resize(fusedRingSystem.size(), -1);
  d_fromLocalIndex.resize(numFusedRings, -1);
  unsigned int j = 0;
  for (unsigned int i = 0; i < fusedRingSystem.size(); ++i) {
    if (fusedRingSystem.test(i)) {
      d_toLocalIndex[i] = j;
      d_fromLocalIndex[j++] = i;
    }
  }
  boost::dynamic_bitset<> onStack(1 << numFusedRings);
  std::stack<perm_type> todo;
  auto localRingIdx = d_toLocalIndex.at(ringIdx);
  CHECK_INVARIANT(localRingIdx != -1, "bad localRingIdx");
  auto perm = static_cast<perm_type>(1 << localRingIdx);
  onStack.set(perm);
  todo.push(perm);
  d_permutations.push_back(perm);
  while (!todo.empty()) {
    perm = todo.top();
    todo.pop();
    localRingIdx = 0;
    auto permTmp = perm;
    while (permTmp) {
      if (permTmp & 1) {
        auto ringIdx = d_fromLocalIndex.at(localRingIdx);
        CHECK_INVARIANT(ringIdx != -1, "bad ringIdx");
        const auto &fusedRings = fusedRingInfo->d_fusedRings.at(ringIdx);
        for (unsigned int ringNbrIdx = 0; ringNbrIdx < fusedRings.size();
             ++ringNbrIdx) {
          if (fusedRings.test(ringNbrIdx)) {
            auto localRingNbrIdx = d_toLocalIndex.at(ringNbrIdx);
            CHECK_INVARIANT(localRingNbrIdx != -1, "bad localRingNbrIdx");
            auto newPerm = perm | static_cast<perm_type>(1 << localRingNbrIdx);
            if (!onStack.test(newPerm)) {
              onStack.set(newPerm);
              todo.push(newPerm);
              d_permutations.push_back(newPerm);
            }
          }
        }
      }
      ++localRingIdx;
      permTmp >>= 1;
    }
  }
}

void RingInfo::FusedRingInfo::init(RingInfo *ringInfo) {
  d_ringInfo = ringInfo;
  const auto &bondRings = d_ringInfo->d_bondRings;
  auto numRings = bondRings.size();
  d_fusedRings.resize(numRings);
  for (auto &fusedRing : d_fusedRings) {
    fusedRing.resize(numRings);
  }
  for (unsigned int bi = 0; bi < d_ringInfo->d_bondMembers.size(); ++bi) {
    if (d_ringInfo->numBondRings(bi) <= 1) {
      continue;
    }
    const auto &ringIndices = d_ringInfo->d_bondMembers.at(bi);
    for (unsigned int i = 0; i < ringIndices.size() - 1; ++i) {
      unsigned int ringIdx1 = ringIndices[i];
      for (unsigned int j = i + 1; j < ringIndices.size(); ++j) {
        unsigned int ringIdx2 = ringIndices[j];
        d_fusedRings[ringIdx1].set(ringIdx2);
        d_fusedRings[ringIdx2].set(ringIdx1);
      }
    }
  }
  d_numFusedBonds.resize(numRings, 0);
  for (unsigned int ri = 0; ri < numRings; ++ri) {
    d_numFusedBonds[ri] += std::count_if(
        bondRings[ri].begin(), bondRings[ri].end(),
        [this](unsigned int bi) { return d_ringInfo->numBondRings(bi) > 1; });
  }
}

void RingInfo::FusedRingInfo::initFusedRingSystems() {
  const auto &bondRings = d_ringInfo->d_bondRings;
  auto numRings = bondRings.size();
  auto maxBondIdx = d_ringInfo->d_bondMembers.size();
  d_fusedRingSystems.resize(numRings);
  boost::dynamic_bitset<> visited(numRings);
  for (unsigned int i = 0; i < numRings; ++i) {
    visited.reset();
    auto &fusedRingSystem = d_fusedRingSystems[i];
    fusedRingSystem.resize(numRings);
    addNbrRings(fusedRingSystem, visited, i);
  }
  d_bondRingsAsBitset.resize(numRings);
  for (unsigned int i = 0; i < numRings; ++i) {
    auto &bondRingAsBitset = d_bondRingsAsBitset[i];
    bondRingAsBitset.resize(maxBondIdx);
    for (unsigned int bi : bondRings.at(i)) {
      bondRingAsBitset.set(bi);
    }
  }
  d_fusedRingSizesVec.reserve(numRings);
  for (unsigned int i = 0; i < numRings; ++i) {
    d_fusedRingSizesVec.emplace_back(this, i);
  }
  d_bondFusedRingSizes.resize(maxBondIdx);
}

const boost::dynamic_bitset<> &RingInfo::FusedRingInfo::ringSizesForBond(
    unsigned int idx) {
  checkInitialized();
  if (!hasRingFusionInfoForBond(idx)) {
    throw FusedSystemTooLarge(
        "The fused system the bond belongs to is too large");
  }
  auto &bondFusedRingSizes = d_bondFusedRingSizes.at(idx);
  if (bondFusedRingSizes.empty()) {
    bondFusedRingSizes.resize(d_ringInfo->d_bondMembers.size());
    const auto &rings = d_ringInfo->d_bondMembers.at(idx);
    auto ringIdxMask = getRingIdxMask(rings);
    for (auto ringIdx : rings) {
      auto ringIdxMaskTmp = ringIdxMask;
      ringIdxMaskTmp.reset(ringIdx);
      auto &fusedRingSizes = d_fusedRingSizesVec[ringIdx];
      fusedRingSizes.setMask(ringIdxMaskTmp);
      fusedRingSizes.reset();
      while (1) {
        auto ringSize = fusedRingSizes.next();
        if (ringSize == -1) {
          break;
        }
        bondFusedRingSizes.set(ringSize);
      }
    }
  }
  return bondFusedRingSizes;
}

void RingInfo::FusedRingInfo::addNbrRings(boost::dynamic_bitset<> &ringNbrs,
                                          boost::dynamic_bitset<> &visited,
                                          unsigned int ringIdx) {
  if (!visited.test(ringIdx)) {
    visited.set(ringIdx);
    ringNbrs.set(ringIdx);
    const auto &fusedRings = d_fusedRings.at(ringIdx);
    for (unsigned int i = 0; i < fusedRings.size(); ++i) {
      if (fusedRings.test(i)) {
        addNbrRings(ringNbrs, visited, i);
      }
    }
  }
}
boost::dynamic_bitset<> RingInfo::FusedRingInfo::getRingIdxMask(
    const RingInfo::INT_VECT &bondRings) {
  boost::dynamic_bitset<> ringIdxMask(d_fusedRingSystems.size());
  std::for_each(
      bondRings.begin(), bondRings.end(),
      [&ringIdxMask](const auto ringIdx) { ringIdxMask.set(ringIdx); });
  return ringIdxMask;
}

void RingInfo::RingNbrPermutations::setMask(
    const boost::dynamic_bitset<> &mask) {
  d_mask = 0;
  for (unsigned int i = 0; i < mask.size(); ++i) {
    if (mask.test(i) && d_toLocalIndex.at(i) != -1) {
      d_mask |= (1 << d_toLocalIndex[i]);
    }
  }
}

int RingInfo::RingNbrPermutations::next() {
  if (d_permutationIdx == d_permutations.size()) {
    return -1;
  }
  auto perm = d_permutations[d_permutationIdx++];
  if (perm & d_mask) {
    return next();
  }
  unsigned int localRingIdx = 0;
  boost::dynamic_bitset<> envelope(
      d_fusedRingInfo->d_ringInfo->d_bondMembers.size());
  while (perm) {
    if (perm & 1) {
      auto ringIdx = d_fromLocalIndex.at(localRingIdx);
      CHECK_INVARIANT(ringIdx != -1, "bad localRingIdx");
      envelope ^= d_fusedRingInfo->d_bondRingsAsBitset.at(ringIdx);
    }
    ++localRingIdx;
    perm >>= 1;
  }
  return static_cast<int>(envelope.count());
}

#ifdef RDK_USE_URF
unsigned int RingInfo::numRingFamilies() const {
  checkInitialized();
  return d_atomRingFamilies.size();
};

unsigned int RingInfo::numRelevantCycles() const {
  checkInitialized();
  return rdcast<unsigned int>(RDL_getNofRC(dp_urfData.get()));
};

unsigned int RingInfo::addRingFamily(const INT_VECT &atomIndices,
                                     const INT_VECT &bondIndices) {
  checkInitialized();
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
