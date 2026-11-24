//
// Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDKitBase.h"
#include "Bridgeheads.h"

#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>

namespace RDKit {
namespace MolOps {

namespace details {
// returns indices for all pairs of rings that share at least two bonds along
// with the shared bond indices
std::vector<
    std::pair<std::pair<unsigned int, unsigned int>, std::vector<unsigned int>>>
findRingPairs(const ROMol &mol, unsigned int minCommonBonds) {
  const auto &ringInfo = mol.getRingInfo();
  const auto &rings = ringInfo->bondRings();
  std::vector<boost::dynamic_bitset<>> ringMasks(rings.size());
  for (auto i = 0u; i < rings.size(); ++i) {
    ringMasks[i].resize(mol.getNumBonds());
    for (auto b : rings[i]) {
      ringMasks[i].set(b);
    }
  }
  std::vector<std::pair<std::pair<unsigned int, unsigned int>,
                        std::vector<unsigned int>>>
      ringPairs;
  for (size_t i = 0; i < rings.size(); ++i) {
    for (size_t j = i + 1; j < rings.size(); ++j) {
      auto commonBonds = ringMasks[i] & ringMasks[j];
      if (commonBonds.count() < minCommonBonds) {
        continue;
      }
      std::vector<unsigned int> sharedBonds;
      size_t pos = commonBonds.find_first();
      while (pos != boost::dynamic_bitset<>::npos) {
        sharedBonds.push_back(pos);
        pos = commonBonds.find_next(pos);
      }
      ringPairs.emplace_back(std::make_pair(i, j), sharedBonds);
    }
  }
  return ringPairs;
}
}  // namespace details

std::vector<Bridgehead> findBridgeheads(const ROMol &mol) {
  // FIX: verify that SSSR is actually sufficient for bridgehead detection
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(mol);
  }

  // find pairs of rings that share at least two bonds:
  const auto ringPairs = details::findRingPairs(mol, 2);

  // mark all bonds that are part of the shared rings
  const auto ringinfo = mol.getRingInfo();
  boost::dynamic_bitset<> bondsAreInSharedRings(mol.getNumBonds());
  for (const auto &[pr, sharedBonds] : ringPairs) {
    for (const auto pri : {pr.first, pr.second}) {
      for (const auto bi : ringinfo->bondRings()[pri]) {
        bondsAreInSharedRings.set(bi);
      }
    }
  }

  std::unordered_map<unsigned int, unsigned int> atomSharedRingBondCount;
  // lambda to count the number of shared ring bonds for at atom
  auto getAtomSharedRingBondCount = [&atomSharedRingBondCount,
                                     &bondsAreInSharedRings,
                                     &mol](unsigned int ai) {
    if (atomSharedRingBondCount.find(ai) == atomSharedRingBondCount.end()) {
      unsigned int ringBondCount = 0;
      for (const auto bnd : mol.atomBonds(mol.getAtomWithIdx(ai))) {
        if (bondsAreInSharedRings.test(bnd->getIdx())) {
          ++ringBondCount;
        }
      }
      atomSharedRingBondCount[ai] = ringBondCount;
    }
    return atomSharedRingBondCount[ai];
  };

  const auto atomRings = mol.getRingInfo()->atomRings();
  std::vector<Bridgehead> bridgeheads;
  for (const auto &[pr, sharedBonds] : ringPairs) {
    const auto &[i, j] = pr;
    const auto &ringi = atomRings[i];
    const auto &ringj = atomRings[j];
    // Find atoms with at least three ring bonds that are common to both rings.
    // We know that these can only be involved in the shared bonds, so we can
    // limit the atoms we consider to the atoms involved in those bonds.
    std::set<unsigned int> commonAtoms;
    for (const auto bi : sharedBonds) {
      auto bondi = mol.getBondWithIdx(bi);
      for (const auto ai : {bondi->getBeginAtomIdx(), bondi->getEndAtomIdx()}) {
        if (getAtomSharedRingBondCount(ai) < 3) {
          continue;
        }
        if (std::ranges::find(ringi, ai) != ringi.end() &&
            std::ranges::find(ringj, ai) != ringj.end()) {
          commonAtoms.insert(ai);
        }
      }
    }

    if (commonAtoms.size() == 2) {
      auto it = commonAtoms.begin();
      unsigned int atom1 = *it;
      ++it;
      unsigned int atom2 = *it;
      bridgeheads.emplace_back(
          atom1, atom2, std::vector<unsigned int>(ringi.begin(), ringi.end()),
          std::vector<unsigned int>(ringj.begin(), ringj.end()));
    }
  }

  // Implementation for finding bridgeheads goes here
  return bridgeheads;
}

std::vector<unsigned int> findBridgeheadAtoms(const ROMol &mol) {
  auto bridgeheads = findBridgeheads(mol);
  std::set<unsigned int> bridgeheadAtoms;
  for (const auto &bh : bridgeheads) {
    bridgeheadAtoms.insert(bh.atom1);
    bridgeheadAtoms.insert(bh.atom2);
  }
  std::vector<unsigned int> result(bridgeheadAtoms.begin(),
                                   bridgeheadAtoms.end());
  // Implementation for finding bridgehead atoms goes here
  return result;
}

}  // namespace MolOps
}  // namespace RDKit

std::ostream &operator<<(std::ostream &target,
                         const RDKit::MolOps::Bridgehead &bh) {
  target << "Bridgehead(atom1=" << bh.atom1 << ", atom2=" << bh.atom2
         << ", ring1=[";
  std::copy(bh.ring1.begin(), bh.ring1.end(),
            std::ostream_iterator<unsigned int>(target, ", "));
  target << "], ring2=[";
  std::copy(bh.ring2.begin(), bh.ring2.end(),
            std::ostream_iterator<unsigned int>(target, ", "));
  target << "])";
  return target;
}