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

namespace RDKit {
namespace MolOps {

namespace details {
// returns indices for all pairs of rings that share at least two bonds
std::vector<std::pair<unsigned int, unsigned int>> findRingPairs(
    const ROMol &mol, unsigned int minCommonBonds) {
  const auto &ringInfo = mol.getRingInfo();
  const auto &rings = ringInfo->bondRings();
  std::vector<boost::dynamic_bitset<>> ringMasks(rings.size());
  for (auto i = 0u; i < rings.size(); ++i) {
    ringMasks[i].resize(mol.getNumBonds());
    for (auto b : rings[i]) {
      ringMasks[i].set(b);
    }
  }
  std::vector<std::pair<unsigned int, unsigned int>> ringPairs;
  for (size_t i = 0; i < rings.size(); ++i) {
    for (size_t j = i + 1; j < rings.size(); ++j) {
      auto commonBonds = ringMasks[i] & ringMasks[j];
      if (commonBonds.count() >= minCommonBonds) {
        ringPairs.emplace_back(i, j);
      }
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
  std::vector<Bridgehead> bridgeheads;
  bridgeheads.resize(mol.getNumAtoms());
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