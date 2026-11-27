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
findRingPairs(const ROMol &mol, unsigned int minCommonBonds,
              const std::vector<unsigned int> &onlyWithAtoms) {
  const auto &ringInfo = mol.getRingInfo();
  const auto &rings = ringInfo->bondRings();

  std::vector<boost::dynamic_bitset<>> ringMasks(rings.size());
  for (auto i = 0u; i < rings.size(); ++i) {
    ringMasks[i].resize(mol.getNumBonds());
    bool hasAtom = onlyWithAtoms.empty();
    if (!onlyWithAtoms.empty()) {
      hasAtom = std::ranges::any_of(rings[i], [&mol, &onlyWithAtoms](auto b) {
        const auto bnd = mol.getBondWithIdx(b);
        const auto begIdx = bnd->getBeginAtomIdx();
        const auto endIdx = bnd->getEndAtomIdx();
        return std::ranges::any_of(onlyWithAtoms, [&begIdx, &endIdx](auto ai) {
          return begIdx == ai || endIdx == ai;
        });
      });
    }
    if (!hasAtom) {
      continue;
    }
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

namespace {

std::vector<Bridgehead> bridgeheadFinder(
    const ROMol &mol, const std::vector<unsigned int> &onlyWithAtoms) {
  // find pairs of rings that share at least two bonds:
  const auto ringPairs = details::findRingPairs(mol, 2, onlyWithAtoms);

  // mark all bonds that are part of the shared rings
  const auto ringinfo = mol.getRingInfo();
  boost::dynamic_bitset<> bondsAreInSharedRings(mol.getNumBonds());
  const auto bondRings = ringinfo->bondRings();
  for (const auto &[pr, sharedBonds] : ringPairs) {
    for (const auto pri : {pr.first, pr.second}) {
      for (const auto bi : bondRings[pri]) {
        bondsAreInSharedRings.set(bi);
      }
    }
  }

  std::unordered_map<unsigned int, unsigned int> atomSharedRingBondCount;
  // lambda to count the number of shared ring bonds for at atom
  auto getAtomSharedRingBondCount = [&bondsAreInSharedRings, &mol](
                                        unsigned int ai, unsigned int i,
                                        unsigned int j) {
    const auto bondRings = mol.getRingInfo()->bondRings();
    const auto &ringi = bondRings[i];
    const auto &ringj = bondRings[j];
    unsigned int ringBondCount = 0;
    for (const auto bnd : mol.atomBonds(mol.getAtomWithIdx(ai))) {
      const auto bi = bnd->getIdx();
      if (bondsAreInSharedRings.test(bi) &&
          (std::ranges::find(ringi, bi) != ringi.end() ||
           std::ranges::find(ringj, bi) != ringj.end())) {
        ++ringBondCount;
      }
    }
    return ringBondCount;
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
        if (getAtomSharedRingBondCount(ai, i, j) < 3) {
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
      if (!onlyWithAtoms.empty()) {
        if (std::ranges::find(onlyWithAtoms, atom1) == onlyWithAtoms.end() &&
            std::ranges::find(onlyWithAtoms, atom2) == onlyWithAtoms.end()) {
          continue;
        }
      }
      bridgeheads.emplace_back(
          atom1, atom2, std::vector<unsigned int>(ringi.begin(), ringi.end()),
          std::vector<unsigned int>(ringj.begin(), ringj.end()));
    }
  }

  // Implementation for finding bridgeheads goes here
  return bridgeheads;
}
}  // namespace

std::vector<Bridgehead> findBridgeheads(
    ROMol &mol, const std::vector<unsigned int> &onlyWithAtoms) {
  // We don't get a complete list of all of the bridghead pairs unless we have
  // the symmetric set of smallest rings.
  if (!mol.getRingInfo()->isSymmSssr()) {
    MolOps::symmetrizeSSSR(mol);
  }
  return bridgeheadFinder(mol, onlyWithAtoms);
}

std::vector<unsigned int> findBridgeheadAtoms(
    const ROMol &mol, const std::vector<unsigned int> &onlyWithAtoms) {
  // We can find all of the bridgehead atoms with only the SSSR rings, we don't
  // need the symmetric set
  if (!mol.getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(mol);
  }
  auto bridgeheads = bridgeheadFinder(mol, onlyWithAtoms);
  std::set<unsigned int> bridgeheadAtoms;
  for (const auto &bh : bridgeheads) {
    if (!onlyWithAtoms.empty()) {
      if (std::ranges::find(onlyWithAtoms, bh.atom1) != onlyWithAtoms.end()) {
        bridgeheadAtoms.insert(bh.atom1);
      }
      if (std::ranges::find(onlyWithAtoms, bh.atom2) != onlyWithAtoms.end()) {
        bridgeheadAtoms.insert(bh.atom2);
      }
    } else {
      bridgeheadAtoms.insert(bh.atom1);
      bridgeheadAtoms.insert(bh.atom2);
    }
  }
  std::vector<unsigned int> result(bridgeheadAtoms.begin(),
                                   bridgeheadAtoms.end());
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