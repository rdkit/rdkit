//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <unordered_map>

#include <boost/container/small_vector.hpp>

#include "RingSystemFilter.h"

#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <RDGeneral/utils.h>

#include "RingSystemFilter.h"

// Enable this for additional checks when adding
// new patterns
// #define DEBUG_NEW_PATTERNS 1

using namespace RDKit;

namespace {

// Lazily initialize the patterns
const std::vector<ROMol> &getPatterns() {
  // Note: the "1" and "-1" values for the "parity" atom property are arbitrary:
  // We just care whether the parities for each pair of atoms are the same or
  // different. They don't represent any specific CW/CCW value, since depending
  // on construction of the mol being matched, they can be one or the other as
  // long as the same/opposite relationships are preserved.
  static std::vector<ROMol> patterns = {

      // Norbornane
      //      0 - 1
      //    /       \
      //   5 -  6  - 2
      //    \       /
      //      4 - 3

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]1[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:2.parity.-1:5.parity.1|)SMARTS"_smarts),

      // Adamantane
      //     0  -   1
      //           /
      //   /      8    \
      //         /
      //  5 -6- 7        2
      //         \
      //   \      9    /
      //           \
      //     4  -   3

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]1[C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]3[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]1[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]([#6-0,#7-0,#7+1,#8-0,#16-0]2)[#6-0,#7-0,#7+1,#8-0,#16-0]3)SMARTS"
           R"SMARTS(
|atomProp:1.parity.1:3.parity.1:5.parity.-1:7.parity.1|)SMARTS"_smarts),

      // C5_O
      //     0 - 1
      //     |     \
      //     | 3  - 2
      //     |  \  /
      //     5 - 4

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~1)SMARTS"
           R"SMARTS( |atomProp:2.parity.-1:4.parity.-1|)SMARTS"_smarts),

      // 331-bicyclononane
      //       0
      //    /     \
      //  7         1
      //  |         |
      //  6 -  8  - 2
      //  |         |
      //  5         3
      //    \     /
      //       4

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]([#6-0,#7-0,#7+1,#8-0,#16-0]~1)[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:2.parity.1:6.parity.1|)SMARTS"_smarts),

      // 321-bicyclooctane
      //       0
      //    /     \
      //  6         1
      //  |         |
      //  5 -  7  - 2
      //  |         |
      //  4 -  -  - 3

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]([#6-0,#7-0,#7+1,#8-0,#16-0]~1)[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:2.parity.-1:5.parity.-1|)SMARTS"_smarts),

      // 14-bicyclohept-1,4-dione
      //     0 - 1
      //   /       \
      //  6    4    2
      //   \  / \  /
      //     5 - 3

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]-1~[#6-0,#7-0,#7+1,#8-0,#16-0]-[C-0X3]-[C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2-[C-0X3]-1)SMARTS"
           R"SMARTS( |atomProp:3.parity.1:5.parity.1|)SMARTS"_smarts),

      // misc_1
      //       N1
      //     /   \
      //  =C0     C2=
      //    \     /
      //    C9 - C3
      //   /       \
      // C6 -  5  - C4
      //   \       /
      //     7 - 8

      std::move(
          *R"SMARTS([C-0X3]-1-[N-0X3]-[C-0X3]-[C-0X4]-2-[C-0X4]-3-[#6-0,#7-0,#7+1,#8-0,#16-0]-[C-0X4](-[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]-3)-[C-0X4]-1-2)SMARTS"
           R"SMARTS( |atomProp:3.parity.1:9.parity.1|)SMARTS"_smarts),

      // 2,2,2-bicyclooctane
      //     0 - 1
      //   /       \
      //  5 -6 - 7- 2
      //   \       /
      //     4 - 3
      //

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]1[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:2.parity.1:5.parity.-1|)SMARTS"_smarts),

      // 2,3,3-bicyclodecane
      //       0
      //      / \
      //     7   1
      //   /       \
      //  6 -8 -9 - 2
      //   \       /
      //     5   3
      //      \ /
      //       4

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]([#6-0,#7-0,#7+1,#8-0,#16-0]1)[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:2.parity.-1:6.parity.1|)SMARTS"_smarts),

      // 2,2,3-bicyclononane
      //       0
      //      / \
      //     6   1
      //   /       \
      //  5 -7 - 8- 2
      //   \       /
      //     4 - 3

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]([#6-0,#7-0,#7+1,#8-0,#16-0]~1)[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:2.parity.1:5.parity.1|)SMARTS"_smarts),

      // 1,1,3-bicycloheptane
      //     0 - 1
      //   /       \
      //  5    5  - 2
      //   \  /    /
      //     4 - 3

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]([#6-0,#7-0,#7+1,#8-0,#16-0]~1)[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:2.parity.-1:4.parity.-1|)SMARTS"_smarts),

      // misc_2
      //       N8
      //     /   \
      //  =C7     C9=
      //    \     /
      //    C6 - C10
      //   /       \
      // C5- 4 - 3 -C2
      //   \       /
      //     0 - 1

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]1[C-0X4,N+1X4]1[C-0X3,N-OX3]~[C-OX3,N-0X3]~[C-0X3,N-OX3][C-0X4,N+1X4]12)SMARTS"
           R"SMARTS( |atomProp:6.parity.-1:10.parity.-1|)SMARTS"_smarts),

      // misc_3
      //       1
      //     / | \
      //    0  |  2
      //    \  |  /
      //    11 - 3
      //   /   |   \
      // 10    |    4
      //   \   |   /
      //     9 - 5
      //    /  |  \
      //   8   |   6
      //    \  |  /
      //       7

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]1[C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]3[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]4[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]4[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]13)SMARTS"
           R"SMARTS(
|atomProp:1.parity.1:3.parity.-1:5.parity.-1:7.parity.1:9.parity.-1:11.parity.-1|)SMARTS"_smarts),

      // misc_4
      //       4
      //    / /   \
      //  5  9 - 8  3
      //  |     /   |
      //  6 -  7  - 2
      //  |         |
      //  0 -  -  - 1

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]1~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]3[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4,N+1X4]1[C-0X4,N+1X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]3)SMARTS"
           R"SMARTS(
|atomProp:2.parity.1:4.parity.1:6.parity.-1:7.parity.-1|)SMARTS"_smarts),

      // fused_5-5_membered_rings

      std::move(*R"SMARTS([*]~1~[*]~[C-0X4]-2~[*]~[*]~[*]~[C-0X4]-2[*]~1)SMARTS"
                 R"SMARTS( |atomProp:2.parity.1:6.parity.1|)SMARTS"_smarts),

      // cyclopropacyclohexane
      //     0 - 1
      //   /       \
      //  6    4    2
      //   \  / \  /
      //     5 - 3

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]~[C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~1)SMARTS"
           R"SMARTS( |atomProp:3.parity.1:5.parity.1|)SMARTS"_smarts),

      // cyclopropacyclooctane
      //     0 - 1
      //   /       \
      //  8         2
      //  |         |
      //  7    5    3
      //   \  / \  /
      //     6 - 4

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0]~[C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~1)SMARTS"
           R"SMARTS( |atomProp:5.parity.1:7.parity.1|)SMARTS"_smarts),

      // 421-bicyclononane
      //    0 - 1
      //   /     \
      //  7       2
      //  |       |
      //  6 - 8 - 3
      //  |       |
      //  5   -   4

      std::move(
          *R"SMARTS([#6-0,#7-0,#7+1,#8-0,#16-0]~1~[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]2[#6-0,#7-0,#7+1,#8-0,#16-0]~[#6-0,#7-0,#7+1,#8-0,#16-0][C-0X4]([#6-0,#7-0,#7+1,#8-0,#16-0]~1)[#6-0,#7-0,#7+1,#8-0,#16-0]2)SMARTS"
           R"SMARTS( |atomProp:3.parity.1:6.parity.1|)SMARTS"_smarts),

  };
  return patterns;
}

// Check whether an atom's parity is flipped between the pattern and the actual
// molecule
bool isParitySwapped(const RDKit::ROMol &mol, const RDKit::ROMol &pattern,
                     const std::vector<int> &forwardMapping,
                     const std::vector<int> &reverseMapping,
                     std::unordered_map<int, bool> &swappedParitiesCache,
                     int patternAtomIdx) {
  using BondIdxVector = boost::container::small_vector<int, 4>;

  // Check the cache first
  if (auto cached = swappedParitiesCache.find(patternAtomIdx);
      cached != swappedParitiesCache.end()) {
    return cached->second;
  }

  // Build the pattern reference
  BondIdxVector patternBndIndices;
  auto patternAtom = pattern.getAtomWithIdx(patternAtomIdx);
  for (const auto &bnd : pattern.atomBonds(patternAtom)) {
    patternBndIndices.push_back(bnd->getIdx());
  }

  // the relevant (pseudo)chiral atoms in the query have to be
  // specified enough in the pattern to make sense of the parity.
  // This means that we should, at most, have one unmapped atom,
  // which should be ok.
  if (patternBndIndices.size() != 4) {
    patternBndIndices.push_back(-1);
  }

#ifdef DEBUG_NEW_PATTERNS
  if (patternBndIndices.size() != 4) {
    throw std::logic_error("unexpected number of neighbors");
  }
#endif

  // Now build the mol reference
  BondIdxVector molBndIndices;
  auto molAtomIdx = forwardMapping[patternAtomIdx];
  auto molAtom = mol.getAtomWithIdx(molAtomIdx);
  for (const auto &nbr : mol.atomNeighbors(molAtom)) {
    auto ptnNbrIdx = reverseMapping[nbr->getIdx()];
    if (ptnNbrIdx == -1) {
      molBndIndices.push_back(-1);
    } else {
      auto bnd = pattern.getBondBetweenAtoms(patternAtomIdx, ptnNbrIdx);
      if (bnd == nullptr) {
        molBndIndices.push_back(-1);
      } else {
        molBndIndices.push_back(bnd->getIdx());
      }
    }
  }

  if (molBndIndices.size() != 4) {
    molBndIndices.push_back(-1);
  }

#ifdef DEBUG_NEW_PATTERNS
  if (molBndIndices.size() != 4) {
    throw std::logic_error("unexpected number of neighbors");
  }
#endif

  int nSwaps =
      RDKit::countSwapsToInterconvert(patternBndIndices, molBndIndices);

  bool isSwapped = nSwaps % 2;
  swappedParitiesCache[patternAtomIdx] = isSwapped;

  return isSwapped;
}
}  // namespace

void getRingPatternsParityRelations(
    const RDKit::ROMol &mol,
    std::set<std::pair<unsigned int, unsigned int>> &same,
    std::set<std::pair<unsigned int, unsigned int>> &opposite) {
  // Do not do the extra work of deduplicating the matches, since
  // we already deduplicate the parities
  static SubstructMatchParameters p;
  p.uniquify = false;

  same.clear();
  opposite.clear();

  std::vector<int> forwardMapping;
  std::vector<int> reverseMapping;
  for (const auto &pattern : getPatterns()) {
    // The SubstructMatch is limited to 1000 matches, so we might miss
    // some matches. If this happens, and we miss some parity paris,
    // we might see some stereoisomers that would have been filtered
    // out, but no isomers that shouldn't be discarded will be lost
    auto matches = SubstructMatch(mol, pattern, p);
    if (matches.empty()) {
      continue;
    }

    forwardMapping.resize(pattern.getNumAtoms());  // pattern -> mol
    reverseMapping.resize(mol.getNumAtoms());      // mol -> pattern

    for (const auto &match : matches) {
      // The match is a substructure, so it is not unlikely
      // that the reverse map is left with some -1 values
      std::ranges::fill(reverseMapping, -1);

      // No need to reset forwardMapping: all the pattern atoms
      // must be matched, so the mapping will be completely overwritten
      for (auto [qIdx, molIdx] : match) {
        forwardMapping[qIdx] = molIdx;
        reverseMapping[molIdx] = qIdx;
      }

      // For patterns with more than 2 parities, these will
      // be visited multiple times, so cache them here.
      std::unordered_map<int, bool> swappedParitiesCache;

      int iParity = 0;
      int jParity = 0;
      for (auto i = 0u; i + 1 < pattern.getNumAtoms(); ++i) {
        auto iQAtom = pattern.getAtomWithIdx(i);
        if (!iQAtom->getPropIfPresent("parity", iParity)) {
          continue;
        }

        for (auto j = i + 1; j < pattern.getNumAtoms(); ++j) {
          auto jQAtom = pattern.getAtomWithIdx(j);
          if (!jQAtom->getPropIfPresent("parity", jParity)) {
            continue;
          }

          bool paritiesSame = iParity == jParity;
          if (isParitySwapped(mol, pattern, forwardMapping, reverseMapping,
                              swappedParitiesCache, i)) {
            paritiesSame = !paritiesSame;
          }
          if (isParitySwapped(mol, pattern, forwardMapping, reverseMapping,
                              swappedParitiesCache, j)) {
            paritiesSame = !paritiesSame;
          }

          // ensure the lowest index is the first element of the pair:
          // this helps avoid duplicates in the sets
          auto molIdxI = forwardMapping[i];
          auto molIdxJ = forwardMapping[j];
          auto pair = (molIdxI < molIdxJ ? std::make_pair(molIdxI, molIdxJ)
                                         : std::make_pair(molIdxJ, molIdxI));
          if (paritiesSame) {
            same.insert(std::move(pair));
          } else {
            opposite.insert(std::move(pair));
          }
        }
      }
    }
  }
}
