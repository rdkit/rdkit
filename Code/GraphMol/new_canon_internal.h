//
//  Copyright (C) 2026 Kevin Boyd
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Internal helpers shared by the GraphMol canonicalization implementation.
// This header is intentionally not installed.

#ifndef RD_NEW_CANON_INTERNAL_H
#define RD_NEW_CANON_INTERNAL_H

#include "new_canon.h"

#include <compare>
#include <cstdint>
#include <vector>

namespace RDKit {
namespace Canon {
namespace detail {

struct StereoAtomCompareCode {
  std::uint64_t primary{0};
  std::uint64_t secondary{0};

  // The defaulted comparison orders primary first, then secondary.
  auto operator<=>(const StereoAtomCompareCode &) const = default;
};

class StereoAtomCompareFunctor {
 public:
  Canon::canon_atom *dp_atoms{nullptr};
  const ROMol *dp_mol{nullptr};
  const StereoAtomCompareCode *dp_compareCodes{nullptr};
  const boost::dynamic_bitset<> *dp_atomsInPlay{nullptr};
  const boost::dynamic_bitset<> *dp_bondsInPlay{nullptr};
  bool df_useNbrs{false};

  StereoAtomCompareFunctor(
      Canon::canon_atom *atoms, const ROMol &mol,
      const std::vector<StereoAtomCompareCode> &compareCodes,
      const boost::dynamic_bitset<> *atomsInPlay,
      const boost::dynamic_bitset<> *bondsInPlay)
      : dp_atoms(atoms),
        dp_mol(&mol),
        dp_compareCodes(compareCodes.data()),
        dp_atomsInPlay(atomsInPlay),
        dp_bondsInPlay(bondsInPlay) {}

  int operator()(int i, int j) const {
    if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
      return 0;
    }

    if (dp_atoms[i].index != dp_atoms[j].index) {
      return dp_atoms[i].index < dp_atoms[j].index ? -1 : 1;
    }

    int atomMapI = 0;
    int atomMapJ = 0;
    if (dp_atoms[i].atom->getAtomicNum() == 0) {
      dp_atoms[i].atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                         atomMapI);
    }
    if (dp_atoms[j].atom->getAtomicNum() == 0) {
      dp_atoms[j].atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                         atomMapJ);
    }
    if (atomMapI != atomMapJ) {
      return atomMapI < atomMapJ ? -1 : 1;
    }

    if (dp_atoms[i].degree != dp_atoms[j].degree) {
      return dp_atoms[i].degree < dp_atoms[j].degree ? -1 : 1;
    }

    const auto &lhs = dp_compareCodes[i];
    const auto &rhs = dp_compareCodes[j];
    if (const auto cmp = lhs <=> rhs; cmp != 0) {
      return cmp < 0 ? -1 : 1;
    }

    if (df_useNbrs) {
      if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
        updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
      }
      if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
        updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
      }
      for (unsigned int ii = 0;
           ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size();
           ++ii) {
        const auto cmp =
            bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
        if (cmp) {
          return cmp;
        }
      }
      if (dp_atoms[i].bonds.size() != dp_atoms[j].bonds.size()) {
        return dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size() ? -1 : 1;
      }
    }
    return 0;
  }
};

void rankStereoAtoms(StereoAtomCompareFunctor &ftor, std::vector<int> &order,
                     const boost::dynamic_bitset<> *atomsInPlay,
                     const boost::dynamic_bitset<> *bondsInPlay);

}  // namespace detail
}  // namespace Canon
}  // namespace RDKit

#endif  // RD_NEW_CANON_INTERNAL_H
