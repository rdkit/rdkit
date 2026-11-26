//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "JointRing.h"
#include "RDGeneral/types.h"
#include "ROMol.h"

namespace RDKit {
std::vector<JointRingConstraint> collectJRConstraints(const ROMol &mol) {
  std::vector<JointRingConstraint> jrToAtomIdx;
  auto transfer = [&jrToAtomIdx](const auto &b, const auto value,
                                 const bool exclude = false) {
    for (const auto idx : b) {
      if (idx >= jrToAtomIdx.size()) {
        ASSERT_INVARIANT(idx < UINT_MAX, "Corruped ring indices");
        jrToAtomIdx.resize(idx + 1u);
      }
      if (exclude) {
        jrToAtomIdx[idx].excluded.emplace_back(value);
      } else {
        jrToAtomIdx[idx].joint.emplace_back(value);
      }
    }
  };
  for (const auto *atom : mol.atoms()) {
    if (!atom->hasProp(common_properties::_CXJointRingIDs)) {
      continue;
    }
    const auto &jrIndices = rdvalue_cast<const JointRingIds &>(
        atom->getDict().getRDValue(common_properties::_CXJointRingIDs));
    transfer(jrIndices.joint, atom->getIdx());
    transfer(jrIndices.excluded, atom->getIdx(), true);
  }
  return jrToAtomIdx;
}
}  // namespace RDKit