//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/SynthonSpaceSearch/SynthonSpaceHitSet.h>

namespace RDKit::SynthonSpaceSearch {
SynthonSpaceHitSet::SynthonSpaceHitSet(
    const SynthonSet &reaction, const std::vector<std::vector<size_t>> &stu,
    const std::vector<std::shared_ptr<ROMol>> &fragSet)
    : d_reaction(&reaction) {
  synthonsToUse.reserve(stu.size());
  const auto &synthons = reaction.getSynthons();
  for (size_t i = 0; i < stu.size(); ++i) {
    synthonsToUse.emplace_back();
    synthonsToUse[i].reserve(stu[i].size());
    for (size_t j = 0; j < stu[i].size(); ++j) {
      synthonsToUse[i].emplace_back(synthons[i][stu[i][j]]);
    }
  }

  frags = fragSet;
  numHits = std::accumulate(
      stu.begin(), stu.end(), size_t(1),
      [](const int prevRes, const std::vector<size_t> &s2) -> size_t {
        return prevRes * s2.size();
      });
}

SynthonSpaceFPHitSet::SynthonSpaceFPHitSet(
    const SynthonSet &reaction, const std::vector<std::vector<size_t>> &stu,
    const std::vector<std::shared_ptr<ROMol>> &fragSet)
    : SynthonSpaceHitSet(reaction, stu, fragSet) {
  synthonFPs.reserve(stu.size());
  for (size_t i = 0; i < stu.size(); ++i) {
    synthonFPs.emplace_back();
    synthonFPs[i].reserve(stu[i].size());
    const auto &sfps = reaction.getSynthons();
    for (size_t j = 0; j < stu[i].size(); ++j) {
      synthonFPs[i].emplace_back(sfps[i][stu[i][j]].second->getFP().get());
    }
  }
  addFP = reaction.getAddFP().get();
  subtractFP = reaction.getSubtractFP().get();
}

SynthonSpaceShapeHitSet::SynthonSpaceShapeHitSet(
    const SynthonSet &reaction, const std::vector<std::vector<size_t>> &stu,
    const std::vector<std::shared_ptr<ROMol>> &fragSet,
    const std::vector<SynthonShapeInput *> &fShapes,
    const std::vector<unsigned int> &sSetOrder)
    : SynthonSpaceHitSet(reaction, stu, fragSet),
      fragShapes(fShapes),
      synthonSetOrder(sSetOrder) {
  // it may be that there are fewer entries in synthonSetOrder than
  // there are synthon sets in the SynthonSet.  That occurs if the
  // fragment set was smaller than the SynthonSet.  Pad it out if so.
  for (size_t i = 0; i < reaction.getSynthons().size(); ++i) {
    if (auto it = std::ranges::find(synthonSetOrder, i);
        it == synthonSetOrder.end()) {
      synthonSetOrder.emplace_back(i);
    }
  }
}

}  // namespace RDKit::SynthonSpaceSearch
