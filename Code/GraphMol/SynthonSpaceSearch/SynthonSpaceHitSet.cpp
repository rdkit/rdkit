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
    const std::vector<std::unique_ptr<ROMol>> &fragSet)
    : d_reaction(&reaction) {
  synthonsToUse.reserve(stu.size());
  const auto &synthons = reaction.getSynthons();
  for (size_t i = 0; i < stu.size(); ++i) {
    synthonsToUse.emplace_back();
    synthonsToUse[i].reserve(stu[i].size());
    for (size_t j = 0; j < stu[i].size(); ++j) {
      synthonsToUse[i].emplace_back(
          std::make_pair(synthons[i][stu[i][j]].first,
                         synthons[i][stu[i][j]].second->getOrigMol().get()));
    }
  }

  frags.reserve(fragSet.size());
  std::transform(
      fragSet.begin(), fragSet.end(), std::back_inserter(frags),
      [](const std::unique_ptr<ROMol> &f) -> ROMol * { return f.get(); });

  numHits = std::accumulate(
      stu.begin(), stu.end(), size_t(1),
      [](const int prevRes, const std::vector<size_t> &s2) -> size_t {
        return prevRes * s2.size();
      });
}

SynthonSpaceFPHitSet::SynthonSpaceFPHitSet(
    const SynthonSet &reaction, const std::vector<std::vector<size_t>> &stu,
    const std::vector<std::unique_ptr<ROMol>> &fragSet)
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

}  // namespace RDKit::SynthonSpaceSearch
