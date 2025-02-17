//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef SYNTHONSPACEHITSET_H
#define SYNTHONSPACEHITSET_H

#include <numeric>
#include <string>
#include <vector>

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>

namespace RDKit::SynthonSpaceSearch {
// Holds the information about a set of hits.  The molecules can be built
// by making all combinations of synthons, one taken from each synthon set.
// The synthons are copied from the SynthonSet so it is not needed to
// build the hits.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceHitSet {
  SynthonSpaceHitSet() = delete;
  SynthonSpaceHitSet(const SynthonSet &reaction,
                     const std::vector<std::vector<size_t>> &stu,
                     const std::vector<std::shared_ptr<ROMol>> &fragSet)
      : reactionId(reaction.getId()) {
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
    for (size_t i = 0; i < fragSet.size(); ++i) {
      frags.push_back(fragSet[i].get());
    }
    numHits = std::accumulate(
        stu.begin(), stu.end(), size_t(1),
        [](const int prevRes, const std::vector<size_t> &s2) -> size_t {
          return prevRes * s2.size();
        });
  }
  virtual ~SynthonSpaceHitSet() = default;

  std::string reactionId;
  // The fragments that this hitset is derived from, useful for debugging.
  std::vector<const ROMol *> frags;
  std::vector<std::vector<std::pair<std::string, const ROMol *>>> synthonsToUse;
  size_t numHits{0};
};

// This sub-class holds results from a SynthonSpaceFingerprintSearch.
// It needs the synthon fingerprints as well.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceFPHitSet
    : SynthonSpaceHitSet {
  SynthonSpaceFPHitSet() = delete;
  SynthonSpaceFPHitSet(const SynthonSet &reaction,
                       const std::vector<std::vector<size_t>> &stu,
                       const std::vector<std::shared_ptr<ROMol>> &fragSet)
      : SynthonSpaceHitSet(reaction, stu, fragSet) {
    synthonFPs.reserve(stu.size());
    for (size_t i = 0; i < stu.size(); ++i) {
      synthonFPs.emplace_back();
      synthonFPs[i].reserve(stu[i].size());
      const auto &sfps = reaction.getSynthons();
      for (size_t j = 0; j < stu[i].size(); ++j) {
        synthonFPs[i].emplace_back(
            new ExplicitBitVect(*sfps[i][stu[i][j]].second->getFP()));
      }
    }
    addFP.reset(new ExplicitBitVect(*reaction.getAddFP()));
    subtractFP.reset(new ExplicitBitVect(*reaction.getSubtractFP()));
  }

  std::vector<std::vector<std::unique_ptr<ExplicitBitVect>>> synthonFPs;
  std::unique_ptr<ExplicitBitVect> addFP;
  std::unique_ptr<ExplicitBitVect> subtractFP;
};
}  // namespace RDKit::SynthonSpaceSearch
#endif  // SYNTHONSPACEHITSET_H
