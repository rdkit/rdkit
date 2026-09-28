//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/SynthonSpaceSearch/SearchResults.h>

namespace RDKit::SynthonSpaceSearch {
SearchResults::SearchResults(std::vector<std::unique_ptr<ROMol>> &&mols,
                             std::unique_ptr<ROMol> &&bestHitFound,
                             const std::uint64_t maxNumRes, bool timedOut,
                             bool cancelled)
    : d_maxNumResults(maxNumRes), d_timedOut(timedOut), d_cancelled(cancelled) {
  d_hitMolecules = std::move(mols);
  mols.clear();
  d_bestHitFound = std::move(bestHitFound);
  bestHitFound.reset();
}

SearchResults::SearchResults(const SearchResults &other)
    : d_maxNumResults(other.d_maxNumResults),
      d_timedOut(other.d_timedOut),
      d_cancelled(other.d_cancelled) {
  for (const auto &hm : other.d_hitMolecules) {
    d_hitMolecules.emplace_back(new ROMol(*hm));
  }
  if (other.d_bestHitFound) {
    d_bestHitFound.reset(new ROMol(*other.d_bestHitFound));
  }
}

void SearchResults::mergeResults(SearchResults &other) {
  d_maxNumResults += other.d_maxNumResults;
  if (other.d_timedOut) {
    d_timedOut = true;
  }
  if (other.d_cancelled) {
    d_cancelled = true;
  }
  if (d_molNames.empty()) {
    for (size_t i = 0; i < d_hitMolecules.size(); i++) {
      d_molNames.insert(std::make_pair(
          d_hitMolecules[i]->getProp<std::string>(common_properties::_Name),
          i));
    }
  }
  d_hitMolecules.reserve(d_hitMolecules.size() + other.d_hitMolecules.size());
  for (auto &mol : other.d_hitMolecules) {
    if (auto it = d_molNames.find(
            mol->getProp<std::string>(common_properties::_Name));
        it == d_molNames.end()) {
      d_molNames.insert(
          std::make_pair(mol->getProp<std::string>(common_properties::_Name),
                         d_hitMolecules.size()));
      d_hitMolecules.emplace_back(std::move(mol));
    } else {
      auto &currHit = d_hitMolecules[it->second];
      double currSim, newSim;
      if (currHit->getPropIfPresent<double>("Similarity", currSim) &&
          mol->getPropIfPresent<double>("Similarity", newSim)) {
        if (newSim > currSim) {
          currHit = std::move(mol);
        }
      }
    }
  }
  if (d_bestHitFound && other.d_bestHitFound) {
    if (d_bestHitFound->getProp<double>("Similarity") <
        other.d_bestHitFound->getProp<double>("Similarity")) {
      d_bestHitFound.reset(new ROMol(*other.d_bestHitFound));
    }
  } else if (!d_bestHitFound && other.d_bestHitFound) {
    d_bestHitFound.reset(new ROMol(*other.d_bestHitFound));
  }
  other.d_hitMolecules.clear();
  other.d_maxNumResults = 0;
}
}  // namespace RDKit::SynthonSpaceSearch