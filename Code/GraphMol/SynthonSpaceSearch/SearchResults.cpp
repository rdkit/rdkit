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
                             const std::uint64_t maxNumRes, bool timedOut,
                             bool cancelled)
    : d_maxNumResults(maxNumRes), d_timedOut(timedOut), d_cancelled(cancelled) {
  d_hitMolecules = std::move(mols);
  mols.clear();
}

SearchResults::SearchResults(const SearchResults &other)
    : d_maxNumResults(other.d_maxNumResults),
      d_timedOut(other.d_timedOut),
      d_cancelled(other.d_cancelled) {
  for (const auto &hm : other.d_hitMolecules) {
    d_hitMolecules.emplace_back(new ROMol(*hm));
  }
}

void SearchResults::mergeResults(SearchResults &other) {
  d_maxNumResults += other.d_maxNumResults;
  if (other.d_timedOut) d_timedOut = true;
  if (other.d_cancelled) d_cancelled = true;
  if (d_molNames.empty()) {
    for (const auto &mol : d_hitMolecules) {
      d_molNames.insert(mol->getProp<std::string>(common_properties::_Name));
    }
  }
  d_hitMolecules.reserve(d_hitMolecules.size() + other.d_hitMolecules.size());
  for (auto &mol : other.d_hitMolecules) {
    if (auto it = d_molNames.find(
            mol->getProp<std::string>(common_properties::_Name));
        it == d_molNames.end()) {
      d_molNames.insert(mol->getProp<std::string>(common_properties::_Name));
      d_hitMolecules.emplace_back(std::move(mol));
    }
  }
  other.d_hitMolecules.clear();
  other.d_maxNumResults = 0;
}
}  // namespace RDKit::SynthonSpaceSearch