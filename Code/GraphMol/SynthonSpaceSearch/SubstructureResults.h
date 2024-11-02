//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_SYNTHONSPACE_SUBSTRUCTURERESULTS_H
#define RDKIT_SYNTHONSPACE_SUBSTRUCTURERESULTS_H

#include <GraphMol/ROMol.h>

namespace RDKit::SynthonSpaceSearch {
class RDKIT_SYNTHONSPACESEARCH_EXPORT SubstructureResults {
 public:
  explicit SubstructureResults() : d_maxNumResults(0) {}
  SubstructureResults(std::vector<std::unique_ptr<ROMol>> &&mols,
                      size_t maxNumRes);
  SubstructureResults(const SubstructureResults &other) = delete;
  SubstructureResults(SubstructureResults &&other) = delete;
  ~SubstructureResults() = default;

  SubstructureResults &operator=(const SubstructureResults &other);
  SubstructureResults &operator=(SubstructureResults &&other) = default;

  /*!
   * Returns the upper bound on the number of results the search might return.
   * There may be fewer than this in practice for several reasons such as
   * duplicate reagent sets being removed or the final product not matching the
   * query even though the synthons suggested it would.
   *
   * @return int
   */
  size_t getMaxNumResults() const { return d_maxNumResults; }
  /*!
   * Returns the hits from the search. Not necessarily all those possible,
   * just the maximum number requested.
   *
   * @return std::vector<std::unique_ptr<ROMol>>
   */
  const std::vector<std::unique_ptr<ROMol>> &getHitMolecules() const {
    return d_hitMolecules;
  }

  // Return a copy of this object.
  std::unique_ptr<SubstructureResults> clone() const;

 private:
  std::vector<std::unique_ptr<ROMol>> d_hitMolecules;
  size_t d_maxNumResults;
};

inline SubstructureResults::SubstructureResults(
    std::vector<std::unique_ptr<ROMol>> &&mols, size_t maxNumRes)
    : d_maxNumResults(maxNumRes) {
  d_hitMolecules = std::move(mols);
  mols.clear();
}

std::unique_ptr<SubstructureResults> SubstructureResults::clone() const {
  std::unique_ptr<SubstructureResults> retVal =
      std::make_unique<SubstructureResults>();

  retVal->d_maxNumResults = d_maxNumResults;
  retVal->d_hitMolecules.resize(d_hitMolecules.size());
  for (const auto &hm : d_hitMolecules) {
    retVal->d_hitMolecules.emplace_back(new ROMol(*hm));
  }
  return retVal;
}
}  // namespace RDKit::SynthonSpaceSearch

#endif  // RDKIT_SUBSTRUCTURERESULTS_H
