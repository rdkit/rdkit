//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_SYNTHONSPACE_SEARCHRESULTS_H
#define RDKIT_SYNTHONSPACE_SEARCHRESULTS_H

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>

namespace RDKit::SynthonSpaceSearch {
class RDKIT_SYNTHONSPACESEARCH_EXPORT SearchResults {
 public:
  explicit SearchResults() : d_maxNumResults(0) {}
  SearchResults(std::vector<std::unique_ptr<ROMol>> &&mols, size_t maxNumRes,
                bool timedOut, bool cancelled);
  SearchResults(const SearchResults &other);
  SearchResults(SearchResults &&other) = default;
  ~SearchResults() = default;

  SearchResults &operator=(const SearchResults &other);
  SearchResults &operator=(SearchResults &&other) = default;

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

  /*!
   * Returns whether the search timed out or not,
   * @return bool
   */
  bool getTimedOut() const { return d_timedOut; }
  /*!
   * Returns whether the search was cancelled or not,
   * @return bool
   */
  bool getCancelled() const { return d_cancelled; }

 private:
  std::vector<std::unique_ptr<ROMol>> d_hitMolecules;
  size_t d_maxNumResults;
  bool d_timedOut{false};
  bool d_cancelled{false};
};

inline SearchResults::SearchResults(std::vector<std::unique_ptr<ROMol>> &&mols,
                                    const size_t maxNumRes, bool timedOut,
                                    bool cancelled)
    : d_maxNumResults(maxNumRes), d_timedOut(timedOut), d_cancelled(cancelled) {
  d_hitMolecules = std::move(mols);
  mols.clear();
}

inline SearchResults::SearchResults(const SearchResults &other)
    : d_maxNumResults(other.d_maxNumResults),
      d_timedOut(other.d_timedOut),
      d_cancelled(other.d_cancelled) {
  for (const auto &hm : other.d_hitMolecules) {
    d_hitMolecules.emplace_back(new ROMol(*hm));
  }
}
}  // namespace RDKit::SynthonSpaceSearch

#endif  // RDKIT_SYNTHONSPACE_SEARCHRESULTS_H
