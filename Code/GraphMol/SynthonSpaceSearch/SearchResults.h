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

#include <unordered_map>
#include <functional>

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>

namespace RDKit::SynthonSpaceSearch {

// takes vector of search results; returns true if enough hits have been
// returned, false if the search should continue.
// Invoking the callback transfers ownership of the molecules to the
// callee, which avoids an extra copy of the molecule.
using SearchResultCallback =
    std::function<bool(std::vector<std::unique_ptr<ROMol>> &)>;

// A class holding a set of results from a search.  Contains the hit
// molecules and information about how the search progressed, whether
// it timed out etc.
class RDKIT_SYNTHONSPACESEARCH_EXPORT SearchResults {
 public:
  explicit SearchResults() : d_maxNumResults(0) {}
  SearchResults(std::vector<std::unique_ptr<ROMol>> &&mols,
                std::unique_ptr<ROMol> &&bestHit, std::uint64_t maxNumRes,
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
  std::uint64_t getMaxNumResults() const { return d_maxNumResults; }
  /*!
   * Returns the hit molecules from the search. Not necessarily all
   * those possible, just the maximum number requested.
   *
   * @return const std::vector<std::unique_ptr<ROMol>> &
   */
  const std::vector<std::unique_ptr<ROMol>> &getHitMolecules() const {
    return d_hitMolecules;
  }

  /*!
   * Returns the best hit found in a similarity search, even when none
   * were under the search threshold.  May be empty, for example if it
   * was a substructure search.  Also, may be empty if none of the synthons
   * matched within the fragment similarity threshold i.e. nothing came close
   * to a match.
   *
   * @return const std::unique_ptr<ROMol> &
   */
  const std::unique_ptr<ROMol> &getBestHit() const { return d_bestHitFound; }

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

  // Merge other into this, keeping only molecules with unique
  // names and destroying contents of other on exit.  If an
  // incoming molecule is the same name as one already present
  // and both have "Similarity" prop, then the one with the
  // higher value is kept.
  void mergeResults(SearchResults &other);

 private:
  std::vector<std::unique_ptr<ROMol>> d_hitMolecules;
  // Only used when merging in another set, so will be
  // filled in then if needed, empty otherwise.
  std::unordered_map<std::string, size_t> d_molNames;

  std::uint64_t d_maxNumResults;
  bool d_timedOut{false};
  bool d_cancelled{false};

  // Even if there were no hits within the cutoff, the similarity
  // searches will but the best hit it found in this molecule, with
  // the similarity value in the "Similarity" property as for the
  // genuine hits.
  std::unique_ptr<ROMol> d_bestHitFound;
};

}  // namespace RDKit::SynthonSpaceSearch

#endif  // RDKIT_SYNTHONSPACE_SEARCHRESULTS_H
