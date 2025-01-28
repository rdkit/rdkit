//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// This file declares an abstract base class for searching a synthon
// space.  Concrete base classes include SynthonSpaceSubstructureSearcher
// and SynthonSpaceFingerprintSearcher.

#ifndef SYNTHONSPACESEARCHER_H
#define SYNTHONSPACESEARCHER_H

#include <chrono>

#include <boost/random.hpp>

#include <RDGeneral/export.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>

using Clock = std::chrono::steady_clock;
using TimePoint = std::chrono::time_point<Clock>;

namespace RDKit {
class ROMol;

namespace SynthonSpaceSearch {

class SynthonSpaceSearcher {
 public:
  SynthonSpaceSearcher() = delete;
  SynthonSpaceSearcher(const ROMol &query,
                       const SynthonSpaceSearchParams &params,
                       SynthonSpace &space);
  SynthonSpaceSearcher(const SynthonSpaceSearcher &other) = delete;
  SynthonSpaceSearcher(SynthonSpaceSearcher &&other) = delete;
  SynthonSpaceSearcher &operator=(const SynthonSpaceSearcher &other) = delete;
  SynthonSpaceSearcher &operator=(SynthonSpaceSearcher &&other) = delete;

  virtual ~SynthonSpaceSearcher() = default;

  SearchResults search();

  SynthonSpace &getSpace() const { return d_space; }
  const ROMol &getQuery() const { return d_query; }
  const SynthonSpaceSearchParams &getParams() const { return d_params; }

 private:
  std::unique_ptr<boost::mt19937> d_randGen;

  const ROMol &d_query;
  const SynthonSpaceSearchParams &d_params;
  SynthonSpace &d_space;

  // Do the search of this fragSet against the SynthonSpace in the
  // appropriate way, for example by substructure or fingerprint
  // similarity.
  virtual std::vector<SynthonSpaceHitSet> searchFragSet(
      std::vector<std::unique_ptr<ROMol>> &fragSet) const = 0;
  // Make the hit, constructed from a specific combination of
  // synthons in the SynthonSet, and verify that it matches the
  // query in the appropriate way.  There'll be 1 entry in synthNums
  // for each synthon list in the reaction.  Returns an empty pointer
  // if the hit isn't accepted for whatever reason.
  std::unique_ptr<ROMol> buildAndVerifyHit(
      const std::unique_ptr<SynthonSet> &reaction,
      const std::vector<size_t> &synthNums,
      std::set<std::string> &resultsNames) const;
  // Some of the search methods (Rascal, for example) can do a quick
  // check on whether this set of synthons can match the query without having to
  // build the full molecule from the synthons.  They will over-ride this
  // function which by default passes everything.
  virtual bool quickVerify(
      [[maybe_unused]] const std::unique_ptr<SynthonSet> &reaction,
      [[maybe_unused]] const std::vector<size_t> &synthNums) const {
    return true;
  }
  virtual bool verifyHit(const ROMol &mol) const = 0;

  // Build the molecules from the synthons identified in reagentsToUse.
  // There should be bitset in reagentsToUse for each reagent set.
  // If not, it will fail.  Checks that all the results produced match the
  // query.  totHits is the maximum number of hits that are possible from
  // the hitsets, including duplicates.  Duplicates by name are not returned,
  // but duplicate SMILES from different reactions will be.  Hitsets will
  // be re-ordered on exit.
  void buildHits(std::vector<SynthonSpaceHitSet> &hitsets, size_t totHits,
                 const TimePoint *endTime, bool &timedOut,
                 std::vector<std::unique_ptr<ROMol>> &results) const;
  void buildAllHits(const std::vector<SynthonSpaceHitSet> &hitsets,
                    std::set<std::string> &resultsNames,
                    const TimePoint *endTime, bool &timedOut,
                    std::vector<std::unique_ptr<ROMol>> &results) const;
  void buildRandomHits(const std::vector<SynthonSpaceHitSet> &hitsets,
                       size_t totHits, std::set<std::string> &resultsNames,
                       const TimePoint *endTime, bool &timedOut,
                       std::vector<std::unique_ptr<ROMol>> &results) const;
  // get the subset of synthons for the given reaction to use for this
  // enumeration.
  std::vector<std::vector<ROMol *>> getSynthonsToUse(
      const std::vector<boost::dynamic_bitset<>> &synthonsToUse,
      const std::string &reaction_id) const;
};

}  // namespace SynthonSpaceSearch
}  // namespace RDKit
#endif  // SYNTHONSPACESEARCHER_H
