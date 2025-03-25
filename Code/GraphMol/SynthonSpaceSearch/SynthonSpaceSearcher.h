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
#include <random>

#include <RDGeneral/export.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpace.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceHitSet.h>
#include <GraphMol/SynthonSpaceSearch/SearchResults.h>
#include <boost/spirit/home/support/common_terminals.hpp>

using Clock = std::chrono::steady_clock;
using TimePoint = std::chrono::time_point<Clock>;

namespace RDKit {
class ROMol;

namespace SynthonSpaceSearch {

// Abstract base class for searching the SynthonSpace.
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

  // Do the search of this fragSet against the SynthonSet in the
  // appropriate way, for example by substructure or fingerprint
  // similarity.
  virtual std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchFragSet(
      const std::vector<std::unique_ptr<ROMol>> &fragSet,
      const SynthonSet &reaction) const = 0;

  // Make the hit, constructed from a specific combination of
  // synthons in the hitset, and verify that it matches the
  // query in the appropriate way.  There'll be 1 entry in synthNums
  // for each synthon list in the hitset.  Returns an empty pointer
  // if the hit isn't accepted for whatever reason.
  std::unique_ptr<ROMol> buildAndVerifyHit(
      const SynthonSpaceHitSet *hitset,
      const std::vector<size_t> &synthNums) const;

 private:
  std::unique_ptr<std::mt19937> d_randGen;

  const ROMol &d_query;
  const SynthonSpaceSearchParams &d_params;
  SynthonSpace &d_space;

  // Some of the search methods might need extra setup of the fragment
  // sets.  The FingerprintSearcher, for example, needs fingerprints
  // for all the fragments.  The SubstructureSearcher needs connector
  // regions and information about them.
  virtual void extraSearchSetup(
      [[maybe_unused]] std::vector<std::vector<std::unique_ptr<ROMol>>>
          &fragSets) {}

  std::vector<std::unique_ptr<SynthonSpaceHitSet>> doTheSearch(
      std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets,
      const TimePoint *endTime, bool &timedOut, std::uint64_t &totHits);

  // Some of the search methods (fingerprints, for example) can do a quick
  // check on whether this set of synthons can match the query without having to
  // build the full molecule from the synthons.  They will over-ride this
  // function which by default passes everything.
  virtual bool quickVerify(
      [[maybe_unused]] const SynthonSpaceHitSet *hitset,
      [[maybe_unused]] const std::vector<size_t> &synthNums) const {
    return true;
  }
  // Checks that the given molecule is definitely a hit according to
  // the derived class' criteria.
  virtual bool verifyHit(const ROMol &mol) const = 0;

  // Build the molecules from the synthons identified in hitsets.
  // Checks that all the results produced match the
  // query.  Duplicates by name are not returned,
  // but duplicate SMILES from different reactions will be.
  // Hitsets will be re-ordered on exit.
  void buildHits(std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
                 const TimePoint *endTime, bool &timedOut,
                 std::vector<std::unique_ptr<ROMol>> &results) const;
  void buildAllHits(
      const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
      const TimePoint *endTime, bool &timedOut,
      std::vector<std::unique_ptr<ROMol>> &results) const;
  void makeHitsFromToTry(
      const std::vector<
          std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>> &toTry,
      const TimePoint *endTime,
      std::vector<std::unique_ptr<ROMol>> &results) const;
  void processToTrySet(
      std::vector<std::pair<const SynthonSpaceHitSet *, std::vector<size_t>>>
          &toTry,
      const TimePoint *endTime,
      std::vector<std::unique_ptr<ROMol>> &results) const;

  // get the subset of synthons for the given reaction to use for this
  // enumeration.
  std::vector<std::vector<ROMol *>> getSynthonsToUse(
      const std::vector<boost::dynamic_bitset<>> &synthonsToUse,
      const std::string &reaction_id) const;
};

#if 0
  // Build the molecules from the synthons identified in hitsets.
  // Checks that all the results produced match the
  // query.  totHits is the maximum number of hits that are possible from
  // the hitsets, including duplicates.  Duplicates by name are not returned,
  // but duplicate SMILES from different reactions will be.  Hitsets will
  // be re-ordered on exit.
  void buildHits(std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
                 size_t totHits, const TimePoint *endTime, bool &timedOut,
                 std::vector<std::unique_ptr<ROMol>> &results) const;
  void buildAllHits(
      const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
      std::set<std::string> &resultsNames, const TimePoint *endTime,
      bool &timedOut, std::vector<std::unique_ptr<ROMol>> &results) const;
  void buildRandomHits(
      const std::vector<std::unique_ptr<SynthonSpaceHitSet>> &hitsets,
      size_t totHits, std::set<std::string> &resultsNames,
      const TimePoint *endTime, bool &timedOut,
      std::vector<std::unique_ptr<ROMol>> &results) const;
  // get the subset of synthons for the given reaction to use for this
  // enumeration.
  std::vector<std::vector<ROMol *>> getSynthonsToUse(
      const std::vector<boost::dynamic_bitset<>> &synthonsToUse,
      const std::string &reaction_id) const;
};
#endif

}  // namespace SynthonSpaceSearch
}  // namespace RDKit
#endif  // SYNTHONSPACESEARCHER_H
