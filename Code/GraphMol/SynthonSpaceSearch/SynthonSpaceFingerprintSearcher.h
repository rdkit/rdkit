//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// This file declares a concrete class derived from SynthonSpaceSearcher
// that does fingerprint similarity searching of the SynthonSpace.

#ifndef SYNTHONSPACEFINGERPRINTSEARCHER_H
#define SYNTHONSPACEFINGERPRINTSEARCHER_H

#include <RDGeneral/export.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>

namespace RDKit::SynthonSpaceSearch {

class SynthonSpaceFingerprintSearcher : public SynthonSpaceSearcher {
 public:
  SynthonSpaceFingerprintSearcher() = delete;
  SynthonSpaceFingerprintSearcher(
      const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
      const SynthonSpaceSearchParams &params, SynthonSpace &space);

 private:
  std::unique_ptr<ExplicitBitVect> d_queryFP;
  std::unique_ptr<ExplicitBitVect> d_foldedQueryFP;

  const FingerprintGenerator<std::uint64_t> &d_fpGen;
  // These are the fingerprints for the fragments in this search.
  // It's thread-safe because each search creates its own Searcher
  // object so multiple searches in different threads will be in
  // different Searcher objects.  The fingerprints are keyed on
  // the address of the corresponding fragment.
  std::map<void *, std::unique_ptr<ExplicitBitVect>> d_fragFPs;

  void extraSearchSetup(
      std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) override;

  std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchFragSet(
      std::vector<std::unique_ptr<ROMol>> &fragSet,
      const SynthonSet &reaction) const override;
  bool quickVerify(const SynthonSpaceHitSet *hitset,
                   const std::vector<size_t> &synthNums) const override;
  bool verifyHit(const ROMol &hit) const override;
};
}  // namespace RDKit::SynthonSpaceSearch

#endif  // SYNTHONSPACEFINGERPRINTSEARCHER_H
