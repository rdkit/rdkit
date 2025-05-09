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

// Concrete class that does the search by fingerprint similarity.
class SynthonSpaceFingerprintSearcher : public SynthonSpaceSearcher {
 public:
  SynthonSpaceFingerprintSearcher() = delete;
  SynthonSpaceFingerprintSearcher(
      const ROMol &query, const FingerprintGenerator<std::uint64_t> &fpGen,
      const SynthonSpaceSearchParams &params, SynthonSpace &space);

  std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchFragSet(
      const std::vector<std::unique_ptr<ROMol>> &fragSet,
      const SynthonSet &reaction) const override;

 private:
  std::unique_ptr<ExplicitBitVect> d_queryFP;

  const FingerprintGenerator<std::uint64_t> &d_fpGen;
  // These are the fingerprints for the fragments in this search.
  // The fingerprint in d_fragFPs are keyed on the addresses of the
  // corresponding fragment.  There are usually multiple fragments
  // with the same SMILES and this way the fingerprints are
  // generated the minimum number of times.
  // d_fragFPPool is never read, it is just used as a repository
  // for the fragFPs for the lifetime of the search.
  std::vector<std::unique_ptr<ExplicitBitVect>> d_fragFPPool;
  std::vector<std::pair<void *, ExplicitBitVect *>> d_fragFPs;

  void extraSearchSetup(
      std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) override;

  bool quickVerify(const SynthonSpaceHitSet *hitset,
                   const std::vector<size_t> &synthNums) const override;
  bool verifyHit(const ROMol &hit) const override;
};
}  // namespace RDKit::SynthonSpaceSearch

#endif  // SYNTHONSPACEFINGERPRINTSEARCHER_H
