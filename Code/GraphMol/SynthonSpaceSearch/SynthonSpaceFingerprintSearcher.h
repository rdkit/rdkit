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
  const FingerprintGenerator<std::uint64_t> &d_fpGen;

  std::vector<SynthonSpaceHitSet> searchFragSet(
      std::vector<std::unique_ptr<ROMol>> &fragSet) const override;
  bool quickVerify(const std::unique_ptr<SynthonSet> &reaction,
                   const std::vector<size_t> &synthNums) const override;
  bool verifyHit(const ROMol &hit) const override;
};
}  // namespace RDKit::SynthonSpaceSearch

#endif  // SYNTHONSPACEFINGERPRINTSEARCHER_H
