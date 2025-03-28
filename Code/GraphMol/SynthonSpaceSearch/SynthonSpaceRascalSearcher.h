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
// that does RASCAL similarity searching of the SynthonSpace.

#ifndef SYNTHONSPACERASCALSEARCHER_H
#define SYNTHONSPACERASCALSEARCHER_H

#include <RDGeneral/export.h>

#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>

namespace RDKit::RascalMCES {
struct RascalOptions;
}

namespace RDKit::SynthonSpaceSearch {

class SynthonSpaceRascalSearcher : public SynthonSpaceSearcher {
 public:
  SynthonSpaceRascalSearcher() = delete;
  SynthonSpaceRascalSearcher(const ROMol &query,
                             const RascalMCES::RascalOptions &options,
                             const SynthonSpaceSearchParams &params,
                             SynthonSpace &space);

  std::vector<std::unique_ptr<SynthonSpaceHitSet>> searchFragSet(
      const std::vector<std::unique_ptr<ROMol>> &fragSet,
      const SynthonSet &reaction) const override;

 private:
  void extraSearchSetup(
      std::vector<std::vector<std::unique_ptr<ROMol>>> &fragSets) override;

  bool quickVerify(const SynthonSpaceHitSet *hitset,
                   const std::vector<size_t> &synthNums) const override;

  bool verifyHit(const ROMol &hit) const override;

  const RascalMCES::RascalOptions &d_rascalOptions;
  RascalMCES::RascalOptions d_rascalFragOptions;
};

}  // namespace RDKit::SynthonSpaceSearch

#endif  // SYNTHONSPACERASCALSEARCHER_H
