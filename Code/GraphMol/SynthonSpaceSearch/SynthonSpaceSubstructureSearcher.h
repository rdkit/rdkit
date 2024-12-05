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
// that does substructure searching of the SynthonSpace.

#ifndef SYNTHONSPACESUBSTRUCTURESEARCHER_H
#define SYNTHONSPACESUBSTRUCTURESEARCHER_H

#include <RDGeneral/export.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSpaceSearcher.h>

namespace RDKit::SynthonSpaceSearch {

class SynthonSpaceSubstructureSearcher : public SynthonSpaceSearcher {
 public:
  SynthonSpaceSubstructureSearcher() = delete;
  SynthonSpaceSubstructureSearcher(const ROMol &query,
                                   const SynthonSpaceSearchParams &params,
                                   SynthonSpace &space)
      : SynthonSpaceSearcher(query, params, space) {}

 private:
  std::vector<SynthonSpaceHitSet> searchFragSet(
      std::vector<std::unique_ptr<ROMol>> &fragSet) const override;
  bool verifyHit(const ROMol &hit) const override;
};

}  // namespace RDKit::SynthonSpaceSearch

#endif  // SYNTHONSPACESUBSTRUCTURESEARCHER_H
