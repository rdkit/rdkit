//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef SYNTHONSPACEHITSET_H
#define SYNTHONSPACEHITSET_H

#include <numeric>
#include <string>
#include <vector>

#include <RDGeneral/export.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SynthonSpaceSearch/SynthonSet.h>

namespace RDKit::SynthonSpaceSearch {
// Holds the information about a set of hits.  The molecules can be built
// by making all combinations of synthons, one taken from each synthon set.
// The synthons are copied from the SynthonSet so it is not needed to
// build the hits.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceHitSet {
  SynthonSpaceHitSet() = delete;
  SynthonSpaceHitSet(const SynthonSet &reaction,
                     const std::vector<std::vector<size_t>> &stu,
                     const std::vector<std::unique_ptr<ROMol>> &fragSet);
  SynthonSpaceHitSet(const SynthonSpaceHitSet &lhs) = delete;
  SynthonSpaceHitSet(SynthonSpaceHitSet &&lhs) = delete;
  virtual ~SynthonSpaceHitSet() = default;
  SynthonSpaceHitSet &operator=(const SynthonSpaceHitSet &rhs) = delete;
  SynthonSpaceHitSet &operator=(SynthonSpaceHitSet &&rhs) = delete;

  const SynthonSet *d_reaction{nullptr};
  std::vector<std::vector<std::pair<std::string, const ROMol *>>> synthonsToUse;
  size_t numHits{0};
  // The fragments that this hitset is derived from, useful for debugging.
  std::vector<const ROMol *> frags;
};

// This sub-class holds results from a SynthonSpaceFingerprintSearch.
// It needs pointers to the synthon fingerprints as well.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceFPHitSet
    : SynthonSpaceHitSet {
  SynthonSpaceFPHitSet() = delete;
  SynthonSpaceFPHitSet(const SynthonSet &reaction,
                       const std::vector<std::vector<size_t>> &stu,
                       const std::vector<std::unique_ptr<ROMol>> &fragSet);
  SynthonSpaceFPHitSet(const SynthonSpaceFPHitSet &lhs) = delete;
  SynthonSpaceFPHitSet(SynthonSpaceFPHitSet &&lhs) = delete;
  ~SynthonSpaceFPHitSet() = default;
  SynthonSpaceFPHitSet &operator=(const SynthonSpaceFPHitSet &rhs) = delete;
  SynthonSpaceFPHitSet &operator=(SynthonSpaceFPHitSet &&rhs) = delete;

  std::vector<std::vector<ExplicitBitVect *>> synthonFPs;
  ExplicitBitVect *addFP;
  ExplicitBitVect *subtractFP;
};

}  // namespace RDKit::SynthonSpaceSearch
#endif  // SYNTHONSPACEHITSET_H
