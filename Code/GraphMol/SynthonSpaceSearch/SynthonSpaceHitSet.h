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
  // The outer vector of synthonsToUse is the size of the synthon lists
  // in the reaction.  So synthonsToUse[0] is all the synthons selected
  // from reaction->synthons[0], likewise for synthonsToUse[1] etc.
  // There should be at least one entry in each outer vector, but there
  // may be different numbers in each.
  std::vector<std::vector<std::pair<std::string, const Synthon *>>>
      synthonsToUse;
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
  ~SynthonSpaceFPHitSet() override = default;
  SynthonSpaceFPHitSet &operator=(const SynthonSpaceFPHitSet &rhs) = delete;
  SynthonSpaceFPHitSet &operator=(SynthonSpaceFPHitSet &&rhs) = delete;

  std::vector<std::vector<ExplicitBitVect *>> synthonFPs;
  ExplicitBitVect *addFP;
  ExplicitBitVect *subtractFP;
};

// This sub-class holds results from a SynthonSpaceShapeSearch.
// It needs pointers to the fragment shapes and the order of the
// synthon sets as they matched the fragments.
struct RDKIT_SYNTHONSPACESEARCH_EXPORT SynthonSpaceShapeHitSet
    : SynthonSpaceHitSet {
  SynthonSpaceShapeHitSet() = delete;
  SynthonSpaceShapeHitSet(const SynthonSet &reaction,
                          const std::vector<std::vector<size_t>> &stu,
                          const std::vector<std::unique_ptr<ROMol>> &fragSet,
                          const std::vector<SearchShapeInput *> &fShapes,
                          const std::vector<unsigned int> &sSetOrder);
  SynthonSpaceShapeHitSet(const SynthonSpaceShapeHitSet &lhs) = delete;
  SynthonSpaceShapeHitSet(SynthonSpaceShapeHitSet &&lhs) = delete;
  ~SynthonSpaceShapeHitSet() override = default;
  SynthonSpaceShapeHitSet &operator=(const SynthonSpaceShapeHitSet &rhs) =
      delete;
  SynthonSpaceShapeHitSet &operator=(SynthonSpaceShapeHitSet &&rhs) = delete;

  // The shapes corresponding to the frags, for the approximate similarity
  // calculation.
  std::vector<SearchShapeInput *> fragShapes;
  // The order of the synthon sets in the SynthonSet as they matched the
  // fragments.
  std::vector<unsigned int> synthonSetOrder;
};

}  // namespace RDKit::SynthonSpaceSearch
#endif  // SYNTHONSPACEHITSET_H
