//
//  Copyright (C) 2017-2022 Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_DECOMP_DATA
#define RGROUP_DECOMP_DATA

#include "RGroupCore.h"
#include "RGroupScore.h"
#include "RGroupFingerprintScore.h"
#include <vector>
#include <map>

namespace RDKit {

extern const std::string _rgroupInputDummy;

struct RGroupDecompData {
  // matches[mol_idx] == vector of potential matches
  std::map<int, RCore> cores;
  std::map<std::string, int> newCores;  // new "cores" found along the way
  int newCoreLabel = EMPTY_CORE_LABEL;
  // this caches the running product of permutations
  // across calls to process()
  size_t permutationProduct = 1;
  // this caches the size of the previous matches vector
  // such that the size of the current chunk can be inferred
  size_t previousMatchSize = 0;
  // the default for Greedy/GreedyChunks is keeping only the best
  // permutation after each call to process()
  bool prunePermutations = true;
  RGroupDecompositionParameters params;

  std::vector<std::vector<RGroupMatch>> matches;
  std::set<int> labels;
  std::vector<size_t> permutation;
  unsigned int pruneLength = 0U;
  FingerprintVarianceScoreData prunedFingerprintVarianceScoreData;
  std::map<int, std::vector<int>> userLabels;

  std::vector<int> processedRlabels;

  std::map<int, int> finalRlabelMapping;
  mutable RGroupScorer rGroupScorer;

  RGroupDecompData(const RWMol &inputCore,
                   RGroupDecompositionParameters inputParams);

  RGroupDecompData(const std::vector<ROMOL_SPTR> &inputCores,
                   RGroupDecompositionParameters inputParams);

  void addCore(const ROMol &inputCore);

  void prepareCores();

  void setRlabel(Atom *atom, int rlabel);

  int getRlabel(Atom *atom) const;

  double scoreFromPrunedData(const std::vector<size_t> &permutation,
                             bool reset = true);

  void prune();

  // Return the RGroups with the current "best" permutation
  //  of matches.
  std::vector<RGroupMatch> GetCurrentBestPermutation() const;

  class UsedLabels {
   public:
    std::set<int> labels_used;

    bool add(int rlabel);

    int next();
  };

  void addCoreUserLabels(const RWMol &core, std::set<int> &userLabels);

  void addAtoms(RWMol &mol,
                const std::vector<std::pair<Atom *, Atom *>> &atomsToAdd);

  bool replaceHydrogenCoreDummy(const RGroupMatch &match, RWMol &core,
                                const Atom &atom, const int currentLabel,
                                const int rLabel);

  void relabelCore(RWMol &core, std::map<int, int> &mappings,
                   UsedLabels &used_labels, const std::set<int> &indexLabels,
                   const std::map<int, std::vector<int>> &extraAtomRLabels,
                   const RGroupMatch *const match = nullptr);

  void relabelRGroup(RGroupData &rgroup, const std::map<int, int> &mappings);

  // relabel the core and sidechains using the specified user labels
  //  if matches exist for non labelled atoms, these are added as well
  void relabel();

  double score(const std::vector<size_t> &permutation,
               FingerprintVarianceScoreData *fingerprintVarianceScoreData =
                   nullptr) const;

  RGroupDecompositionProcessResult process(bool pruneMatches,
                                           bool finalize = false);

private:
  void addInputCore(const ROMol &inputCore);
};
}  // namespace RDKit

#endif
