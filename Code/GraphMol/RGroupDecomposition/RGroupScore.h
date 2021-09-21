//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_SCORE_H
#define RGROUP_SCORE_H

#include "RGroupMatch.h"
#include <vector>
#include <deque>
#include <set>
namespace RDKit {

//! iterate through all possible permutations of the rgroups
struct CartesianProduct {
  std::vector<size_t> permutation;
  std::vector<size_t> sizes;
  std::deque<size_t> bases;
  size_t maxPermutations;
  size_t permutationCount;
  CartesianProduct(const std::vector<size_t> &inputSizes)
      : permutation(inputSizes.size(), 0),
        sizes(inputSizes),
        permutationCount(0) {
    maxPermutations = 1;
    for (unsigned long size : sizes) {
      bases.push_front(maxPermutations);
      maxPermutations *= size;  // may overflow....
    }
  }

  bool next() {
    ++permutationCount;
    if (permutationCount == 1) {
      return true;
    }

    return increment(0);
  }

  size_t value(const std::vector<size_t> &p) const {
    size_t v = 0;
    for (size_t i = 0; i < p.size(); ++i) {
      v += bases[i] * p[i];
    }
    return v;
  }

  size_t value() { return value(permutation); }

  bool increment(size_t rowToIncrement) {
    if (permutationCount > maxPermutations) {
      return false;
    }

    permutation[rowToIncrement] += 1;
    size_t max_index_of_row = sizes[rowToIncrement] - 1;
    if (permutation[rowToIncrement] > max_index_of_row) {
      permutation[rowToIncrement] = 0;
      return increment(rowToIncrement + 1);
    }
    return true;
  }
};

class RDKIT_RGROUPDECOMPOSITION_EXPORT RGroupScorer {
 public:
  RGroupScorer(){};
  RGroupScorer(const std::vector<std::vector<size_t>> &permutations,
               double score);
  double matchScore(const std::vector<size_t> &permutation,
                    const std::vector<std::vector<RGroupMatch>> &matches,
                    const std::set<int> &labels);
  void setBestPermutation(const std::vector<size_t> &permutation, double score);
  const std::vector<size_t> &getBestPermutation() const {
    return d_saved.permutation;
  }
  void startProcessing();
  void pushTieToStore(const std::vector<size_t> &permutation);
  void breakTies(const std::vector<std::vector<RGroupMatch>> &matches,
                 const std::set<int> &labels,
                 const std::unique_ptr<CartesianProduct> &iterator,
                 const std::chrono::steady_clock::time_point &t0,
                 double timeout);
  void clearTieStore();
  size_t tieStoreSize() const { return d_store.size(); }
  double getBestScore() const { return d_bestScore; }

 private:
  void restoreInitialState() { d_current = d_initial; }
  struct RLabelData {
    int numRGroups = 0;
    std::vector<std::map<std::string, unsigned int>> matchSetVect;
    std::map<std::set<int>, size_t> linkerMatchSet;
  };
  struct State {
    void computeTieBreakingCriteria(
        const std::vector<std::vector<RGroupMatch>> &matches,
        const std::vector<int> &orderedLabels, std::vector<int> &heavyCounts) {
      // heavyCounts is a vector which has the same size of labels
      // for each label we add an increment if a molecule
      // bears an R-group at that label
      size_t offset = matches.size() - permutation.size();
      // numMatchedUserRGroups counts the total number of user labelled r
      // groups filled in this permutation.  We want to maximize this number
      size_t i = 0;
      for (int label : orderedLabels) {
        for (size_t m = 0; m < permutation.size(); ++m) {  // for each molecule
          // for each molecule, check if we add an R-group at this negative
          // label if we do, count it once. So we know how many different
          // negative labels we have filled: we prefer permutations which fill
          // less, as it means we have added less groups on different positions
          const auto &match = matches[m + offset][permutation[m]];
          auto rg = match.rgroups.find(label);
          if (rg != match.rgroups.end() && !rg->second->is_hydrogen) {
            if (label < 0 && heavyCounts.at(i) == 0) {
              ++numAddedRGroups;
            } else if (label > 0) {
              ++numMatchedUserRGroups;
            }
            ++heavyCounts[i];
          }
        }
        ++i;
      }
    }

    int N = 0;
    int numAddedRGroups = 0;
    int numMatchedUserRGroups = 0;
    std::map<int, int> heavyCountPerLabel;
    std::map<int, RLabelData> labelDataMap;
    std::vector<size_t> permutation;
  };
  double d_bestScore = 0.0;
  State d_current;
  State d_initial;
  State d_saved;
  std::deque<State> d_store;
};

}  // namespace RDKit
#endif
