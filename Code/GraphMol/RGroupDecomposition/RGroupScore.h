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

  size_t value() {
    size_t v = 0;
    for (size_t i = 0; i < permutation.size(); ++i) {
      v += bases[i] * permutation[i];
    }
    return v;
  }

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
  
double score(const std::vector<size_t> &permutation,
             const std::vector<std::vector<RGroupMatch>> &matches,
             const std::set<int> &labels);
}
#endif
