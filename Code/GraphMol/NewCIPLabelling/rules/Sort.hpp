//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include "Priority.hpp"
#include <vector>

namespace RDKit {
namespace NewCIPLabelling {

template <typename A, typename B>
class SequenceRule;

/**
 * A simple insertion sort for ligands. The number of ligands is not likely to
 * be very larger as such doing a merge sort would have little benefit.
 *
 */
template <typename A, typename B>
class Sort {
 private:
  unsigned ruleMax = 0;

  const std::vector<const SequenceRule<A, B>*> rules;

 public:
  Sort(const SequenceRule<A, B>* comparator) : rules{{comparator}} {}

  Sort(const std::vector<const SequenceRule<A, B>*>& comparators)
      : rules{comparators} {}

  const std::vector<const SequenceRule<A, B>*>& getRules() const {
    return rules;
  }

  Priority prioritise(const Node<A, B>* node,
                      std::vector<Edge<A, B>*>& edges) const {
    return prioritise(node, edges, true);
  }

  Priority prioritise(const Node<A, B>* node, std::vector<Edge<A, B>*>& edges,
                      bool deep) const {
    bool unique = true;
    int numPseudoAsym = 0;

    for (auto i = 0u; i < edges.size(); ++i)
      for (auto j = i; j > 0; --j) {
        int cmp = compareLigands(node, edges[j - 1], edges[j], deep);

        if (cmp < -1 || cmp > +1) {
          ++numPseudoAsym;
        }

        if (cmp < 0) {
          swap(edges, j, j - 1);
        } else {
          if (cmp == 0) {
            unique = false;
          }
          break;
        }
      }

    return Priority(unique, ruleMax, numPseudoAsym == 1);
  }

  int compareLigands(const Node<A, B>* node, const Edge<A, B>* a,
                     const Edge<A, B>* b, bool deep) const {
    // ensure 'up' edges are moved to the front
    if (!a->isBeg(node) && b->isBeg(node)) {
      return +1;
    } else if (a->isBeg(node) && !b->isBeg(node)) {
      return -1;
    }

    for (auto i = 0u; i < rules.size(); ++i) {
      const auto& rule = rules[i];
      int cmp = rule->getComparision(a, b, deep);

      if (cmp != 0) {
        // Is this just statistics ?
        const_cast<Sort<A, B>*>(this)->ruleMax = std::max(ruleMax, i);
        return cmp;
      }
    }
    return 0;
  }

  void swap(std::vector<Edge<A, B>*>& list, int i, int j) const {
    std::swap(list[i], list[j]);
  }

  std::vector<std::vector<Edge<A, B>*>> getGroups(
      const std::vector<Edge<A, B>*>& sorted) const {
    // would be nice to have this integrated whilst sorting - may provide a
    // small speed increase but as most of our lists are small we take use
    // ugly sort then group approach
    auto groups = std::vector<std::vector<Edge<A, B>*>>{};

    Edge<A, B>* prev = nullptr;
    for (auto* edge : sorted) {
      if (prev == nullptr ||
          compareLigands(prev->getBeg(), prev, edge, true) != 0) {
        groups.emplace_back();
      }
      prev = edge;
      groups.back().push_back(edge);
    }

    return groups;
  }
};  // namespace NewCIPLabelling

}  // namespace NewCIPLabelling
}  // namespace RDKit