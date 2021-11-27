//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Sort.h"
#include "rules/SequenceRule.h"

namespace RDKit {
namespace CIPLabeler {

Sort::Sort(const SequenceRule *comparator) : d_rules{comparator} {}

Sort::Sort(std::vector<const SequenceRule *> comparators)
    : d_rules{std::move(comparators)} {}

const std::vector<const SequenceRule *> &Sort::getRules() const {
  return d_rules;
}

Priority Sort::prioritize(const Node *node, std::vector<Edge *> &edges,
                          bool deep) const {
  bool unique = true;
  int numPseudoAsym = 0;

  for (auto i = 0u; i < edges.size(); ++i) {
    for (auto j = i; j > 0; --j) {
      int cmp = compareSubstituents(node, edges[j - 1], edges[j], deep);

      if (cmp < -1 || cmp > +1) {
        ++numPseudoAsym;
      }

      if (cmp < 0) {
        std::swap(edges[j], edges[j - 1]);
      } else {
        if (cmp == 0) {
          unique = false;
        }
        break;
      }
    }
  }

  return {unique, numPseudoAsym == 1};
}

int Sort::compareSubstituents(const Node *node, const Edge *a, const Edge *b,
                              bool deep) const {
  // ensure 'out' edges are moved to the front
  if (!a->isBeg(node) && b->isBeg(node)) {
    return +1;
  } else if (a->isBeg(node) && !b->isBeg(node)) {
    return -1;
  }

  for (const auto &rule : d_rules) {
    int cmp = rule->getComparision(a, b, deep);

    if (cmp != 0) {
      return cmp;
    }
  }
  return 0;
}

std::vector<std::vector<Edge *>> Sort::getGroups(
    const std::vector<Edge *> &sorted) const {
  // would be nice to have this integrated whilst sorting - may provide a
  // small speed increase but as most of our lists are small we take use
  // ugly sort then group approach
  std::vector<std::vector<Edge *>> groups;

  Edge *prev = nullptr;
  for (auto *edge : sorted) {
    if (prev == nullptr ||
        compareSubstituents(prev->getBeg(), prev, edge, true) != 0) {
      groups.emplace_back();
    }
    prev = edge;
    groups.back().push_back(edge);
  }

  return groups;
}

}  // namespace CIPLabeler
}  // namespace RDKit