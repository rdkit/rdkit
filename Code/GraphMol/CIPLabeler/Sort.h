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
#pragma once

#include <vector>
#include "Priority.hpp"

namespace RDKit {
namespace CIPLabeler {

class SequenceRule;
class Edge;
class Node;

/**
 * A simple insertion sort for substituents. The number of substituents is not
 * likely to
 * be very larger as such doing a merge sort would have little benefit.
 *
 */
class Sort {

public:
  Sort(const SequenceRule *comparator);

  Sort(std::vector<const SequenceRule *> comparators);

  const std::vector<const SequenceRule *> &getRules() const;

  Priority prioritise(const Node *node, std::vector<Edge *> &edges,
                      bool deep = true) const;

  std::vector<std::vector<Edge *>>
  getGroups(const std::vector<Edge *> &sorted) const;

private:
  const std::vector<const SequenceRule *> d_rules;

  int compareSubstituents(const Node *node, const Edge *a, const Edge *b,
                          bool deep) const;
}; // namespace CIPLabeler

} // namespace CIPLabeler
} // namespace RDKit