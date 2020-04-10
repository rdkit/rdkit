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
 * A simple insertion sort for ligands. The number of ligands is not likely to
 * be very larger as such doing a merge sort would have little benefit.
 *
 */
class Sort {
private:
  const std::vector<const SequenceRule *> rules;

public:
  Sort(const SequenceRule *comparator);

  Sort(std::vector<const SequenceRule *> comparators);

  const std::vector<const SequenceRule *> &getRules() const;

  Priority prioritise(const Node *node, std::vector<Edge *> &edges) const;

  Priority prioritise(const Node *node, std::vector<Edge *> &edges,
                      bool deep) const;

  int compareLigands(const Node *node, const Edge *a, const Edge *b,
                     bool deep) const;

  void swap(std::vector<Edge *> &list, int i, int j) const;

  std::vector<std::vector<Edge *>>
  getGroups(const std::vector<Edge *> &sorted) const;
}; // namespace CIPLabeler

} // namespace CIPLabeler
} // namespace RDKit