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

#include "SequenceRule.h"

namespace RDKit {
namespace CIPLabeler {

/**
 * A descriptor pair rule. This rule defines that like descriptor pairs have
 * priority over unlike descriptor pairs.
 *
 */
class Rule4b : public SequenceRule {
private:
  const Descriptor ref;

public:
  Rule4b() = delete;

  Rule4b(const CIPMol *mol);

  Rule4b(const CIPMol *mol, Descriptor ref);

  std::vector<Descriptor> getReferenceDescriptors(const Node *node) const;

  int compare(const Edge *a, const Edge *b) const override;

private:
  bool hasDescriptors(const Node *node) const;

  bool getReference(const std::vector<const Node *> &nodes,
                    std::vector<Descriptor> &result) const;

  std::vector<std::vector<const Node *>> initialLevel(const Node *node) const;

  std::vector<std::vector<const Node *>>
  getNextLevel(const std::vector<std::vector<const Node *>> &prevLevel) const;

  std::vector<const Node *>
  toNodeList(const std::vector<Edge *> &eqEdges) const;

  /**
   * Reduce the number of combinations by not including terminal ligands in
   * the permuting. They can't be stereocentres and so won't contribute the
   * the like / unlike std::vector.
   *
   * @param edges a std::vector of edges
   * @return a std::vector of non-terminal ligands
   */
  std::vector<const Edge *>
  getLigandsToSort(const Node *node, const std::vector<Edge *> edges) const;

  std::vector<PairList>
  newPairLists(const std::vector<Descriptor> &descriptors) const;

  void fillPairs(const Node *beg, PairList &plist) const;

  int comparePairs(const Node *a, const Node *b, Descriptor refA,
                   Descriptor refB) const;

  Sort getRefSorter(const SequenceRule *replacement_rule) const;
};

} // namespace CIPLabeler
} // namespace RDKit