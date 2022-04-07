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

#include "SequenceRule.h"

#include "../CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

SequenceRule::SequenceRule() = default;

SequenceRule::~SequenceRule() = default;

Descriptor SequenceRule::getBondLabel(const Edge *edge) const {
  Bond *bond = edge->getBond();
  if (bond == nullptr) {
    return Descriptor::NONE;
  }
  Descriptor label = edge->getAux();
  if (label != Descriptor::NONE) {
    return label;
  }
  return label;
}

int SequenceRule::getComparision(const Edge *a, const Edge *b) const {
  return getComparision(a, b, true);
}

int SequenceRule::getComparision(const Edge *a, const Edge *b,
                                 bool deep) const {
  return deep ? recursiveCompare(a, b) : compare(a, b);
}

const Sort *SequenceRule::getSorter() const {
  if (dp_sorter == nullptr) {
    const_cast<SequenceRule *>(this)->setSorter(new Sort(this));
  }
  return dp_sorter.get();
}

int SequenceRule::recursiveCompare(const Edge *a, const Edge *b) const {
  int cmp = compare(a, b);
  if (cmp != 0) {
    return cmp;
  }

  auto aQueue = std::vector<const Edge *>({a});
  auto bQueue = std::vector<const Edge *>({b});

  for (auto pos = 0u; pos < aQueue.size() && pos < bQueue.size(); ++pos) {
    a = aQueue[pos];
    b = bQueue[pos];
    auto as = a->getEnd()->getEdges();
    auto bs = b->getEnd()->getEdges();

    // shallow sort first of all
    sort(a->getEnd(), as, false);
    sort(b->getEnd(), bs, false);

    int sizediff = three_way_comparison(static_cast<int>(as.size()),
                                        static_cast<int>(bs.size()));

    {
      auto aIt = as.begin();
      auto bIt = bs.begin();
      for (; aIt != as.end() && bIt != bs.end(); ++aIt, ++bIt) {
        Node *aNode = a->getEnd();
        Node *bNode = b->getEnd();
        Edge *aEdge = *aIt;
        Edge *bEdge = *bIt;

        if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
          continue;
        }

        cmp = compare(aEdge, bEdge);
        if (cmp != 0) {
          return cmp;
        }
      }
    }

    if (sizediff != 0) {
      return sizediff;
    }

    sort(a->getEnd(), as);
    sort(b->getEnd(), bs);

    {
      auto aIt = as.begin();
      auto bIt = bs.begin();
      for (; aIt != as.end() && bIt != bs.end(); ++aIt, ++bIt) {
        Node *aNode = a->getEnd();
        Node *bNode = b->getEnd();
        Edge *aEdge = *aIt;
        Edge *bEdge = *bIt;

        if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
          continue;
        }

        cmp = compare(aEdge, bEdge);
        if (cmp != 0) {
          return cmp;
        }

        aQueue.push_back(aEdge);
        bQueue.push_back(bEdge);
      }
    }
  }
  return 0;
}

void SequenceRule::setSorter(const Sort *sorter) { dp_sorter.reset(sorter); }

Priority SequenceRule::sort(const Node *node, std::vector<Edge *> &edges,
                            bool deep) const {
  return getSorter()->prioritize(node, edges, deep);
}

Priority SequenceRule::sort(const Node *node,
                            std::vector<Edge *> &edges) const {
  return sort(node, edges, true);
}

bool SequenceRule::areUpEdges(Node *aNode, Node *bNode, Edge *aEdge,
                              Edge *bEdge) const {
  // step over 'up' edges
  if (aEdge->isEnd(aNode)) {
    // if b is 'down' something's not right!
    if (!bEdge->isEnd(bNode)) {
      throw std::runtime_error("Something unexpected!");
    }
    return true;
  }
  return false;
}

}  // namespace CIPLabeler
}  // namespace RDKit
