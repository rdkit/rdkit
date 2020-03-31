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

#include <memory>
#include <queue>

#include "../BaseMol.hpp"
#include "../Edge.hpp"
#include "../Node.hpp"
#include "Sort.hpp"

namespace RDKit {
namespace NewCIPLabelling {

namespace {
int integer_compare(int x, int y) { return x < y ? -1 : (x == y ? 0 : 1); }
} // namespace

template <typename A, typename B> class SequenceRule {
private:
  const BaseMol<A, B> *mol;

protected:
  std::unique_ptr<const Sort<A, B>> sorter = nullptr;

public:
  SequenceRule() = delete;

  SequenceRule(const BaseMol<A, B> *mol) : mol{mol} {}

  virtual ~SequenceRule() = default;

  virtual const BaseMol<A, B> *getMol() const { return mol; }

  Descriptor getBondLabel(const Edge<A, B> *edge) const {
    B bond = edge->getBond();
    if (bond == nullptr) {
      return Descriptor::NONE;
    }
    Descriptor label = edge->getAux();
    if (label != Descriptor::NONE) {
      return label;
    }
    return label;
  }

  virtual int getNumSubRules() const { return 1; }

  virtual bool isPseudoAsymmetric() const { return false; }

  int getComparision(const Edge<A, B> *a, const Edge<A, B> *b) const {
    return getComparision(a, b, true);
  }

  virtual int getComparision(const Edge<A, B> *a, const Edge<A, B> *b,
                             bool deep) const {
    return deep ? recursiveCompare(a, b) : compare(a, b);
  }

  virtual const Sort<A, B> *getSorter() const {
    if (sorter == nullptr) {
      const_cast<SequenceRule<A, B> *>(this)->setSorter(new Sort<A, B>(this));
    }
    return sorter.get();
  }

  int recursiveCompare(const Edge<A, B> *a, const Edge<A, B> *b) const {
    // pseudo atoms (atomic no. 0) match all
    if (mol->getAtomicNum(a->getEnd()->getAtom()) == 0 ||
        mol->getAtomicNum(b->getEnd()->getAtom()) == 0) {
      return 0;
    }

    int cmp = compare(a, b);
    if (cmp != 0) {
      return cmp;
    }

    auto aQueue = std::queue<const Edge<A, B> *>();
    auto bQueue = std::queue<const Edge<A, B> *>();

    aQueue.push(a);
    bQueue.push(b);

    while (!aQueue.empty() && !bQueue.empty()) {
      a = aQueue.front();
      aQueue.pop();
      b = bQueue.front();
      bQueue.pop();
      auto as = a->getEnd()->getEdges();
      auto bs = b->getEnd()->getEdges();

      // shallow sort first of all
      sort(a->getEnd(), as, false);
      sort(b->getEnd(), bs, false);

      int sizediff = integer_compare(static_cast<int>(as.size()),
                                     static_cast<int>(bs.size()));

      {
        auto aIt = as.begin();
        auto bIt = bs.begin();
        for (; aIt != as.end() && bIt != bs.end(); ++aIt, ++bIt) {
          Node<A, B> *aNode = a->getEnd();
          Node<A, B> *bNode = b->getEnd();
          Edge<A, B> *aEdge = *aIt;
          Edge<A, B> *bEdge = *bIt;

          if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
            continue;
          }

          // pseudo atoms (atomic no. 0) match all
          if (mol->getAtomicNum(aEdge->getEnd()->getAtom()) == 0 ||
              mol->getAtomicNum(bEdge->getEnd()->getAtom()) == 0) {
            return 0;
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
          Node<A, B> *aNode = a->getEnd();
          Node<A, B> *bNode = b->getEnd();
          Edge<A, B> *aEdge = *aIt;
          Edge<A, B> *bEdge = *bIt;

          if (areUpEdges(aNode, bNode, aEdge, bEdge)) {
            continue;
          }

          // pseudo atoms (atomic no. 0) match all
          if (mol->getAtomicNum(aEdge->getEnd()->getAtom()) == 0 ||
              mol->getAtomicNum(bEdge->getEnd()->getAtom()) == 0) {
            return 0;
          }

          cmp = compare(aEdge, bEdge);
          if (cmp != 0) {
            return cmp;
          }

          aQueue.push(aEdge);
          bQueue.push(bEdge);
        }
      }
    }
    return 0;
  }

  void setSorter(const Sort<A, B> *sorter) { this->sorter.reset(sorter); }

  Priority sort(const Node<A, B> *node, std::vector<Edge<A, B> *> &edges,
                bool deep) const {
    return getSorter()->prioritise(node, edges, deep);
  }

  Priority sort(const Node<A, B> *node,
                std::vector<Edge<A, B> *> &edges) const {
    return sort(node, edges, true);
  }

  virtual int compare(const Edge<A, B> *a, const Edge<A, B> *b) const = 0;

private:
  bool areUpEdges(Node<A, B> *aNode, Node<A, B> *bNode, Edge<A, B> *aEdge,
                  Edge<A, B> *bEdge) const {
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
};

} // namespace NewCIPLabelling
} // namespace RDKit