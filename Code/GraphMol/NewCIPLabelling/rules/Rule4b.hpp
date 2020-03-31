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

#include <cassert>
#include <vector>

#include "../Edge.hpp"
#include "../Node.hpp"
#include "Pairlist.hpp"
#include "SequenceRule.hpp"

namespace RDKit {
namespace NewCIPLabelling {

/**
 * A descriptor pair rule. This rule defines that like descriptor pairs have
 * priority over unlike descriptor pairs.
 *
 */
template <typename A, typename B> class Rule4b : public SequenceRule<A, B> {
private:
  const Descriptor ref;

public:
  Rule4b() = delete;

  Rule4b(const BaseMol<A, B> *mol)
      : SequenceRule<A, B>(mol), ref{Descriptor::NONE} {}

  Rule4b(const BaseMol<A, B> *mol, Descriptor ref)
      : SequenceRule<A, B>(mol), ref{ref} {}

  std::vector<Descriptor>
  getReferenceDescriptors(const Node<A, B> *node) const {
    auto result = std::vector<Descriptor>();
    auto prev = initialLevel(node);
    while (!prev.empty()) {
      for (const auto &nodes : prev) {
        if (getReference(nodes, result)) {
          return result;
        }
      }
      prev = getNextLevel(prev);
    }
    return {};
  }

  int compare(const Edge<A, B> *a, const Edge<A, B> *b) const override {
    if (a->getBeg()->getDigraph()->getCurrRoot() != a->getBeg() ||
        b->getBeg()->getDigraph()->getCurrRoot() != b->getBeg()) {
      if (ref == Descriptor::NONE) {
        return 0;
      }
      Descriptor aDesc = a->getEnd()->getAux();
      Descriptor bDesc = b->getEnd()->getAux();
      if (aDesc != Descriptor::NONE && bDesc != Descriptor::NONE &&
          aDesc != Descriptor::ns && bDesc != Descriptor::ns) {
        bool alike = PairList::ref(ref) == PairList::ref(aDesc);
        bool blike = PairList::ref(ref) == PairList::ref(bDesc);
        if (alike && !blike) {
          return +1;
        }
        if (blike && !alike) {
          return -1;
        }
      }
      return 0;
    } else {
      auto list1 = newPairLists(getReferenceDescriptors(a->getEnd()));

      auto list2 = newPairLists(getReferenceDescriptors(b->getEnd()));

      if (list1.empty() != list2.empty()) {
        throw std::runtime_error("Ligands should be topologically equivalent!");
      }
      if (list1.size() == 1) {
        return comparePairs(a->getEnd(), b->getEnd(),
                            list1[0].getRefDescriptor(),
                            list2[0].getRefDescriptor());
      } else if (list1.size() > 1) {
        for (auto &plist : list1) {
          fillPairs(a->getEnd(), plist);
        }
        for (auto &plist : list2) {
          fillPairs(b->getEnd(), plist);
        }

        std::sort(list1.rbegin(), list1.rend());
        std::sort(list2.rbegin(), list2.rend());

        for (auto i = 0u; i < list1.size(); ++i) {
          int cmp = list1[0].compareTo(list2[0]);
          if (cmp != 0) {
            return cmp;
          }
        }
      }
      return 0;
    }
  }

private:
  bool hasDescriptors(const Node<A, B> *node) const {
    auto queue = std::vector<const Node<A, B> *>({node});

    for (auto pos = 0u; pos < queue.size(); ++pos) {
      const auto n = queue[pos];

      if (n->getAux() != nullptr) {
        return true;
      }
      for (const auto &e : n->getEdges()) {
        if (e->getEnd() == n) {
          continue;
        }
        if (this->getBondLabel(e) != Descriptor::NONE) {
          return true;
        }
        queue.push_back(e->getEnd());
      }
    }
    return false;
  }

  bool getReference(const std::vector<const Node<A, B> *> &nodes,
                    std::vector<Descriptor> &result) const {
    int right = 0;
    int left = 0;
    for (const auto &node : nodes) {
      auto desc = node->getAux();
      switch (desc) {
      case Descriptor::NONE:
        continue;
      case Descriptor::R:
      case Descriptor::M:
      case Descriptor::seqCis:
        ++right;
        break;
      case Descriptor::S:
      case Descriptor::P:
      case Descriptor::seqTrans:
        ++left;
        break;
      default:
        break;
      }
    }
    if (right + left == 0) {
      return false;
    } else if (right > left) {
      result.push_back(Descriptor::R);
      return true;
    } else if (right < left) {
      result.push_back(Descriptor::S);
      return true;
    } else {
      result.push_back(Descriptor::R);
      result.push_back(Descriptor::S);
      return true;
    }
  }

  std::vector<std::vector<const Node<A, B> *>>
  initialLevel(const Node<A, B> *node) const {
    return {{{node}}};
  }

  std::vector<std::vector<const Node<A, B> *>> getNextLevel(
      const std::vector<std::vector<const Node<A, B> *>> &prevLevel) const {
    auto nextLevel = std::vector<std::vector<const Node<A, B> *>>();
    nextLevel.reserve(4 * prevLevel.size());

    for (const auto &prev : prevLevel) {
      auto tmp = std::vector<std::vector<std::vector<Edge<A, B> *>>>();
      for (const auto &node : prev) {
        auto edges = node->getNonTerminalOutEdges();
        this->sort(node, edges);
        tmp.push_back(this->getSorter()->getGroups(edges));
      }

      // check sizes
      int size = -1;
      for (auto i = 0u; i < tmp.size(); ++i) {
        int localSize = tmp[0].size();
        if (size < 0) {
          size = localSize;
        } else if (size != localSize) {
          throw std::runtime_error("Something unexpected!");
        }
      }

      for (int i = 0; i < size; ++i) {
        auto eq = std::vector<const Node<A, B> *>();
        for (const auto &aTmp : tmp) {
          auto tmpNodes = toNodeList(aTmp[i]);
          eq.insert(eq.end(), tmpNodes.begin(), tmpNodes.end());
        }
        if (!eq.empty()) {
          nextLevel.push_back(eq);
        }
      }
    }
    return nextLevel;
  }

  std::vector<const Node<A, B> *>
  toNodeList(const std::vector<Edge<A, B> *> &eqEdges) const {
    auto eqNodes = std::vector<const Node<A, B> *>();
    eqNodes.reserve(eqEdges.size());

    for (const auto &edge : eqEdges) {
      eqNodes.push_back(edge->getEnd());
    }
    return eqNodes;
  }

  void visit(std::vector<PairList> &plists, const Node<A, B> *beg) const {
    auto queue = std::vector<Node<A, B> *>({beg});

    for (auto pos = 0u; pos < queue.sizd(); ++pos) {
      const auto node = queue[pos];

      // append any descriptors to the std::vector
      for (const auto &plist : plists) {
        plist.add(node->getAux());
      }

      for (const auto &e : node->getEdges()) {
        if (e->isBeg(node)) {
          queue.push_back(e->getEnd());
        }
      }
    }
  }

  /**
   * Reduce the number of combinations by not including terminal ligands in
   * the permuting. They can't be stereocentres and so won't contribute the
   * the like / unlike std::vector.
   *
   * @param edges a std::vector of edges
   * @return a std::vector of non-terminal ligands
   */
  std::vector<const Edge<A, B> *>
  getLigandsToSort(const Node<A, B> *node,
                   const std::vector<Edge<A, B> *> edges) const {
    auto filtered = std::vector<const Edge<A, B> *>();
    for (const auto &edge : edges) {
      if (edge->isEnd(node) || edge->getEnd()->isTerminal()) {
        continue;
      }
      if (!hasDescriptors(node)) {
        continue;
      }
      filtered.push_back(edge);
    }
    return filtered;
  }

  std::vector<PairList>
  newPairLists(const std::vector<Descriptor> &descriptors) const {
    auto pairs = std::vector<PairList>();
    for (Descriptor descriptor : descriptors) {
      pairs.emplace_back(descriptor);
    }
    return pairs;
  }

  void fillPairs(const Node<A, B> *beg, PairList &plist) const {
    const Rule4b<A, B> replacement_rule(this->getMol(),
                                        plist.getRefDescriptor());
    const auto &sorter = getRefSorter(&replacement_rule);
    auto queue = std::vector<const Node<A, B> *>({beg});

    for (auto pos = 0u; pos < queue.size(); ++pos) {
      const auto node = queue[pos];

      plist.add(node->getAux());
      auto edges = node->getEdges();
      sorter.prioritise(node, edges);
      for (const auto &edge : edges) {
        if (edge->isBeg(node) && !edge->getEnd()->isTerminal()) {
          queue.push_back(edge->getEnd());
        }
      }
    }
  }

  int comparePairs(const Node<A, B> *a, const Node<A, B> *b, Descriptor refA,
                   Descriptor refB) const {
    const Rule4b<A, B> replacementA(this->getMol(), refA);
    const Rule4b<A, B> replacementB(this->getMol(), refB);
    const auto &aSorter = getRefSorter(&replacementA);
    const auto &bSorter = getRefSorter(&replacementB);
    auto aQueue = std::vector<const Node<A, B> *>({a});
    auto bQueue = std::vector<const Node<A, B> *>({b});

    for (auto pos = 0u; pos < aQueue.size() && pos < bQueue.size(); ++pos) {
      const auto aNode = aQueue[pos];
      const auto bNode = bQueue[pos];

      const auto &desA = PairList::ref(aNode->getAux());
      const auto &desB = PairList::ref(bNode->getAux());

      if (desA == refA && desB != refB) {
        return +1;
      } else if (desA != refA && desB == refB) {
        return -1;
      }

      auto edges = aNode->getEdges();
      aSorter.prioritise(aNode, edges);
      for (const auto &edge : edges) {
        if (edge->isBeg(aNode) && !edge->getEnd()->isTerminal()) {
          aQueue.push_back(edge->getEnd());
        }
      }

      edges = bNode->getEdges();
      bSorter.prioritise(bNode, edges);
      for (const auto &edge : edges) {
        if (edge->isBeg(bNode) && !edge->getEnd()->isTerminal()) {
          bQueue.push_back(edge->getEnd());
        }
      }
    }
    return 0;
  }

  Sort<A, B> getRefSorter(const SequenceRule<A, B> *replacement_rule) const {
    const auto &rules = this->getSorter()->getRules();

    assert(std::find(rules.begin(), rules.end(), this) != rules.end());

    std::vector<const SequenceRule<A, B> *> new_rules;
    new_rules.reserve(rules.size());
    for (const auto &rule : rules) {
      if (this != rule) {
        new_rules.push_back(rule);
      }
    }
    new_rules.push_back(replacement_rule);
    return {new_rules};
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit