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

#include <algorithm>

#include <RDGeneral/Invariant.h>

#include "Rule4b.h"

#include "../Digraph.h"
#include "Pairlist.hpp"

namespace RDKit {
namespace CIPLabeler {

Rule4b::Rule4b(const CIPMol *mol) : SequenceRule(mol), ref{Descriptor::NONE} {}

Rule4b::Rule4b(const CIPMol *mol, Descriptor ref)
    : SequenceRule(mol), ref{ref} {}

std::vector<Descriptor>
Rule4b::getReferenceDescriptors(const Node *node) const {
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

int Rule4b::compare(const Edge *a, const Edge *b) const {
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
      return comparePairs(a->getEnd(), b->getEnd(), list1[0].getRefDescriptor(),
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

bool Rule4b::hasDescriptors(const Node *node) const {
  auto queue = std::vector<const Node *>({node});

  for (auto pos = 0u; pos < queue.size(); ++pos) {
    const auto n = queue[pos];

    if (n->getAux() != Descriptor::NONE) {
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

bool Rule4b::getReference(const std::vector<const Node *> &nodes,
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

std::vector<std::vector<const Node *>>
Rule4b::initialLevel(const Node *node) const {
  return {{node}};
}

std::vector<std::vector<const Node *>> Rule4b::getNextLevel(
    const std::vector<std::vector<const Node *>> &prevLevel) const {
  auto nextLevel = std::vector<std::vector<const Node *>>();
  nextLevel.reserve(4 * prevLevel.size());

  for (const auto &prev : prevLevel) {
    auto tmp = std::vector<std::vector<std::vector<Edge *>>>();
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
      auto eq = std::vector<const Node *>();
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

std::vector<const Node *>
Rule4b::toNodeList(const std::vector<Edge *> &eqEdges) const {
  auto eqNodes = std::vector<const Node *>();
  eqNodes.reserve(eqEdges.size());

  for (const auto &edge : eqEdges) {
    eqNodes.push_back(edge->getEnd());
  }
  return eqNodes;
}

/**
 * Reduce the number of combinations by not including terminal ligands in
 * the permuting. They can't be stereocentres and so won't contribute the
 * the like / unlike std::vector.
 *
 * @param edges a std::vector of edges
 * @return a std::vector of non-terminal ligands
 */
std::vector<const Edge *>
Rule4b::getLigandsToSort(const Node *node,
                         const std::vector<Edge *> edges) const {
  auto filtered = std::vector<const Edge *>();
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
Rule4b::newPairLists(const std::vector<Descriptor> &descriptors) const {
  auto pairs = std::vector<PairList>();
  for (Descriptor descriptor : descriptors) {
    pairs.emplace_back(descriptor);
  }
  return pairs;
}

void Rule4b::fillPairs(const Node *beg, PairList &plist) const {
  const Rule4b replacement_rule(this->getMol(), plist.getRefDescriptor());
  const auto &sorter = getRefSorter(&replacement_rule);
  auto queue = std::vector<const Node *>({beg});

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

int Rule4b::comparePairs(const Node *a, const Node *b, Descriptor refA,
                         Descriptor refB) const {
  const Rule4b replacementA(this->getMol(), refA);
  const Rule4b replacementB(this->getMol(), refB);
  const auto &aSorter = getRefSorter(&replacementA);
  const auto &bSorter = getRefSorter(&replacementB);
  auto aQueue = std::vector<const Node *>({a});
  auto bQueue = std::vector<const Node *>({b});

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

Sort Rule4b::getRefSorter(const SequenceRule *replacement_rule) const {
  const auto &rules = this->getSorter()->getRules();

  CHECK_INVARIANT(std::find(rules.begin(), rules.end(), this) != rules.end(),
                  "Rule4b instance not in rule set");

  std::vector<const SequenceRule *> new_rules;
  new_rules.reserve(rules.size());
  for (const auto &rule : rules) {
    if (this != rule) {
      new_rules.push_back(rule);
    }
  }
  new_rules.push_back(replacement_rule);
  return {new_rules};
}

} // namespace CIPLabeler
} // namespace RDKit