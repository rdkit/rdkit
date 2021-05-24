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
#include <list>

#include <RDGeneral/Invariant.h>

#include "Rule4b.h"

#include "../Digraph.h"
#include "Pairlist.h"

namespace RDKit {
namespace CIPLabeler {

Rule4b::Rule4b() = default;

Rule4b::Rule4b(Descriptor ref) : d_ref{ref} {}

std::vector<Descriptor> Rule4b::getReferenceDescriptors(
    const Node *node) const {
  std::vector<Descriptor> result;
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
  const auto &aBeg = a->getBeg();
  const auto &aEnd = a->getEnd();
  const auto &bBeg = b->getBeg();
  const auto &bEnd = b->getEnd();
  if (aBeg->getDigraph()->getCurrentRoot() != aBeg ||
      bBeg->getDigraph()->getCurrentRoot() != bBeg) {
    if (d_ref == Descriptor::NONE) {
      return 0;
    }
    Descriptor aDesc = aEnd->getAux();
    Descriptor bDesc = bEnd->getAux();
    if (aDesc != Descriptor::NONE && bDesc != Descriptor::NONE &&
        aDesc != Descriptor::ns && bDesc != Descriptor::ns) {
      bool alike = PairList::ref(d_ref) == PairList::ref(aDesc);
      bool blike = PairList::ref(d_ref) == PairList::ref(bDesc);
      if (alike && !blike) {
        return +1;
      }
      if (blike && !alike) {
        return -1;
      }
    }
    return 0;
  } else {
    auto list1 = newPairLists(getReferenceDescriptors(aEnd));

    auto list2 = newPairLists(getReferenceDescriptors(bEnd));

    if (list1.empty() != list2.empty()) {
      throw std::runtime_error(
          "Substituents should be topologically equivalent!");
    }
    if (list1.size() == 1) {
      return comparePairs(aEnd, bEnd, list1[0].getRefDescriptor(),
                          list2[0].getRefDescriptor());
    } else if (list1.size() > 1) {
      for (auto &plist : list1) {
        fillPairs(aEnd, plist);
      }
      for (auto &plist : list2) {
        fillPairs(bEnd, plist);
      }

      std::sort(list1.rbegin(), list1.rend());
      std::sort(list2.rbegin(), list2.rend());

      for (auto i = 0u; i < list1.size(); ++i) {
        int cmp = list1[i].compareTo(list2[i]);
        if (cmp != 0) {
          return cmp;
        }
      }
    }
    return 0;
  }
}

bool Rule4b::hasDescriptors(const Node *node) const {
  auto queue = std::list<const Node *>({node});

  for (const auto &node : queue) {
    if (node->getAux() != Descriptor::NONE) {
      return true;
    }
    for (const auto &e : node->getEdges()) {
      if (e->getEnd() == node) {
        continue;
      }
      if (getBondLabel(e) != Descriptor::NONE) {
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

std::vector<std::vector<const Node *>> Rule4b::initialLevel(
    const Node *node) const {
  return {{node}};
}

std::vector<std::vector<const Node *>> Rule4b::getNextLevel(
    const std::vector<std::vector<const Node *>> &prevLevel) const {
  std::vector<std::vector<const Node *>> nextLevel;
  nextLevel.reserve(4 * prevLevel.size());

  for (const auto &prev : prevLevel) {
    std::vector<std::vector<std::vector<Edge *>>> tmp;
    for (const auto &node : prev) {
      auto edges = node->getNonTerminalOutEdges();
      sort(node, edges);
      tmp.push_back(getSorter()->getGroups(edges));
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
      std::vector<const Node *> eq;
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

std::vector<const Node *> Rule4b::toNodeList(
    const std::vector<Edge *> &eqEdges) const {
  std::vector<const Node *> eqNodes;
  eqNodes.reserve(eqEdges.size());
  for (const auto &edge : eqEdges) {
    eqNodes.push_back(edge->getEnd());
  }
  return eqNodes;
}

std::vector<PairList> Rule4b::newPairLists(
    const std::vector<Descriptor> &descriptors) const {
  std::vector<PairList> pairs;
  pairs.reserve(descriptors.size());
  for (Descriptor descriptor : descriptors) {
    pairs.emplace_back(descriptor);
  }
  return pairs;
}

void Rule4b::fillPairs(const Node *beg, PairList &plist) const {
  const Rule4b replacement_rule(plist.getRefDescriptor());
  const auto &sorter = getRefSorter(&replacement_rule);
  auto queue = std::list<const Node *>({beg});

  for (const auto &node : queue) {
    plist.add(node->getAux());
    auto edges = node->getEdges();
    sorter.prioritize(node, edges);
    for (const auto &edge : edges) {
      if (edge->isBeg(node) && !edge->getEnd()->isTerminal()) {
        queue.push_back(edge->getEnd());
      }
    }
  }
}

int Rule4b::comparePairs(const Node *a, const Node *b, Descriptor refA,
                         Descriptor refB) const {
  const Rule4b replacementA(refA);
  const Rule4b replacementB(refB);
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
    aSorter.prioritize(aNode, edges);
    for (const auto &edge : edges) {
      if (edge->isBeg(aNode) && !edge->getEnd()->isTerminal()) {
        aQueue.push_back(edge->getEnd());
      }
    }

    edges = bNode->getEdges();
    bSorter.prioritize(bNode, edges);
    for (const auto &edge : edges) {
      if (edge->isBeg(bNode) && !edge->getEnd()->isTerminal()) {
        bQueue.push_back(edge->getEnd());
      }
    }
  }
  return 0;
}

Sort Rule4b::getRefSorter(const SequenceRule *replacement_rule) const {
  const auto &rules = getSorter()->getRules();

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

}  // namespace CIPLabeler
}  // namespace RDKit