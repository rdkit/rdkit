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

#include <RDGeneral/Invariant.h>

#include "Rule5New.h"

#include "../Digraph.h"

namespace RDKit {
namespace CIPLabeler {

Rule5New::Rule5New(const CIPMol *mol)
    : SequenceRule(mol), ref{Descriptor::NONE} {}

Rule5New::Rule5New(const CIPMol *mol, Descriptor ref)
    : SequenceRule(mol), ref{ref} {}

std::vector<Descriptor> Rule5New::getReferenceDescriptors(const Node *node) {
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

int Rule5New::compare(const Edge *a, const Edge *b) const {
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
    auto listRA = PairList(Descriptor::R);
    auto listRB = PairList(Descriptor::R);
    auto listSA = PairList(Descriptor::S);
    auto listSB = PairList(Descriptor::S);
    fillPairs(a->getEnd(), listRA);
    fillPairs(a->getEnd(), listSA);
    fillPairs(b->getEnd(), listRB);
    fillPairs(b->getEnd(), listSB);
    int cmpR = listRA.compareTo(listRB);
    int cmpS = listSA.compareTo(listSB);
    // -2/+2 for psuedo-asymetric
    // -1/+1 if not (e.g. the R > R and S > S lists)
    if (cmpR < 0) {
      return cmpS < 0 ? -1 : -2;
    } else if (cmpR > 0) {
      return cmpS > 0 ? +1 : +2;
    } else {
      return 0;
    }
  }
}

bool Rule5New::hasDescriptors(const Node *node) const {
  auto queue = std::vector<const Node *>({node});

  for (auto pos = 0u; pos < queue.size(); ++pos) {
    const auto n = queue[0];

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

bool Rule5New::getReference(const std::vector<const Node *> &nodes,
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

bool Rule5New::isPseudoAsymmetric() const { return true; }

std::vector<std::vector<const Node *>>
Rule5New::initialLevel(const Node *node) const {
  return {{node}};
}

std::vector<std::vector<const Node *>> Rule5New::getNextLevel(
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
Rule5New::toNodeList(const std::vector<Edge *> &eqEdges) const {
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
Rule5New::getLigandsToSort(const Node *node,
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
Rule5New::newPairLists(const std::vector<Descriptor> &descriptors) const {
  auto pairs = std::vector<PairList>();
  for (Descriptor descriptor : descriptors) {
    pairs.emplace_back(descriptor);
  }
  return pairs;
}

void Rule5New::fillPairs(const Node *beg, PairList &plist) const {
  const Rule5New replacement_rule(this->getMol(), plist.getRefDescriptor());
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

int Rule5New::comparePairs(const Node *a, const Node *b, Descriptor refA,
                           Descriptor refB) const {
  const Rule5New replacementA(this->getMol(), refA);
  const Rule5New replacementB(this->getMol(), refB);
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

Sort Rule5New::getRefSorter(const SequenceRule *replacement_rule) const {
  const auto &rules = this->getSorter()->getRules();

  CHECK_INVARIANT(std::find(rules.begin(), rules.end(), this) != rules.end(),
                  "Rule5New instance not in rule set");

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