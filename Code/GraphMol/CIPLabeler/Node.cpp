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
#include <vector>

#include "Digraph.h"
#include "Edge.h"
#include "Node.h"
#include "CIPMol.h"

namespace RDKit {
namespace CIPLabeler {

Node *Node::newTerminalChild(int idx, Atom *atom, int flags) const {
  int new_dist = flags & DUPLICATE ? visit[idx] : dist + 1;
  auto new_visit = std::vector<char>{};

  if (flags & BOND_DUPLICATE) {
    auto frac = g->getMol()->getFractionalAtomicNum(this->atom);
    if (frac.denominator() > 1) {
      return new Node(g, std::move(new_visit), atom, std::move(frac), new_dist,
                      flags);
    }
  }

  auto frac = Fraction(atom ? atom->getAtomicNum() : 1);
  return new Node(g, std::move(new_visit), atom, std::move(frac), new_dist,
                  flags);
}

Node::Node(Digraph *g, std::vector<char> &&visit, Atom *atom, Fraction &&frac,
           int dist, int flags)
    : g{g}, atom{atom}, dist{dist}, d_atomic_num{std::move(frac)}, flags{flags},
      visit{std::move(visit)} {
  if (flags & DUPLICATE) {
    this->edges.reserve(4);
  };
  if (this->visit.empty() || flags & DUPLICATE) {
    this->flags |= EXPANDED;
  }
}

Node::~Node() {
  for (auto &edge : edges) {
    if (edge->isBeg(this)) {
      delete edge->getEnd();
      delete edge;
    }
  }
}

Digraph *Node::getDigraph() const { return g; }

Atom *Node::getAtom() const { return atom; }

int Node::getDistance() const { return dist; }

Fraction Node::getAtomicNumFraction() const { return d_atomic_num; }

int Node::getAtomicNum() const {
  if (atom == nullptr) {
    return 1;
  }
  return atom->getAtomicNum();
};

double Node::getMassNum() const {
  if (atom == nullptr) {
    return 0;
  }
  return atom->getIsotope();
}

Descriptor Node::getAux() const { return aux; }

bool Node::isSet(int mask) const { return mask & flags; }

bool Node::isDuplicate() const { return flags & DUPLICATE; }

bool Node::isTerminal() const {
  return visit.empty() || (isExpanded() && edges.size() == 1);
}

bool Node::isExpanded() const { return flags & EXPANDED; }

Node *Node::newChild(int idx, Atom *atom) const {
  auto new_visit = visit;
  new_visit[idx] = static_cast<char>(dist + 1);
  auto frac = Fraction(atom ? atom->getAtomicNum() : 1);
  return new Node(g, std::move(new_visit), atom, std::move(frac), dist + 1, 0);
}

Node *Node::newBondDuplicateChild(int idx, Atom *atom) const {
  return newTerminalChild(idx, atom, BOND_DUPLICATE);
}

Node *Node::newRingDuplicateChild(int idx, Atom *atom) const {
  return newTerminalChild(idx, atom, RING_DUPLICATE);
}

Node *Node::newImplicitHydrogenChild() const {
  return newTerminalChild(-1, nullptr, IMPL_HYDROGEN);
}

void Node::add(Edge *e) { edges.push_back(e); }

void Node::setAux(Descriptor desc) { aux = desc; }

std::vector<Edge *> Node::getEdges() const {
  if (!isExpanded()) {
    auto non_const_this = const_cast<Node *>(this);
    non_const_this->flags |= EXPANDED;
    g->expand(non_const_this);
  }
  return edges;
}

std::vector<Edge *> Node::getEdges(Atom *end) const {
  auto res = std::vector<Edge *>();
  for (auto &edge : getEdges()) {
    if (edge->getEnd()->isDuplicate()) {
      continue;
    };
    if (end == edge->getBeg()->getAtom() || end == edge->getEnd()->getAtom()) {
      res.push_back(edge);
    }
  }
  return res;
}

std::vector<Edge *> Node::getOutEdges() const {
  auto edges = std::vector<Edge *>();
  for (auto &edge : getEdges()) {
    if (edge->isBeg(this)) {
      edges.push_back(edge);
    }
  }
  return edges;
}

std::vector<Edge *> Node::getNonTerminalOutEdges() const {
  auto edges = std::vector<Edge *>();
  for (auto &edge : getEdges()) {
    if (edge->isBeg(this) && !edge->getEnd()->isTerminal())
      edges.push_back(edge);
  }
  return edges;
}

} // namespace CIPLabeler
} // namespace RDKit
