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
  int new_dist = flags & DUPLICATE ? d_visit[idx] : d_dist + 1;
  std::vector<char> new_visit;

  if (flags & BOND_DUPLICATE) {
    auto frac = dp_g->getMol().getFractionalAtomicNum(dp_atom);
    if (frac.denominator() > 1) {
      return &dp_g->addNode(std::move(new_visit), atom, std::move(frac),
                            new_dist, flags);
    }
  }

  auto atomic_num = atom ? atom->getAtomicNum() : 1;
  return &dp_g->addNode(std::move(new_visit), atom, atomic_num, new_dist,
                        flags);
}

Node::Node(Digraph *g, std::vector<char> &&visit, Atom *atom,
           boost::rational<int> &&frac, int dist, int flags)
    : dp_g{g},
      dp_atom{atom},
      d_dist{dist},
      d_atomic_num{std::move(frac)},
      d_flags{flags},
      d_visit{std::move(visit)} {
  if (d_flags & DUPLICATE) {
    d_edges.reserve(4);
    d_atomic_mass = 0.;
  } else {
    const auto &table = RDKit::PeriodicTable::getTable();
    auto atomic_number = getAtomicNum();
    auto isotope = getMassNum();
    if (isotope == 0u) {
      d_atomic_mass = table->getAtomicWeight(atomic_number);
    } else {
      d_atomic_mass = table->getMassForIsotope(atomic_number, isotope);
    }
  }
  if (d_visit.empty() || d_flags & DUPLICATE) {
    d_flags |= EXPANDED;
  }
}

Digraph *Node::getDigraph() const { return dp_g; }

Atom *Node::getAtom() const { return dp_atom; }

int Node::getDistance() const { return d_dist; }

boost::rational<int> Node::getAtomicNumFraction() const { return d_atomic_num; }

int Node::getAtomicNum() const {
  if (dp_atom == nullptr) {
    return 1;
  }
  return dp_atom->getAtomicNum();
};

unsigned Node::getMassNum() const {
  if (dp_atom == nullptr || isDuplicate()) {
    return 0u;
  }
  return dp_atom->getIsotope();
}

double Node::getAtomicMass() const { return d_atomic_mass; }

Descriptor Node::getAux() const { return d_aux; }

bool Node::isSet(int mask) const { return mask & d_flags; }

bool Node::isDuplicate() const { return d_flags & DUPLICATE; }

bool Node::isTerminal() const {
  return d_visit.empty() || (isExpanded() && d_edges.size() == 1);
}

bool Node::isExpanded() const { return d_flags & EXPANDED; }

bool Node::isVisited(int idx) const { return d_visit[idx] != 0; }

Node *Node::newChild(int idx, Atom *atom) const {
  auto new_visit = d_visit;
  new_visit[idx] = static_cast<char>(d_dist + 1);
  auto atomic_num = atom ? atom->getAtomicNum() : 1;
  return &dp_g->addNode(std::move(new_visit), atom, atomic_num, d_dist + 1, 0);
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

void Node::add(Edge *e) { d_edges.push_back(e); }

void Node::setAux(Descriptor desc) { d_aux = desc; }

const std::vector<Edge *> &Node::getEdges() const {
  if (!isExpanded()) {
    auto non_const_this = const_cast<Node *>(this);
    non_const_this->d_flags |= EXPANDED;
    dp_g->expand(non_const_this);
  }
  return d_edges;
}

std::vector<Edge *> Node::getEdges(Atom *end) const {
  std::vector<Edge *> res;
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

std::vector<Edge *> Node::getNonTerminalOutEdges() const {
  std::vector<Edge *> edges;
  for (auto &edge : getEdges()) {
    if (edge->isBeg(this) && !edge->getEnd()->isTerminal()) {
      edges.push_back(edge);
    }
  }
  return edges;
}

}  // namespace CIPLabeler
}  // namespace RDKit
