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

#include "Configuration.h"

#include "../Priority.hpp"

namespace RDKit {
namespace CIPLabeler {

Edge *Configuration::findInternalEdge(const std::vector<Edge *> &edges,
                                      Atom *f1, Atom *f2) {
  for (const auto &edge : edges) {
    if (edge->getBeg()->isDuplicate() || edge->getEnd()->isDuplicate()) {
      continue;
    }
    if (isInternalEdge(edge, f1, f2)) {
      return edge;
    }
  }
  return nullptr;
}

bool Configuration::isInternalEdge(const Edge *edge, Atom *f1, Atom *f2) {
  const auto &beg = edge->getBeg();
  const auto &end = edge->getEnd();
  if (f1 == beg->getAtom() && f2 == end->getAtom()) {
    return true;
  } else if (f1 == end->getAtom() && f2 == beg->getAtom()) {
    return true;
  }
  return false;
}

void Configuration::removeInternalEdges(std::vector<Edge *> &edges, Atom *f1,
                                        Atom *f2) {
  std::vector<Edge *> new_edges;
  for (auto &&e : edges) {
    if (!isInternalEdge(e, f1, f2)) {
      new_edges.push_back(std::move(e));
    }
  }
  std::swap(edges, new_edges);
}

void Configuration::removeDuplicatedEdges(std::vector<Edge *> &&edges) {
  std::vector<Edge *> new_edges;
  for (const auto &e : edges) {
    if (!e->getEnd()->isDuplicate()) {
      new_edges.push_back(e);
    }
  }
  std::swap(edges, new_edges);
}

Configuration::Configuration() = default;

Configuration::Configuration(Atom *focus, std::vector<Atom *> &&carriers,
                             int cfg)
    : foci{focus}, carriers{std::move(carriers)}, cfg{cfg} {};

Configuration::Configuration(std::vector<Atom *> &&foci,
                             std::vector<Atom *> &&carriers, int cfg)
    : foci{std::move(foci)}, carriers{std::move(carriers)}, cfg{cfg} {}

Configuration::~Configuration() = default;

Atom *Configuration::getFocus() const { return foci[0]; }

const std::vector<Atom *> &Configuration::getFoci() const { return foci; }

int Configuration::getConfig() const { return cfg; }

const std::vector<Atom *> &Configuration::getCarriers() const {
  return carriers;
}

std::shared_ptr<Digraph> Configuration::getDigraph() const {
  if (digraph == nullptr) {
    throw std::runtime_error("Digraph has not been set.");
  }
  return digraph;
}

void Configuration::setDigraph(Digraph *digraph) {
  this->digraph.reset(digraph);
}

Descriptor Configuration::label(Node *node, Digraph *digraph,
                                const SequenceRule *comp) {
  (void)node;
  (void)digraph;
  (void)comp;

  return Descriptor::UNKNOWN;
}

} // namespace CIPLabeler
} // namespace RDKit