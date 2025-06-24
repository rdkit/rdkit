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

#include "../Priority.h"
#include "../rules/Rules.h"

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

bool Configuration::isDuplicateOrHydrogenEdge(const Edge *edge) {
  return edge->getBeg()->isDuplicateOrH() || edge->getEnd()->isDuplicateOrH();
}

void Configuration::removeDuplicatesAndHs(std::vector<Edge *> &edges) {
  std::vector<Edge *> new_edges;
  for (auto &&e : edges) {
    if (!isDuplicateOrHydrogenEdge(e)) {
      new_edges.push_back(std::move(e));
    }
  }
  edges = std::move(new_edges);
}

void Configuration::setCarriers(std::vector<Atom *> &&carriers) {
  d_carriers = std::move(carriers);
}

Configuration::Configuration(const CIPMol &mol, Atom *focus)
    : d_foci{focus}, d_digraph{mol, focus} {};

Configuration::Configuration(const CIPMol &mol, std::vector<Atom *> &&foci,
                             bool atropisomerMode)
    : d_foci{std::move(foci)}, d_digraph{mol, d_foci[0], atropisomerMode} {}

Configuration::~Configuration() = default;

Atom *Configuration::getFocus() const { return d_foci[0]; }

const std::vector<Atom *> &Configuration::getFoci() const { return d_foci; }

const std::vector<Atom *> &Configuration::getCarriers() const {
  return d_carriers;
}

Digraph &Configuration::getDigraph() { return d_digraph; }

Descriptor Configuration::label(Node *node, Digraph &digraph,
                                const Rules &comp) {
  (void)node;
  (void)digraph;
  (void)comp;

  return Descriptor::UNKNOWN;
}

}  // namespace CIPLabeler
}  // namespace RDKit