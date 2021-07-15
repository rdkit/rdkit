//
// Digraph is the core data structure for determining
// Cahn–Ingold–Prelog (CIP) chirality of a molecule.
//
// It's a "directed graph" - meaning that each bond
// has a start and an end. For CIP determination,
// the start points back towards the atom that is
// being labelled.
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

#include <list>
#include <vector>

#include <boost/rational.hpp>

#include "TooManyNodesException.h"

namespace RDKit {

class Atom;
class Bond;

namespace CIPLabeler {

class Node;
class Edge;
class CIPMol;

/**
 * A class to hold directed acyclic graphs representing the molecule.
 *
 * The root of the DAG is one of the foci of the configuration for
 * which the label is being calculated. The tmproot may be set to
 * other nodes that may become relevant in the calculation.
 *
 */
class Digraph {
 public:
  Digraph() = delete;
  Digraph(const Digraph &) = delete;
  Digraph &operator=(const Digraph &) = delete;

  Digraph(const CIPMol &mol, Atom *atom);

  const CIPMol &getMol() const;

  Node *getOriginalRoot() const;

  Node *getCurrentRoot() const;

  int getNumNodes() const;

  /**
   * Get all nodes which refer to `atom` in order of
   * distance from the root.
   */
  std::vector<Node *> getNodes(Atom *atom) const;

  /**
   * Access the reference atom for Rule 6 (if one is set).
   */
  Atom *getRule6Ref() const;

  /**
   * Used exclusively for Rule 6, we set one atom as the reference.
   * @param ref reference atom
   */
  void setRule6Ref(Atom *ref);

  /**
   * Sets the root node of this digraph by flipping the directions
   * of edges as required.
   *
   * This is more efficient than building a new Digraph, but is
   * only valid for neighboring Nodes.
   *
   * @param newroot the new root
   */
  void changeRoot(Node *newroot);

  void expand(Node *beg);

  Node &addNode(std::vector<char> &&visit, Atom *atom,
                boost::rational<int> &&frac, int dist, int flags);

 private:
  const CIPMol &d_mol;

  // The node from which the Digraph is first initialized.
  // It matches the atom that is being labeled.
  Node *dp_origin = nullptr;

  // The current root of the Digraph
  Node *dp_root = nullptr;

  Atom *dp_rule6Ref = nullptr;

  // We can't store these in a vector, as adding new items will
  // cause it to reallocate and invalidate the references
  std::list<Node> d_nodes;
  std::list<Edge> d_edges;

  void addEdge(Node *beg, Bond *bond, Node *end);
};

}  // namespace CIPLabeler
}  // namespace RDKit
