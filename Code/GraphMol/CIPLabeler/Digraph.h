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

#include <vector>

#include "TooManyNodesException.hpp"

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
private:
  CIPMol *mol = nullptr;
  Node *root = nullptr;
  Node *tmproot = nullptr;
  int numNodes = 0;
  Atom *rule6Ref = nullptr;

  void addEdge(Node *beg, Bond *bond, Node *end);

public:
  Digraph();
  Digraph(CIPMol *mol);
  Digraph(CIPMol *mol, Atom *atom);

  ~Digraph();

  CIPMol *getMol() const;

  Node *getRoot() const;

  Node *getCurrRoot() const;

  int getNumNodes() const;

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

  Node *init(Atom *atom);

  /**
   * Sets the root node of this digraph by flipping the directions
   * of edges as required.
   *
   * @param newroot the new root
   */
  void changeRoot(Node *newroot);

  void expand(Node *beg);

  void expandAll();
};

} // namespace CIPLabeler
} // namespace RDKit
