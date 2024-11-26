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

#include "Descriptor.h"
#include "Mancude.h"
#include "Edge.h"

namespace RDKit {

class Atom;

namespace CIPLabeler {

class Digraph;

class Node {
 public:
  /**
   * Flag indicates whether the node has been expanded.
   */
  static const int EXPANDED = 0x1;

  /**
   * Flag indicates whether the node was duplicated
   * at a ring closure.
   */
  static const int RING_DUPLICATE = 0x2;

  /**
   * Flag indicates whether the node was duplicated
   * at a bond with order &gt; 1.
   */
  static const int BOND_DUPLICATE = 0x4;

  /**
   * Mask to check if a node is duplicated.
   */

  static const int DUPLICATE = RING_DUPLICATE | BOND_DUPLICATE;

  /**
   * Node was created for an implicit hydrogen,
   * the 'atom' value will be null.
   */
  static const int IMPL_HYDROGEN = 0x8;

  /**
   * Mask to check if a node is duplicated or created for an implicit H (not a
   * primary node).
   */
  static const int DUPLICATE_OR_H =
      RING_DUPLICATE | BOND_DUPLICATE | IMPL_HYDROGEN;

  Node() = delete;
  Node(const Node &) = delete;
  Node &operator=(const Node &) = delete;

  Node(Digraph *g, std::vector<char> &&visit, Atom *atom,
       boost::rational<int> &&frac, int dist, int flags);

  Digraph *getDigraph() const;

  Atom *getAtom() const;

  int getDistance() const;

  boost::rational<int> getAtomicNumFraction() const;

  int getAtomicNum() const;

  unsigned getMassNum() const;

  double getAtomicMass() const;

  Descriptor getAux() const;

  bool isSet(int mask) const;

  bool isDuplicate() const;

  bool isDuplicateOrH() const;

  bool isTerminal() const;

  bool isExpanded() const;

  bool isVisited(int idx) const;

  Node *newChild(int idx, Atom *atom) const;

  Node *newBondDuplicateChild(int idx, Atom *atom) const;

  Node *newRingDuplicateChild(int idx, Atom *atom) const;

  Node *newImplicitHydrogenChild() const;

  void add(Edge *e);

  void setAux(Descriptor desc);

  const std::vector<Edge *> &getEdges() const;

  std::vector<Edge *> getEdges(Atom *end) const;

  std::vector<Edge *> getNonTerminalOutEdges() const;

 private:
  Digraph *dp_g;
  Atom *dp_atom;
  int d_dist;
  boost::rational<int> d_atomic_num;
  double d_atomic_mass;
  Descriptor d_aux = Descriptor::NONE;
  int d_flags = 0x0;

  std::vector<Edge *> d_edges;

  std::vector<char> d_visit;

  Node *newTerminalChild(int idx, Atom *atom, int flags) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
