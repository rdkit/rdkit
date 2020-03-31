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

#include "Descriptor.hpp"
#include "Digraph.hpp"

namespace RDKit {
namespace NewCIPLabelling {

template <typename A, typename B> class Edge;

template <typename A, typename B> class Node {
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
  static const int DUPLICATE = 0x6;

  /**
   * Node was created for an implicit hydrogen,
   * the 'atom' value will be null.
   */
  static const int IMPL_HYDROGEN = 0x8;

private:
  Digraph<A, B> *g;
  A atom;
  int dist;
  short atomic_num_numerator;
  short atomic_num_denominator;
  Descriptor aux = Descriptor::NONE;
  int flags = 0x0;

  std::vector<Edge<A, B> *> edges;

  Node<A, B> *newTerminalChild(int idx, A atom, int flags) const {
    int new_dist =
        static_cast<char>(((flags & DUPLICATE) != 0 ? visit[idx] : dist + 1));
    int anum;
    int aden;
    if ((flags & BOND_DUPLICATE) != 0) {
      auto frac = g->getMol()->getFractionalAtomicNum(this->atom);
      if (frac.denominator() > 1) {
        anum = frac.numerator();
        aden = frac.denominator();
      } else {
        anum = g->getMol()->getAtomicNum(atom);
        aden = 1;
      }
    } else {
      anum = g->getMol()->getAtomicNum(atom);
      aden = 1;
    }

    auto new_visit = std::vector<char>();
    return new Node(g, std::move(new_visit), atom, anum, aden, new_dist, flags);
  }

public:
  std::vector<char> visit;

  Node() = delete;
  Node(const Node &) = delete;
  Node &operator=(const Node &) = delete;

  Node(Digraph<A, B> *g, std::vector<char> &&visit, A atom, short anum,
       short aden, int dist, int flags)
      : g{g}, atom{atom}, dist{dist}, atomic_num_numerator{anum},
        atomic_num_denominator{aden}, flags{flags}, visit{std::move(visit)} {
    if ((flags & DUPLICATE) != 0) {
      this->edges.reserve(4);
    };
    if (this->visit.empty() || (flags & DUPLICATE) != 0) {
      this->flags |= EXPANDED;
    }
  }

  ~Node() {
    for (auto &edge : edges) {
      if (edge->isBeg(this)) {
        delete edge->getEnd();
        delete edge;
      }
    }
  }

  Digraph<A, B> *getDigraph() const { return g; }

  A getAtom() const { return atom; }

  int getDistance() const { return dist; }

  int getAtomicNumNumerator() const { return atomic_num_numerator; }

  int getAtomicNumDenominator() const { return atomic_num_denominator; }

  Descriptor getAux() const { return aux; }

  bool isSet(int mask) const { return (mask & flags) != 0; }

  bool isDuplicate() const { return (flags & DUPLICATE) != 0; }

  bool isTerminal() const {
    return visit.empty() || ((flags & EXPANDED) != 0 && edges.size() == 1);
  }

  bool isExpanded() const { return (flags & EXPANDED) != 0; }

  Node<A, B> *newChild(int idx, A atom) const {
    auto new_visit = visit;
    new_visit[idx] = static_cast<char>(dist + 1);
    return new Node<A, B>(g, std::move(new_visit), atom,
                          g->getMol()->getAtomicNum(atom), 1, dist + 1, 0);
  }

  Node<A, B> *newBondDuplicateChild(int idx, A atom) const {
    return newTerminalChild(idx, atom, BOND_DUPLICATE);
  }

  Node<A, B> *newRingDuplicateChild(int idx, A atom) const {
    return newTerminalChild(idx, atom, RING_DUPLICATE);
  }

  Node<A, B> *newImplicitHydrogenChild() const {
    return newTerminalChild(-1, nullptr, IMPL_HYDROGEN);
  }

  void add(Edge<A, B> *e) { edges.push_back(e); }

  void setAux(Descriptor desc) { aux = desc; }

  std::vector<Edge<A, B> *> getEdges() const {
    if ((flags & EXPANDED) == 0) {
      auto non_const_this = const_cast<Node<A, B> *>(this);
      non_const_this->flags |= EXPANDED;
      g->expand(non_const_this);
    }
    return edges;
  }

  std::vector<Edge<A, B> *> getEdges(A end) const {
    auto res = std::vector<Edge<A, B> *>();
    for (auto &edge : getEdges()) {
      if (edge->getEnd()->isDuplicate()) {
        continue;
      };
      if (end == edge->getBeg()->getAtom() ||
          end == edge->getEnd()->getAtom()) {
        res.push_back(edge);
      }
    }
    return res;
  }

  std::vector<Edge<A, B> *> getOutEdges() const {
    auto edges = std::vector<Edge<A, B> *>();
    for (auto &edge : getEdges()) {
      if (edge->isBeg(this)) {
        edges.push_back(edge);
      }
    }
    return edges;
  }

  std::vector<Edge<A, B> *> getNonTerminalOutEdges() const {
    auto edges = std::vector<Edge<A, B> *>();
    for (auto &edge : getEdges()) {
      if (edge->isBeg(this) && !edge->getEnd()->isTerminal())
        edges.push_back(edge);
    }
    return edges;
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit
