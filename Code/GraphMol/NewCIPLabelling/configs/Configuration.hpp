//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include "../Descriptor.hpp"
#include "../rules/Priority.hpp"
#include <memory>
#include <vector>

namespace RDKit {
namespace NewCIPLabelling {

template <typename A, typename B> class BaseMol;

template <typename A, typename B> class Digraph;

template <typename A, typename B> class Edge;

template <typename A, typename B> class Node;

template <typename A, typename B> class SequenceRule;

template <typename A, typename B> class Configuration {
private:
  /**
   * Foci are the atoms on which the configuration is based,
   * and which will carry the label. E.g., the chiral atom in
   * a tetrahedral chirality, or the bond ends in a double bond.
   */
  const std::vector<A> foci;

  /**
   * Carriers are the atoms neighboring the foci that define the
   * configuration. E.g., for a chiral atom, its four neighbors
   * define a parity; for a double bond, one neighbor on each
   * side of the bond defines it as Cis or Trans.
   */
  const std::vector<A> carriers;

  int cfg;

  std::shared_ptr<Digraph<A, B>> digraph = nullptr;

protected:
  Edge<A, B> *findInternalEdge(const std::vector<Edge<A, B> *> &edges, A f1,
                               A f2) {
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

  bool isInternalEdge(const Edge<A, B> *edge, A f1, A f2) {
    const auto &beg = edge->getBeg();
    const auto &end = edge->getEnd();
    if (f1 == beg->getAtom() && f2 == end->getAtom()) {
      return true;
    } else if (f1 == end->getAtom() && f2 == beg->getAtom()) {
      return true;
    }
    return false;
  }

  void removeInternalEdges(std::vector<Edge<A, B> *> &edges, A f1, A f2) {
    std::vector<Edge<A, B> *> new_edges;
    for (auto &&e : edges) {
      if (!isInternalEdge(e, f1, f2)) {
        new_edges.push_back(std::move(e));
      }
    }
    std::swap(edges, new_edges);
  }

  void removeDuplicatedEdges(std::vector<Edge<A, B> *> &&edges) {
    std::vector<Edge<A, B> *> new_edges;
    for (const auto &e : edges) {
      if (!e->getEnd()->isDuplicate()) {
        new_edges.push_back(e);
      }
    }
    std::swap(edges, new_edges);
  }

public:
  template <typename T>
  static int parity4(const std::vector<T> &trg, const std::vector<T> &ref) {
    if (ref.size() != 4 || trg.size() != ref.size()) {
      throw std::runtime_error("Parity vectors must have size 4.");
    }

    if (ref[0] == trg[0]) {
      if (ref[1] == trg[1]) {
        // a,b,c,d -> a,b,c,d
        if (ref[2] == trg[2] && ref[3] == trg[3])
          return 2;
        // a,b,c,d -> a,b,d,c
        if (ref[2] == trg[3] && ref[3] == trg[2])
          return 1;
      } else if (ref[1] == trg[2]) {
        // a,b,c,d -> a,c,b,d
        if (ref[2] == trg[1] && ref[3] == trg[3])
          return 1;
        // a,b,c,d -> a,c,d,b
        if (ref[2] == trg[3] && ref[3] == trg[1])
          return 2;
      } else if (ref[1] == trg[3]) {
        // a,b,c,d -> a,d,c,b
        if (ref[2] == trg[2] && ref[3] == trg[1])
          return 1;
        // a,b,c,d -> a,d,b,c
        if (ref[2] == trg[1] && ref[3] == trg[2])
          return 2;
      }
    } else if (ref[0] == trg[1]) {
      if (ref[1] == trg[0]) {
        // a,b,c,d -> b,a,c,d
        if (ref[2] == trg[2] && ref[3] == trg[3])
          return 1;
        // a,b,c,d -> b,a,d,c
        if (ref[2] == trg[3] && ref[3] == trg[2])
          return 2;
      } else if (ref[1] == trg[2]) {
        // a,b,c,d -> b,c,a,d
        if (ref[2] == trg[0] && ref[3] == trg[3])
          return 2;
        // a,b,c,d -> b,c,d,a
        if (ref[2] == trg[3] && ref[3] == trg[0])
          return 1;
      } else if (ref[1] == trg[3]) {
        // a,b,c,d -> b,d,c,a
        if (ref[2] == trg[2] && ref[3] == trg[0])
          return 2;
        // a,b,c,d -> b,d,a,c
        if (ref[2] == trg[0] && ref[3] == trg[2])
          return 1;
      }
    } else if (ref[0] == trg[2]) {
      if (ref[1] == trg[1]) {
        // a,b,c,d -> c,b,a,d
        if (ref[2] == trg[0] && ref[3] == trg[3])
          return 1;
        // a,b,c,d -> c,b,d,a
        if (ref[2] == trg[3] && ref[3] == trg[0])
          return 2;
      } else if (ref[1] == trg[0]) {
        // a,b,c,d -> c,a,b,d
        if (ref[2] == trg[1] && ref[3] == trg[3])
          return 2;
        // a,b,c,d -> c,a,d,b
        if (ref[2] == trg[3] && ref[3] == trg[1])
          return 1;
      } else if (ref[1] == trg[3]) {
        // a,b,c,d -> c,d,a,b
        if (ref[2] == trg[0] && ref[3] == trg[1])
          return 2;
        // a,b,c,d -> c,d,b,a
        if (ref[2] == trg[1] && ref[3] == trg[0])
          return 1;
      }
    } else if (ref[0] == trg[3]) {
      if (ref[1] == trg[1]) {
        // a,b,c,d -> d,b,c,a
        if (ref[2] == trg[2] && ref[3] == trg[0])
          return 1;
        // a,b,c,d -> d,b,a,c
        if (ref[2] == trg[0] && ref[3] == trg[2])
          return 2;
      } else if (ref[1] == trg[2]) {
        // a,b,c,d -> d,c,b,a
        if (ref[2] == trg[1] && ref[3] == trg[0])
          return 2;
        // a,b,c,d -> d,c,a,b
        if (ref[2] == trg[0] && ref[3] == trg[1])
          return 1;
      } else if (ref[1] == trg[0]) {
        // a,b,c,d -> d,a,c,b
        if (ref[2] == trg[2] && ref[3] == trg[1])
          return 2;
        // a,b,c,d -> d,a,b,c
        if (ref[2] == trg[1] && ref[3] == trg[2])
          return 1;
      }
    }

    // We should never hit this, but the compiler still complains
    // about a missing return statement.
    return 0;
  }

  Configuration() = default;

  Configuration(A focus, std::vector<A> &&carriers, int cfg)
      : foci{{focus}}, carriers{std::move(carriers)}, cfg{cfg} {};

  Configuration(std::vector<A> &&foci, std::vector<A> &&carriers, int cfg)
      : foci{std::move(foci)}, carriers{std::move(carriers)}, cfg{cfg} {}

  virtual ~Configuration() = default;

  A getFocus() const { return foci[0]; }

  const std::vector<A> &getFoci() const { return foci; }

  int getConfig() const { return cfg; }

  const std::vector<A> &getCarriers() const { return carriers; }

  std::shared_ptr<Digraph<A, B>> getDigraph() const {
    if (digraph == nullptr) {
      throw std::runtime_error("Digraph has not been set.");
    }
    return digraph;
  }

  void setDigraph(Digraph<A, B> *digraph) { this->digraph.reset(digraph); }

  virtual Descriptor label(Node<A, B> *node, Digraph<A, B> *digraph,
                           const SequenceRule<A, B> *comp) {
    (void)node;
    (void)digraph;
    (void)comp;

    return Descriptor::UNKNOWN;
  }

  virtual Descriptor label(const SequenceRule<A, B> *comp) = 0;

  virtual void setPrimaryLabel(BaseMol<A, B> *mol, Descriptor desc) = 0;
}; // namespace NewCIPLabelling

} // namespace NewCIPLabelling
} // namespace RDKit