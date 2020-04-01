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

namespace RDKit {
namespace CIPLabelling {

template <typename A, typename B> class Node;

template <typename A, typename B> class Edge {
private:
  Node<A, B> *beg;
  Node<A, B> *end;
  B bond;
  Descriptor aux = Descriptor::NONE;

public:
  Edge() = delete;
  Edge(const Edge &) = delete;
  Edge &operator=(const Edge &) = delete;

  Edge(Node<A, B> *beg, Node<A, B> *end, B bond)
      : beg{beg}, end{end}, bond{bond} {}

  Node<A, B> *getOther(Node<A, B> *node) const {
    if (isBeg(node)) {
      return getEnd();
    } else if (isEnd(node)) {
      return getBeg();
    } else {
      throw std::runtime_error("Not an end-point of this edge!");
    }
  }

  Node<A, B> *getBeg() const { return beg; }

  Node<A, B> *getEnd() const { return end; }

  B getBond() const { return bond; }

  Descriptor getAux() const { return aux; }

  bool isBeg(const Node<A, B> *node) const { return node == beg; }

  bool isEnd(const Node<A, B> *node) const { return node == end; }

  void setAux(Descriptor aux) { this->aux = aux; }

  void flip() { std::swap(beg, end); }
};

} // namespace CIPLabelling
} // namespace RDKit
