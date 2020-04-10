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
#include <stdexcept>

#include "Edge.h"
#include "Node.h"

namespace RDKit {
namespace CIPLabeler {

Edge::Edge(Node *beg, Node *end, Bond *bond) : beg{beg}, end{end}, bond{bond} {}

Node *Edge::getOther(Node *node) const {
  if (isBeg(node)) {
    return getEnd();
  } else if (isEnd(node)) {
    return getBeg();
  } else {
    throw std::runtime_error("Not an end-point of this edge!");
  }
}

Node *Edge::getBeg() const { return beg; }

Node *Edge::getEnd() const { return end; }

Bond *Edge::getBond() const { return bond; }

Descriptor Edge::getAux() const { return aux; }

bool Edge::isBeg(const Node *node) const { return node == beg; }

bool Edge::isEnd(const Node *node) const { return node == end; }

void Edge::setAux(Descriptor aux) { this->aux = aux; }

void Edge::flip() { std::swap(beg, end); }

} // namespace CIPLabeler
} // namespace RDKit
