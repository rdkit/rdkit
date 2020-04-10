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

#include <stdexcept>
#include <memory>
#include <vector>

#include "../Descriptor.hpp"
#include "../Edge.h"
#include "../Node.h"
#include "../Sort.h"
#include "Pairlist.hpp"

namespace RDKit {
namespace CIPLabeler {

class CIPMol;

namespace {
inline int integer_compare(int x, int y) {
  return x < y ? -1 : (x == y ? 0 : 1);
}
} // namespace

class SequenceRule {
private:
  const CIPMol *mol;

protected:
  std::unique_ptr<const Sort> sorter = nullptr;

public:
  SequenceRule() = delete;

  SequenceRule(const CIPMol *mol);

  virtual ~SequenceRule();

  virtual const CIPMol *getMol() const;

  Descriptor getBondLabel(const Edge *edge) const;

  virtual int getNumSubRules() const;

  virtual bool isPseudoAsymmetric() const;

  int getComparision(const Edge *a, const Edge *b) const;

  virtual int getComparision(const Edge *a, const Edge *b, bool deep) const;

  virtual const Sort *getSorter() const;

  int recursiveCompare(const Edge *a, const Edge *b) const;

  void setSorter(const Sort *sorter);

  Priority sort(const Node *node, std::vector<Edge *> &edges, bool deep) const;

  Priority sort(const Node *node, std::vector<Edge *> &edges) const;

  virtual int compare(const Edge *a, const Edge *b) const = 0;

private:
  bool areUpEdges(Node *aNode, Node *bNode, Edge *aEdge, Edge *bEdge) const;
};

} // namespace CIPLabeler
} // namespace RDKit