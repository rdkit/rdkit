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

#include "Rule6.h"

#include "../Digraph.h"

namespace RDKit {
namespace CIPLabeler {

Rule6::Rule6() = default;

int Rule6::compare(const Edge *a, const Edge *b) const {
  const auto &digraph = a->getBeg()->getDigraph();
  const auto &ref = digraph->getRule6Ref();
  if (ref == nullptr) {
    return 0;
  }
  const auto &aAtom = a->getEnd()->getAtom();
  const auto &bAtom = b->getEnd()->getAtom();
  if (ref == aAtom && ref != bAtom) {
    return +1;  // a is ref (has priority)
  } else if (ref != aAtom && ref == bAtom) {
    return -1;  // b is ref (has priority)
  }
  return 0;
}

}  // namespace CIPLabeler
}  // namespace RDKit