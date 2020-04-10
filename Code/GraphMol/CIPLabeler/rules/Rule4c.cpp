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

#include "Rule4c.h"

namespace RDKit {
namespace CIPLabeler {

int Rule4c::ord(Descriptor lab) {
  switch (lab) {
  case Descriptor::m:
  case Descriptor::r:
    return 2;
  case Descriptor::p:
  case Descriptor::s:
    return 1;
  default:
    return 0;
  }
}

Rule4c::Rule4c(const CIPMol *mol) : SequenceRule(mol) {}

int Rule4c::compare(const Edge *a, const Edge *b) const {
  // m vs p
  int aOrdinal = ord(this->getBondLabel(a));
  int bOrdinal = ord(this->getBondLabel(b));
  int cmp = integer_compare(aOrdinal, bOrdinal);
  if (cmp != 0) {
    return cmp;
  }
  // r vs s
  aOrdinal = ord(a->getEnd()->getAux());
  bOrdinal = ord(b->getEnd()->getAux());
  return integer_compare(aOrdinal, bOrdinal);
}

} // namespace CIPLabeler
} // namespace RDKit