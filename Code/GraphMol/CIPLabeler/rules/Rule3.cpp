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

#include "Rule3.h"

namespace RDKit {
namespace CIPLabeler {

int Rule3::ord(Descriptor lab) {
  switch (lab) {
  case Descriptor::E:
    return 1;
  case Descriptor::Z:
    return 2;
  default:
    return 0;
  }
}

Rule3::Rule3(const CIPMol *mol) : SequenceRule(mol) {}

int Rule3::compare(const Edge *a, const Edge *b) const {
  return integer_compare(ord(a->getEnd()->getAux()),
                         ord(b->getEnd()->getAux()));
}

} // namespace CIPLabeler
} // namespace RDKit