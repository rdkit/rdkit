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

#include "Rule4a.h"

namespace RDKit {
namespace CIPLabeler {

namespace {
int ord(Descriptor lab) {
  switch (lab) {
    case Descriptor::UNKNOWN:
    case Descriptor::ns:
    case Descriptor::NONE:
      return 0;
    case Descriptor::r:
    case Descriptor::s:
    case Descriptor::m:
    case Descriptor::p:
    case Descriptor::E:
    case Descriptor::Z:
      return 1;
    case Descriptor::R:
    case Descriptor::S:
    case Descriptor::M:
    case Descriptor::P:
    case Descriptor::seqTrans:
    case Descriptor::seqCis:
      return 2;
    default:
      throw std::logic_error("Invalid stereo descriptor");
  }
}
}  // namespace

Rule4a::Rule4a() = default;

int Rule4a::compare(const Edge *a, const Edge *b) const {
  int aOrdinal = ord(getBondLabel(a));
  int bOrdinal = ord(getBondLabel(b));
  int cmp = three_way_comparison(aOrdinal, bOrdinal);
  if (cmp != 0) {
    return cmp;
  }
  aOrdinal = ord(a->getEnd()->getAux());
  bOrdinal = ord(b->getEnd()->getAux());
  return three_way_comparison(aOrdinal, bOrdinal);
}

}  // namespace CIPLabeler
}  // namespace RDKit