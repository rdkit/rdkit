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

#include "Rule5.h"

namespace RDKit {
namespace CIPLabeler {

namespace {
int ord(Descriptor lab) {
  switch (lab) {
    case Descriptor::M:
    case Descriptor::R:
    case Descriptor::seqCis:
      return 2;
    case Descriptor::P:
    case Descriptor::S:
    case Descriptor::seqTrans:
      return 1;
    default:
      return 0;
  }
}
}  // namespace

Rule5::Rule5() = default;

int Rule5::compare(const Edge *a, const Edge *b) const {
  int aOrdinal = ord(getBondLabel(a));
  int bOrdinal = ord(getBondLabel(b));
  int cmp = three_way_comparison(aOrdinal, bOrdinal);
  if (cmp != 0) {
    return cmp;
  }
  aOrdinal = ord(a->getEnd()->getAux());
  bOrdinal = ord(b->getEnd()->getAux());
  return 2 * three_way_comparison(aOrdinal, bOrdinal);
}

}  // namespace CIPLabeler
}  // namespace RDKit