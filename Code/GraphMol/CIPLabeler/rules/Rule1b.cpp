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

#include "Rule1b.h"

namespace RDKit {
namespace CIPLabeler {

Rule1b::Rule1b() = default;

int Rule1b::compare(const Edge *a, const Edge *b) const {
  if (IUPAC_2013) {
    return -three_way_comparison(a->getEnd()->getDistance(),
                                 b->getEnd()->getDistance());
  } else {
    if (a->getEnd()->isSet(Node::RING_DUPLICATE) &&
        b->getEnd()->isSet(Node::RING_DUPLICATE)) {
      return -three_way_comparison(a->getEnd()->getDistance(),
                                   b->getEnd()->getDistance());
    } else {
      if (a->getEnd()->isSet(Node::RING_DUPLICATE) &&
          !b->getEnd()->isSet(Node::RING_DUPLICATE)) {
        return +1;
      }
      if (!a->getEnd()->isSet(Node::RING_DUPLICATE) &&
          b->getEnd()->isSet(Node::RING_DUPLICATE)) {
        return -1;
      }
      return 0;
    }
  }
}

}  // namespace CIPLabeler
}  // namespace RDKit