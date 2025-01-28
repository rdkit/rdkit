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

#include <GraphMol/PeriodicTable.h>

#include "Rule2.h"

#include "../CIPMol.h"

namespace RDKit {

namespace CIPLabeler {

Rule2::Rule2() = default;

int Rule2::compare(const Edge *a, const Edge *b) const {
  auto a_end = a->getEnd();
  auto b_end = b->getEnd();

  auto aAtomNum = a_end->getAtomicNum();
  auto bAtomNum = b_end->getAtomicNum();

  if (aAtomNum == 0 && bAtomNum == 0) {
    return 0;
  } else if (aAtomNum == 0 || bAtomNum == 0) {
    // This should be caught by Rule 1a, but just in case
    return three_way_comparison(aAtomNum, bAtomNum);
  }

  auto aMassNum = a_end->getMassNum();
  auto bMassNum = b_end->getMassNum();
  if (aMassNum == 0u && bMassNum == 0u) {
    return 0;
  }

  auto aweight = a_end->getAtomicMass();
  auto bweight = b_end->getAtomicMass();

  return three_way_comparison(aweight, bweight);
}

}  // namespace CIPLabeler
}  // namespace RDKit
