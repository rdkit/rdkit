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

#include "GraphMol/PeriodicTable.h"

#include "Rule2.h"

#include "../CIPMol.h"

namespace RDKit {

namespace CIPLabeler {

Rule2::Rule2() = default;

int Rule2::compare(const Edge *a, const Edge *b) const {
  auto aAtomNum = a->getEnd()->getAtomicNum();
  auto bAtomNum = b->getEnd()->getAtomicNum();
  if (aAtomNum == 0 || bAtomNum == 0) {
    return 0;
  }

  int aMassNum = a->getEnd()->isDuplicate() ? 0 : a->getEnd()->getMassNum();
  int bMassNum = b->getEnd()->isDuplicate() ? 0 : b->getEnd()->getMassNum();
  if (aMassNum == 0 && bMassNum == 0) {
    return 0;
  }

  const auto &table = RDKit::PeriodicTable::getTable();

  auto aweight = aMassNum ? table->getMassForIsotope(aAtomNum, aMassNum)
                          : table->getAtomicWeight(aAtomNum);
  auto bweight = bMassNum ? table->getMassForIsotope(bAtomNum, bMassNum)
                          : table->getAtomicWeight(bAtomNum);

  return three_way_comparison(aweight, bweight);
}

} // namespace CIPLabeler
} // namespace RDKit