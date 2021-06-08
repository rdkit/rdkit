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

#include "SequenceRule.h"

namespace RDKit {

namespace CIPLabeler {

/**
 * <b>Sequence Rule 2</b>
 * <i>"Higher atomic mass number precedes lower."</i>
 *
 * To resolve the ambiguity of what the "atomic mass"
 * is in case of duplicate nodes, isotpes, etc, this is
 * implmemented as the "proposed" rule 2 from the original
 * paper:
 *
 * <i>"Higher mass precedes lower mass, where mass is defined
 * in the case of a duplicate node as 0, an atom with isotope
 * indicated as its exact isotopic mass, and in all other
 * cases as the element’s atomic weight."</i>
 */
class Rule2 : public SequenceRule {
 public:
  Rule2();

  int compare(const Edge *a, const Edge *b) const override;
};

}  // namespace CIPLabeler
}  // namespace RDKit
