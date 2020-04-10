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
 * <b>Sequence Rule 3</b>
 * <i>"‘seqcis’ (‘Z’) precedes ‘seqtrans’ (‘E’) and this order precedes
 * nonstereogenic double bonds"</i>
 *
 */
class Rule3 : public SequenceRule {
private:
  static int ord(Descriptor lab);

public:
  Rule3() = delete;

  Rule3(const CIPMol *mol);

  int compare(const Edge *a, const Edge *b) const override;
};

} // namespace CIPLabeler
} // namespace RDKit