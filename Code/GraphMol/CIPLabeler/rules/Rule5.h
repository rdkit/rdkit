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
 * <b>Sequence Rule 5</b>
 * <i>An atom or group with descriptor 'R','M' and 'seqCis'
 * has priority over its enantiomorph 'S', 'P' or 'seqTrans'.
 * </i>
 */
class Rule5 : public SequenceRule {
 public:
  Rule5();

  int compare(const Edge *a, const Edge *b) const override;
};

}  // namespace CIPLabeler
}  // namespace RDKit
