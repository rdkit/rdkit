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
 * <b>Sequence Rule 1b</b>
 * <i>"A duplicate atom node whose corresponding nonduplicated atom
 * node is the root or is closer to the root ranks higher than
 * a duplicate atom node whose corresponding nonduplicated atom
 * node is farther from the root."</i>
 *
 */
class Rule1b : public SequenceRule {
 public:
  Rule1b();

  int compare(const Edge *a, const Edge *b) const override;

  /**
   * Flag indicates whether to match the problematic
   * IUPAC 2013 recommendations for Rule 1B.
   */
 private:
  static const bool IUPAC_2013 = false;
};

}  // namespace CIPLabeler
}  // namespace RDKit
