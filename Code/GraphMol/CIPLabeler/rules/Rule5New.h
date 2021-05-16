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
 * A descriptor pair rule. This rule defines that like descriptor pairs have
 * priority over unlike descriptor pairs.
 *
 */
class Rule5New : public SequenceRule {
 public:
  Rule5New();

  Rule5New(Descriptor ref);

  int compare(const Edge *a, const Edge *b) const override;

 private:
  const Descriptor d_ref = Descriptor::NONE;

  void fillPairs(const Node *beg, PairList &plist) const;

  Sort getRefSorter(const SequenceRule *replacement_rule) const;
};

}  // namespace CIPLabeler
}  // namespace RDKit
