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

#include "../Mancude.hpp"
#include "SequenceRule.hpp"

namespace RDKit {
namespace NewCIPLabelling {

template <typename A, typename B> class Rule1a : public SequenceRule<A, B> {
public:
  Rule1a() = delete;

  Rule1a(const BaseMol<A, B> *mol) : SequenceRule<A, B>(mol) {}

  int compare(const Edge<A, B> *a, const Edge<A, B> *b) const override {
    const int anum = a->getEnd()->getAtomicNumNumerator();
    const int aden = a->getEnd()->getAtomicNumDenominator();
    const int bnum = b->getEnd()->getAtomicNumNumerator();
    const int bden = b->getEnd()->getAtomicNumDenominator();
    if (anum == 0 || bnum == 0) {
      return 0;
    }
    if (aden == 1 && bden == 1) {
      return integer_compare(anum, bnum);
    }
    return Mancude::Fraction::compare(anum, aden, bnum, bden);
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit