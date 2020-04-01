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

#include "../Isotope.hpp"
#include "SequenceRule.hpp"

namespace RDKit {
namespace CIPLabelling {

/**
 * <b>Sequence Rule 1b</b>
 * <i>"A duplicate atom node whose corresponding nonduplicated atom
 * node is the root or is closer to the root ranks higher than
 * a duplicate atom node whose corresponding nonduplicated atom
 * node is farther from the root."</i>
 *
 * @param <A> generic atom class
 */
template <typename A, typename B> class Rule2 : public SequenceRule<A, B> {
private:
  const BaseMol<A, B> *mol;

public:
  Rule2() = delete;

  Rule2(const BaseMol<A, B> *mol) : SequenceRule<A, B>(mol), mol{mol} {}

  int compare(const Edge<A, B> *a, const Edge<A, B> *b) const override {
    int aAtomNum = mol->getAtomicNum(a->getEnd()->getAtom());
    int bAtomNum = mol->getAtomicNum(b->getEnd()->getAtom());
    if (aAtomNum == 0 || bAtomNum == 0) {
      return 0;
    }
    int aMassNum = a->getEnd()->isDuplicate()
                       ? 0
                       : mol->getMassNum(a->getEnd()->getAtom());
    int bMassNum = b->getEnd()->isDuplicate()
                       ? 0
                       : mol->getMassNum(b->getEnd()->getAtom());
    if (aMassNum == 0 && bMassNum == 0) {
      return 0;
    }
    auto aweight = Isotope::findAtomicWeight(aAtomNum, aMassNum);
    auto bweight = Isotope::findAtomicWeight(bAtomNum, bMassNum);

    if (aweight == bweight) {
      return 0;
    } else if (aweight < bweight) {
      return -1;
    } else {
      return 1;
    }
  }
};

} // namespace CIPLabelling
} // namespace RDKit