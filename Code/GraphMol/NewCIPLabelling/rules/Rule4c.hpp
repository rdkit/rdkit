//
//
//  Copyright (C) 2020 Greg Landrum and T5 Informatics GmbH
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include "SequenceRule.hpp"

namespace RDKit {
namespace NewCIPLabelling {

/**
 * <b>Sequence Rule 4c</b>
 * <i>‘r’ precedes ‘s’ and ‘m’ precedes ‘p’</i>
 *
 * @param <A> generic atom class
 */
template <typename A, typename B>
class Rule4c : public SequenceRule<A, B> {
 private:
  static int ord(Descriptor lab) {
    switch (lab) {
      case Descriptor::m:
      case Descriptor::r:
        return 2;
      case Descriptor::p:
      case Descriptor::s:
        return 1;
      default:
        return 0;
    }
  }

 public:
  Rule4c() = delete;

  Rule4c(const BaseMol<A, B>* mol) : SequenceRule<A, B>(mol) {}

  int compare(const Edge<A, B>* a, const Edge<A, B>* b) const override {
    // m vs p
    int aOrdinal = ord(this->getBondLabel(a));
    int bOrdinal = ord(this->getBondLabel(b));
    int cmp = integer_compare(aOrdinal, bOrdinal);
    if (cmp != 0) {
      return cmp;
    }
    // r vs s
    aOrdinal = ord(a->getEnd()->getAux());
    bOrdinal = ord(b->getEnd()->getAux());
    return integer_compare(aOrdinal, bOrdinal);
  }
};

}  // namespace NewCIPLabelling
}  // namespace RDKit