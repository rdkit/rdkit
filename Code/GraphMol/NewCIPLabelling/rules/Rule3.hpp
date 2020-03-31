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
 * <b>Sequence Rule 3</b>
 * <i>"‘seqcis’ (‘Z’) precedes ‘seqtrans’ (‘E’) and this order precedes
 * nonstereogenic double bonds"</i>
 *
 * @param <A> generic atom class
 */
template <typename A, typename B> class Rule3 : public SequenceRule<A, B> {
private:
  static int ord(Descriptor lab) {
    switch (lab) {
    case Descriptor::E:
      return 1;
    case Descriptor::Z:
      return 2;
    default:
      return 0;
    }
  }

public:
  Rule3() = delete;

  Rule3(const BaseMol<A, B> *mol) : SequenceRule<A, B>(mol) {}

  int compare(const Edge<A, B> *a, const Edge<A, B> *b) const override {
    return integer_compare(ord(a->getEnd()->getAux()),
                           ord(b->getEnd()->getAux()));
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit