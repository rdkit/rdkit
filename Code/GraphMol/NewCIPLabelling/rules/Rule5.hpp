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
 * <b>Sequence Rule 5</b>
 * <i>An atom or group with descriptor 'R','M' and 'seqCis'
 * has priority over its enantiomorph 'S', 'P' or 'seqTrans'.
 * </i>
 * @param <A> generic atom class
 */
template <typename A, typename B> class Rule5 : public SequenceRule<A, B> {
public:
  Rule5() = delete;

  Rule5(BaseMol<A, B> *mol) : SequenceRule<A, B>(mol) {}

  bool isPseudoAsymmetric() const override { return true; }

  int compare(const Edge<A, B> *a, const Edge<A, B> *b) const override {
    int aOrdinal = ord(this->getBondLabel(a));
    int bOrdinal = ord(this->getBondLabel(b));
    int cmp = integer_compare(aOrdinal, bOrdinal);
    if (cmp != 0) {
      return cmp;
    }
    aOrdinal = ord(a->getEnd()->getAux());
    bOrdinal = ord(b->getEnd()->getAux());
    return 2 * integer_compare(aOrdinal, bOrdinal);
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit