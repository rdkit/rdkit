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
 * <b>Sequence Rule 4a</b>
 * <i>Chiral stereogenic units precede pseudoasymmetric stereogenic
 *    units and these precede nonstereogenic units.</i>
 *
 * @param <A> generic atom class
 */
template <typename A, typename B>
class Rule4a : public SequenceRule<A, B> {
 private:
  static int ord(Descriptor lab) {
    switch (lab) {
      case Descriptor::UNKNOWN:
      case Descriptor::ns:
      case Descriptor::NONE:
        return 0;
      case Descriptor::r:
      case Descriptor::s:
      case Descriptor::m:
      case Descriptor::p:
      case Descriptor::E:
      case Descriptor::Z:
        return 1;
      case Descriptor::R:
      case Descriptor::S:
      case Descriptor::M:
      case Descriptor::P:
      case Descriptor::seqTrans:
      case Descriptor::seqCis:
        return 2;
      default:
        throw std::logic_error("Invalid stereo descriptor");
    }
  }

 public:
  Rule4a() = delete;

  Rule4a(const BaseMol<A, B>* mol) : SequenceRule<A, B>(mol) {}

  int compare(const Edge<A, B>* a, const Edge<A, B>* b) const override {
    int aOrdinal = ord(this->getBondLabel(a));
    int bOrdinal = ord(this->getBondLabel(b));
    int cmp = integer_compare(aOrdinal, bOrdinal);
    if (cmp != 0) {
      return cmp;
    }
    aOrdinal = ord(a->getEnd()->getAux());
    bOrdinal = ord(b->getEnd()->getAux());
    return integer_compare(aOrdinal, bOrdinal);
  }
};

}  // namespace NewCIPLabelling
}  // namespace RDKit