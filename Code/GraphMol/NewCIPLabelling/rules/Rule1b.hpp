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

#include "../Node.hpp"
#include "SequenceRule.hpp"

namespace RDKit {
namespace NewCIPLabelling {

/**
 * <b>Sequence Rule 1b</b>
 * <i>"A duplicate atom node whose corresponding nonduplicated atom
 * node is the root or is closer to the root ranks higher than
 * a duplicate atom node whose corresponding nonduplicated atom
 * node is farther from the root."</i>
 *
 * @param <A> generic atom class
 */
template <typename A, typename B> class Rule1b : public SequenceRule<A, B> {
  /**
   * Flag indicates whether to match the problematic
   * IUPAC 2013 recommendations for Rule 1B.
   */
private:
  static const bool IUPAC_2013 = false;

public:
  Rule1b() = delete;

  Rule1b(const BaseMol<A, B> *mol) : SequenceRule<A, B>(mol) {}

  int compare(const Edge<A, B> *a, const Edge<A, B> *b) const override {
    if (IUPAC_2013) {
      return -integer_compare(a->getEnd()->getDistance(),
                              b->getEnd()->getDistance());
    } else {
      if (a->getEnd()->isSet(Node<A, B>::RING_DUPLICATE) &&
          b->getEnd()->isSet(Node<A, B>::RING_DUPLICATE)) {
        return -integer_compare(a->getEnd()->getDistance(),
                                b->getEnd()->getDistance());
      } else {
        if (a->getEnd()->isSet(Node<A, B>::RING_DUPLICATE) &&
            !b->getEnd()->isSet(Node<A, B>::RING_DUPLICATE)) {
          return +1;
        }
        if (!a->getEnd()->isSet(Node<A, B>::RING_DUPLICATE) &&
            b->getEnd()->isSet(Node<A, B>::RING_DUPLICATE)) {
          return -1;
        }
        return 0;
      }
    }
  }
};

} // namespace NewCIPLabelling
} // namespace RDKit