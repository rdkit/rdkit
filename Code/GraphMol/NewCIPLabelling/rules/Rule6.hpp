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

namespace RDKit
{
namespace NewCIPLabelling
{

/**
 * <b>Sequence Rule 6 (proposed)</b>
 * @param <A> generic atom class
 */
template <typename A, typename B> class Rule6 : public SequenceRule<A, B>
{
  public:
    Rule6() = delete;

    Rule6(const BaseMol<A, B>* mol) : SequenceRule<A, B>(mol) {}

    bool isPseudoAsymmetric() const override
    {
        return true; // comes after Rule 5 so must be true
    }

    int compare(const Edge<A, B>* a, const Edge<A, B>* b) const override
    {
        const auto& digraph = a->getBeg()->getDigraph();
        A ref = digraph->getRule6Ref();
        if (ref == nullptr) {
            return 0;
        }
        A aAtom = a->getEnd()->getAtom();
        A bAtom = b->getEnd()->getAtom();
        if (ref == aAtom && ref != bAtom) {
            return +1; // a is ref (has priority)
        } else if (ref != aAtom && ref == bAtom) {
            return -1; // b is ref (has priority)
        }
        return 0;
    }
};

} // namespace NewCIPLabelling
} // namespace RDKit