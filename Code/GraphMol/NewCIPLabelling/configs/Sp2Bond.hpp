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

#include "Configuration.hpp"

namespace RDKit
{
namespace NewCIPLabelling
{

template <typename A, typename B> class Sp2Bond : public Configuration<A, B>
{
  private:
    B bond;

  public:
    static const int OPPOSITE = 0x1;
    static const int TOGETHER = 0x2;

    Sp2Bond() = default;

    Sp2Bond(B bond, std::vector<A>&& foci, std::vector<A>&& carriers, int cfg)
        : Configuration<A, B>(std::move(foci), std::move(carriers), cfg),
          bond{bond}
    {
    }

    void setPrimaryLabel(BaseMol<A, B>* mol, Descriptor desc) override
    {
        mol->setBondDescriptor(bond, CIP_LABEL_KEY, desc);
    }

    Descriptor label(const SequenceRule<A, B>* comp) override
    {
        const auto& digraph = this->getDigraph();

        const A focus1 = this->getFoci()[0];
        const A focus2 = this->getFoci()[1];

        auto root1 = digraph->getRoot();
        if (root1 == nullptr) {
            root1 = digraph->init(focus1);
        } else {
            digraph->changeRoot(root1);
        }

        const auto& root1_edges = root1->getEdges();
        const auto& internal =
            this->findInternalEdge(root1_edges, focus1, focus2);
        auto filter = [&internal](const Edge<A, B>* e) {
            return e != internal;
        };

        std::vector<Edge<A, B>*> edges1;
        std::copy_if(root1_edges.begin(), root1_edges.end(),
                     std::back_inserter(edges1), filter);

        const auto& priority1 = comp->sort(root1, edges1);
        if (!priority1.isUnique()) {
            return Descriptor::UNKNOWN;
        }

        const auto& root2 = internal->getOther(root1);
        digraph->changeRoot(root2);

        std::vector<Edge<A, B>*> edges2;
        const auto& root2_edges = root2->getEdges();
        std::copy_if(root2_edges.begin(), root2_edges.end(),
                     std::back_inserter(edges2), filter);

        const auto& priority2 = comp->sort(root2, edges2);
        if (!priority2.isUnique()) {
            return Descriptor::UNKNOWN;
        }

        const auto& carriers = this->getCarriers();
        int config = this->getConfig();

        // swap
        if (edges1.size() > 1u &&
            carriers[0] == edges1[1]->getEnd()->getAtom()) {
            config ^= 0x3;
        }
        // swap
        if (edges2.size() > 1u &&
            carriers[1] == edges2[1]->getEnd()->getAtom()) {
            config ^= 0x3;
        }

#if 0
    Stats.INSTANCE.countRule(
        Math.max(priority1.getRuleIdx(), priority2.getRuleIdx()));
#endif

        if (config == TOGETHER) {
            if (priority1.isPseduoAsymettric() !=
                priority2.isPseduoAsymettric()) {
                return Descriptor::seqCis;
            } else {
                return Descriptor::Z;
            }
        } else if (config == OPPOSITE) {
            if (priority1.isPseduoAsymettric() !=
                priority2.isPseduoAsymettric()) {
                return Descriptor::seqTrans;
            } else {
                return Descriptor::E;
            }
        }

        return Descriptor::UNKNOWN;
    }

    Descriptor label(Node<A, B>* root1, Digraph<A, B>* digraph,
                     const SequenceRule<A, B>* rules) override
    {
        const A focus1 = this->getFoci()[0];
        const A focus2 = this->getFoci()[1];

        const auto& internal =
            this->findInternalEdge(root1->getEdges(), focus1, focus2);
        if (internal == nullptr) {
            return Descriptor::UNKNOWN;
        }
        const auto& root2 = internal->getOther(root1);

        auto edges1 = root1->getEdges();
        auto edges2 = root2->getEdges();
        this->removeInternalEdges(edges1, focus1, focus2);
        this->removeInternalEdges(edges2, focus1, focus2);

        auto carriers = std::vector<A>(this->getCarriers());
        int config = this->getConfig();

        if (root1->getAtom() == focus2) {
            std::swap(carriers[0], carriers[1]);
        }

        digraph->changeRoot(root1);
        const auto& priority1 = rules->sort(root1, edges1);
        if (!priority1.isUnique()) {
            return Descriptor::UNKNOWN;
        }
        // swap
        if (edges1.size() > 1 &&
            carriers[0] == edges1[1]->getEnd()->getAtom()) {
            config ^= 0x3;
        }
        digraph->changeRoot(root2);
        const auto& priority2 = rules->sort(root2, edges2);
        if (!priority2.isUnique()) {
            return Descriptor::UNKNOWN;
        }
        // swap
        if (edges2.size() > 1 &&
            carriers[1] == edges2[1]->getEnd()->getAtom()) {
            config ^= 0x3;
        }

        if (config == TOGETHER) {
            if (priority1.isPseduoAsymettric() !=
                priority2.isPseduoAsymettric()) {
                return Descriptor::seqCis;
            } else {
                return Descriptor::Z;
            }
        } else if (config == OPPOSITE) {
            if (priority1.isPseduoAsymettric() !=
                priority2.isPseduoAsymettric()) {
                return Descriptor::seqTrans;
            } else {
                return Descriptor::E;
            }
        }
        return Descriptor::UNKNOWN;
    }
}; // namespace NewCIPLabelling

} // namespace NewCIPLabelling
} // namespace RDKit