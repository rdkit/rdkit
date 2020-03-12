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

#include <vector>

#include "Configuration.hpp"

namespace RDKit
{
namespace NewCIPLabelling
{

template <typename A, typename B> class Tetrahedral : public Configuration<A, B>
{
  public:
    static const int LEFT = 0x1;
    static const int RIGHT = 0x2;

    Tetrahedral() = default;

    Tetrahedral(A focus, std::vector<A>&& carriers, int cfg)
        : Configuration<A, B>(focus, std::move(carriers), cfg){};

    void setPrimaryLabel(BaseMol<A, B>* mol, Descriptor desc) override
    {
        mol->setAtomDescriptor(this->getFocus(), CIP_LABEL_KEY, desc);
    }

    Descriptor label(const SequenceRule<A, B>* comp) override
    {
        auto digraph = this->getDigraph();
        auto root = digraph->getRoot();

        if (root == nullptr) {
            root = digraph->init(this->getFocus());
        } else {
            digraph->changeRoot(root);
        }

        return label(root, comp);
    }

    Descriptor label(Node<A, B>* node, Digraph<A, B>* digraph,
                     const SequenceRule<A, B>* comp) override
    {
        digraph->changeRoot(node);
        return label(node, comp);
    }

  private:
    Descriptor label(Node<A, B>* node, const SequenceRule<A, B>* comp) const
    {
        A focus = this->getFocus();
        auto edges = node->getEdges();

        // something not right!?! bad creation
        if (edges.size() < 3) {
            return Descriptor::ns;
        }

        auto priority = comp->sort(node, edges);

        bool isUnique = priority.isUnique();
        if (!isUnique && edges.size() == 4) {
            if (comp->getNumSubRules() == 3) {
                return Descriptor::UNKNOWN;
            }
            auto partition = comp->getSorter().getGroups(edges);
            if (partition.size() == 2) {
                // a a' b b' and a a' a'' b
                node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
                priority = comp->sort(node, edges);
                node->getDigraph()->setRule6Ref(nullptr);
            } else if (partition.size() == 1) {
                // S4 symmetric case
                node->getDigraph()->setRule6Ref(edges[0]->getEnd()->getAtom());
                comp->sort(node, edges);
                auto nbrs1 =
                    std::vector<Edge<A, B>*>(edges.begin(), edges.end());
                node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
                priority = comp->sort(node, edges);
                auto nbrs2 =
                    std::vector<Edge<A, B>*>(edges.begin(), edges.end());
                if (this->parity4(nbrs1, nbrs2) == 1) {
                    return Descriptor::UNKNOWN;
                }
                node->getDigraph()->setRule6Ref(nullptr);
            }
            if (!priority.isUnique()) {
                return Descriptor::UNKNOWN;
            }
        } else if (!isUnique) {
            return Descriptor::UNKNOWN;
        }

        auto ordered = std::vector<A>(4);
        int idx = 0;
        for (const auto& edge : edges) {
            if (edge->getEnd()->isSet(Node<A, B>::BOND_DUPLICATE) ||
                edge->getEnd()->isSet(Node<A, B>::IMPL_HYDROGEN)) {
                continue;
            }
            ordered[idx] = edge->getEnd()->getAtom();
            ++idx;
        }
        if (idx < 4) {
            ordered[idx] = focus;
        }

#if 0
****** Statistics? ******
    if (node->getDigraph()->getRoot() == node) {
      Stats.INSTANCE.countRule(priority.getRuleIdx());
    }
#endif

        int parity = this->parity4(ordered, this->getCarriers());

        if (parity == 0) {
            throw std::runtime_error(
                "Could not calculate parity! Carrier mismatch");
        }

        int config = this->getConfig();
        if (parity == 1) {
            config ^= 0x3;
        }

        if (config == 0x1) {
            if (priority.isPseduoAsymettric()) {
                return Descriptor::s;
            } else {
                return Descriptor::S;
            }
        } else if (config == 0x2) {
            if (priority.isPseduoAsymettric()) {
                return Descriptor::r;
            } else {
                return Descriptor::R;
            }
        }

        return Descriptor::UNKNOWN;
    }
};

} // namespace NewCIPLabelling
} // namespace RDKit