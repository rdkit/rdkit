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

#include "Tetrahedral.h"

namespace RDKit {
namespace CIPLabeler {

Tetrahedral::Tetrahedral() = default;

Tetrahedral::Tetrahedral(Atom *focus, std::vector<Atom *> &&carriers, int cfg)
    : Configuration(focus, std::move(carriers), cfg){};

void Tetrahedral::setPrimaryLabel(CIPMol *mol, Descriptor desc) {
  mol->setAtomDescriptor(this->getFocus(), CIP_LABEL_KEY, desc);
}

Descriptor Tetrahedral::label(const SequenceRule *comp) {
  auto digraph = this->getDigraph();
  auto root = digraph->getRoot();

  if (root == nullptr) {
    root = digraph->init(this->getFocus());
  } else {
    digraph->changeRoot(root);
  }

  return label(root, comp);
}

Descriptor Tetrahedral::label(Node *node, Digraph *digraph,
                              const SequenceRule *comp) {
  digraph->changeRoot(node);
  return label(node, comp);
}

Descriptor Tetrahedral::label(Node *node, const SequenceRule *comp) const {
  auto focus = this->getFocus();
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
    auto partition = comp->getSorter()->getGroups(edges);
    if (partition.size() == 2) {
      // a a' b b' and a a' a'' b
      node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
      priority = comp->sort(node, edges);
      node->getDigraph()->setRule6Ref(nullptr);
    } else if (partition.size() == 1) {
      // S4 symmetric case
      node->getDigraph()->setRule6Ref(edges[0]->getEnd()->getAtom());
      comp->sort(node, edges);
      auto nbrs1 = std::vector<Edge *>(edges.begin(), edges.end());
      node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
      priority = comp->sort(node, edges);
      auto nbrs2 = std::vector<Edge *>(edges.begin(), edges.end());
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

  auto ordered = std::vector<Atom *>(4);
  int idx = 0;
  for (const auto &edge : edges) {
    if (edge->getEnd()->isSet(Node::BOND_DUPLICATE) ||
        edge->getEnd()->isSet(Node::IMPL_HYDROGEN)) {
      continue;
    }
    ordered[idx] = edge->getEnd()->getAtom();
    ++idx;
  }
  if (idx < 4) {
    ordered[idx] = focus;
  }

  int parity = this->parity4(ordered, this->getCarriers());

  if (parity == 0) {
    throw std::runtime_error("Could not calculate parity! Carrier mismatch");
  }

  int config = this->getConfig();
  if (parity == 1) {
    config ^= 0x3;
  }

  if (config == 0x1) {
    if (priority.isPseudoAsymetric()) {
      return Descriptor::s;
    } else {
      return Descriptor::S;
    }
  } else if (config == 0x2) {
    if (priority.isPseudoAsymetric()) {
      return Descriptor::r;
    } else {
      return Descriptor::R;
    }
  }

  return Descriptor::UNKNOWN;
}

} // namespace CIPLabeler
} // namespace RDKit