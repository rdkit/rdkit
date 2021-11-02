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
#include "../rules/Rules.h"

namespace RDKit {
namespace CIPLabeler {

Tetrahedral::Tetrahedral(const CIPMol &mol, Atom *focus)
    : Configuration(mol, focus) {
  CHECK_INVARIANT(focus, "bad atom")
  CHECK_INVARIANT(focus->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW ||
                      focus->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW,
                  "bad config")

  std::vector<Atom *> carriers;
  carriers.reserve(4);
  for (auto &nbr : mol.getNeighbors(focus)) {
    carriers.push_back(nbr);
  }
  if (carriers.size() < 4) {
    // Implicit H -- use the central atom instead of a dummy H
    carriers.push_back(focus);
  }
  if (carriers.size() < 4) {
    // Trigonal pyramid centers with an implicit H need a phantom
    // atom as fourth carrier. This one must be represented differently
    // than the implicit H.
    carriers.push_back(nullptr);
  }
  POSTCONDITION(carriers.size() == 4, "configurtion must have 4 carriers");

  setCarriers(std::move(carriers));
};

void Tetrahedral::setPrimaryLabel(Descriptor desc) {
  switch (desc) {
    case Descriptor::R:
    case Descriptor::S:
    case Descriptor::r:
    case Descriptor::s:
      getFocus()->setProp(common_properties::_CIPCode, to_string(desc));
      return;
    case Descriptor::seqTrans:
    case Descriptor::seqCis:
    case Descriptor::E:
    case Descriptor::Z:
    case Descriptor::M:
    case Descriptor::P:
    case Descriptor::m:
    case Descriptor::p:
    case Descriptor::SP_4:
    case Descriptor::TBPY_5:
    case Descriptor::OC_6:
      throw std::runtime_error(
          "Received a Descriptor that is not supported for atoms");
    default:
      throw std::runtime_error("Received an invalid Atom Descriptor");
  }
}

Descriptor Tetrahedral::label(const Rules &comp) {
  auto &digraph = getDigraph();

  auto root = digraph.getOriginalRoot();
  if (digraph.getCurrentRoot() != root) {
    digraph.changeRoot(root);
  }

  return label(root, comp);
}

Descriptor Tetrahedral::label(Node *node, Digraph &digraph, const Rules &comp) {
  digraph.changeRoot(node);
  return label(node, comp);
}

Descriptor Tetrahedral::label(Node *node, const Rules &comp) const {
  auto focus = getFocus();
  auto edges = node->getEdges();

  // something not right!?! bad creation
  if (edges.size() < 3) {
    return Descriptor::ns;
  }

  auto priority = comp.sort(node, edges);

  bool isUnique = priority.isUnique();
  if (!isUnique && edges.size() == 4) {
    if (comp.getNumSubRules() == 3) {
      return Descriptor::UNKNOWN;
    }
    auto partition = comp.getSorter()->getGroups(edges);
    if (partition.size() == 2) {
      node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
      priority = comp.sort(node, edges);
      node->getDigraph()->setRule6Ref(nullptr);
    } else if (partition.size() == 1) {
      // S4 symmetric case
      node->getDigraph()->setRule6Ref(edges[0]->getEnd()->getAtom());
      comp.sort(node, edges);
      auto nbrs1 = std::vector<Edge *>(edges.begin(), edges.end());

      node->getDigraph()->setRule6Ref(edges[1]->getEnd()->getAtom());
      priority = comp.sort(node, edges);

      node->getDigraph()->setRule6Ref(nullptr);

      if (parity4(nbrs1, edges) == 1) {
        return Descriptor::UNKNOWN;
      }
    }
    if (!priority.isUnique()) {
      return Descriptor::UNKNOWN;
    }
  } else if (!isUnique) {
    return Descriptor::UNKNOWN;
  }

  // if we are resolving a trigonal pyramid with an implicit H,
  // the 4th carrier will be a nullptr: we need to add a phantom
  // atom, which will always have the lowest priority, so that
  // it must be different than the representation of the implicit H.
  auto ordered = std::vector<Atom *>(4, nullptr);
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

  int parity = parity4(ordered, getCarriers());

  if (parity == 0) {
    throw std::runtime_error("Could not calculate parity! Carrier mismatch");
  }

  auto config = focus->getChiralTag();
  if (parity == 1) {
    if (config == Atom::CHI_TETRAHEDRAL_CCW) {
      config = Atom::CHI_TETRAHEDRAL_CW;
    } else {
      config = Atom::CHI_TETRAHEDRAL_CCW;
    }
  }

  if (config == Atom::CHI_TETRAHEDRAL_CCW) {
    if (priority.isPseudoAsymetric()) {
      return Descriptor::s;
    } else {
      return Descriptor::S;
    }
  } else if (config == Atom::CHI_TETRAHEDRAL_CW) {
    if (priority.isPseudoAsymetric()) {
      return Descriptor::r;
    } else {
      return Descriptor::R;
    }
  }

  return Descriptor::UNKNOWN;
}

}  // namespace CIPLabeler
}  // namespace RDKit