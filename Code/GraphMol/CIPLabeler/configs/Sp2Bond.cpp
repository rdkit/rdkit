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
#include <GraphMol/Chirality.h>

#include "Sp2Bond.h"
#include "../Sort.h"
#include "../rules/Rules.h"

namespace RDKit {
namespace CIPLabeler {

Sp2Bond::Sp2Bond(const CIPMol &mol, Bond *bond, Atom *startAtom, Atom *endAtom,
                 Bond::BondStereo cfg)
    : Configuration(mol, {startAtom, endAtom}), dp_bond{bond}, d_cfg{cfg} {
  CHECK_INVARIANT(startAtom && endAtom, "bad foci")
  CHECK_INVARIANT(d_cfg == Bond::STEREOTRANS || d_cfg == Bond::STEREOCIS,
                  "bad config")

  auto stereo_atoms = Chirality::findStereoAtoms(bond);
  CHECK_INVARIANT(stereo_atoms.size() == 2, "incorrect number of stereo atoms")

  std::vector<Atom *> anchors{
      {mol.getAtom(stereo_atoms[0]), mol.getAtom(stereo_atoms[1])}};

  setCarriers(std::move(anchors));
}

void Sp2Bond::setPrimaryLabel(Descriptor desc) {
  switch (desc) {
    case Descriptor::seqTrans:
    case Descriptor::E:
    case Descriptor::seqCis:
    case Descriptor::Z: {
      auto carriers = getCarriers();
      dp_bond->setStereoAtoms(carriers[0]->getIdx(), carriers[1]->getIdx());
      dp_bond->setStereo(d_cfg);
      dp_bond->setProp(common_properties::_CIPCode, to_string(desc));
      return;
    }
    case Descriptor::R:
    case Descriptor::S:
    case Descriptor::r:
    case Descriptor::s:
    case Descriptor::M:
    case Descriptor::P:
    case Descriptor::m:
    case Descriptor::p:
    case Descriptor::SP_4:
    case Descriptor::TBPY_5:
    case Descriptor::OC_6:
      throw std::runtime_error(
          "Received a Descriptor that is not supported for double bonds");
    default:
      throw std::runtime_error("Received an invalid Bond Descriptor");
  }
}

Descriptor Sp2Bond::label(const Rules &comp) {
  auto &digraph = getDigraph();
  auto root1 = digraph.getOriginalRoot();
  if (digraph.getCurrentRoot() != root1) {
    digraph.changeRoot(root1);
  }

  return label(root1, digraph, comp);
}

Descriptor Sp2Bond::label(Node *root1, Digraph &digraph, const Rules &comp) {
  const auto &focus1 = getFoci()[0];
  const auto &focus2 = getFoci()[1];

  const auto &internal = findInternalEdge(root1->getEdges(), focus1, focus2);
  if (internal == nullptr) {
    return Descriptor::UNKNOWN;
  }
  const auto &root2 = internal->getOther(root1);

  auto edges1 = root1->getEdges();
  auto edges2 = root2->getEdges();
  removeInternalEdges(edges1, focus1, focus2);
  removeInternalEdges(edges2, focus1, focus2);

  auto carriers = std::vector<Atom *>(getCarriers());
  auto config = d_cfg;

  if (root1->getAtom() == focus2) {
    std::swap(carriers[0], carriers[1]);
  }

  digraph.changeRoot(root1);
  const auto &priority1 = comp.sort(root1, edges1);
  if (!priority1.isUnique()) {
    return Descriptor::UNKNOWN;
  }
  // swap
  if (edges1.size() > 1 && carriers[0] == edges1[1]->getEnd()->getAtom()) {
    if (config == Bond::STEREOCIS) {
      config = Bond::STEREOTRANS;
    } else {
      config = Bond::STEREOCIS;
    }
  }
  digraph.changeRoot(root2);
  const auto &priority2 = comp.sort(root2, edges2);
  if (!priority2.isUnique()) {
    return Descriptor::UNKNOWN;
  }
  // swap
  if (edges2.size() > 1 && carriers[1] == edges2[1]->getEnd()->getAtom()) {
    if (config == Bond::STEREOCIS) {
      config = Bond::STEREOTRANS;
    } else {
      config = Bond::STEREOCIS;
    }
  }

  if (config == Bond::STEREOCIS) {
    if (priority1.isPseudoAsymetric() != priority2.isPseudoAsymetric()) {
      return Descriptor::seqCis;
    } else {
      return Descriptor::Z;
    }
  } else if (config == Bond::STEREOTRANS) {
    if (priority1.isPseudoAsymetric() != priority2.isPseudoAsymetric()) {
      return Descriptor::seqTrans;
    } else {
      return Descriptor::E;
    }
  }
  return Descriptor::UNKNOWN;
}

}  // namespace CIPLabeler
}  // namespace RDKit