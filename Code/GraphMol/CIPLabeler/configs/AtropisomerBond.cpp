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
#include <GraphMol/Atropisomers.h>

#include "AtropisomerBond.h"
#include "../Sort.h"
#include "../rules/Rules.h"

namespace RDKit {
namespace CIPLabeler {
AtropisomerBond::AtropisomerBond(const CIPMol &mol, Bond *bond, Atom *startAtom,
                                 Atom *endAtom, Bond::BondStereo cfg)
    : Configuration(mol, {startAtom, endAtom}, true),
      dp_bond{bond},
      d_cfg{cfg} {
  CHECK_INVARIANT(startAtom && endAtom, "bad foci")
  CHECK_INVARIANT(d_cfg == Bond::STEREOATROPCW || d_cfg == Bond::STEREOATROPCCW,
                  "bad config")

  Atropisomers::AtropAtomAndBondVec atomAndBondVecs[2];
  if (!Atropisomers::getAtropisomerAtomsAndBonds(bond, atomAndBondVecs,
                                                 bond->getOwningMol())) {
    return;  // not an atropisomer
  }
  auto atom1 = mol.getAtom(atomAndBondVecs[0].second[0]->getOtherAtomIdx(
      atomAndBondVecs[0].first->getIdx()));
  auto atom2 = mol.getAtom(atomAndBondVecs[1].second[0]->getOtherAtomIdx(
      atomAndBondVecs[1].first->getIdx()));

  std::vector<Atom *> anchors{atom1, atom2};

  setCarriers(std::move(anchors));
}

void AtropisomerBond::setPrimaryLabel(Descriptor desc) {
  switch (desc) {
    case Descriptor::M:
    case Descriptor::P:
    case Descriptor::m:
    case Descriptor::p: {
      auto carriers = getCarriers();
      // dp_bond->setStereoAtoms(carriers[0]->getIdx(), carriers[1]->getIdx());
      // dp_bond->setStereo(d_cfg);
      dp_bond->setProp(common_properties::_CIPCode, to_string(desc));
      return;
    }
    case Descriptor::R:
    case Descriptor::S:
    case Descriptor::r:
    case Descriptor::s:
    case Descriptor::SP_4:
    case Descriptor::TBPY_5:
    case Descriptor::OC_6:
    case Descriptor::seqTrans:
    case Descriptor::E:
    case Descriptor::seqCis:
    case Descriptor::Z:
      throw std::runtime_error(
          "Received a Descriptor that is not supported for atropisomer bonds");
    default:
      throw std::runtime_error("Received an invalid Bond Descriptor");
  }
}

Descriptor AtropisomerBond::label(const Rules &comp) {
  auto &digraph = getDigraph();
  auto root1 = digraph.getOriginalRoot();
  if (digraph.getCurrentRoot() != root1) {
    digraph.changeRoot(root1);
  }

  return label(root1, digraph, comp);
}

Descriptor AtropisomerBond::label(Node *root1, Digraph &digraph,
                                  const Rules &comp) {
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

  removeDuplicatesAndHs(edges1);
  removeDuplicatesAndHs(edges2);

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
    if (config == Bond::STEREOATROPCCW) {
      config = Bond::STEREOATROPCW;
    } else {
      config = Bond::STEREOATROPCCW;
    }
  }
  digraph.changeRoot(root2);
  const auto &priority2 = comp.sort(root2, edges2);
  if (!priority2.isUnique()) {
    return Descriptor::UNKNOWN;
  }
  // swap
  if (edges2.size() > 1 && carriers[1] == edges2[1]->getEnd()->getAtom()) {
    if (config == Bond::STEREOATROPCCW) {
      config = Bond::STEREOATROPCW;
    } else {
      config = Bond::STEREOATROPCCW;
    }
  }

  if (config == Bond::STEREOATROPCCW) {
    if (priority1.isPseudoAsymetric() || priority2.isPseudoAsymetric()) {
      return Descriptor::m;
    } else {
      return Descriptor::M;
    }
  } else if (config == Bond::STEREOATROPCW) {
    if (priority1.isPseudoAsymetric() || priority2.isPseudoAsymetric()) {
      return Descriptor::p;
    } else {
      return Descriptor::P;
    }
  }
  return Descriptor::UNKNOWN;
}

}  // namespace CIPLabeler
}  // namespace RDKit