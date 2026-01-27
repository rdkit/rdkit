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
#include <RDGeneral/types.h>

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
      dp_bond->setProp(common_properties::_CIPCode, to_string(desc));
      dp_bond->setProp(common_properties::_CIPPrioritizedAnchors,
                       d_ranked_anchors, true);
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

bool AtropisomerBond::hasPrimaryLabel() const {
  return dp_bond->hasProp(common_properties::_CIPCode);
}

void AtropisomerBond::resetPrimaryLabel() const {
  dp_bond->clearProp(common_properties::_CIPCode);
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

  d_ranked_anchors.clear();

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

  {
    // This is mostly the same as in Sp2Bonds, but I doubt the anchors will be
    // implicit Hs in this case.

    // At this point, edges1 and edges2 are sorted by priority starting from
    // this node. Record that now! - they may be resorted after processing
    // other nodes.

    // As weird as it seems, these may actually be implicit Hs: Rule 2
    // in the paper on which this code is based states that,
    // in CIP ranks, H > 1H, so implicit H actually has a higher
    // priority than 1H (!!!). getAtomIdx() returns Atom::NOATOM
    // if that is the case.
    auto carrier1_idx = edges1[0]->getEnd()->getAtomIdx();
    auto carrier2_idx = edges2[0]->getEnd()->getAtomIdx();

    // Make sure the stereo atoms are in the right order
    if (edges1[0]->getBeg()->getAtom() == focus1) {
      d_ranked_anchors.assign({carrier1_idx, carrier2_idx});
    } else if (edges2[0]->getBeg()->getAtom() == focus1) {
      d_ranked_anchors.assign({carrier2_idx, carrier1_idx});
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