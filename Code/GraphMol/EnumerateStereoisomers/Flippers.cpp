//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/Chirality.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/EnumerateStereoisomers/Flippers.h>

namespace RDKit::EnumerateStereoisomers::details {

AtomFlipper::AtomFlipper(RWMol &mol, const Chirality::StereoInfo &si)
    : Flipper() {
  dp_atom = mol.getAtomWithIdx(si.centeredOn);
}

void AtomFlipper::flip(bool flag) {
  if (flag) {
    dp_atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  } else {
    dp_atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }
}

BondFlipper::BondFlipper(RWMol &mol, const Chirality::StereoInfo &si)
    : Flipper() {
  dp_bond = mol.getBondWithIdx((si.centeredOn));
  auto stereoAtoms = dp_bond->getStereoAtoms();
  if (stereoAtoms.empty()) {
    if (si.controllingAtoms[0] != Chirality::StereoInfo::NOATOM &&
        si.controllingAtoms[2] != Chirality::StereoInfo::NOATOM) {
      dp_bond->setStereoAtoms(si.controllingAtoms[0], si.controllingAtoms[2]);
    } else {
      dp_bond = nullptr;
    }
  }
}

void BondFlipper::flip(bool flag) {
  if (flag) {
    dp_bond->setStereo(Bond::STEREOCIS);
  } else {
    dp_bond->setStereo(Bond::STEREOTRANS);
  }
}

StereoGroupFlipper::StereoGroupFlipper(const StereoGroup &sg) : Flipper() {
  for (auto a : sg.getAtoms()) {
    d_original_parities.emplace_back(std::make_pair(a, a->getChiralTag()));
  }
}

void StereoGroupFlipper::flip(bool flag) {
  if (flag) {
    for (auto &a : d_original_parities) {
      a.first->setChiralTag(a.second);
    }
  } else {
    for (auto &[atom, chiralType] : d_original_parities) {
      if (chiralType == Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
        atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
      } else if (chiralType == Atom::ChiralType::CHI_TETRAHEDRAL_CCW) {
        atom->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
      }
    }
  }
}

}  // namespace RDKit::EnumerateStereoisomers::details
