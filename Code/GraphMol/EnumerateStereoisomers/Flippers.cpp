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
#include <GraphMol/ROMol.h>
#include <GraphMol/EnumerateStereoisomers/Flippers.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace RDKit {
namespace EnumerateStereoisomers {
namespace details {
namespace {
void clearUnknownStereo(Bond *bond) {
  // Clear both spellings only when assigning the configuration they describe.
  bond->clearProp(common_properties::_UnknownStereo);
  if (bond->getBondDir() == Bond::UNKNOWN ||
      bond->getBondDir() == Bond::EITHERDOUBLE) {
    bond->setBondDir(Bond::NONE);
  }
}
}  // namespace

AtomFlipper::AtomFlipper(ROMol &mol, const Chirality::StereoInfo &si)
    : Flipper(), d_atomIdx(si.centeredOn) {
  for (auto bond : mol.atomBonds(mol.getAtomWithIdx(d_atomIdx))) {
    if (bond->getBeginAtomIdx() == d_atomIdx) {
      clearUnknownStereo(bond);
    }
  }
}

void AtomFlipper::flip(ROMol &mol, bool flag) {
  auto atom = mol.getAtomWithIdx(d_atomIdx);
  if (flag) {
    atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CW);
  } else {
    atom->setChiralTag(Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
  }
}

BondFlipper::BondFlipper(ROMol &mol, const Chirality::StereoInfo &si)
    : Flipper(), d_bondIdx(si.centeredOn) {
  auto bond = mol.getBondWithIdx(d_bondIdx);
  auto stereoAtoms = bond->getStereoAtoms();
  if (stereoAtoms.empty()) {
    if (si.controllingAtoms[0] != Atom::NOATOM &&
        si.controllingAtoms[2] != Atom::NOATOM) {
      bond->setStereoAtoms(si.controllingAtoms[0], si.controllingAtoms[2]);
    } else {
      d_valid = false;
    }
  }
  if (d_valid) {
    clearUnknownStereo(bond);
    for (auto atom : {bond->getBeginAtom(), bond->getEndAtom()}) {
      for (auto neighbor : mol.atomBonds(atom)) {
        if (neighbor->getBondType() == Bond::SINGLE) {
          clearUnknownStereo(neighbor);
        }
      }
    }
  }
}

void BondFlipper::flip(ROMol &mol, bool flag) {
  auto bond = mol.getBondWithIdx(d_bondIdx);
  if (flag) {
    bond->setStereo(Bond::STEREOCIS);
  } else {
    bond->setStereo(Bond::STEREOTRANS);
  }
}

bool BondFlipper::isValid() const { return d_valid; }

StereoGroupFlipper::StereoGroupFlipper(const StereoGroup &sg) : Flipper() {
  for (auto a : sg.getAtoms()) {
    d_original_parities.emplace_back(
        std::make_pair(a->getIdx(), a->getChiralTag()));
  }
}

void StereoGroupFlipper::flip(ROMol &mol, bool flag) {
  if (flag) {
    for (auto &a : d_original_parities) {
      mol.getAtomWithIdx(a.first)->setChiralTag(a.second);
    }
  } else {
    for (auto &[atom, chiralType] : d_original_parities) {
      if (chiralType == Atom::ChiralType::CHI_TETRAHEDRAL_CW) {
        mol.getAtomWithIdx(atom)->setChiralTag(
            Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
      } else if (chiralType == Atom::ChiralType::CHI_TETRAHEDRAL_CCW) {
        mol.getAtomWithIdx(atom)->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
      }
    }
  }
}

AtropisomerFlipper::AtropisomerFlipper(const Chirality::StereoInfo &si)
    : d_bondIdx(si.centeredOn) {
  d_ctrlAtoms =
      std::vector<unsigned int>{si.controllingAtoms[0], si.controllingAtoms[1],
                                si.controllingAtoms[2], si.controllingAtoms[3]};
}

void AtropisomerFlipper::flip(ROMol &mol, bool flag) {
  auto bond = mol.getBondWithIdx(d_bondIdx);
  auto b1 = mol.getBondBetweenAtoms(bond->getBeginAtomIdx(), d_ctrlAtoms[0]);
  auto b2 = mol.getBondBetweenAtoms(bond->getBeginAtomIdx(), d_ctrlAtoms[1]);
  auto b3 = mol.getBondBetweenAtoms(bond->getEndAtomIdx(), d_ctrlAtoms[2]);
  auto b4 = mol.getBondBetweenAtoms(bond->getEndAtomIdx(), d_ctrlAtoms[3]);
  b3->setBondDir(Bond::BondDir::NONE);
  b4->setBondDir(Bond::BondDir::NONE);
  if (!flag) {
    bond->setStereo(Bond::STEREOATROPCW);
    if (d_ctrlAtoms[0] < d_ctrlAtoms[1]) {
      b1->setBondDir(Bond::BondDir::BEGINDASH);
      b2->setBondDir(Bond::BondDir::NONE);
    } else {
      b2->setBondDir(Bond::BondDir::BEGINDASH);
      b1->setBondDir(Bond::BondDir::NONE);
    }
  } else {
    bond->setStereo(Bond::STEREOATROPCCW);
    if (d_ctrlAtoms[0] < d_ctrlAtoms[1]) {
      b1->setBondDir(Bond::BondDir::BEGINWEDGE);
      b2->setBondDir(Bond::BondDir::NONE);
    } else {
      b2->setBondDir(Bond::BondDir::BEGINWEDGE);
      b1->setBondDir(Bond::BondDir::NONE);
    }
  }
}

}  // namespace details
}  // namespace EnumerateStereoisomers
}  // namespace RDKit
