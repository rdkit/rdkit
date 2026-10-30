//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RWMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Chirality.h>
#include "SpacialScore.h"

#include <vector>

namespace RDKit {
namespace Descriptors {

double calcSPS(const ROMol &mol, bool normalize) {
  RWMol molCp(mol);
  MolOps::findPotentialStereoBonds(molCp, false);

  const auto nAtoms = molCp.getNumAtoms();
  std::vector<bool> isStereo(nAtoms, false);
  // (pseudo)stereo centers, including unassigned ones
  bool cleanIt = false;
  bool flagPossible = true;
  for (const auto &si :
       Chirality::findPotentialStereo(molCp, cleanIt, flagPossible)) {
    if (si.type == Chirality::StereoType::Atom_Tetrahedral) {
      isStereo[si.centeredOn] = true;
    }
  }
  // atoms in stereo double bonds
  for (const auto bond : molCp.bonds()) {
    if (bond->getBondType() == Bond::BondType::DOUBLE &&
        bond->getStereo() != Bond::BondStereo::STEREONONE) {
      isStereo[bond->getBeginAtomIdx()] = true;
      isStereo[bond->getEndAtomIdx()] = true;
    }
  }

  if (!molCp.getRingInfo()->isSssrOrBetter()) {
    MolOps::findSSSR(molCp);
  }
  int score = 0;
  for (const auto atom : molCp.atoms()) {
    int hyb;
    switch (atom->getHybridization()) {
      case Atom::HybridizationType::SP:
        hyb = 1;
        break;
      case Atom::HybridizationType::SP2:
        hyb = 2;
        break;
      case Atom::HybridizationType::SP3:
        hyb = 3;
        break;
      default:
        hyb = 4;
    }
    int stereo = isStereo[atom->getIdx()] ? 2 : 1;
    int ring = 1;
    if (!atom->getIsAromatic() &&
        molCp.getRingInfo()->numAtomRings(atom->getIdx())) {
      ring = 2;
    }
    int deg = atom->getDegree();
    score += hyb * stereo * ring * deg * deg;
  }
  if (!normalize) {
    return score;
  }
  auto numHeavy = molCp.getNumHeavyAtoms();
  if (!numHeavy) {
    return 0.0;
  }
  return static_cast<double>(score) / numHeavy;
}

}  // namespace Descriptors
}  // namespace RDKit
