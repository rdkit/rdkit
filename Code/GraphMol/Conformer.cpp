//
//  Copyright (C) 2001-2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Conformer.h"
#include "ROMol.h"

namespace RDKit {

void Conformer::setOwningMol(ROMol *mol) {
  PRECONDITION(mol, "");
  dp_mol = mol;
}

void Conformer::setOwningMol(ROMol &mol) { setOwningMol(&mol); }

const RDGeom::POINT3D_VECT &Conformer::getPositions() const {
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
  }
  return d_positions;
}

RDGeom::POINT3D_VECT &Conformer::getPositions() { return d_positions; }

const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const {
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
  }
  URANGE_CHECK(atomId, d_positions.size());
  return d_positions.at(atomId);
}

RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) {
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
  }
  URANGE_CHECK(atomId, d_positions.size());
  return d_positions.at(atomId);
}
}  // namespace RDKit
