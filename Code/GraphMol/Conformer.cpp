// $Id$
//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
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

void Conformer::initFromOther(const Conformer &conf) {
  RDProps::operator=(conf);
  dp_mol = conf.dp_mol;
  auto nat = conf.getNumAtoms();
  d_positions.resize(nat);
  for (unsigned i = 0; i < nat; i++) {
    d_positions[i] = conf.getAtomPos(i);
  }
  d_id = conf.getId();
  df_is3D = conf.is3D();
}

Conformer::Conformer(const Conformer &conf) : RDProps() { initFromOther(conf); }
Conformer &Conformer::operator=(const Conformer &other) {
  if (this == &other) {
    return *this;
  }
  initFromOther(other);
  return *this;
}

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

RDGeom::POINT3D_VECT &Conformer::getPositions() {
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
  }
  return d_positions;
}

const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const {
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
  }
  URANGE_CHECK(atomId, d_positions.size());
  return d_positions[atomId];
}

RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) {
  if (dp_mol) {
    PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
  }
  URANGE_CHECK(atomId, d_positions.size());
  return d_positions[atomId];
}
}  // namespace RDKit
