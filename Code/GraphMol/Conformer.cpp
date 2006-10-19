// $Id$
//
//  Copyright (C) 2001-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include "Conformer.h"
#include "ROMol.h"

namespace RDKit {
  Conformer::Conformer() {
    dp_mol = 0;
    d_positions.clear();
    d_id = 0;
  }

  Conformer::Conformer(unsigned int numAtoms) {
    dp_mol = 0;
    if(numAtoms){
      d_positions.resize(numAtoms, RDGeom::Point3D(0.0, 0.0, 0.0));
    } else {
      d_positions.resize(0);
      d_positions.clear();
    }
  }

  Conformer::Conformer(const Conformer &conf) {
    dp_mol = 0;
    int i, nat = conf.getNumAtoms();
    d_positions.reserve(nat);
    
    for (i = 0; i < nat; i++) {
      d_positions.push_back(conf.getAtomPos(i));
    }
    d_id = conf.getId();
  }

  void Conformer::setOwningMol(ROMol *mol) {
    PRECONDITION(mol, "");
    dp_mol = mol;
  }

  void Conformer::setOwningMol(ROMol &mol) {setOwningMol(&mol);}

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
    RANGE_CHECK(0,atomId,d_positions.size()-1);
    return d_positions[atomId];
  } 
  
  RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) {
    if (dp_mol) {
      PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    }
    RANGE_CHECK(0,atomId,d_positions.size()-1);
    return d_positions[atomId];
  }
  
}
