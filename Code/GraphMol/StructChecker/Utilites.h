//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include "../RDKitBase.h"
#include "StructChecker.h"

namespace RDKit {
namespace StructureCheck {

struct Neighbourhood {          // a set of an atom neighbours
  std::vector<unsigned> Atoms;  // indices of atoms
  std::vector<unsigned> Bonds;  // indices of bonds
};

void SetupNeighbourhood(const ROMol &mol,
                        std::vector<Neighbourhood> &neighbour_array);
bool getMolAtomPoints(const ROMol &mol,
                      std::vector<RDGeom::Point3D> &atomPoint, bool twod=false);

std::string LogNeighbourhood(
    const ROMol &mol, unsigned int idx,
    const std::vector<Neighbourhood> &bneighbour_array);
}
}
