//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Utilites.h"

namespace RDKit {
 namespace StructureCheck {

void SetupNeighbourhood(const ROMol &mol, std::vector<Neighbourhood> &neighbours) {
    neighbours.clear();
    neighbours.resize(mol.getNumAtoms());

    for (unsigned i = 0; i < mol.getNumBonds(); i++) {
        const Bond* bond = mol.getBondWithIdx(i);
        unsigned a1 = bond->getBeginAtomIdx();
        unsigned a2 = bond->getEndAtomIdx();

        neighbours[a1].Atoms.push_back(a2);
        neighbours[a1].Bonds.push_back(i);

        neighbours[a2].Atoms.push_back(a1);
        neighbours[a2].Bonds.push_back(i);
    }
}

 }// namespace StructureCheck
} // namespace RDKit
