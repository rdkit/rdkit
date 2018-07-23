//
//  Copyright (C) 2017 Sereina Riniker
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_TORSIONPREFERENCES_H_
#define _RD_TORSIONPREFERENCES_H_
#include <vector>

namespace RDKit {
class ROMol;
}

namespace ForceFields {
namespace CrystalFF {
//! Get the experimental torsional angles in a molecule
RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &mol, std::vector<std::vector<int> > &expTorsionAtoms,
    std::vector<std::pair<std::vector<int>, std::vector<double> > > &
        expTorsionAngles,
    std::vector<std::vector<int> > &improperAtoms, bool useExpTorsions = false,
    bool useBasicKnowledge = false, unsigned int version = 1, bool verbose = false);
}
}

#endif
