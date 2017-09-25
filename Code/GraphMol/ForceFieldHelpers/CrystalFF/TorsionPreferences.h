//
//  Copyright (C) 2015 Sereina Riniker
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_TORSIONPREFERENCES_H_
#define _RD_TORSIONPREFERENCES_H_

namespace RDKit {
class ROMol;
}

namespace ForceFields {
namespace CrystalFF {
//! Get the experimental torsional angles in a molecule
void getExperimentalTorsions(
    const RDKit::ROMol &mol, std::vector<std::vector<int> > &expTorsionAtoms,
    std::vector<std::pair<std::vector<int>, std::vector<double> > > &
        expTorsionAngles,
    std::vector<std::vector<int> > &improperAtoms, bool useExpTorsions = false,
    bool useBasicKnowledge = false, bool verbose = false);
}
}

#endif
