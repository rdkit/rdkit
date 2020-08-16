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
}  // namespace RDKit

namespace ForceFields {
namespace CrystalFF {
struct CrystalFFDetails {
  std::vector<std::vector<int>> expTorsionAtoms;
  std::vector<std::pair<std::vector<int>, std::vector<double>>>
      expTorsionAngles;
  std::vector<std::vector<int>> improperAtoms;
  std::vector<std::pair<int, int>> bonds;
  std::vector<std::vector<int>> angles;
  std::vector<int> atomNums;
};

//! Get the experimental torsional angles in a molecule
RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &mol, CrystalFFDetails &details,
    bool useExpTorsions = false, 
    bool useSmallRingTorsions = false, bool useMacrocycleTorsions = false,
    bool useBasicKnowledge = false,
    unsigned int version = 1, bool verbose = false);
}  // namespace CrystalFF
}  // namespace ForceFields

#endif
