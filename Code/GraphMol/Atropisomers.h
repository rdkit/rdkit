//
//  Copyright (C) 2004-203 Tad hurst/CDD and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>

#ifndef RD_ATROPISOMERS_H
#define RD_ATROPISOMERS_H

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Chirality.h>
#include <string>
#include <stdexcept>

namespace RDKit {
namespace Atropisomers {
using AtropAtomAndBondVec = std::pair<Atom *, std::vector<Bond *>>;
RDKIT_GRAPHMOL_EXPORT void detectAtropisomerChirality(ROMol &mol,
                                                      const Conformer *conf);
RDKIT_GRAPHMOL_EXPORT void wedgeBondsFromAtropisomers(
    const ROMol &mol, const Conformer *conf,
    std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds);

RDKIT_GRAPHMOL_EXPORT bool doesMolHaveAtropisomers(const ROMol &mol);

RDKIT_GRAPHMOL_EXPORT bool getAtropisomerAtomsAndBonds(
    const Bond *bond, AtropAtomAndBondVec atomAndBonds[2], const ROMol &mol);

RDKIT_GRAPHMOL_EXPORT void getAllAtomIdsForStereoGroup(
    const ROMol &mol, const StereoGroup &group,
    std::vector<unsigned int> &atomIds,
    const std::map<int, std::unique_ptr<RDKit::Chirality::WedgeInfoBase>>
        &wedgeBonds);
}  // namespace Atropisomers
}  // namespace RDKit
#endif
