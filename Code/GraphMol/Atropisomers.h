//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
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
#include <string>
#include <stdexcept>

namespace RDKit {
RDKIT_FILEPARSERS_EXPORT void DetectAtropisomerChirality(ROMol &mol,
                                                         const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT void WedgeBondsFromAtropisomers(
    const ROMol &mol, const Conformer *conf, const INT_MAP_INT &wedgeBonds);

RDKIT_FILEPARSERS_EXPORT bool doesMolHaveAtropisomers(const ROMol &mol);

RDKIT_FILEPARSERS_EXPORT bool GetAtropisomerAtomsAndBonds(
    const Bond *bond, Atom *atoms[2], std::vector<Bond *> bonds[2],
    const ROMol &mol);

RDKIT_FILEPARSERS_EXPORT void getAllAtomIdsForStereoGroup(
    const ROMol &mol, const StereoGroup &group,
    std::vector<unsigned int> &atomIds);
}  // namespace RDKit
#endif
