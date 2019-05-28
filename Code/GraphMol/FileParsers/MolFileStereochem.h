//
//  Copyright (C) 2004-2017 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_MOL_FILE_STEREOCHEM_H
#define RD_MOL_FILE_STEREOCHEM_H

#include <GraphMol/RDKitBase.h>

namespace RDKit {
//! deprecated, please use MolOps::assignChiralTypesFromBondDirs instead
RDKIT_FILEPARSERS_EXPORT void DetectAtomStereoChemistry(RWMol &mol,
                                                        const Conformer *conf);
//! deprecated, please use MolOps::detectBondStereoChemistry instead
RDKIT_FILEPARSERS_EXPORT void DetectBondStereoChemistry(ROMol &mol,
                                                        const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT void WedgeMolBonds(ROMol &mol, const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT void WedgeBond(Bond *bond, unsigned int fromAtomIdx,
                                        const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT INT_MAP_INT pickBondsToWedge(const ROMol &mol);
RDKIT_FILEPARSERS_EXPORT void ClearSingleBondDirFlags(ROMol &mol);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, unsigned int fromAtomIdx, const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, const INT_MAP_INT &wedgeBonds, const Conformer *conf);
}  // namespace RDKit
#endif
