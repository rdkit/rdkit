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
//! set wavy bonds around double bonds with STEREOANY stereo
// if \c addWhenImpossible is nonzero, a neighboring single will be made wavy
// even if it connects to a chiral center or double bond with specified stereo.
//   one example of this would be the middle double bond in C/C=C/C=C/C=C/C
//   (if that's set to STEREOANY after constructing the molecule)
// Otherwise, no wavy bond will be set
//  thresholds:
//     1000: double bond with no stereochemistry indicated
//    10000: double bond with specified stereochemistry
//   100000: chiral atom
//  1000000: direction already set
//
struct RDKIT_FILEPARSERS_EXPORT StereoBondThresholds {
  const static unsigned DBL_BOND_NO_STEREO = 1000;
  const static unsigned DBL_BOND_SPECIFIED_STEREO = 10000;
  const static unsigned CHIRAL_ATOM = 100000;
  const static unsigned DIRECTION_SET = 1000000;
};
RDKIT_FILEPARSERS_EXPORT void addWavyBondsForStereoAny(
    ROMol &mol,
    unsigned addWhenImpossible = StereoBondThresholds::DBL_BOND_NO_STEREO);

//! picks the bonds which should be wedged
// \returns a map from bond idx -> controlling atom idx
RDKIT_FILEPARSERS_EXPORT INT_MAP_INT pickBondsToWedge(const ROMol &mol);
RDKIT_FILEPARSERS_EXPORT void ClearSingleBondDirFlags(ROMol &mol);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, unsigned int fromAtomIdx, const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, const INT_MAP_INT &wedgeBonds, const Conformer *conf);
}  // namespace RDKit
#endif
