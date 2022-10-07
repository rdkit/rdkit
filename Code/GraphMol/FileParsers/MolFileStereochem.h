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
struct RDKIT_FILEPARSERS_EXPORT StereoBondThresholds {
  const static unsigned DBL_BOND_NO_STEREO =
      1000;  //!< neighboring double bond without stereo info
  const static unsigned DBL_BOND_SPECIFIED_STEREO =
      10000;  //!< neighboring double bond with stereo specified
  const static unsigned CHIRAL_ATOM =
      100000;  //!< atom with specified chirality
  const static unsigned DIRECTION_SET =
      1000000;  //!< single bond with the direction already set
};
//! set wavy bonds around double bonds with STEREOANY stereo
/*!
 \param mol molecule to be modified
 \param clearDoubleBondFlags when this is true flags for unknown double bond
   stereo will also be removed.
 \param addWhenImpossible if nonzero a neighboring single bond will be made wavy
   even if it connects to a chiral center or double bond with specified stereo.
   one example of this would be the middle double bond in C/C=C/C=C/C=C/C (if
   that's set to STEREOANY after constructing the molecule) Otherwise, no wavy
   bond will be set
*/
RDKIT_FILEPARSERS_EXPORT void addWavyBondsForStereoAny(
    ROMol &mol, bool clearDoubleBondFlags = true,
    unsigned addWhenImpossible = StereoBondThresholds::DBL_BOND_NO_STEREO);

//! picks the bonds which should be wedged
/// \returns a map from bond idx -> controlling atom idx
RDKIT_FILEPARSERS_EXPORT INT_MAP_INT pickBondsToWedge(const ROMol &mol);
//! deprecated, please use MolOps::clearSingleBondDirFlags instead
RDKIT_FILEPARSERS_EXPORT void ClearSingleBondDirFlags(ROMol &mol);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, unsigned int fromAtomIdx, const Conformer *conf);
RDKIT_FILEPARSERS_EXPORT Bond::BondDir DetermineBondWedgeState(
    const Bond *bond, const INT_MAP_INT &wedgeBonds, const Conformer *conf);
//! Forces use of atom wedging from MolBlock if present.
/*
 \param mol: molecule to have its wedges altered
 */
RDKIT_FILEPARSERS_EXPORT void reapplyMolBlockWedging(ROMol &mol);

//! Set double bonds with unspecified stereo to STEREOANY and add wavy bonds to
///  potential stereocenters with unspecified chirality
RDKIT_FILEPARSERS_EXPORT void markUnspecifiedStereoAsUnknown(ROMol &mol,
                                                             int confId = -1);

}  // namespace RDKit
#endif
