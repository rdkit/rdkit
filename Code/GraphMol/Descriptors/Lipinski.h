//
//  Copyright (C) 2007-2011 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file Lipinski.h

  \brief Contains Lipinski and Lipinski-like descriptors. Use MolDescriptors.h
  in client code.

*/
#include <RDGeneral/export.h>
#ifndef __RD_LIPINSKI_H__
#define __RD_LIPINSKI_H__
#include "RegisterDescriptor.h"

namespace RDKit {
class ROMol;
namespace Descriptors {

const std::string lipinskiHBAVersion = "1.0.0";
//! calculates the standard Lipinski HBA definition (number of Ns and Os)
RDKIT_DESCRIPTORS_EXPORT unsigned int calcLipinskiHBA(const ROMol &mol);

const std::string lipinskiHBDVersion = "2.0.0";
//! calculates the standard Lipinski HBA definition (number of N-H and O-H
// bonds)
RDKIT_DESCRIPTORS_EXPORT unsigned int calcLipinskiHBD(const ROMol &mol);

enum NumRotatableBondsOptions {
  Default = -1,
  NonStrict = 0,
  Strict = 1,
  StrictLinkages = 2,
};

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumRotatableBondsVersion;
//! calculates the number of rotatable bonds
/*!
  \param mol           the molecule of interest
  \param strict        if Strict, a stricter definition of rotable bonds is used
                           this excludes amides, esters, etc.
                       if StrictLinkages, a much stricter definition that
                           handles rotatable bonds between rings as well.
                       if Default - uses the default choice (normally Strict)
*/
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumRotatableBonds(
    const ROMol &mol, NumRotatableBondsOptions useStrictDefinition = Default);

//! calculates the number of rotatable bonds ( backwards compatibility function,
//!  deprecated, please use calcNumRotatableBonds(const ROMol&, int)
/*!
  \param mol           the molecule of interest
  \param strict        if Strict == true, uses NumRotatableBondsOptions::Strict
*/
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumRotatableBonds(const ROMol &mol,
                                                            bool strict);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumHBDVersion;
//! calculates the number of H-bond donors
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumHBD(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumHBAVersion;
//! calculates the number of H-bond acceptors
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumHBA(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumHeteroatomsVersion;
//! calculates the number of heteroatoms
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumHeteroatoms(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumAmideBondsVersion;
//! calculates the number of amide bonds
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumAmideBonds(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string FractionCSP3Version;
//! calculates the fraction of carbons that are SP3 hybridized
RDKIT_DESCRIPTORS_EXPORT double calcFractionCSP3(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumRingsVersion;
//! calculates the number of SSSR rings
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumRings(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumAromaticRingsVersion;
//! calculates the number of aromatic SSSR rings
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumAromaticRings(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumAliphaticRingsVersion;
//! calculates the number of aliphatic (at least one non-aromatic bond) SSSR
// rings
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumAliphaticRings(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumSaturatedRingsVersion;
//! calculates the number of saturated SSSR rings
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumSaturatedRings(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumHeterocyclesVersion;
//! calculates the number of SSSR heterocycles
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumHeterocycles(const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string
    NumAromaticHeterocyclesVersion;
//! calculates the number of aromatic SSSR heterocycles
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumAromaticHeterocycles(
    const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumAromaticCarbocyclesVersion;
//! calculates the number of aromatic SSSR carbocycles
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumAromaticCarbocycles(
    const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string
    NumSaturatedHeterocyclesVersion;
//! calculates the number of saturated SSSR heterocycles
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumSaturatedHeterocycles(
    const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string
    NumSaturatedCarbocyclesVersion;
//! calculates the number of saturated SSSR carbocycles
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumSaturatedCarbocycles(
    const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string
    NumAliphaticHeterocyclesVersion;
//! calculates the number of aliphatic (at least one non-aromatic bond) SSSR
// heterocycles
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumAliphaticHeterocycles(
    const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string
    NumAliphaticCarbocyclesVersion;
//! calculates the number of aliphatic (at least one non-aromatic bond) SSSR
// carbocycles
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumAliphaticCarbocycles(
    const ROMol &mol);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumSpiroAtomsVersion;
//! calculates the number of spiro atoms (atoms shared between rings that share
// exactly one atom)
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumSpiroAtoms(
    const ROMol &mol, std::vector<unsigned int> *atoms = nullptr);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumBridgeheadAtomsVersion;
//! calculates the number of bridgehead atoms (atoms shared between rings that
// share at least two bonds)
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumBridgeheadAtoms(
    const ROMol &mol, std::vector<unsigned int> *atoms = nullptr);

RDKIT_DESCRIPTORS_EXPORT extern const std::string NumAtomStereoCentersVersion;
//! calculates the total number of atom stereo centers
RDKIT_DESCRIPTORS_EXPORT unsigned numAtomStereoCenters(const ROMol &mol);

//! calculates the number of unspecified stereo atom stereo centers
RDKIT_DESCRIPTORS_EXPORT extern const std::string
    NumUnspecifiedAtomStereoCentersVersion;
RDKIT_DESCRIPTORS_EXPORT unsigned numUnspecifiedAtomStereoCenters(
    const ROMol &mol);

//! Helper function to register the descriptors with the descriptor service
RDKIT_DESCRIPTORS_EXPORT void registerDescriptors();
}  // end of namespace Descriptors
}  // end of namespace RDKit

#endif
