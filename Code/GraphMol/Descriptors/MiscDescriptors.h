//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file MiscDescriptors.h

  \brief C++ ports of simple descriptors from rdkit.Chem.Descriptors

*/
#include <RDGeneral/export.h>
#ifndef RD_MISCDESCRIPTORS_H
#define RD_MISCDESCRIPTORS_H

#include <string>
#include <vector>

namespace RDKit {
class ROMol;
namespace Descriptors {

//! The average molecular weight of the molecule ignoring hydrogens
RDKIT_DESCRIPTORS_EXPORT double calcHeavyAtomMolWt(const ROMol &mol);
const std::string HeavyAtomMolWtVersion = "1.0.0";

//! The number of valence electrons the molecule has
RDKIT_DESCRIPTORS_EXPORT int calcNumValenceElectrons(const ROMol &mol);
const std::string NumValenceElectronsVersion = "1.1.0";

//! The number of radical electrons the molecule has (says nothing about spin
//! state)
RDKIT_DESCRIPTORS_EXPORT unsigned int calcNumRadicalElectrons(const ROMol &mol);
const std::string NumRadicalElectronsVersion = "1.1.0";

struct PartialCharges {
  std::vector<double> charges;
  double maxPartialCharge;
  double minPartialCharge;
  double maxAbsPartialCharge;
  double minAbsPartialCharge;
};

//! The largest Gasteiger partial charge in the molecule
/*!
  A molecule with no atoms returns 500 and -500 for max and min partial charges
*/
RDKIT_DESCRIPTORS_EXPORT PartialCharges calcPartialCharges(const ROMol &mol);
const std::string PartialChargeVersion = "1.0.0";
  
//! The largest Gasteiger partial charge in the molecule
/*!
  a molecule with no atoms returns -500
*/
RDKIT_DESCRIPTORS_EXPORT double calcMaxPartialCharge(const ROMol &mol);
const std::string MaxPartialChargeVersion = "1.0.0";
//! The smallest Gasteiger partial charge in the molecule
/*!
  a molecule with no atoms returns 500
*/
RDKIT_DESCRIPTORS_EXPORT double calcMinPartialCharge(const ROMol &mol);
const std::string MinPartialChargeVersion = "1.0.0";
//! max(abs(MinPartialCharge), abs(MaxPartialCharge))
RDKIT_DESCRIPTORS_EXPORT double calcMaxAbsPartialCharge(const ROMol &mol);
const std::string MaxAbsPartialChargeVersion = "1.0.0";
//! min(abs(MinPartialCharge), abs(MaxPartialCharge))
RDKIT_DESCRIPTORS_EXPORT double calcMinAbsPartialCharge(const ROMol &mol);
const std::string MinAbsPartialChargeVersion = "1.0.0";

//! Number of unique Morgan environments (unfolded count fingerprint with the
//! default Morgan generator settings) divided by the number of heavy atoms
/*!
  FpDensityMorgan1, FpDensityMorgan2 and FpDensityMorgan3 correspond to
  radius 1, 2 and 3. Molecules with no heavy atoms return 0.
*/
RDKIT_DESCRIPTORS_EXPORT double calcFpDensityMorgan(const ROMol &mol,
                                                    unsigned int radius);
const std::string FpDensityMorganVersion = "1.0.0";

}  // namespace Descriptors
}  // namespace RDKit
#endif
