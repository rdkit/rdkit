//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file EState.h

  \brief C++ ports of the EState descriptors from rdkit.Chem.EState
  (EState indices and the hybrid EState-VSA descriptors).

*/
#include <RDGeneral/export.h>
#ifndef RD_ESTATE_DESCRIPTORS_H
#define RD_ESTATE_DESCRIPTORS_H

#include <string>
#include <vector>

namespace RDKit {
class ROMol;
namespace Descriptors {

//! returns the EState index for each atom in the molecule
/*!
  Reference: Hall, Mohney and Kier. JCICS _31_ 76-81 (1991)
*/
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEStateIndices(
    const ROMol &mol);
const std::string EStateIndicesVersion = "1.0.0";

//! the maximum EState index (0 for molecules with no atoms)
RDKIT_DESCRIPTORS_EXPORT double calcMaxEStateIndex(const ROMol &mol);
const std::string MaxEStateIndexVersion = "1.0.0";
//! the minimum EState index (0 for molecules with no atoms)
RDKIT_DESCRIPTORS_EXPORT double calcMinEStateIndex(const ROMol &mol);
const std::string MinEStateIndexVersion = "1.0.0";
//! the maximum absolute EState index (0 for molecules with no atoms)
RDKIT_DESCRIPTORS_EXPORT double calcMaxAbsEStateIndex(const ROMol &mol);
const std::string MaxAbsEStateIndexVersion = "1.0.0";
//! the minimum absolute EState index (0 for molecules with no atoms)
RDKIT_DESCRIPTORS_EXPORT double calcMinAbsEStateIndex(const ROMol &mol);
const std::string MinAbsEStateIndexVersion = "1.0.0";

//! the default bins used for the EState_VSA descriptors
RDKIT_DESCRIPTORS_EXPORT extern const std::vector<double> estateBins;
//! the default bins used for the VSA_EState descriptors
RDKIT_DESCRIPTORS_EXPORT extern const std::vector<double> vsaEStateBins;

//! EState_VSA1..EState_VSA11: Labute ASA contributions binned by EState index
/*!
  \param mol   the molecule of interest
  \param bins  bin boundaries, defaults to \c estateBins

  \return a vector of size bins.size()+1
*/
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcEState_VSA(
    const ROMol &mol, const std::vector<double> *bins = nullptr);
const std::string EState_VSAVersion = "1.0.1";

//! VSA_EState1..VSA_EState10: EState indices binned by Labute ASA
//! contribution
/*!
  \param mol   the molecule of interest
  \param bins  bin boundaries, defaults to \c vsaEStateBins

  \return a vector of size bins.size()+1
*/
RDKIT_DESCRIPTORS_EXPORT std::vector<double> calcVSA_EState(
    const ROMol &mol, const std::vector<double> *bins = nullptr);
const std::string VSA_EStateVersion = "1.0.0";

}  // namespace Descriptors
}  // namespace RDKit
#endif
