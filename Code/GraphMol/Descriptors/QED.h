//
//  Copyright (C) 2026 Glysade LLC and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

/*! \file QED.h

  \brief C++ port of rdkit.Chem.QED

  QED :  quantitative estimation of drug-likeness:
  Bickerton, G.R.; Paolini, G.V.; Besnard, J.; Muresan, S.; Hopkins, A.L.
  (2012) 'Quantifying the chemical beauty of drugs', Nature Chemistry, 4,
  90-98 [https://doi.org/10.1038/nchem.1243]

*/
#include <RDGeneral/export.h>
#ifndef RD_QED_H
#define RD_QED_H

#include <string>

namespace RDKit {
class ROMol;
namespace Descriptors {

//! the properties used by QED (also used to hold the weights)
struct RDKIT_DESCRIPTORS_EXPORT QEDProperties {
  double MW = 0.0;
  double ALOGP = 0.0;
  double HBA = 0.0;
  double HBD = 0.0;
  double PSA = 0.0;
  double ROTB = 0.0;
  double AROM = 0.0;
  double ALERTS = 0.0;
};

RDKIT_DESCRIPTORS_EXPORT extern const QEDProperties QEDWeightsMax;
RDKIT_DESCRIPTORS_EXPORT extern const QEDProperties QEDWeightsMean;
RDKIT_DESCRIPTORS_EXPORT extern const QEDProperties QEDWeightsNone;

//! calculates the properties that are required for QED
/*!
  Hs are removed from a copy of the molecule before the properties are
  calculated.
*/
RDKIT_DESCRIPTORS_EXPORT QEDProperties calcQEDProperties(const ROMol &mol);

//! calculates the weighted QED score from precomputed properties
RDKIT_DESCRIPTORS_EXPORT double calcQED(
    const QEDProperties &props, const QEDProperties &weights = QEDWeightsMean);

//! calculates the QED score (with the mean weights by default)
RDKIT_DESCRIPTORS_EXPORT double calcQED(
    const ROMol &mol, const QEDProperties &weights = QEDWeightsMean);
const std::string QEDVersion = "1.1.0";

}  // namespace Descriptors
}  // namespace RDKit
#endif
