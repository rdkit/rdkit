//
//  Copyright (C) 2013-2024 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_UFFUTILS_H
#define RD_UFFUTILS_H
#include <tuple>
#include <Geometry/point.h>

namespace ForceFields {
namespace UFF {
namespace Utils {
//! calculates and returns the cosine of the Y angle in an improper torsion
//! (see UFF paper, equation 19)
RDKIT_FORCEFIELD_EXPORT double calculateCosY(const RDGeom::Point3D &iPoint,
                                             const RDGeom::Point3D &jPoint,
                                             const RDGeom::Point3D &kPoint,
                                             const RDGeom::Point3D &lPoint);

//! calculates and returns the UFF force constant for an improper torsion
/*!

  \param at2AtomicNum   atomic number for atom 2
  \param isCBoundToO    boolean flag; true if atom 2 is sp2 carbon bound to sp2
  oxygen

  \return the force constant

*/
RDKIT_FORCEFIELD_EXPORT std::tuple<double, double, double, double>
calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO);
}  // namespace Utils
}  // namespace UFF
}  // namespace ForceFields

#endif