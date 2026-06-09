//
//  Copyright (C) 2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_FINITEDIFFERENCE_H
#define RD_FINITEDIFFERENCE_H

namespace ForceFields {
class ForceField;

/// \brief Compares analytic and central-difference gradients for a force field.
///
/// \param ff        an initialized ForceField
/// \param stepSize  the displacement used for the central difference
///
/// \return the maximum absolute deviation between analytic and
///         finite-difference gradient components
RDKIT_FORCEFIELD_EXPORT double calcFiniteDifference(ForceField &ff,
                                                    double stepSize = 1e-5);

}  // namespace ForceFields

#endif
