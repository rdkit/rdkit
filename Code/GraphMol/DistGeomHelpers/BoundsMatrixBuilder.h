//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_BOUNDS_MATRIX_BUILDER_H
#define RD_BOUNDS_MATRIX_BUILDER_H

#include <DistGeom/BoundsMatrix.h>
#include "Embedder.h"

namespace RDKit {
class ROMol;
namespace DGeomHelpers {
//! Set default upper and lower distance bounds in a distance matrix
/*!
  \param mmat        pointer to the bounds matrix to be altered
  \param defaultMin  default value for the lower distance bounds
  \param defaultMax  default value for the upper distance bounds
*/
RDKIT_DISTGEOMHELPERS_EXPORT void initBoundsMat(DistGeom::BoundsMatrix *mmat,
                                                double defaultMin = 0.0,
                                                double defaultMax = 1000.0);
/*! \overload
 */
RDKIT_DISTGEOMHELPERS_EXPORT void initBoundsMat(DistGeom::BoundsMatPtr mmat,
                                                double defaultMin = 0.0,
                                                double defaultMax = 1000.0);

RDKIT_DISTGEOMHELPERS_EXPORT void setTopolBounds(
    const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    const EmbedParameters &params, bool scaleVDW = false,
    bool set15bounds = true, bool set14bounds = true, bool set13bounds = true);

//! Set upper and lower distance bounds between atoms in a molecule based on
/// topology
/*!
  This consists of setting 1-2, 1-3 and 1-4 distance based on bond lengths,
  bond angles and torsion angle ranges. Optionally 1-5 bounds can also be set,
  in particular, for path that contain rigid 1-4 paths.
  The final step involves setting lower bound to the sum of the vdW radii for
  the remaining atom pairs.
  \param mol          The molecule of interest
  \param mmat         Bounds matrix to the bounds are written
  \param set15bounds  If true try to set 1-5 bounds also based on topology
  \param scaleVDW     Ignored.
  \param useMacrocycle14config  If 1-4 distances bound heuristics for
  macrocycles is used <b>Note</b> For some strained systems the bounds matrix
  resulting from setting 1-5 bounds may fail triangle smoothing. In these cases
  it is recommended to back out and recompute the bounds matrix with no 1-5
  bounds and with vdW scaling.
  \param forceTransAmides  If true, amide bonds are enforced to be trans
  \param set14bounds  If true, set 1-4 distance bounds based on topology
  \param set13bounds  If true, set 1-3 distance bounds based on topology
*/
inline void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
                           bool set15bounds = true, bool scaleVDW = false,
                           bool useMacrocycle14config = false,
                           bool forceTransAmides = true,
                           bool set14bounds = true, bool set13bounds = true) {
  EmbedParameters params{.useMacrocycle14config = useMacrocycle14config,
                         .forceTransAmides = forceTransAmides};
  setTopolBounds(mol, mmat, params, scaleVDW, set15bounds, set14bounds,
                 set13bounds);
}

/* ! \overload */
RDKIT_DISTGEOMHELPERS_EXPORT void setTopolBounds(
    const ROMol &mol, DistGeom::BoundsMatPtr mmat,
    std::vector<std::pair<int, int>> &bonds,
    std::vector<std::vector<int>> &angles, const EmbedParameters &params,
    bool scaleVDW = false, bool set15bounds = true, bool set14bounds = true,
    bool set13bounds = true);
/*! \overload for experimental torsion angle preferences
 */
inline void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
                           std::vector<std::pair<int, int>> &bonds,
                           std::vector<std::vector<int>> &angles,
                           bool set15bounds = true, bool scaleVDW = false,
                           bool useMacrocycle14config = false,
                           bool forceTransAmides = true,
                           bool set14bounds = true, bool set13bounds = true) {
  EmbedParameters params{.useMacrocycle14config = useMacrocycle14config,
                         .forceTransAmides = forceTransAmides};
  setTopolBounds(mol, mmat, bonds, angles, params, scaleVDW, set15bounds,
                 set14bounds, set13bounds);
}

//! generate the vectors of bonds and angles used by (ET)KDG
RDKIT_DISTGEOMHELPERS_EXPORT void collectBondsAndAngles(
    const ROMol &mol, std::vector<std::pair<int, int>> &bonds,
    std::vector<std::vector<int>> &angles);

}  // namespace DGeomHelpers
}  // namespace RDKit
#endif
