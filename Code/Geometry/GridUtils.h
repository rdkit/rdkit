//
//   Copyright (C) 2003-2007 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _GRIDUTILS_H_20050126
#define _GRIDUTILS_H_20050126

#include <vector>

namespace RDGeom {
class UniformGrid3D;
class Point3D;

//! calculate the tversky index between the shapes encoded on two grids
/*!

   tversky(S1,S2) =  | S1&S2 | / ( alpha * ( | S1 | - | S1&S2 | ) + beta * ( |
   S2 | - | S1&S2 | )  + | S1&S2 | )

*/

template <class GRIDTYPE>
RDKIT_RDGEOMETRYLIB_EXPORT double tverskyIndex(const GRIDTYPE &grid1,
                                               const GRIDTYPE &grid2,
                                               double alpha, double beta);

//! calculate the tanimoto distance between the shapes encoded on two grids
/*!

   tanimoto(S1,S2) =  1 - ( | S1&S2 | / | S1|S2 | )

*/

template <class GRIDTYPE>
RDKIT_RDGEOMETRYLIB_EXPORT double tanimotoDistance(const GRIDTYPE &grid1,
                                                   const GRIDTYPE &grid2);
//! calculate the protrude distance between the shapes encoded on two grids
/*!

   protrude(S1,S2) = ( | S1|S2 | - | S1&S2 | ) / | S1 |

*/
template <class GRIDTYPE>
RDKIT_RDGEOMETRYLIB_EXPORT double protrudeDistance(const GRIDTYPE &grid1,
                                                   const GRIDTYPE &grid2);

//! calculate the grid centroid within a window of a point
RDKIT_RDGEOMETRYLIB_EXPORT Point3D
computeGridCentroid(const UniformGrid3D &grid, const Point3D &pt,
                    double windowRadius, double &weightSum);

//! find terminal points of a shape encoded on a grid
//!  this is part of the subshape implementation
RDKIT_RDGEOMETRYLIB_EXPORT std::vector<Point3D> findGridTerminalPoints(
    const UniformGrid3D &grid, double windowRadius, double inclusionFraction);
}  // namespace RDGeom

#endif
