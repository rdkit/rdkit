//
//  Copyright (C) 2004-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_DISTGEOMUTILS_H_
#define _RD_DISTGEOMUTILS_H_

#include "BoundsMatrix.h"
#include <Numerics/SymmMatrix.h>
#include <map>
#include <Geometry/point.h>
#include "ChiralSet.h"
#include <RDGeneral/utils.h>
#include <boost/dynamic_bitset.hpp>

namespace ForceFields {
class ForceField;
namespace CrystalFF {
struct CrystalFFDetails;
}
}  // namespace ForceFields

namespace DistGeom {

//! Pick a distance matrix at random such that the
//!  distance satisfy the bounds in the BoundsMatrix
/*!
  \param mmat     Bounds matrix
  \param distmat  Storage for randomly chosen distances
  \param seed     the random number seed to use

  \return the largest element of the distance matrix
 */
RDKIT_DISTGEOMETRY_EXPORT double pickRandomDistMat(
    const BoundsMatrix &mmat, RDNumeric::SymmMatrix<double> &distmat,
    int seed = -1);
//! \overload
RDKIT_DISTGEOMETRY_EXPORT double pickRandomDistMat(
    const BoundsMatrix &mmat, RDNumeric::SymmMatrix<double> &distmat,
    RDKit::double_source_type &rng);

//! Compute an initial embedded in 3D based on a distance matrix
/*!
  This function follows the embed algorithm mentioned in
  "Distance Geometry and Molecular Conformation" by G.M.Crippen and T.F.Havel
  (pages 312-313)

  \param distmat     Distance matrix
  \param positions     A vector of pointers to Points to write out the resulting
  coordinates
  \param randNegEig  If set to true and if any of the eigen values are negative,
  we will
                     pick the corresponding components of the coordinates at
  random
  \param numZeroFail Fail embedding is more this many (or more) eigen values are
  zero
  \param seed        the random number seed to use

  \return true if the embedding was successful
*/
RDKIT_DISTGEOMETRY_EXPORT bool computeInitialCoords(
    const RDNumeric::SymmMatrix<double> &distmat,
    RDGeom::PointPtrVect &positions, bool randNegEig = false,
    unsigned int numZeroFail = 2, int seed = -1);
//! \overload
RDKIT_DISTGEOMETRY_EXPORT bool computeInitialCoords(
    const RDNumeric::SymmMatrix<double> &distmat,
    RDGeom::PointPtrVect &positions, RDKit::double_source_type &rng,
    bool randNegEig = false, unsigned int numZeroFail = 2);

//! places atoms randomly in a box
/*!
  \param positions     A vector of pointers to Points to write out the resulting
  coordinates
  \param boxSize     the side-length of the cubic box
  \param seed        the random number seed to use

  \return true if the coordinate generation was successful
*/
RDKIT_DISTGEOMETRY_EXPORT bool computeRandomCoords(
    RDGeom::PointPtrVect &positions, double boxSize, int seed = -1);
//! \overload
RDKIT_DISTGEOMETRY_EXPORT bool computeRandomCoords(
    RDGeom::PointPtrVect &positions, double boxSize,
    RDKit::double_source_type &rng);

//! Setup the error function for violation of distance bounds as a forcefield
/*!
  This is based on function E3 on page 311 of "Distance Geometry in Molecular
  Modeling" Jeffrey M.Blaney and J.Scott Dixon, Review in Computational
  Chemistry,
  Volume V

  \param mmat            Distance bounds matrix
  \param positions       A vector of pointers to Points to write out the
  resulting coordinates
  \param csets           The vector of chiral points (type: ChiralSet)
  \param weightChiral    weight to be used to enforce chirality
  \param weightFourthDim another chiral weight
  \param extraWeights    an optional set of weights for distance bounds
  violations
  \param basinSizeTol  Optional: any distance bound with a basin (distance
  between max and
                       min bounds) larger than this value will not be included
  in the force
                       field used to cleanup the structure.

  \return a pointer to a ForceField suitable for cleaning up the violations.
    <b>NOTE:</b> the caller is responsible for deleting this force field.

*/
RDKIT_DISTGEOMETRY_EXPORT ForceFields::ForceField *constructForceField(
    const BoundsMatrix &mmat, RDGeom::PointPtrVect &positions,
    const VECT_CHIRALSET &csets, double weightChiral = 1.0,
    double weightFourthDim = 0.1,
    std::map<std::pair<int, int>, double> *extraWeights = nullptr,
    double basinSizeTol = 5.0, boost::dynamic_bitset<> *fixedPts = nullptr);

//! Force field with experimental torsion angle preferences and 1-2/1-3 distance
/// constraints
/*!

  \param mmat            Distance bounds matrix
  \param positions       A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param etkdgDetails    Contains information about the ETKDG force field

  <b>NOTE:</b> the caller is responsible for deleting this force field.

*/
RDKIT_DISTGEOMETRY_EXPORT ForceFields::ForceField *construct3DForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails);
//! Force field with experimental torsion angle preferences and 1-2/1-3 distance
/// constraints, as well as atom pairwise Columbic interactions
/*!

  \param mmat            Distance bounds matrix
  \param positions       A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param etkdgDetails    Contains information about the ETKDG force field
  \param CPCI            Contains which atom pair(s) have what strength of
  attractive/repulsive electrostatic interaction(s)

  <b>NOTE:</b> the caller is responsible for deleting this force field.

*/
RDKIT_DISTGEOMETRY_EXPORT ForceFields::ForceField *construct3DForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails,
    const std::map<std::pair<unsigned int, unsigned int>, double> &CPCI);
//! Force field with experimental torsion angle preferences and 1-2/1-3 distance
/// constraints
/*!

  \param mmat            Distance bounds matrix
  \param positions       A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param etkdgDetails    Contains information about the ETKDG force field

  <b>NOTE:</b> the caller is responsible for deleting this force field.

*/
RDKIT_DISTGEOMETRY_EXPORT ForceFields::ForceField *constructPlain3DForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const ForceFields::CrystalFF::CrystalFFDetails &etkdgDetails);

//! Force field with only improper terms
/*!

  \param mmat            Distance bounds matrix
  \param positions       A vector of pointers to 3D Points to write out the
  resulting coordinates
  \param improperAtoms   A list of groups of 4 atom indices for inversion terms
  \param angles          List of lists with the three angle indices and whether
  the center atom in the angle is SP hybridized for every angle in the molecule.
  \param atomNums        A list of atomic numbers for all atoms in the molecule

  \return a pointer to a ForceField with improper terms
    <b>NOTE:</b> the caller is responsible for deleting this force field.

*/
RDKIT_DISTGEOMETRY_EXPORT ForceFields::ForceField *
construct3DImproperForceField(
    const BoundsMatrix &mmat, RDGeom::Point3DPtrVect &positions,
    const std::vector<std::vector<int>> &improperAtoms,
    const std::vector<std::vector<int>> &angles,
    const std::vector<int> &atomNums);

}  // namespace DistGeom

#endif
