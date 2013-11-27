//
//  Copyright (C) 2004-2007 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_DISTGEOMUTILS_H_
#define _RD_DISTGEOMUTILS_H_

#include "BoundsMatrix.h"
#include <Numerics/SymmMatrix.h>
#include <map>
#include <Geometry/point.h>
#include "ChiralSet.h"
#include <RDGeneral/utils.h>

namespace ForceFields {
  class ForceField;
}


namespace DistGeom {

  //! Pick a distance matrix at random such that the
  //!  distance satisfy the bounds in the BoundsMatrix
  /*!
    \param mmat     Bounds matrix
    \param distmat  Storage for randomly chosen distances
    \param seed     the random number seed to use

    \return the largest element of the distance matrix
   */
  double pickRandomDistMat(const BoundsMatrix &mmat, 
                           RDNumeric::SymmMatrix<double> &distmat,
                           int seed=-1);
  //! \overload
  double pickRandomDistMat(const BoundsMatrix &mmat, 
                           RDNumeric::SymmMatrix<double> &distmat,
                           RDKit::double_source_type &rng);

  //! Compute an initial embedded in 3D based on a distance matrix
  /*! 
    This function follows the embed algorithm mentioned in 
    "Distance Geometry and Molecular Conformation" by G.M.Crippen and T.F.Havel
    (pages 312-313) 

    \param distmat     Distance matrix
    \param positions     A vector of pointers to Points to write out the resulting coordinates
    \param randNegEig  If set to true and if any of the eigen values are negative, we will
                       pick the corresponding components of the coordinates at random
    \param numZeroFail Fail embedding is more this many (or more) eigen values are zero
    \param seed        the random number seed to use

    \return true if the embedding was successful
  */
  bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distmat,  
                            RDGeom::PointPtrVect &positions, bool randNegEig=false, 
                            unsigned int numZeroFail=2,
                            int seed=-1);
  //! \overload
  bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distmat,  
                            RDGeom::PointPtrVect &positions,
                            RDKit::double_source_type &rng,
                            bool randNegEig=false, 
                            unsigned int numZeroFail=2
                            );

  //! places atoms randomly in a box
  /*! 
    \param positions     A vector of pointers to Points to write out the resulting coordinates
    \param boxSize     the side-length of the cubic box
    \param seed        the random number seed to use

    \return true if the coordinate generation was successful
  */
  bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
                           int seed=-1);
  //! \overload
  bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
                           RDKit::double_source_type &rng);

  //! Setup the error function for violation of distance bounds as a forcefield
  /*! 
    This is based on function E3 on page 311 of "Distance Geometry in Molecular
    Modeling" Jeffrey M.Blaney and J.Scott Dixon, Review in Computational Chemistry,
    Volume V

    \param mmat            Distance bounds matrix
    \param positions       A vector of pointers to Points to write out the resulting coordinates
    \param csets           The vector of chiral points (type: ChiralSet)
    \param weightChiral    weight to be used to enforce chirality
    \param weightFourthDim another chiral weight
    \param extraWeights    an optional set of weights for distance bounds violations
    \param basinSizeTol  Optional: any distance bound with a basin (distance between max and
                         min bounds) larger than this value will not be included in the force
                         field used to cleanup the structure.

    \return a pointer to a ForceField suitable for cleaning up the violations.
      <b>NOTE:</b> the caller is responsible for deleting this force field.

  */
  ForceFields::ForceField *constructForceField(const BoundsMatrix &mmat,
                                               RDGeom::PointPtrVect &positions, const VECT_CHIRALSET &csets,
                                               double weightChiral=1.0,
                                               double weightFourthDim=0.1,
                                               std::map< std::pair<int,int>,double> *extraWeights=0,
                                               double basinSizeTol=5.0);

}
    
#endif
