//
//  Copyright (C) 2004-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_ALIGN_POINTS_H__
#define __RD_ALIGN_POINTS_H__

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>

namespace RDNumeric {
  
  namespace Alignments {
    
    //! \brief Compute an optimal alignment (minimum sum of squared distance) between
    //! two sets of points in 3D 
    /*!
      \param refPoints      A vector of pointers to the reference points
      \param probePoints    A vector of pointers to the points to be aligned to the refPoints
      \param trans          A RDGeom::Transform3D object to capture the necessary transformation
      \param weights        A vector of weights for each of the points
      \param reflect        Add reflection is true
      \param maxIterations  Maximum number of iterations

      \return The sum of squared distances between the points

      <b>Note</b> 
      This function returns the sum of squared distance (SSR) not the RMSD
      RMSD = sqrt(SSR/numPoints)
    */
    double AlignPoints(const RDGeom::Point3DConstPtrVect &refPoints, 
                       const RDGeom::Point3DConstPtrVect &probePoints, 
                       RDGeom::Transform3D &trans,
                       const DoubleVector *weights=0, bool reflect=false, 
                       unsigned int maxIterations=50);
  }
}

#endif
      
