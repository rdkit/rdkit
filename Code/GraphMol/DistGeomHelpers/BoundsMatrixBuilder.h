//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_BOUNDS_MATRIX_BUILDER_H_
#define _RD_BOUNDS_MATRIX_BUILDER_H_

#include <DistGeom/BoundsMatrix.h>

namespace RDKit {
  class ROMol;
  namespace DGeomHelpers {
    //! Set default upper and lower distance bounds in a distance matrix
    /*!   
      \param mmat        pointer to the bounds matrix to be altered
      \param defaultMin  default value for the lower distance bounds
      \param defaultMax  default value for the upper distance bounds

    */
    void initBoundsMat(DistGeom::BoundsMatrix *mmat,double defaultMin=0.0,
		       double defaultMax=1000.0);
    void initBoundsMat(DistGeom::BoundsMatPtr mmat,double defaultMin=0.0,
		       double defaultMax=1000.0);

    //! Set upper and lower distance bounds between atoms in a molecule based on topology
    /*!   
      This consists of setting 1-2, 1-3 and 1-4 distance based on bond lengths,
      bond angles and torsion angle ranges. Optionally 1-5 bounds can also be set,
      in particular, for path that contain rigid 1-4 paths. 

      The final step involves setting lower bound to the sum of the vdW radii for
      the remaining atom pairs.
      
      \param mol          The molecule of interest
      \param mmat         Bounds matrix to the bounds are written
      \param set15bounds  If true try to set 1-5 bounds also based on topology
      \param scaleVDW     If true scale the sum of the vdW radii while setting lower bounds
                          so that a smaller value (0.7*(vdw1 + vdw2) ) is used for paths
			  that are less five bonds apart.

      <b>Note</b>
      For some strained systems the bounds matrix resulting from setting 1-5 bounds may
      fail triangle smoothing. In these cases it is recommended to back out and
      recompute the bounds matrix with no 1-5 bounds and with vdW scaling. 
    */
    void setTopolBounds(const ROMol &mol, DistGeom::BoundsMatPtr mmat,
			bool set15bounds=true, bool scaleVDW=false);
  }
}
#endif
