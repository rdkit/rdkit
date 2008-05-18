//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef __RD_TRIANGLE_SMOOTH_H__
#define __RD_TRIANGLE_SMOOTH_H__

#include "BoundsMatrix.h"


namespace DistGeom {
  //! Smooth the upper and lower bound in a metric matrix so that triangle 
  //! inequality is not violated
  /*!
    This an implementation of the O(N^3) algorithm given on pages 252-253 of 
    "Distance Geometry and Molecular Conformation" by G.M.Crippen and T.F.Havel
    Research Studies Press, 1988. There are other (slightly) more implementations
    (see pages 301-302 in the above book), but that is for later

    \param boundsMat  A pointer to the distance bounds matrix

  */
  bool triangleSmoothBounds(BoundsMatrix *boundsMat);
  //! \overload
  bool triangleSmoothBounds(BoundsMatPtr boundsMat);
}

#endif

