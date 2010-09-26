//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef __RD_DIST_UTILS_H__
#define __RD_DIST_UTILS_H__

#include <math.h>
#include "point.h"
#include "Transform3D.h"
#include "Transform.h"

namespace RDGeom {

  /*! \brief Compute the 13 distance between points give the 12 distances
   *  and the angle between the axes.
   */
  inline double compute13Dist(double d1, double d2, double angle) {
    double res = d1*d1 + d2*d2 - 2*d1*d2*cos(angle);
    return sqrt(res);
  }

  /*! \brief Compute the 14 distances give the 12 distance and the angles
   *
   *   This is computed by aligning the d2 axis with the x-axis (with atom 2 at 
   *   the origin. Atom 1 is made to lie int he xy-plane with a +ve y-coordinate
   *   and finally the coordinates for atom 4 are computed. 
   *
   * ARGUMENTS:
   *   d1 - distance between atoms 1 and 2
   *   d2 - distance between atoms 2 and 3
   *   d3 - distance between atoms 3 and 4
   *   ang12 - angle between the axes d1 and d2
   *   ang23 - angle between the axes d2 and d3
   *   torAng - torsional agnle of the axis d2
   *
   * NOTE:
   *   we are specifically calling this function compute14Dist3D because
   *   we assume the torsional angle can take any value including 0 and 180 deg.
   *   However, if using either 0 or 180 as the torsional angle (which is often 
   *   the case) the user is recommended to use the specialized functions below
   *   instead of this function; they will be speedier.
   */
  inline double compute14Dist3D(double d1, double d2, double d3, 
                                double ang12, double ang23, double torAng) {
    // location of atom1
    Point3D p1(d1*cos(ang12), d1*sin(ang12), 0.0);

    // location of atom 4 if the rosion angle was 0
    Point3D p4(d2-d3*cos(ang23), d3*sin(ang23), 0.0);
    
    // now we will rotate p4 about the x-axis by the desired torsion angle
    Transform3D trans;
    trans.SetRotation(torAng, X_Axis);
    trans.TransformPoint(p4);

    // find the distance
    p4 -= p1;
    return p4.length();
  }

  /*! \brief Compute the 14 distances give the 12 distance and bond angle
   *  for cis configuration
   *
   *  This is simply a special case of the above function compute14Dist3D;
   *  with torsion angle set to 0. However, this function should be speedier
   */
  inline double compute14DistCis(double d1, double d2, double d3, 
                                double ang12, double ang23) {
    double dx = d2 - d3*cos(ang23) - d1*cos(ang12);
    double dy = d3*sin(ang23) - d1*sin(ang12);
    double res = dx*dx + dy*dy;
    return sqrt(res);
  }

  
   /*! \brief Compute the 14 distances give the 12 distance and bond angle
   *  for trans configuration
   *
   *  This is simply a special case of the above function compute14Dist3D;
   *  with torsion angle set to 180. However, this function should be speedier
   */
  inline double compute14DistTrans(double d1, double d2, double d3, 
                                double ang12, double ang23) {
    double dx = d2 - d3*cos(ang23) - d1*cos(ang12);
    double dy = d3*sin(ang23) + d1*sin(ang12);
    double res = dx*dx + dy*dy;
    return sqrt(res);
  } 
}

#endif
