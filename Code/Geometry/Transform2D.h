//
//   Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_TRANSFORM2D_H__
#define __RD_TRANSFORM2D_H__

#include "Transform.h"
#include <Numerics/SquareMatrix.h>

namespace RDGeom {
class Point2D;
const unsigned int DIM_2D = 3;

class RDKIT_RDGEOMETRYLIB_EXPORT Transform2D
    : public RDNumeric::SquareMatrix<double> {
 public:
  //! \brief Constructor
  /*!
    Initialize to an identity matrix transformation
    This is a 3x3 matrix that includes the rotation and translation parts
    see Foley's "Introduction to Computer Graphics" for the representation

    Operator *= and = are provided by the parent class square matrix.
    Operator *= needs some explanation, since the order matters. This transform
    gets set to
    the combination other and the current state of this transform
    If this_old and this_new are the states of this object before and after this
    function
    we have
            this_new(point) = this_old(other(point))
  */
  Transform2D() : RDNumeric::SquareMatrix<double>(DIM_2D, 0.0) {
    for (unsigned int i = 0; i < DIM_2D; i++) {
      unsigned int id = i * (DIM_2D + 1);
      d_data[id] = 1.0;
    }
  }

  void setToIdentity();

  void TransformPoint(Point2D &pt) const;

  void SetTranslation(const Point2D &pt);

  /*! \brief Set the transform so that the specified points are aligned
   *
   * The resulting transformation will align pt1 with ref1, and rotation
   * pt2 such that the line betweem (pt1, pt2) will align with
   * with the line (ref1, ref2)
   */
  void SetTransform(const Point2D &ref1, const Point2D &ref2,
                    const Point2D &pt1, const Point2D &pt2);

  /*! \brief Set the trans form to counterclock wise rotation by the specified
   *value around point
   *
   * ARGUMENTS:
   *   - pt : point about which to rotate
   *   - angle : the angle by which to rotate, in radians
   */
  void SetTransform(const Point2D &pt, double angle);

 private:
};
}  // namespace RDGeom

/*! \brief Combine two transforms and return the results as a new transform
 *
 * The order is important here, on two transforms t1 and t2
 * t3 = t1*t2
 * The resulting transform t3 has the folliwng effect
 *  t3(point) = t1(t2(point))
 */
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Transform2D operator*(
    const RDGeom::Transform2D &t1, const RDGeom::Transform2D &t2);

#endif

