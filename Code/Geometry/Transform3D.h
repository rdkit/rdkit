//
//   Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef __RD_TRANSFORM3D_H__
#define __RD_TRANSFORM3D_H__

#include "Transform.h"

#include <Numerics/SquareMatrix.h>

namespace RDGeom {
class Point3D;
const unsigned int DIM_3D = 4;

class RDKIT_RDGEOMETRYLIB_EXPORT Transform3D
    : public RDNumeric::SquareMatrix<double> {
 public:
  //!  Constructor
  /*!
    Initialize to an identity matrix transformation.
    This is a 4x4 matrix that includes the rotation and translation parts
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

  Transform3D() : RDNumeric::SquareMatrix<double>(DIM_3D, 0.0) {
    for (unsigned int i = 0; i < DIM_3D; i++) {
      unsigned int id = i * (DIM_3D + 1);
      d_data[id] = 1.0;
    }
  }

  void setToIdentity();

  void TransformPoint(Point3D &pt) const;

  /*! \brief Set the translation vector
   */
  void SetTranslation(const Point3D &move);

  /*! \brief set the rotation matrix
   *
   * The rotation matrix is set to rotation by the specified angle
   * about the specified axis
   */
  void SetRotation(double angle, AxisType axis);

  /*! \brief set the rotation matrix
   *
   * The rotation matrix is set to rotation by the specified angle
   * about an arbitrary axis.
   * Note: if the axis is not a unit vector scaling will also occur.
   * This can be ensured by a call to Point3D#normalize() prior to calling
   * this method
   */
  void SetRotation(double angle, const Point3D &axis);
  void SetRotation(double cosT, double sinT, const Point3D &axis);

  //! Set the rotation matrix from a quaternion
  void SetRotationFromQuaternion(double quaternion[4]);

  //! Reflect the rotation
  void Reflect();

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
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Transform3D operator*(
    const RDGeom::Transform3D &t1, const RDGeom::Transform3D &t2);

/*! \brief Transform a point:
 *
 */
RDKIT_RDGEOMETRYLIB_EXPORT RDGeom::Point3D operator*(
    const RDGeom::Transform3D &t, const RDGeom::Point3D &pt);

#endif

