//  $Id$
//
//   Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <cstring>
#include "Transform.h"
#include "Transform2D.h"
#include <cmath>
#include "point.h"

namespace RDGeom {

void Transform2D::setToIdentity() {
  double *data = d_data.get();
  memset(static_cast<void *>(data), 0, d_dataSize * sizeof(double));
  for (unsigned int i = 0; i < DIM_2D; i++) {
    unsigned int id = i * (DIM_2D + 1);
    data[id] = 1.0;
  }
}

void Transform2D::TransformPoint(Point2D &pt) const {
  double *data = d_data.get();
  double x = data[0] * pt.x + data[1] * pt.y + data[2];
  double y = data[3] * pt.x + data[4] * pt.y + data[5];

  pt.x = x;
  pt.y = y;
}

void Transform2D::SetTranslation(const Point2D &pt) {
  unsigned int i = DIM_2D - 1;
  double *data = d_data.get();
  data[i] = pt.x;
  i += DIM_2D;
  data[i] = pt.y;
  i += DIM_2D;
  data[i] = 1.0;
}

void Transform2D::SetTransform(const Point2D &pt, double angle) {
  this->setToIdentity();

  Transform2D trans1;
  trans1.SetTranslation(-pt);
  double *data = d_data.get();
  // set the rotation
  data[0] = cos(angle);
  data[1] = -sin(angle);
  data[3] = sin(angle);
  data[4] = cos(angle);

  (*this) *= trans1;

  // translation back to the original coordinate
  Transform2D trans2;
  trans2.SetTranslation(pt);
  trans2 *= (*this);

  // now combine them
  this->assign(trans2);
}

void Transform2D::SetTransform(const Point2D &ref1, const Point2D &ref2,
                               const Point2D &pt1, const Point2D &pt2) {
  // compute the angle between the two vectors
  Point2D rvec = ref2 - ref1;
  Point2D pvec = pt2 - pt1;

  double dp = rvec.dotProduct(pvec);
  double lp = (rvec.length()) * (pvec.length());
  if (lp <= 0.0) {
    this->setToIdentity();
    return;
  }

  double cval = dp / lp;
  if (cval < -1.0) {
    cval = -1.0;
  } else if (cval > 1.0) {
    cval = 1.0;
  }

  double ang = acos(cval);

  // figure out if we have to do clock wise or anti clock wise rotation
  double cross = (pvec.x) * (rvec.y) - (pvec.y) * (rvec.x);

  if (cross < 0.0) {
    ang *= -1.0;
  }
  this->setToIdentity();
  // set the rotation
  double *data = d_data.get();
  data[0] = cos(ang);
  data[1] = -sin(ang);
  data[3] = sin(ang);
  data[4] = cos(ang);

  // apply this rotation to pt1 and compute the translation
  Point2D npt1 = pt1;
  this->TransformPoint(npt1);
  data[DIM_2D - 1] = ref1.x - npt1.x;
  data[2 * DIM_2D - 1] = ref1.y - npt1.y;
}
}  // namespace RDGeom

RDGeom::Transform2D operator*(const RDGeom::Transform2D &t1,
                              const RDGeom::Transform2D &t2) {
  RDGeom::Transform2D res;
  RDNumeric::multiply(t1, t2, res);
  return res;
};
