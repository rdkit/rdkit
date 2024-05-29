// $Id$
//
//  Copyright (C) 2013 Paolo Tosco
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "PositionConstraint.h"
#include <algorithm>
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace UFF {
PositionConstraintContrib::PositionConstraintContrib(ForceField *owner,
                                                     unsigned int idx,
                                                     double maxDispl,
                                                     double forceConst) {
  PRECONDITION(owner, "bad owner");
  const RDGeom::PointPtrVect &pos = owner->positions();
  URANGE_CHECK(idx, pos.size());

  dp_forceField = owner;
  d_atIdx = idx;
  d_maxDispl = maxDispl;
  d_pos0 = *((RDGeom::Point3D *)pos[idx]);
  d_forceConstant = forceConst;
}

double PositionConstraintContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  RDGeom::Point3D p(pos[3 * d_atIdx], pos[3 * d_atIdx + 1],
                    pos[3 * d_atIdx + 2]);
  double dist = (p - d_pos0).length();
  double distTerm = std::max(dist - d_maxDispl, 0.0);
  double res = 0.5 * d_forceConstant * distTerm * distTerm;

  return res;
}

void PositionConstraintContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  RDGeom::Point3D p(pos[3 * d_atIdx], pos[3 * d_atIdx + 1],
                    pos[3 * d_atIdx + 2]);
  double dist = (p - d_pos0).length();

  double preFactor = 0.0;
  if (dist > d_maxDispl) {
    preFactor = dist - d_maxDispl;
  } else {
    return;
  }
  preFactor *= d_forceConstant;

  for (unsigned int i = 0; i < 3; ++i) {
    double dGrad = preFactor * (p[i] - d_pos0[i]) / std::max(dist, 1.0e-8);
    grad[3 * d_atIdx + i] += dGrad;
  }
}
}  // namespace UFF
}  // namespace ForceFields
