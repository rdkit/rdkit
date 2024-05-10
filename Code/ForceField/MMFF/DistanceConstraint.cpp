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
#include "DistanceConstraint.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace MMFF {
DistanceConstraintContrib::DistanceConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, double minLen,
    double maxLen, double forceConst) {
  PRECONDITION(owner, "bad owner");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  PRECONDITION(maxLen >= minLen, "bad bounds");

  dp_forceField = owner;
  d_end1Idx = idx1;
  d_end2Idx = idx2;
  d_minLen = minLen;
  d_maxLen = maxLen;
  d_forceConstant = forceConst;
}

DistanceConstraintContrib::DistanceConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, bool relative,
    double minLen, double maxLen, double forceConst) {
  PRECONDITION(owner, "bad owner");
  const RDGeom::PointPtrVect &pos = owner->positions();
  URANGE_CHECK(idx1, pos.size());
  URANGE_CHECK(idx2, pos.size());
  PRECONDITION(maxLen >= minLen, "bad bounds");

  double dist = 0.0;
  if (relative) {
    RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
    RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
    dist = (p1 - p2).length();
  }
  dp_forceField = owner;
  d_end1Idx = idx1;
  d_end2Idx = idx2;
  d_minLen = std::max(dist + minLen, 0.0);
  d_maxLen = std::max(dist + maxLen, 0.0);
  d_forceConstant = forceConst;
}

double DistanceConstraintContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);
  double distTerm = 0.0;
  if (dist < d_minLen) {
    distTerm = d_minLen - dist;
  } else if (dist > d_maxLen) {
    distTerm = dist - d_maxLen;
  }
  double res = 0.5 * d_forceConstant * distTerm * distTerm;

  return res;
}

void DistanceConstraintContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  double dist = dp_forceField->distance(d_end1Idx, d_end2Idx, pos);

  double preFactor = 0.0;
  if (dist < d_minLen) {
    preFactor = dist - d_minLen;
  } else if (dist > d_maxLen) {
    preFactor = dist - d_maxLen;
  } else {
    return;
  }
  preFactor *= d_forceConstant;

  double *end1Coords = &(pos[3 * d_end1Idx]);
  double *end2Coords = &(pos[3 * d_end2Idx]);
  for (unsigned int i = 0; i < 3; ++i) {
    double dGrad =
        preFactor * (end1Coords[i] - end2Coords[i]) / std::max(dist, 1.0e-8);
    grad[3 * d_end1Idx + i] += dGrad;
    grad[3 * d_end2Idx + i] -= dGrad;
  }
}
}  // namespace MMFF
}  // namespace ForceFields
