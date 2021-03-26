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
#include "AngleConstraint.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace UFF {
AngleConstraintContrib::AngleConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    double minAngleDeg, double maxAngleDeg, double forceConst) {
  PRECONDITION(owner, "bad owner");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());
  PRECONDITION(!(minAngleDeg > maxAngleDeg),
               "minAngleDeg must be <= maxAngleDeg");
  RDKit::ForceFieldsHelper::normalizeAngleDeg(minAngleDeg);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(maxAngleDeg);

  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  d_minAngleDeg = minAngleDeg;
  d_maxAngleDeg = maxAngleDeg;
  d_forceConstant = forceConst;
}

AngleConstraintContrib::AngleConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    bool relative, double minAngleDeg, double maxAngleDeg, double forceConst) {
  PRECONDITION(owner, "bad owner");
  const RDGeom::PointPtrVect &pos = owner->positions();
  URANGE_CHECK(idx1, pos.size());
  URANGE_CHECK(idx2, pos.size());
  URANGE_CHECK(idx3, pos.size());
  PRECONDITION(!(minAngleDeg > maxAngleDeg),
               "minAngleDeg must be <= maxAngleDeg");

  double angle = 0.0;
  if (relative) {
    RDGeom::Point3D p1 = *((RDGeom::Point3D *)pos[idx1]);
    RDGeom::Point3D p2 = *((RDGeom::Point3D *)pos[idx2]);
    RDGeom::Point3D p3 = *((RDGeom::Point3D *)pos[idx3]);
    double dist1 = (p1 - p2).length();
    double dist2 = (p3 - p2).length();
    RDGeom::Point3D p12 = (p1 - p2) / dist1;
    RDGeom::Point3D p32 = (p3 - p2) / dist2;
    double cosTheta = p12.dotProduct(p32);
    clipToOne(cosTheta);
    angle = RAD2DEG * acos(cosTheta);
  }
  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  minAngleDeg += angle;
  maxAngleDeg += angle;
  RDKit::ForceFieldsHelper::normalizeAngleDeg(minAngleDeg);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(maxAngleDeg);
  d_minAngleDeg = minAngleDeg;
  d_maxAngleDeg = maxAngleDeg;
  d_forceConstant = forceConst;
}

double AngleConstraintContrib::computeAngleTerm(double angle) const {
  double angleTerm = 0.0;
  if (angle < d_minAngleDeg) {
    angleTerm = angle - d_minAngleDeg;
  } else if (angle > d_maxAngleDeg) {
    angleTerm = angle - d_maxAngleDeg;
  }
  return angleTerm;
}

double AngleConstraintContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);

  RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
  double rLengthSq[2] = {(std::max)(1.0e-5, r[0].lengthSq()),
                         (std::max)(1.0e-5, r[1].lengthSq())};
  double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
  clipToOne(cosTheta);
  double angle = RAD2DEG * acos(cosTheta);
  double angleTerm = computeAngleTerm(angle);
  double res = d_forceConstant * angleTerm * angleTerm;

  return res;
}

void AngleConstraintContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);

  double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
                  &(grad[3 * d_at3Idx])};
  RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
  double rLengthSq[2] = {(std::max)(1.0e-5, r[0].lengthSq()),
                         (std::max)(1.0e-5, r[1].lengthSq())};
  double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
  clipToOne(cosTheta);

  double angle = RAD2DEG * acos(cosTheta);
  double angleTerm = computeAngleTerm(angle);

  double dE_dTheta = 2.0 * RAD2DEG * d_forceConstant * angleTerm;

  RDGeom::Point3D rp = r[1].crossProduct(r[0]);
  double prefactor = dE_dTheta / (std::max)(1.0e-5, rp.length());
  double t[2] = {-prefactor / rLengthSq[0], prefactor / rLengthSq[1]};
  RDGeom::Point3D dedp[3];
  dedp[0] = r[0].crossProduct(rp) * t[0];
  dedp[2] = r[1].crossProduct(rp) * t[1];
  dedp[1] = -dedp[0] - dedp[2];
  for (unsigned int i = 0; i < 3; ++i) {
    g[i][0] += dedp[i].x;
    g[i][1] += dedp[i].y;
    g[i][2] += dedp[i].z;
  }
}
}  // namespace UFF
}  // namespace ForceFields
