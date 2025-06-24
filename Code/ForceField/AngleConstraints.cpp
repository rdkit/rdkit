//
//  Copyright (C) 2024 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "AngleConstraints.h"
#include "ForceField.h"
#include <cmath>
#include <RDGeneral/Invariant.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ForceFields {
constexpr double RAD2DEG = 180.0 / M_PI;

AngleConstraintContribs::AngleConstraintContribs(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void AngleConstraintContribs::addContrib(unsigned int idx1, unsigned int idx2,
                                         unsigned int idx3, double minAngleDeg,
                                         double maxAngleDeg,
                                         double forceConst) {
  RANGE_CHECK(0.0, minAngleDeg, 180.0);
  RANGE_CHECK(0.0, maxAngleDeg, 180.0);
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  PRECONDITION(maxAngleDeg >= minAngleDeg,
               "minAngleDeg must be <= maxAngleDeg");
  d_contribs.emplace_back(idx1, idx2, idx3, minAngleDeg, maxAngleDeg,
                          forceConst);
}
void AngleConstraintContribs::addContrib(unsigned int idx1, unsigned int idx2,
                                         unsigned int idx3, bool relative,
                                         double minAngleDeg, double maxAngleDeg,
                                         double forceConst) {
  const RDGeom::PointPtrVect &pos = dp_forceField->positions();
  URANGE_CHECK(idx1, pos.size());
  URANGE_CHECK(idx2, pos.size());
  URANGE_CHECK(idx3, pos.size());
  PRECONDITION(maxAngleDeg >= minAngleDeg,
               "minAngleDeg must be <= maxAngleDeg");
  if (relative) {
    const RDGeom::Point3D &p1 = *((RDGeom::Point3D *)pos[idx1]);
    const RDGeom::Point3D &p2 = *((RDGeom::Point3D *)pos[idx2]);
    const RDGeom::Point3D &p3 = *((RDGeom::Point3D *)pos[idx3]);
    const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
    const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
                                 std::max(1.0e-5, r[1].lengthSq())};
    double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
    cosTheta = std::clamp(cosTheta, -1.0, 1.0);
    const double angle = RAD2DEG * acos(cosTheta);
    minAngleDeg += angle;
    maxAngleDeg += angle;
  }
  RANGE_CHECK(0.0, minAngleDeg, 180.0);
  RANGE_CHECK(0.0, maxAngleDeg, 180.0);
  d_contribs.emplace_back(idx1, idx2, idx3, minAngleDeg, maxAngleDeg,
                          forceConst);
}

double AngleConstraintContribs::computeAngleTerm(
    const double &angle, const AngleConstraintContribsParams &contrib) const {
  double angleTerm = 0.0;
  if (angle < contrib.minAngle) {
    angleTerm = angle - contrib.minAngle;
  } else if (angle > contrib.maxAngle) {
    angleTerm = angle - contrib.maxAngle;
  }
  return angleTerm;
}

double AngleConstraintContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double accum = 0.0;
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
                             pos[3 * contrib.idx1 + 2]);
    const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
                             pos[3 * contrib.idx2 + 2]);
    const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
                             pos[3 * contrib.idx3 + 2]);
    const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
    const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
                                 std::max(1.0e-5, r[1].lengthSq())};
    double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
    cosTheta = std::clamp(cosTheta, -1.0, 1.0);
    const double angle = RAD2DEG * acos(cosTheta);
    const double angleTerm = computeAngleTerm(angle, contrib);
    accum += contrib.forceConstant * angleTerm * angleTerm;
  }
  return accum;
}

void AngleConstraintContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
                             pos[3 * contrib.idx1 + 2]);
    const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
                             pos[3 * contrib.idx2 + 2]);
    const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
                             pos[3 * contrib.idx3 + 2]);
    const RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
    const double rLengthSq[2] = {std::max(1.0e-5, r[0].lengthSq()),
                                 std::max(1.0e-5, r[1].lengthSq())};
    double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
    cosTheta = std::clamp(cosTheta, -1.0, 1.0);
    const double angle = RAD2DEG * acos(cosTheta);
    const double angleTerm = computeAngleTerm(angle, contrib);

    const double dE_dTheta = 2.0 * RAD2DEG * contrib.forceConstant * angleTerm;

    const RDGeom::Point3D rp = r[1].crossProduct(r[0]);
    const double prefactor = dE_dTheta / std::max(1.0e-5, rp.length());
    const double t[2] = {-prefactor / rLengthSq[0], prefactor / rLengthSq[1]};
    RDGeom::Point3D dedp[3];
    dedp[0] = r[0].crossProduct(rp) * t[0];
    dedp[2] = r[1].crossProduct(rp) * t[1];
    dedp[1] = -dedp[0] - dedp[2];
    double *g[3] = {&(grad[3 * contrib.idx1]), &(grad[3 * contrib.idx2]),
                    &(grad[3 * contrib.idx3])};
    for (unsigned int i = 0; i < 3; ++i) {
      g[i][0] += dedp[i].x;
      g[i][1] += dedp[i].y;
      g[i][2] += dedp[i].z;
    }
  }
}

}  // namespace ForceFields