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
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace UFF {
AngleConstraintContribs::AngleConstraintContribs(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void AngleConstraintContribs::addContrib(unsigned int idx1, unsigned int idx2,
                                         unsigned int idx3, double minAngleDeg,
                                         double maxAngleDeg,
                                         double forceConst) {
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  PRECONDITION(maxAngleDeg >= minAngleDeg,
               "minAngleDeg must be <= maxAngleDeg");
  RDKit::ForceFieldsHelper::normalizeAngleDeg(minAngleDeg);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(maxAngleDeg);
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
  minAngleDeg += angle;
  maxAngleDeg += angle;
  RDKit::ForceFieldsHelper::normalizeAngleDeg(minAngleDeg);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(maxAngleDeg);
  d_contribs.emplace_back(idx1, idx2, idx3, minAngleDeg, maxAngleDeg,
                          forceConst);
}
std::tuple<RDGeom::Point3D, RDGeom::Point3D, RDGeom::Point3D>
AngleConstraintContribs::getPoints(
    double *pos, const AngleConstraintContribsParams &contrib) const {
  RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
                     pos[3 * contrib.idx1 + 2]);
  RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
                     pos[3 * contrib.idx2 + 2]);
  RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
                     pos[3 * contrib.idx3 + 2]);
  return std::make_tuple(p1, p2, p3);
}

double AngleConstraintContribs::computeAngleTerm(
    double angle, const AngleConstraintContribsParams &contrib) const {
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
    auto [p1, p2, p3] = getPoints(pos, contrib);
    RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
    double rLengthSq[2] = {(std::max)(1.0e-5, r[0].lengthSq()),
                           (std::max)(1.0e-5, r[1].lengthSq())};
    double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
    clipToOne(cosTheta);
    double angle = RAD2DEG * acos(cosTheta);
    double angleTerm = computeAngleTerm(angle, contrib);
    accum += contrib.forceConstant * angleTerm * angleTerm;
  }
  return accum;
}

void AngleConstraintContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  for (const auto &contrib : d_contribs) {
    auto [p1, p2, p3] = getPoints(pos, contrib);
    double *g[3] = {&(grad[3 * contrib.idx1]), &(grad[3 * contrib.idx2]),
                    &(grad[3 * contrib.idx3])};
    RDGeom::Point3D r[2] = {p1 - p2, p3 - p2};
    double rLengthSq[2] = {(std::max)(1.0e-5, r[0].lengthSq()),
                           (std::max)(1.0e-5, r[1].lengthSq())};
    double cosTheta = r[0].dotProduct(r[1]) / sqrt(rLengthSq[0] * rLengthSq[1]);
    clipToOne(cosTheta);

    double angle = RAD2DEG * acos(cosTheta);
    double angleTerm = computeAngleTerm(angle, contrib);

    double dE_dTheta = 2.0 * RAD2DEG * contrib.forceConstant * angleTerm;

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
}

}  // namespace UFF
}  // namespace ForceFields