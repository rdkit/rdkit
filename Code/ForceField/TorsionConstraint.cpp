//
//  Copyright (C) 2013-2024 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ForceField.h"
#include "TorsionConstraint.h"
#include <RDGeneral/BoostStartInclude.h>
#include <RDGeneral/BoostEndInclude.h>
#include <RDGeneral/Invariant.h>

#include <boost/math/special_functions/round.hpp>

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace ForceFields {
constexpr double RAD2DEG = 180.0 / M_PI;

inline void checkPrecondition(const ForceField *owner, unsigned int idx1,
                              unsigned int idx2, unsigned int idx3,
                              unsigned int idx4, double minDihedralDeg,
                              double maxDihedralDeg) {
  PRECONDITION(owner, "bad owner");
  PRECONDITION(!(minDihedralDeg > maxDihedralDeg),
               "minDihedralDeg must be <= maxDihedralDeg");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());
  URANGE_CHECK(idx4, owner->positions().size());
}

double TorsionConstraintContrib::computeDihedralTerm(double dihedral) const {
  double dihedralTarget = dihedral;
  if (!(dihedral > d_minDihedralDeg && dihedral < d_maxDihedralDeg) &&
      !(dihedral > d_minDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg) &&
      !(dihedral < d_maxDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg)) {
    double dihedralMinTarget = dihedral - d_minDihedralDeg;
    RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralMinTarget);
    double dihedralMaxTarget = dihedral - d_maxDihedralDeg;
    RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralMaxTarget);
    if (fabs(dihedralMinTarget) < fabs(dihedralMaxTarget)) {
      dihedralTarget = d_minDihedralDeg;
    } else {
      dihedralTarget = d_maxDihedralDeg;
    }
  }
  double dihedralTerm = dihedral - dihedralTarget;
  RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralTerm);
  return dihedralTerm;
}

void TorsionConstraintContrib::setParameters(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    unsigned int idx4, double minDihedralDeg, double maxDihedralDeg,
    double forceConst) {
  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  d_at4Idx = idx4;
  RDKit::ForceFieldsHelper::normalizeAngleDeg(minDihedralDeg);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(maxDihedralDeg);
  d_minDihedralDeg = minDihedralDeg;
  d_maxDihedralDeg = maxDihedralDeg;
  d_forceConstant = forceConst;
}

TorsionConstraintContrib::TorsionConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    unsigned int idx4, double minDihedralDeg, double maxDihedralDeg,
    double forceConst) {
  checkPrecondition(owner, idx1, idx2, idx3, idx4, minDihedralDeg,
                    maxDihedralDeg);
  setParameters(owner, idx1, idx2, idx3, idx4, minDihedralDeg, maxDihedralDeg,
                forceConst);
}

TorsionConstraintContrib::TorsionConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    unsigned int idx4, bool relative, double minDihedralDeg,
    double maxDihedralDeg, double forceConst) {
  checkPrecondition(owner, idx1, idx2, idx3, idx4, minDihedralDeg,
                    maxDihedralDeg);
  if (relative) {
    double dihedral;
    RDKit::ForceFieldsHelper::computeDihedral(owner->positions(), idx1, idx2,
                                              idx3, idx4, &dihedral);
    dihedral *= RAD2DEG;
    minDihedralDeg += dihedral;
    maxDihedralDeg += dihedral;
  }
  setParameters(owner, idx1, idx2, idx3, idx4, minDihedralDeg, maxDihedralDeg,
                forceConst);
}

double TorsionConstraintContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double dihedral;
  RDKit::ForceFieldsHelper::computeDihedral(pos, d_at1Idx, d_at2Idx, d_at3Idx,
                                            d_at4Idx, &dihedral);
  dihedral *= RAD2DEG;
  double dihedralTerm = computeDihedralTerm(dihedral);
  double res = d_forceConstant * dihedralTerm * dihedralTerm;

  return res;
}

void TorsionConstraintContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
                  &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};

  RDGeom::Point3D r[4];
  RDGeom::Point3D t[2];
  double d[2];
  double dihedral;
  RDKit::ForceFieldsHelper::computeDihedral(
      pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, &dihedral, nullptr, r, t, d);
  dihedral *= RAD2DEG;
  double dihedralTerm = computeDihedralTerm(dihedral);
  double dE_dPhi = 2.0 * RAD2DEG * d_forceConstant * dihedralTerm;

  double d23 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);
  RDGeom::Point3D r31(pos[3 * d_at3Idx] - pos[3 * d_at1Idx],
                      pos[3 * d_at3Idx + 1] - pos[3 * d_at1Idx + 1],
                      pos[3 * d_at3Idx + 2] - pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D r42(pos[3 * d_at4Idx] - pos[3 * d_at2Idx],
                      pos[3 * d_at4Idx + 1] - pos[3 * d_at2Idx + 1],
                      pos[3 * d_at4Idx + 2] - pos[3 * d_at2Idx + 2]);
  double prefactor = dE_dPhi / d23;
  RDGeom::Point3D tt[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
  RDGeom::Point3D dedt[2] = {
      tt[0].crossProduct(r[2]) / tt[0].lengthSq() * prefactor,
      tt[1].crossProduct(r[1]) / tt[1].lengthSq() * prefactor};
  RDGeom::Point3D dedp[4] = {
      r[2].crossProduct(dedt[0]),
      r31.crossProduct(dedt[0]) - r[3].crossProduct(dedt[1]),
      r[0].crossProduct(dedt[0]) + r42.crossProduct(dedt[1]),
      r[2].crossProduct(dedt[1])};
  for (unsigned int i = 0; i < 4; ++i) {
    g[i][0] += dedp[i].x;
    g[i][1] += dedp[i].y;
    g[i][2] += dedp[i].z;
  }
}
}  // namespace ForceFields
