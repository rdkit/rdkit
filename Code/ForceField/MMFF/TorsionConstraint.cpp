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
#include "TorsionAngle.h"
#include "TorsionConstraint.h"
#include "Params.h"
#include <cmath>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/math/special_functions/round.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace MMFF {
inline void checkPrecondition(const ForceField *owner, unsigned int idx1,
    unsigned int idx2, unsigned int idx3, unsigned int idx4) {
  PRECONDITION(owner, "bad owner");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());
  URANGE_CHECK(idx4, owner->positions().size());
}

double TorsionConstraintContrib::computeDihedralTerm(double dihedral) const {
  double dihedralTarget = dihedral;
  if (!(dihedral > d_minDihedralDeg && dihedral < d_maxDihedralDeg)
      && !(dihedral > d_minDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg)
      && !(dihedral < d_maxDihedralDeg && d_minDihedralDeg > d_maxDihedralDeg)) {
    double minDihedralTarget = dihedral - d_minDihedralDeg;
    double maxDihedralTarget = dihedral - d_maxDihedralDeg;
    RDKit::ForceFieldsHelper::normalizeAngleDeg(minDihedralTarget);
    RDKit::ForceFieldsHelper::normalizeAngleDeg(maxDihedralTarget);
    dihedralTarget = (fabs(minDihedralTarget) < fabs(maxDihedralTarget)
      ? d_minDihedralDeg : d_maxDihedralDeg);
  }
  double dihedralTerm = dihedral - dihedralTarget;
  RDKit::ForceFieldsHelper::normalizeAngleDeg(dihedralTerm);
  return dihedralTerm;
}

void TorsionConstraintContrib::setParameters(ForceField *owner, unsigned int idx1,
    unsigned int idx2, unsigned int idx3, unsigned int idx4, double minDihedralDeg,
    double maxDihedralDeg, double forceConst) {
  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  d_at4Idx = idx4;
  d_minDihedralDeg = minDihedralDeg;
  d_maxDihedralDeg = maxDihedralDeg;
  d_forceConstant = forceConst;
}

TorsionConstraintContrib::TorsionConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    unsigned int idx4, double minDihedralDeg, double maxDihedralDeg,
    double forceConst) {
  checkPrecondition(owner, idx1, idx2, idx3, idx4);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(minDihedralDeg);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(maxDihedralDeg);
  setParameters(owner, idx1, idx2, idx3, idx4,
    minDihedralDeg, maxDihedralDeg, forceConst);
}

TorsionConstraintContrib::TorsionConstraintContrib(
    ForceField *owner, unsigned int idx1, unsigned int idx2, unsigned int idx3,
    unsigned int idx4, bool relative, double minDihedralDeg,
    double maxDihedralDeg, double forceConst) {
  checkPrecondition(owner, idx1, idx2, idx3, idx4);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(minDihedralDeg);
  RDKit::ForceFieldsHelper::normalizeAngleDeg(maxDihedralDeg);
  if (relative) {
    double dihedral;
    RDKit::ForceFieldsHelper::computeDihedral(
        owner->positions(), idx1, idx2, idx3, idx4, &dihedral);
    dihedral *= RAD2DEG;
    minDihedralDeg += dihedral;
    maxDihedralDeg += dihedral;
    RDKit::ForceFieldsHelper::normalizeAngleDeg(minDihedralDeg);
    RDKit::ForceFieldsHelper::normalizeAngleDeg(maxDihedralDeg);
  }
  setParameters(owner, idx1, idx2, idx3, idx4,
    minDihedralDeg, maxDihedralDeg, forceConst);
}

double TorsionConstraintContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double dihedral;
  RDKit::ForceFieldsHelper::computeDihedral(
    pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, &dihedral);
  dihedral *= RAD2DEG;
  double dihedralTerm = computeDihedralTerm(dihedral);
  static double const c = 0.5 * DEG2RAD * DEG2RAD;
  double res = c * d_forceConstant * dihedralTerm * dihedralTerm;

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
  double cosPhi;
  double dihedral;
  RDKit::ForceFieldsHelper::computeDihedral(
    pos, d_at1Idx, d_at2Idx, d_at3Idx, d_at4Idx, &dihedral, &cosPhi, r, t, d);
  dihedral *= RAD2DEG;
  double sinPhiSq = 1.0 - cosPhi * cosPhi;
  double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
  // dE/dPhi is independent of cartesians:
  double dihedralTerm = computeDihedralTerm(dihedral);
  double dE_dPhi = DEG2RAD * d_forceConstant * dihedralTerm;

  // FIX: use a tolerance here
  // this is hacky, but it's per the
  // recommendation from Niketic and Rasmussen:
  double sinTerm =
      -dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));
  Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
}
}
}
