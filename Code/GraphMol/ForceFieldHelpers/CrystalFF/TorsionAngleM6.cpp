// $Id$
//
//  Copyright (C) 2015 Sereina Riniker
//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "TorsionAngleM6.h"
#include <cmath>
#include <utility>
#include <vector>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/TorsionAngle.h>

namespace ForceFields {
namespace CrystalFF {
using namespace MMFF;

double calcTorsionEnergyM6(const std::vector<double> &V,
                           const std::vector<int> &signs, const double cosPhi) {
  double cosPhi2 = cosPhi * cosPhi;
  double cosPhi3 = cosPhi * cosPhi2;
  double cosPhi4 = cosPhi * cosPhi3;
  double cosPhi5 = cosPhi * cosPhi4;
  double cosPhi6 = cosPhi * cosPhi5;

  double cos2Phi = 2.0 * cosPhi2 - 1.0;
  double cos3Phi = 4.0 * cosPhi3 - 3.0 * cosPhi;
  double cos4Phi = 8.0 * cosPhi4 - 8.0 * cosPhi2 + 1.0;
  double cos5Phi = 16.0 * cosPhi5 - 20.0 * cosPhi3 + 5.0 * cosPhi;
  double cos6Phi = 32.0 * cosPhi6 - 48.0 * cosPhi4 + 18.0 * cosPhi2 - 1.0;

  return (
      V[0] * (1.0 + signs[0] * cosPhi) + V[1] * (1.0 + signs[1] * cos2Phi) +
      V[2] * (1.0 + signs[2] * cos3Phi) + V[3] * (1.0 + signs[3] * cos4Phi) +
      V[4] * (1.0 + signs[4] * cos5Phi) + V[5] * (1.0 + signs[5] * cos6Phi));
}

TorsionAngleContribM6::TorsionAngleContribM6(
    ForceFields::ForceField *owner, unsigned int idx1, unsigned int idx2,
    unsigned int idx3, unsigned int idx4, std::vector<double> V,
    std::vector<int> signs)
    : ForceFieldContrib(owner),
      d_at1Idx(idx1),
      d_at2Idx(idx2),
      d_at3Idx(idx3),
      d_at4Idx(idx4),
      d_V(std::move(V)),
      d_sign(std::move(signs)) {
  PRECONDITION(owner, "bad owner");
  PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
                   (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
               "degenerate points");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());
  URANGE_CHECK(idx4, owner->positions().size());
};

double TorsionAngleContribM6::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                         pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                         pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                         pos[3 * d_at3Idx + 2]);
  RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
                         pos[3 * d_at4Idx + 2]);

  return calcTorsionEnergyM6(
      d_V, d_sign, Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
}

void TorsionAngleContribM6::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                         pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                         pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                         pos[3 * d_at3Idx + 2]);
  RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
                         pos[3 * d_at4Idx + 2]);
  double *g[4] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
                  &(grad[3 * d_at3Idx]), &(grad[3 * d_at4Idx])};

  RDGeom::Point3D r[4] = {iPoint - jPoint, kPoint - jPoint, jPoint - kPoint,
                          lPoint - kPoint};
  RDGeom::Point3D t[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
  double d[2] = {t[0].length(), t[1].length()};
  if (isDoubleZero(d[0]) || isDoubleZero(d[1])) {
    return;
  }
  t[0] /= d[0];
  t[1] /= d[1];
  double cosPhi = t[0].dotProduct(t[1]);
  clipToOne(cosPhi);
  double sinPhiSq = 1.0 - cosPhi * cosPhi;
  double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
  double cosPhi2 = cosPhi * cosPhi;
  double cosPhi3 = cosPhi * cosPhi2;
  double cosPhi4 = cosPhi * cosPhi3;
  double cosPhi5 = cosPhi * cosPhi4;
  // dE/dPhi is independent of cartesians:
  double dE_dPhi =
      (-d_V[0] * d_sign[0] * sinPhi -
       2.0 * d_V[1] * d_sign[1] * (2.0 * cosPhi * sinPhi) -
       3.0 * d_V[2] * d_sign[2] * (4.0 * cosPhi2 * sinPhi - sinPhi) -
       4.0 * d_V[3] * d_sign[3] *
           (8.0 * cosPhi3 * sinPhi - 4.0 * cosPhi * sinPhi) -
       5.0 * d_V[4] * d_sign[4] *
           (16.0 * cosPhi4 * sinPhi - 12.0 * cosPhi2 * sinPhi + sinPhi) -
       6.0 * d_V[4] * d_sign[4] *
           (32.0 * cosPhi5 * sinPhi - 32.0 * cosPhi3 * sinPhi + 6.0 * sinPhi));

  // FIX: use a tolerance here
  // this is hacky, but it's per the
  // recommendation from Niketic and Rasmussen:
  double sinTerm =
      -dE_dPhi * (isDoubleZero(sinPhi) ? (1.0 / cosPhi) : (1.0 / sinPhi));

  Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
}
}  // namespace CrystalFF
}  // namespace ForceFields
