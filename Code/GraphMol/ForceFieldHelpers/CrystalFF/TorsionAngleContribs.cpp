//
//  Copyright (C) 2024 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "TorsionAngleContribs.h"
#include <cmath>
#include <utility>
#include <vector>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/TorsionAngle.h>

namespace ForceFields {
namespace CrystalFF {
namespace {
double calcTorsionEnergyM6(const std::vector<double> &forceConstants,
                           const std::vector<int> &signs, const double cosPhi) {
  const double cosPhi2 = cosPhi * cosPhi;
  const double cosPhi3 = cosPhi * cosPhi2;
  const double cosPhi4 = cosPhi * cosPhi3;
  const double cosPhi5 = cosPhi * cosPhi4;
  const double cosPhi6 = cosPhi * cosPhi5;

  const double cos2Phi = 2.0 * cosPhi2 - 1.0;
  const double cos3Phi = 4.0 * cosPhi3 - 3.0 * cosPhi;
  const double cos4Phi = 8.0 * cosPhi4 - 8.0 * cosPhi2 + 1.0;
  const double cos5Phi = 16.0 * cosPhi5 - 20.0 * cosPhi3 + 5.0 * cosPhi;
  const double cos6Phi = 32.0 * cosPhi6 - 48.0 * cosPhi4 + 18.0 * cosPhi2 - 1.0;

  return (forceConstants[0] * (1.0 + signs[0] * cosPhi) +
          forceConstants[1] * (1.0 + signs[1] * cos2Phi) +
          forceConstants[2] * (1.0 + signs[2] * cos3Phi) +
          forceConstants[3] * (1.0 + signs[3] * cos4Phi) +
          forceConstants[4] * (1.0 + signs[4] * cos5Phi) +
          forceConstants[5] * (1.0 + signs[5] * cos6Phi));
}
}  // namespace

TorsionAngleContribs::TorsionAngleContribs(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void TorsionAngleContribs::addContrib(unsigned int idx1, unsigned int idx2,
                                      unsigned int idx3, unsigned int idx4,
                                      std::vector<double> forceConstants,
                                      std::vector<int> signs) {
  PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
                   (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
               "degenerate points");
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  URANGE_CHECK(idx4, dp_forceField->positions().size());
  d_contribs.emplace_back(idx1, idx2, idx3, idx4, std::move(forceConstants),
                          std::move(signs));
}

double TorsionAngleContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double accum = 0.0;
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
                                 pos[3 * contrib.idx1 + 1],
                                 pos[3 * contrib.idx1 + 2]);
    const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
                                 pos[3 * contrib.idx2 + 1],
                                 pos[3 * contrib.idx2 + 2]);
    const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
                                 pos[3 * contrib.idx3 + 1],
                                 pos[3 * contrib.idx3 + 2]);
    const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
                                 pos[3 * contrib.idx4 + 1],
                                 pos[3 * contrib.idx4 + 2]);
    accum += calcTorsionEnergyM6(
        contrib.forceConstants, contrib.signs,
        MMFF::Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint));
  }
  return accum;
}

void TorsionAngleContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D iPoint(pos[3 * contrib.idx1],
                                 pos[3 * contrib.idx1 + 1],
                                 pos[3 * contrib.idx1 + 2]);
    const RDGeom::Point3D jPoint(pos[3 * contrib.idx2],
                                 pos[3 * contrib.idx2 + 1],
                                 pos[3 * contrib.idx2 + 2]);
    const RDGeom::Point3D kPoint(pos[3 * contrib.idx3],
                                 pos[3 * contrib.idx3 + 1],
                                 pos[3 * contrib.idx3 + 2]);
    const RDGeom::Point3D lPoint(pos[3 * contrib.idx4],
                                 pos[3 * contrib.idx4 + 1],
                                 pos[3 * contrib.idx4 + 2]);
    double *g[4] = {&(grad[3 * contrib.idx1]), &(grad[3 * contrib.idx2]),
                    &(grad[3 * contrib.idx3]), &(grad[3 * contrib.idx4])};

    RDGeom::Point3D r[4] = {iPoint - jPoint, kPoint - jPoint, jPoint - kPoint,
                            lPoint - kPoint};
    RDGeom::Point3D t[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
    double d[2] = {t[0].length(), t[1].length()};
    if (MMFF::isDoubleZero(d[0]) || MMFF::isDoubleZero(d[1])) {
      return;
    }
    t[0] /= d[0];
    t[1] /= d[1];
    double cosPhi = t[0].dotProduct(t[1]);
    cosPhi = std::clamp(cosPhi, -1.0, 1.0);
    const double sinPhiSq = 1.0 - cosPhi * cosPhi;
    const double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
    const double cosPhi2 = cosPhi * cosPhi;
    const double cosPhi3 = cosPhi * cosPhi2;
    const double cosPhi4 = cosPhi * cosPhi3;
    const double cosPhi5 = cosPhi * cosPhi4;
    // dE/dPhi is independent of cartesians:
    const double dE_dPhi =
        (-contrib.forceConstants[0] * contrib.signs[0] * sinPhi -
         2.0 * contrib.forceConstants[1] * contrib.signs[1] *
             (2.0 * cosPhi * sinPhi) -
         3.0 * contrib.forceConstants[2] * contrib.signs[2] *
             (4.0 * cosPhi2 * sinPhi - sinPhi) -
         4.0 * contrib.forceConstants[3] * contrib.signs[3] *
             (8.0 * cosPhi3 * sinPhi - 4.0 * cosPhi * sinPhi) -
         5.0 * contrib.forceConstants[4] * contrib.signs[4] *
             (16.0 * cosPhi4 * sinPhi - 12.0 * cosPhi2 * sinPhi + sinPhi) -
         6.0 * contrib.forceConstants[4] * contrib.signs[4] *
             (32.0 * cosPhi5 * sinPhi - 32.0 * cosPhi3 * sinPhi +
              6.0 * sinPhi));

    // FIX: use a tolerance here
    // this is hacky, but it's per the
    // recommendation from Niketic and Rasmussen:
    double sinTerm = -dE_dPhi * (MMFF::isDoubleZero(sinPhi) ? (1.0 / cosPhi)
                                                            : (1.0 / sinPhi));

    MMFF::Utils::calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
  }
}

}  // namespace CrystalFF
}  // namespace ForceFields