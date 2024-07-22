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

double calcTorsionEnergy(const std::vector<double> &forceConstants,
                         const std::vector<int> &signs, const double cosPhi) {
  PRECONDITION(forceConstants.size() == signs.size(),
               "number of force constants need to match number of signs");
  PRECONDITION(forceConstants.size(),
               "number of force constants and signs need to be > 0");
  auto multiplicity = forceConstants.size();
  if (multiplicity == 1) {
    return forceConstants.front() * (1 + signs.front() * cosPhi);
  }
  auto cosNPhis = std::vector<double>(multiplicity + 1);
  double res = 0.0;
  cosNPhis[0] = 1.0;
  cosNPhis[1] = cosPhi;
  for (unsigned int i = 2; i < multiplicity + 1; ++i) {
    cosNPhis[i] = 2.0 * cosPhi * cosNPhis[i - 1] - cosNPhis[i - 2];
  }
  res = 0.0;
  for (unsigned int i = 0; i < multiplicity; i++) {
    res += forceConstants[i] * (1.0 + signs[i] * cosNPhis[i + 1]);
  }
  return res;
}

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
    RDGeom::Point3D iPoint(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
                           pos[3 * contrib.idx1 + 2]);
    RDGeom::Point3D jPoint(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
                           pos[3 * contrib.idx2 + 2]);
    RDGeom::Point3D kPoint(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
                           pos[3 * contrib.idx3 + 2]);
    RDGeom::Point3D lPoint(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
                           pos[3 * contrib.idx4 + 2]);
    accum += calcTorsionEnergy(
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
    RDGeom::Point3D iPoint(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
                           pos[3 * contrib.idx1 + 2]);
    RDGeom::Point3D jPoint(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
                           pos[3 * contrib.idx2 + 2]);
    RDGeom::Point3D kPoint(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
                           pos[3 * contrib.idx3 + 2]);
    RDGeom::Point3D lPoint(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
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
    MMFF::clipToOne(cosPhi);
    double sinPhiSq = 1.0 - cosPhi * cosPhi;
    double sinPhi = ((sinPhiSq > 0.0) ? sqrt(sinPhiSq) : 0.0);
    double cosPhi2 = cosPhi * cosPhi;
    double cosPhi3 = cosPhi * cosPhi2;
    double cosPhi4 = cosPhi * cosPhi3;
    double cosPhi5 = cosPhi * cosPhi4;
    // dE/dPhi is independent of cartesians:
    double dE_dPhi =
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