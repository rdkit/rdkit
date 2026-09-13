//
//  Copyright (C) 2026 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "GaussianTorsionAngleContribs.h"
#include <cmath>
#include <utility>
#include <vector>
#include <numbers>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <ForceField/MMFF/Params.h>
#include <ForceField/MMFF/TorsionAngle.h>

namespace ForceFields {
namespace CrystalFF {
constexpr double twoPI = 2 * std::numbers::pi;
constexpr double KT = 2.479;

inline double interpolate(const std::vector<double> &table, double phi) {
  const double norm = phi * (table.size() - 1) / std::numbers::pi;
  const std::size_t idx =
      std::min(static_cast<std::size_t>(norm), table.size() - 2);
  const double t = norm - static_cast<double>(idx);
  return table[idx] * (1.0 - t) + table[idx + 1] * t;
}
double getEnergy(const std::vector<double> &heights,
                 const std::vector<double> &positions,
                 const std::vector<double> &widths, const double phi) {
  double result = 0.0;
  for (size_t i = 0; i < heights.size(); ++i) {
    const double a = heights[i];
    const double b = positions[i];
    const double c = widths[i];
    const double cinv = 1 / (c * c);
    const double d1 = phi - b;
    const double d2 = twoPI - (phi + b);
    const double d3 = twoPI + d1;
    const double d4 = twoPI + d2;
    const double t1 = std::exp(-d1 * d1 * cinv);
    const double t2 = std::exp(-d2 * d2 * cinv);
    const double t3 = std::exp(-d3 * d3 * cinv);
    const double t4 = std::exp(-d4 * d4 * cinv);
    result += a * (t1 + t2 + t3 + t4);
  }
  return -KT * std::log(result);
}

double getdEdPhi(const std::vector<double> &heights,
                 const std::vector<double> &positions,
                 const std::vector<double> &widths, double phi) {
  double innerDerivs = 0.0;
  double denominator = 0.0;
  for (size_t i = 0; i < heights.size(); ++i) {
    const double a = heights[i];
    const double b = positions[i];
    const double c = widths[i];
    const double cinv = 1 / (c * c);
    const double d1 = phi - b;
    const double d2 = twoPI - (phi + b);
    const double d3 = twoPI + d1;
    const double d4 = twoPI + d2;
    const double t1 = std::exp(-d1 * d1 * cinv);
    const double t2 = std::exp(-d2 * d2 * cinv);
    const double t3 = std::exp(-d3 * d3 * cinv);
    const double t4 = std::exp(-d4 * d4 * cinv);
    innerDerivs += a * cinv * (-t1 * d1 + t2 * d2 - t3 * d3 + t4 * d4);
    denominator += a * (t1 + t2 + t3 + t4);
  }
  return -2 * KT * innerDerivs / denominator;
}

GaussianTorsionAngleContribs::GaussianTorsionAngleContribs(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void GaussianTorsionAngleContribs::addContrib(
    std::size_t idx1, std::size_t idx2, std::size_t idx3, std::size_t idx4,
    std::vector<double> energies, std::vector<double> gradients,
    double scaling) {
  PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
                   (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
               "degenerate points");
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  URANGE_CHECK(idx4, dp_forceField->positions().size());
  d_contribs.emplace_back(idx1, idx2, idx3, idx4, std::move(energies),
                          std::move(gradients), scaling);
}

double GaussianTorsionAngleContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double accum = 0.0;
  const unsigned int dim = dp_forceField->dimension();
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D iPoint(pos[dim * contrib.idx1],
                                 pos[dim * contrib.idx1 + 1],
                                 pos[dim * contrib.idx1 + 2]);
    const RDGeom::Point3D jPoint(pos[dim * contrib.idx2],
                                 pos[dim * contrib.idx2 + 1],
                                 pos[dim * contrib.idx2 + 2]);
    const RDGeom::Point3D kPoint(pos[dim * contrib.idx3],
                                 pos[dim * contrib.idx3 + 1],
                                 pos[dim * contrib.idx3 + 2]);
    const RDGeom::Point3D lPoint(pos[dim * contrib.idx4],
                                 pos[dim * contrib.idx4 + 1],
                                 pos[dim * contrib.idx4 + 2]);
    const double cosphi =
        MMFF::Utils::calcTorsionCosPhi(iPoint, jPoint, kPoint, lPoint);
    const double Phi = std::acos(cosphi);
    accum += interpolate(contrib.energies, Phi) * contrib.scaling;
  }
  return accum;
}

void calcTorsionGrad(const RDGeom::Point3D *const r,
                     const RDGeom::Point3D *const t, const double *d,
                     double **g, const double &sinTerm, const double &cosPhi) {
  const double dCos_dT[6] = {1.0 / d[0] * (t[1].x - cosPhi * t[0].x),
                             1.0 / d[0] * (t[1].y - cosPhi * t[0].y),
                             1.0 / d[0] * (t[1].z - cosPhi * t[0].z),
                             1.0 / d[1] * (t[0].x - cosPhi * t[1].x),
                             1.0 / d[1] * (t[0].y - cosPhi * t[1].y),
                             1.0 / d[1] * (t[0].z - cosPhi * t[1].z)};

  g[0][0] += sinTerm * (dCos_dT[2] * r[1].y - dCos_dT[1] * r[1].z);
  g[0][1] += sinTerm * (dCos_dT[0] * r[1].z - dCos_dT[2] * r[1].x);
  g[0][2] += sinTerm * (dCos_dT[1] * r[1].x - dCos_dT[0] * r[1].y);

  g[1][0] += sinTerm *
             (dCos_dT[1] * (r[1].z - r[0].z) + dCos_dT[2] * (r[0].y - r[1].y) +
              dCos_dT[4] * (-r[3].z) + dCos_dT[5] * (r[3].y));
  g[1][1] += sinTerm *
             (dCos_dT[0] * (r[0].z - r[1].z) + dCos_dT[2] * (r[1].x - r[0].x) +
              dCos_dT[3] * (r[3].z) + dCos_dT[5] * (-r[3].x));
  g[1][2] += sinTerm *
             (dCos_dT[0] * (r[1].y - r[0].y) + dCos_dT[1] * (r[0].x - r[1].x) +
              dCos_dT[3] * (-r[3].y) + dCos_dT[4] * (r[3].x));

  g[2][0] += sinTerm *
             (dCos_dT[1] * (r[0].z) + dCos_dT[2] * (-r[0].y) +
              dCos_dT[4] * (r[3].z - r[2].z) + dCos_dT[5] * (r[2].y - r[3].y));
  g[2][1] += sinTerm *
             (dCos_dT[0] * (-r[0].z) + dCos_dT[2] * (r[0].x) +
              dCos_dT[3] * (r[2].z - r[3].z) + dCos_dT[5] * (r[3].x - r[2].x));
  g[2][2] += sinTerm *
             (dCos_dT[0] * (r[0].y) + dCos_dT[1] * (-r[0].x) +
              dCos_dT[3] * (r[3].y - r[2].y) + dCos_dT[4] * (r[2].x - r[3].x));

  g[3][0] += sinTerm * (dCos_dT[4] * r[2].z - dCos_dT[5] * r[2].y);
  g[3][1] += sinTerm * (dCos_dT[5] * r[2].x - dCos_dT[3] * r[2].z);
  g[3][2] += sinTerm * (dCos_dT[3] * r[2].y - dCos_dT[4] * r[2].x);
}

void GaussianTorsionAngleContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  const unsigned int dim = dp_forceField->dimension();
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D iPoint(pos[dim * contrib.idx1],
                                 pos[dim * contrib.idx1 + 1],
                                 pos[dim * contrib.idx1 + 2]);
    const RDGeom::Point3D jPoint(pos[dim * contrib.idx2],
                                 pos[dim * contrib.idx2 + 1],
                                 pos[dim * contrib.idx2 + 2]);
    const RDGeom::Point3D kPoint(pos[dim * contrib.idx3],
                                 pos[dim * contrib.idx3 + 1],
                                 pos[dim * contrib.idx3 + 2]);
    const RDGeom::Point3D lPoint(pos[dim * contrib.idx4],
                                 pos[dim * contrib.idx4 + 1],
                                 pos[dim * contrib.idx4 + 2]);
    double *g[4] = {&(grad[dim * contrib.idx1]), &(grad[dim * contrib.idx2]),
                    &(grad[dim * contrib.idx3]), &(grad[dim * contrib.idx4])};
    const RDGeom::Point3D r[4] = {iPoint - jPoint, kPoint - jPoint,
                                  jPoint - kPoint, lPoint - kPoint};
    RDGeom::Point3D t[2] = {r[0].crossProduct(r[1]), r[2].crossProduct(r[3])};
    double d[2] = {t[0].length(), t[1].length()};
    if (MMFF::isDoubleZero(d[0]) || MMFF::isDoubleZero(d[1])) {
      continue;
    }
    t[0] /= d[0];
    t[1] /= d[1];
    double cosPhi = t[0].dotProduct(t[1]);
    cosPhi = std::clamp(cosPhi, -1.0, 1.0);
    const double Phi = std::acos(cosPhi);
    const double sinPhi = std::sqrt(std::max(0.0, 1.0 - cosPhi * cosPhi));
    constexpr double EPSILON = 1e-8;
    const double safeSinPhi = std::max(sinPhi, EPSILON);
    const double dE_dPhi =
        interpolate(contrib.gradients, Phi) * contrib.scaling;
    const double sinTerm = -dE_dPhi / safeSinPhi;
    calcTorsionGrad(r, t, d, g, sinTerm, cosPhi);
  }
}
}  // namespace CrystalFF
}  // namespace ForceFields
