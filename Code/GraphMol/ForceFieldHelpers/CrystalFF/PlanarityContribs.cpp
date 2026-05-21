//
//  Copyright (C) 2026 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include "PlanarityContribs.h"
#include <numbers>
#include <cmath>

namespace ForceFields {
namespace CrystalFF {
namespace {
double calcOopChi(const RDGeom::Point3D &iPoint, const RDGeom::Point3D &jPoint,
                  const RDGeom::Point3D &kPoint,
                  const RDGeom::Point3D &lPoint) {
  RDGeom::Point3D rJI = iPoint - jPoint;
  RDGeom::Point3D rJK = kPoint - jPoint;
  RDGeom::Point3D rJL = lPoint - jPoint;
  rJI /= rJI.length();
  rJK /= rJK.length();
  rJL /= rJL.length();

  RDGeom::Point3D n = rJI.crossProduct(rJK);
  n /= n.length();
  double sinChi = n.dotProduct(rJL);
  sinChi = std::clamp(sinChi, -1.0, 1.0);

  return std::asin(sinChi) * 180.0 / std::numbers::pi;
}

double calcOopBendEnergy(const double chi, const double forceConst) {
  return (0.5 * forceConst * chi * chi);
}

bool isDoubleZero(const double x) { return ((x < 1.0e-10) && (x > -1.0e-10)); }

}  // namespace

PlanarityContribs::PlanarityContribs(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void PlanarityContribs::addContrib(std::size_t idx1, std::size_t idx2,
                                   std::size_t idx3, std::size_t idx4,
                                   double forceConstant) {
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  URANGE_CHECK(idx4, dp_forceField->positions().size());
  d_contribs.emplace_back(idx1, idx2, idx3, idx4, forceConstant);
}

double PlanarityContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double totalEnergy = 0.0;
  const unsigned int dimension = dp_forceField->dimension();
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D p1(pos[dimension * contrib.idx1],
                             pos[dimension * contrib.idx1 + 1],
                             pos[dimension * contrib.idx1 + 2]);
    const RDGeom::Point3D p2(pos[dimension * contrib.idx2],
                             pos[dimension * contrib.idx2 + 1],
                             pos[dimension * contrib.idx2 + 2]);
    const RDGeom::Point3D p3(pos[dimension * contrib.idx3],
                             pos[dimension * contrib.idx3 + 1],
                             pos[dimension * contrib.idx3 + 2]);
    const RDGeom::Point3D p4(pos[dimension * contrib.idx4],
                             pos[dimension * contrib.idx4 + 1],
                             pos[dimension * contrib.idx4 + 2]);

    totalEnergy +=
        calcOopBendEnergy(calcOopChi(p1, p2, p3, p4), contrib.forceConstant);
  }
  return totalEnergy;
}

void PlanarityContribs::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  const unsigned int dimension = dp_forceField->dimension();
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D iPoint(pos[dimension * contrib.idx1],
                                 pos[dimension * contrib.idx1 + 1],
                                 pos[dimension * contrib.idx1 + 2]);
    const RDGeom::Point3D jPoint(pos[dimension * contrib.idx2],
                                 pos[dimension * contrib.idx2 + 1],
                                 pos[dimension * contrib.idx2 + 2]);
    const RDGeom::Point3D kPoint(pos[dimension * contrib.idx3],
                                 pos[dimension * contrib.idx3 + 1],
                                 pos[dimension * contrib.idx3 + 2]);
    const RDGeom::Point3D lPoint(pos[dimension * contrib.idx4],
                                 pos[dimension * contrib.idx4 + 1],
                                 pos[dimension * contrib.idx4 + 2]);

    RDGeom::Point3D rJI = iPoint - jPoint;
    RDGeom::Point3D rJK = kPoint - jPoint;
    RDGeom::Point3D rJL = lPoint - jPoint;
    const double dJI = rJI.length();
    const double dJK = rJK.length();
    const double dJL = rJL.length();
    if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
      return;
    }
    rJI /= dJI;
    rJK /= dJK;
    rJL /= dJL;

    RDGeom::Point3D n = (-rJI).crossProduct(rJK);
    n /= n.length();
    double sinChi = rJL.dotProduct(n);
    sinChi = std::clamp(sinChi, -1.0, 1.0);
    double cosChiSq = 1.0 - sinChi * sinChi;
    double cosChi = std::max(((cosChiSq > 0.0) ? sqrt(cosChiSq) : 0.0), 1.0e-8);
    double chi = asin(sinChi) * 180 / std::numbers::pi;
    double cosTheta = rJI.dotProduct(rJK);
    cosTheta = std::clamp(cosTheta, -1.0, 1.0);
    double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
    const double sinTheta =
        std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);

    const double dE_dChi =
        contrib.forceConstant * chi * 180.0 / std::numbers::pi;
    const RDGeom::Point3D t1 = rJL.crossProduct(rJK);
    const RDGeom::Point3D t2 = rJI.crossProduct(rJL);
    const RDGeom::Point3D t3 = rJK.crossProduct(rJI);
    const double term1 = cosChi * sinTheta;
    const double term2 = sinChi / (cosChi * sinThetaSq);
    const double tg1[3] = {
        (t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
        (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
        (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
    const double tg3[3] = {
        (t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
        (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
        (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
    const double tg4[3] = {(t3.x / term1 - rJL.x * sinChi / cosChi) / dJL,
                           (t3.y / term1 - rJL.y * sinChi / cosChi) / dJL,
                           (t3.z / term1 - rJL.z * sinChi / cosChi) / dJL};

    double *g1 = &(grad[dimension * contrib.idx1]);
    double *g2 = &(grad[dimension * contrib.idx2]);
    double *g3 = &(grad[dimension * contrib.idx3]);
    double *g4 = &(grad[dimension * contrib.idx4]);
    for (std::size_t i = 0; i < 3; ++i) {
      g1[i] += dE_dChi * tg1[i];
      g2[i] += -dE_dChi * (tg1[i] + tg3[i] + tg4[i]);
      g3[i] += dE_dChi * tg3[i];
      g4[i] += dE_dChi * tg4[i];
    }
  }
}
}  // namespace CrystalFF
}  // namespace ForceFields
