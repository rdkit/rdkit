//
//  Copyright (C) 2024 Niels Maeder and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Inversions.h"
#include "Utils.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
namespace UFF {

InversionContribs::InversionContribs(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void InversionContribs::addContrib(unsigned int idx1, unsigned int idx2,
                                   unsigned int idx3, unsigned int idx4,
                                   int at2AtomicNum, bool isCBoundToO,
                                   double oobForceScalingFactor) {
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  URANGE_CHECK(idx4, dp_forceField->positions().size());
  auto invCoeffForceCon = Utils::calcInversionCoefficientsAndForceConstant(
      at2AtomicNum, isCBoundToO);

  d_contribs.emplace_back(
      idx1, idx2, idx3, idx4, at2AtomicNum, isCBoundToO,
      std::get<1>(invCoeffForceCon), std::get<2>(invCoeffForceCon),
      std::get<3>(invCoeffForceCon),
      std::get<0>(invCoeffForceCon) * oobForceScalingFactor);
}

double InversionContribs::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double accum = 0;
  for (const auto &contrib : d_contribs) {
    const RDGeom::Point3D p1(pos[3 * contrib.idx1], pos[3 * contrib.idx1 + 1],
                             pos[3 * contrib.idx1 + 2]);
    const RDGeom::Point3D p2(pos[3 * contrib.idx2], pos[3 * contrib.idx2 + 1],
                             pos[3 * contrib.idx2 + 2]);
    const RDGeom::Point3D p3(pos[3 * contrib.idx3], pos[3 * contrib.idx3 + 1],
                             pos[3 * contrib.idx3 + 2]);
    const RDGeom::Point3D p4(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
                             pos[3 * contrib.idx4 + 2]);
    const double cosY = Utils::calculateCosY(p1, p2, p3, p4);
    const double sinYSq = 1.0 - cosY * cosY;
    const double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
    // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
    const double cos2W = 2.0 * sinY * sinY - 1.0;
    accum += contrib.forceConstant *
             (contrib.C0 + contrib.C1 * sinY + contrib.C2 * cos2W);
  }
  return accum;
}

void InversionContribs::getGrad(double *pos, double *grad) const {
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
    const RDGeom::Point3D p4(pos[3 * contrib.idx4], pos[3 * contrib.idx4 + 1],
                             pos[3 * contrib.idx4 + 2]);
    double *g1 = &(grad[3 * contrib.idx1]);
    double *g2 = &(grad[3 * contrib.idx2]);
    double *g3 = &(grad[3 * contrib.idx3]);
    double *g4 = &(grad[3 * contrib.idx4]);
    RDGeom::Point3D rJI = p1 - p2;
    RDGeom::Point3D rJK = p3 - p2;
    RDGeom::Point3D rJL = p4 - p2;
    const double dJI = rJI.length();
    const double dJK = rJK.length();
    const double dJL = rJL.length();
    if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
      return;
    }
    rJI.normalize();
    rJK.normalize();
    rJL.normalize();

    RDGeom::Point3D n = (-rJI).crossProduct(rJK);
    n.normalize();
    double cosY = n.dotProduct(rJL);
    cosY = std::clamp(cosY, -1.0, 1.0);
    const double sinYSq = 1.0 - cosY * cosY;
    const double sinY = std::max(sqrt(sinYSq), 1.0e-8);
    double cosTheta = rJI.dotProduct(rJK);
    cosTheta = std::clamp(cosTheta, -1.0, 1.0);
    const double sinThetaSq = 1.0 - cosTheta * cosTheta;
    const double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
    // sin(2 * W) = 2 * sin(W) * cos(W) = 2 * cos(Y) * sin(Y)
    const double dE_dW = -contrib.forceConstant *
                         (contrib.C1 * cosY - 4.0 * contrib.C2 * cosY * sinY);
    const RDGeom::Point3D t1 = rJL.crossProduct(rJK);
    const RDGeom::Point3D t2 = rJI.crossProduct(rJL);
    const RDGeom::Point3D t3 = rJK.crossProduct(rJI);
    const double term1 = sinY * sinTheta;
    const double term2 = cosY / (sinY * sinThetaSq);
    const double tg1[3] = {
        (t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
        (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
        (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
    const double tg3[3] = {
        (t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
        (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
        (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
    const double tg4[3] = {(t3.x / term1 - rJL.x * cosY / sinY) / dJL,
                           (t3.y / term1 - rJL.y * cosY / sinY) / dJL,
                           (t3.z / term1 - rJL.z * cosY / sinY) / dJL};
    for (unsigned int i = 0; i < 3; ++i) {
      g1[i] += dE_dW * tg1[i];
      g2[i] += -dE_dW * (tg1[i] + tg3[i] + tg4[i]);
      g3[i] += dE_dW * tg3[i];
      g4[i] += dE_dW * tg4[i];
    }
  }
}
}  // namespace UFF
}  // namespace ForceFields