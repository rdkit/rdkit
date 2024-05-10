//
//  Copyright (C) 2013 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "Inversion.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
namespace UFF {
namespace Utils {
double calculateCosY(const RDGeom::Point3D &iPoint,
                     const RDGeom::Point3D &jPoint,
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

  return n.dotProduct(rJL);
}

std::tuple<double, double, double, double>
calcInversionCoefficientsAndForceConstant(int at2AtomicNum, bool isCBoundToO) {
  double res = 0.0;
  double C0 = 0.0;
  double C1 = 0.0;
  double C2 = 0.0;
  // if the central atom is sp2 carbon, nitrogen or oxygen
  if ((at2AtomicNum == 6) || (at2AtomicNum == 7) || (at2AtomicNum == 8)) {
    C0 = 1.0;
    C1 = -1.0;
    C2 = 0.0;
    res = (isCBoundToO ? 50.0 : 6.0);
  } else {
    // group 5 elements are not clearly explained in the UFF paper
    // the following code was inspired by MCCCS Towhee's ffuff.F
    double w0 = M_PI / 180.0;
    switch (at2AtomicNum) {
      // if the central atom is phosphorous
      case 15:
        w0 *= 84.4339;
        break;

      // if the central atom is arsenic
      case 33:
        w0 *= 86.9735;
        break;

      // if the central atom is antimonium
      case 51:
        w0 *= 87.7047;
        break;

      // if the central atom is bismuth
      case 83:
        w0 *= 90.0;
        break;
    }
    C2 = 1.0;
    C1 = -4.0 * cos(w0);
    C0 = -(C1 * cos(w0) + C2 * cos(2.0 * w0));
    res = 22.0 / (C0 + C1 + C2);
  }
  res /= 3.0;

  return std::make_tuple(res, C0, C1, C2);
}
}  // end of namespace Utils

InversionContrib::InversionContrib(ForceField *owner, unsigned int idx1,
                                   unsigned int idx2, unsigned int idx3,
                                   unsigned int idx4, int at2AtomicNum,
                                   bool isCBoundToO,
                                   double oobForceScalingFactor) {
  PRECONDITION(owner, "bad owner");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());
  URANGE_CHECK(idx4, owner->positions().size());

  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  d_at4Idx = idx4;

  auto invCoeffForceCon = Utils::calcInversionCoefficientsAndForceConstant(
      at2AtomicNum, isCBoundToO);
  d_forceConstant = oobForceScalingFactor * std::get<0>(invCoeffForceCon);
  d_C0 = std::get<1>(invCoeffForceCon);
  d_C1 = std::get<2>(invCoeffForceCon);
  d_C2 = std::get<3>(invCoeffForceCon);
}

double InversionContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);
  RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
                     pos[3 * d_at4Idx + 2]);

  double cosY = Utils::calculateCosY(p1, p2, p3, p4);
  double sinYSq = 1.0 - cosY * cosY;
  double sinY = ((sinYSq > 0.0) ? sqrt(sinYSq) : 0.0);
  // cos(2 * W) = 2 * cos(W) * cos(W) - 1 = 2 * sin(W) * sin(W) - 1
  double cos2W = 2.0 * sinY * sinY - 1.0;
  double res = d_forceConstant * (d_C0 + d_C1 * sinY + d_C2 * cos2W);
  // std::cout << d_at1Idx + 1 << "," << d_at2Idx + 1 << "," << d_at3Idx + 1 <<
  // "," << d_at4Idx + 1 << " Inversion: " << res << std::endl;

  return res;
}
void InversionContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);
  RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
                     pos[3 * d_at4Idx + 2]);
  double *g1 = &(grad[3 * d_at1Idx]);
  double *g2 = &(grad[3 * d_at2Idx]);
  double *g3 = &(grad[3 * d_at3Idx]);
  double *g4 = &(grad[3 * d_at4Idx]);

  RDGeom::Point3D rJI = p1 - p2;
  RDGeom::Point3D rJK = p3 - p2;
  RDGeom::Point3D rJL = p4 - p2;
  double dJI = rJI.length();
  double dJK = rJK.length();
  double dJL = rJL.length();
  if (isDoubleZero(dJI) || isDoubleZero(dJK) || isDoubleZero(dJL)) {
    return;
  }
  rJI /= dJI;
  rJK /= dJK;
  rJL /= dJL;

  RDGeom::Point3D n = (-rJI).crossProduct(rJK);
  n /= n.length();
  double cosY = n.dotProduct(rJL);
  clipToOne(cosY);
  double sinYSq = 1.0 - cosY * cosY;
  double sinY = std::max(sqrt(sinYSq), 1.0e-8);
  double cosTheta = rJI.dotProduct(rJK);
  clipToOne(cosTheta);
  double sinThetaSq = 1.0 - cosTheta * cosTheta;
  double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
  // sin(2 * W) = 2 * sin(W) * cos(W) = 2 * cos(Y) * sin(Y)
  double dE_dW = -d_forceConstant * (d_C1 * cosY - 4.0 * d_C2 * cosY * sinY);
  RDGeom::Point3D t1 = rJL.crossProduct(rJK);
  RDGeom::Point3D t2 = rJI.crossProduct(rJL);
  RDGeom::Point3D t3 = rJK.crossProduct(rJI);
  double term1 = sinY * sinTheta;
  double term2 = cosY / (sinY * sinThetaSq);
  double tg1[3] = {(t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
                   (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
                   (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
  double tg3[3] = {(t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
                   (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
                   (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
  double tg4[3] = {(t3.x / term1 - rJL.x * cosY / sinY) / dJL,
                   (t3.y / term1 - rJL.y * cosY / sinY) / dJL,
                   (t3.z / term1 - rJL.z * cosY / sinY) / dJL};
  for (unsigned int i = 0; i < 3; ++i) {
    g1[i] += dE_dW * tg1[i];
    g2[i] += -dE_dW * (tg1[i] + tg3[i] + tg4[i]);
    g3[i] += dE_dW * tg3[i];
    g4[i] += dE_dW * tg4[i];
  }
}
}  // namespace UFF
}  // namespace ForceFields
