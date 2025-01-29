//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "OopBend.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>

namespace ForceFields {
namespace MMFF {
namespace Utils {
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
  clipToOne(sinChi);

  return RAD2DEG * asin(sinChi);
}

double calcOopBendForceConstant(const MMFFOop *mmffOopParams) {
  PRECONDITION(mmffOopParams, "no OOP parameters");

  return mmffOopParams->koop;
}

double calcOopBendEnergy(const double chi, const double koop) {
  double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
  return (0.5 * c2 * koop * chi * chi);
}
}  // end of namespace Utils

OopBendContrib::OopBendContrib(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void OopBendContrib::addTerm(unsigned int idx1,
                             unsigned int idx2,
                             unsigned int idx3,
                             unsigned int idx4,
                             const ForceFields::MMFF::MMFFOop *mmffOopParams) {
  PRECONDITION(mmffOopParams, "no OOP parameters");
  PRECONDITION((idx1 != idx2) && (idx1 != idx3) && (idx1 != idx4) &&
                   (idx2 != idx3) && (idx2 != idx4) && (idx3 != idx4),
               "degenerate points");
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  URANGE_CHECK(idx4, dp_forceField->positions().size());

  d_at1Idxs.push_back(idx1);
  d_at2Idxs.push_back(idx2);
  d_at3Idxs.push_back(idx3);
  d_at4Idxs.push_back(idx4);
  d_koop.push_back(mmffOopParams->koop);
}

double OopBendContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  const int numTerms = d_at1Idxs.size();
  double totalEnergy = 0.0;
  for (int i = 0; i < numTerms; ++i) {
    const int d_at1Idx = d_at1Idxs[i];
    const int d_at2Idx = d_at2Idxs[i];
    const int d_at3Idx = d_at3Idxs[i];
    const int d_at4Idx = d_at4Idxs[i];

    RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                       pos[3 * d_at1Idx + 2]);
    RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                       pos[3 * d_at2Idx + 2]);
    RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                       pos[3 * d_at3Idx + 2]);
    RDGeom::Point3D p4(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
                       pos[3 * d_at4Idx + 2]);

    totalEnergy += Utils::calcOopBendEnergy(Utils::calcOopChi(p1, p2, p3, p4), d_koop[i]);
  }
  return totalEnergy;
}

void OopBendContrib::getGrad(double* pos, double* grad) const {
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  PRECONDITION(dp_forceField, "no owner");

  const int numTerms = d_at1Idxs.size();
  for (int i =0; i < numTerms; i++) {
    getSingleGrad(pos, grad, i);
  }
}

void OopBendContrib::getSingleGrad(double *pos, double *grad, unsigned int termIdx) const {

  const int d_at1Idx = d_at1Idxs[termIdx];
  const int d_at2Idx = d_at2Idxs[termIdx];
  const int d_at3Idx = d_at3Idxs[termIdx];
  const int d_at4Idx = d_at4Idxs[termIdx];

  RDGeom::Point3D iPoint(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                         pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D jPoint(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                         pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D kPoint(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                         pos[3 * d_at3Idx + 2]);
  RDGeom::Point3D lPoint(pos[3 * d_at4Idx], pos[3 * d_at4Idx + 1],
                         pos[3 * d_at4Idx + 2]);
  double *g1 = &(grad[3 * d_at1Idx]);
  double *g2 = &(grad[3 * d_at2Idx]);
  double *g3 = &(grad[3 * d_at3Idx]);
  double *g4 = &(grad[3 * d_at4Idx]);

  RDGeom::Point3D rJI = iPoint - jPoint;
  RDGeom::Point3D rJK = kPoint - jPoint;
  RDGeom::Point3D rJL = lPoint - jPoint;
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
  double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
  double sinChi = rJL.dotProduct(n);
  clipToOne(sinChi);
  double cosChiSq = 1.0 - sinChi * sinChi;
  double cosChi = std::max(((cosChiSq > 0.0) ? sqrt(cosChiSq) : 0.0), 1.0e-8);
  double chi = RAD2DEG * asin(sinChi);
  double cosTheta = rJI.dotProduct(rJK);
  clipToOne(cosTheta);
  double sinThetaSq = std::max(1.0 - cosTheta * cosTheta, 1.0e-8);
  double sinTheta =
      std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);

  double dE_dChi = RAD2DEG * c2 * d_koop[termIdx] * chi;
  RDGeom::Point3D t1 = rJL.crossProduct(rJK);
  RDGeom::Point3D t2 = rJI.crossProduct(rJL);
  RDGeom::Point3D t3 = rJK.crossProduct(rJI);
  double term1 = cosChi * sinTheta;
  double term2 = sinChi / (cosChi * sinThetaSq);
  double tg1[3] = {(t1.x / term1 - (rJI.x - rJK.x * cosTheta) * term2) / dJI,
                   (t1.y / term1 - (rJI.y - rJK.y * cosTheta) * term2) / dJI,
                   (t1.z / term1 - (rJI.z - rJK.z * cosTheta) * term2) / dJI};
  double tg3[3] = {(t2.x / term1 - (rJK.x - rJI.x * cosTheta) * term2) / dJK,
                   (t2.y / term1 - (rJK.y - rJI.y * cosTheta) * term2) / dJK,
                   (t2.z / term1 - (rJK.z - rJI.z * cosTheta) * term2) / dJK};
  double tg4[3] = {(t3.x / term1 - rJL.x * sinChi / cosChi) / dJL,
                   (t3.y / term1 - rJL.y * sinChi / cosChi) / dJL,
                   (t3.z / term1 - rJL.z * sinChi / cosChi) / dJL};
  for (unsigned int i = 0; i < 3; ++i) {
    g1[i] += dE_dChi * tg1[i];
    g2[i] += -dE_dChi * (tg1[i] + tg3[i] + tg4[i]);
    g3[i] += dE_dChi * tg3[i];
    g4[i] += dE_dChi * tg4[i];
  }
}
}  // namespace MMFF
}  // namespace ForceFields
