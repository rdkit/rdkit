//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AngleBend.h"
#include "BondStretch.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
namespace MMFF {
namespace Utils {

double calcAngleRestValue(const MMFFAngle *mmffAngleParams) {
  PRECONDITION(mmffAngleParams, "angle parameters not found");

  return mmffAngleParams->theta0;
}

double calcCosTheta(RDGeom::Point3D p1, RDGeom::Point3D p2, RDGeom::Point3D p3,
                    double dist1, double dist2) {
  RDGeom::Point3D p12 = p1 - p2;
  RDGeom::Point3D p32 = p3 - p2;
  double cosTheta = p12.dotProduct(p32) / (dist1 * dist2);
  clipToOne(cosTheta);

  return cosTheta;
}

double calcAngleForceConstant(const MMFFAngle *mmffAngleParams) {
  PRECONDITION(mmffAngleParams, "angle parameters not found");

  return mmffAngleParams->ka;
}

double calcAngleBendEnergy(const double theta0, const double ka, bool isLinear,
                           const double cosTheta) {
  double angle = RAD2DEG * acos(cosTheta) - theta0;
  double const cb = -0.006981317;
  double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;
  double res = 0.0;

  if (isLinear) {
    res = MDYNE_A_TO_KCAL_MOL * ka * (1.0 + cosTheta);
  } else {
    res = 0.5 * c2 * ka * angle * angle * (1.0 + cb * angle);
  }

  return res;
}

void calcAngleBendGrad(RDGeom::Point3D *r, double *dist, double **g,
                       double &dE_dTheta, double &cosTheta, double &sinTheta) {
  // -------
  // dTheta/dx is trickier:
  double dCos_dS[6] = {1.0 / dist[0] * (r[1].x - cosTheta * r[0].x),
                       1.0 / dist[0] * (r[1].y - cosTheta * r[0].y),
                       1.0 / dist[0] * (r[1].z - cosTheta * r[0].z),
                       1.0 / dist[1] * (r[0].x - cosTheta * r[1].x),
                       1.0 / dist[1] * (r[0].y - cosTheta * r[1].y),
                       1.0 / dist[1] * (r[0].z - cosTheta * r[1].z)};

  g[0][0] += dE_dTheta * dCos_dS[0] / (-sinTheta);
  g[0][1] += dE_dTheta * dCos_dS[1] / (-sinTheta);
  g[0][2] += dE_dTheta * dCos_dS[2] / (-sinTheta);

  g[1][0] += dE_dTheta * (-dCos_dS[0] - dCos_dS[3]) / (-sinTheta);
  g[1][1] += dE_dTheta * (-dCos_dS[1] - dCos_dS[4]) / (-sinTheta);
  g[1][2] += dE_dTheta * (-dCos_dS[2] - dCos_dS[5]) / (-sinTheta);

  g[2][0] += dE_dTheta * dCos_dS[3] / (-sinTheta);
  g[2][1] += dE_dTheta * dCos_dS[4] / (-sinTheta);
  g[2][2] += dE_dTheta * dCos_dS[5] / (-sinTheta);
}
}  // end of namespace Utils

AngleBendContrib::AngleBendContrib(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void AngleBendContrib::addTerm(unsigned int idx1,
                               unsigned int idx2,
                               unsigned int idx3,
                               const ForceFields::MMFF::MMFFAngle *mmffAngleParams,
                               const ForceFields::MMFF::MMFFProp *mmffPropParamsCentralAtom) {
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());
  PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),
               "degenerate points");
  d_at1Idxs.push_back(idx1);
  d_at2Idxs.push_back(idx2);
  d_at3Idxs.push_back(idx3);
  d_isLinear.push_back(mmffPropParamsCentralAtom->linh > 0u);
  d_theta0.push_back(mmffAngleParams->theta0);
  d_ka.push_back(mmffAngleParams->ka);
}


double AngleBendContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double res = 0.0;
  const int numTerms = d_at1Idxs.size();

  for (int i = 0; i < numTerms; i++) {
    const int d_at1Idx = d_at1Idxs[i];
    const int d_at2Idx = d_at2Idxs[i];
    const int d_at3Idx = d_at3Idxs[i];

    double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
    double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);

    RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                       pos[3 * d_at1Idx + 2]);
    RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                       pos[3 * d_at2Idx + 2]);
    RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                       pos[3 * d_at3Idx + 2]);

    res += Utils::calcAngleBendEnergy(
        d_theta0[i], d_ka[i], d_isLinear[i],
        Utils::calcCosTheta(p1, p2, p3, dist1, dist2));
  }
  return res;
}

void AngleBendContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  const int numTerms = d_at1Idxs.size();
  for (int i =0; i < numTerms; i++) {
    const int d_at1Idx = d_at1Idxs[i];
    const int d_at2Idx = d_at2Idxs[i];
    const int d_at3Idx = d_at3Idxs[i];

    double dist[2] = {dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
                      dp_forceField->distance(d_at2Idx, d_at3Idx, pos)};

    RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                       pos[3 * d_at1Idx + 2]);
    RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                       pos[3 * d_at2Idx + 2]);
    RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                       pos[3 * d_at3Idx + 2]);
    double *g[3] = {&(grad[3 * d_at1Idx]), &(grad[3 * d_at2Idx]),
                    &(grad[3 * d_at3Idx])};
    RDGeom::Point3D r[2] = {(p1 - p2) / dist[0], (p3 - p2) / dist[1]};
    double cosTheta = r[0].dotProduct(r[1]);
    clipToOne(cosTheta);
    double sinThetaSq = 1.0 - cosTheta * cosTheta;
    double sinTheta =
        std::max(((sinThetaSq > 0.0) ? sqrt(sinThetaSq) : 0.0), 1.0e-8);

    // use the chain rule:
    // dE/dx = dE/dTheta * dTheta/dx

    // dE/dTheta is independent of cartesians:
    double angleTerm = RAD2DEG * acos(cosTheta) - d_theta0[i];
    double const cb = -0.006981317;
    double const c2 = MDYNE_A_TO_KCAL_MOL * DEG2RAD * DEG2RAD;

    double dE_dTheta = (d_isLinear[i] ? -MDYNE_A_TO_KCAL_MOL * d_ka[i] * sinTheta
                                   : RAD2DEG * c2 * d_ka[i] * angleTerm *
                                         (1.0 + 1.5 * cb * angleTerm));

    Utils::calcAngleBendGrad(r, dist, g, dE_dTheta, cosTheta, sinTheta);
  }
}

}  // namespace MMFF
}  // namespace ForceFields
