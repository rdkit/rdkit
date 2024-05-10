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
#include "StretchBend.h"
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

std::pair<double, double> calcStbnForceConstants(
    const MMFFStbn *mmffStbnParams) {
  PRECONDITION(mmffStbnParams, "stretch-bend parameters not found");

  return std::make_pair(mmffStbnParams->kbaIJK, mmffStbnParams->kbaKJI);
}

std::pair<double, double> calcStretchBendEnergy(
    const double deltaDist1, const double deltaDist2, const double deltaTheta,
    const std::pair<double, double> forceConstants) {
  double factor = MDYNE_A_TO_KCAL_MOL * DEG2RAD * deltaTheta;

  return std::make_pair(factor * forceConstants.first * deltaDist1,
                        factor * forceConstants.second * deltaDist2);
}
}  // end of namespace Utils

StretchBendContrib::StretchBendContrib(
    ForceField *owner, const unsigned int idx1, const unsigned int idx2,
    const unsigned int idx3, const MMFFStbn *mmffStbnParams,
    const MMFFAngle *mmffAngleParams, const MMFFBond *mmffBondParams1,
    const MMFFBond *mmffBondParams2) {
  PRECONDITION(owner, "bad owner");
  PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),
               "degenerate points");
  URANGE_CHECK(idx1, owner->positions().size());
  URANGE_CHECK(idx2, owner->positions().size());
  URANGE_CHECK(idx3, owner->positions().size());

  dp_forceField = owner;
  d_at1Idx = idx1;
  d_at2Idx = idx2;
  d_at3Idx = idx3;
  d_restLen1 = Utils::calcBondRestLength(mmffBondParams1);
  d_restLen2 = Utils::calcBondRestLength(mmffBondParams2);
  d_theta0 = Utils::calcAngleRestValue(mmffAngleParams);
  d_forceConstants = Utils::calcStbnForceConstants(mmffStbnParams);
}

double StretchBendContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
  double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);

  std::pair<double, double> stretchBendEnergies = Utils::calcStretchBendEnergy(
      dist1 - d_restLen1, dist2 - d_restLen2,
      RAD2DEG * acos(Utils::calcCosTheta(p1, p2, p3, dist1, dist2)) - d_theta0,
      d_forceConstants);

  return (stretchBendEnergies.first + stretchBendEnergies.second);
}

void StretchBendContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");
  double dist1 = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
  double dist2 = dp_forceField->distance(d_at2Idx, d_at3Idx, pos);

  RDGeom::Point3D p1(pos[3 * d_at1Idx], pos[3 * d_at1Idx + 1],
                     pos[3 * d_at1Idx + 2]);
  RDGeom::Point3D p2(pos[3 * d_at2Idx], pos[3 * d_at2Idx + 1],
                     pos[3 * d_at2Idx + 2]);
  RDGeom::Point3D p3(pos[3 * d_at3Idx], pos[3 * d_at3Idx + 1],
                     pos[3 * d_at3Idx + 2]);
  double *g1 = &(grad[3 * d_at1Idx]);
  double *g2 = &(grad[3 * d_at2Idx]);
  double *g3 = &(grad[3 * d_at3Idx]);

  RDGeom::Point3D p12 = (p1 - p2) / dist1;
  RDGeom::Point3D p32 = (p3 - p2) / dist2;
  double const c5 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
  double cosTheta = p12.dotProduct(p32);
  clipToOne(cosTheta);
  double sinThetaSq = 1.0 - cosTheta * cosTheta;
  double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
  double angleTerm = RAD2DEG * acos(cosTheta) - d_theta0;
  double distTerm = RAD2DEG * (d_forceConstants.first * (dist1 - d_restLen1) +
                               d_forceConstants.second * (dist2 - d_restLen2));
  double dCos_dS1 = 1.0 / dist1 * (p32.x - cosTheta * p12.x);
  double dCos_dS2 = 1.0 / dist1 * (p32.y - cosTheta * p12.y);
  double dCos_dS3 = 1.0 / dist1 * (p32.z - cosTheta * p12.z);

  double dCos_dS4 = 1.0 / dist2 * (p12.x - cosTheta * p32.x);
  double dCos_dS5 = 1.0 / dist2 * (p12.y - cosTheta * p32.y);
  double dCos_dS6 = 1.0 / dist2 * (p12.z - cosTheta * p32.z);

  g1[0] += c5 * (p12.x * d_forceConstants.first * angleTerm +
                 dCos_dS1 / (-sinTheta) * distTerm);
  g1[1] += c5 * (p12.y * d_forceConstants.first * angleTerm +
                 dCos_dS2 / (-sinTheta) * distTerm);
  g1[2] += c5 * (p12.z * d_forceConstants.first * angleTerm +
                 dCos_dS3 / (-sinTheta) * distTerm);

  g2[0] +=
      c5 *
      ((-p12.x * d_forceConstants.first - p32.x * d_forceConstants.second) *
           angleTerm +
       (-dCos_dS1 - dCos_dS4) / (-sinTheta) * distTerm);
  g2[1] +=
      c5 *
      ((-p12.y * d_forceConstants.first - p32.y * d_forceConstants.second) *
           angleTerm +
       (-dCos_dS2 - dCos_dS5) / (-sinTheta) * distTerm);
  g2[2] +=
      c5 *
      ((-p12.z * d_forceConstants.first - p32.z * d_forceConstants.second) *
           angleTerm +
       (-dCos_dS3 - dCos_dS6) / (-sinTheta) * distTerm);

  g3[0] += c5 * (p32.x * d_forceConstants.second * angleTerm +
                 dCos_dS4 / (-sinTheta) * distTerm);
  g3[1] += c5 * (p32.y * d_forceConstants.second * angleTerm +
                 dCos_dS5 / (-sinTheta) * distTerm);
  g3[2] += c5 * (p32.z * d_forceConstants.second * angleTerm +
                 dCos_dS6 / (-sinTheta) * distTerm);
}
}  // namespace MMFF
}  // namespace ForceFields
