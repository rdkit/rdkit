//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
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

StretchBendContrib::StretchBendContrib(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void StretchBendContrib::addTerm(
    const unsigned int idx1, const unsigned int idx2, const unsigned int idx3,
    const ForceFields::MMFF::MMFFStbn *mmffStbnParams,
    const ForceFields::MMFF::MMFFAngle *mmffAngleParams,
    const ForceFields::MMFF::MMFFBond *mmffBondParams1,
    const ForceFields::MMFF::MMFFBond *mmffBondParams2) {
  PRECONDITION(((idx1 != idx2) && (idx2 != idx3) && (idx1 != idx3)),
               "degenerate points");
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  URANGE_CHECK(idx3, dp_forceField->positions().size());

  d_at1Idxs.push_back(idx1);
  d_at2Idxs.push_back(idx2);
  d_at3Idxs.push_back(idx3);
  d_restLen1s.push_back(Utils::calcBondRestLength(mmffBondParams1));
  d_restLen2s.push_back(Utils::calcBondRestLength(mmffBondParams2));
  d_theta0s.push_back(Utils::calcAngleRestValue(mmffAngleParams));
  std::pair<double, double> forceConstants =
      Utils::calcStbnForceConstants(mmffStbnParams);
  d_forceConstants1.push_back(forceConstants.first);
  d_forceConstants2.push_back(forceConstants.second);
}

double StretchBendContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double totalEnergy = 0.0;
  const int numTerms = d_at1Idxs.size();
  for (int i = 0; i < numTerms; i++) {
    const int16_t at1Idx = d_at1Idxs[i];
    const int16_t at2Idx = d_at2Idxs[i];
    const int16_t at3Idx = d_at3Idxs[i];

    double dist1 = dp_forceField->distance(at1Idx, at2Idx, pos);
    double dist2 = dp_forceField->distance(at2Idx, at3Idx, pos);

    RDGeom::Point3D p1(pos[3 * at1Idx], pos[3 * at1Idx + 1],
                       pos[3 * at1Idx + 2]);
    RDGeom::Point3D p2(pos[3 * at2Idx], pos[3 * at2Idx + 1],
                       pos[3 * at2Idx + 2]);
    RDGeom::Point3D p3(pos[3 * at3Idx], pos[3 * at3Idx + 1],
                       pos[3 * at3Idx + 2]);
    std::pair<double, double> forceConstantsPair =
        std::make_pair(d_forceConstants1[i], d_forceConstants2[i]);
    std::pair<double, double> stretchBendEnergies =
        Utils::calcStretchBendEnergy(
            dist1 - d_restLen1s[i], dist2 - d_restLen2s[i],
            RAD2DEG * acos(Utils::calcCosTheta(p1, p2, p3, dist1, dist2)) -
                d_theta0s[i],
            forceConstantsPair);
    totalEnergy += (stretchBendEnergies.first + stretchBendEnergies.second);
  }
  return totalEnergy;
}

void StretchBendContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  const int numTerms = d_at1Idxs.size();
  for (int i = 0; i < numTerms; ++i) {
    const int16_t at1Idx = d_at1Idxs[i];
    const int16_t at2Idx = d_at2Idxs[i];
    const int16_t at3Idx = d_at3Idxs[i];
    const double theta0 = d_theta0s[i];
    const double forceConstant1 = d_forceConstants1[i];
    const double forceConstant2 = d_forceConstants2[i];
    const double restLen1 = d_restLen1s[i];
    const double restLen2 = d_restLen2s[i];

    double dist1 = dp_forceField->distance(at1Idx, at2Idx, pos);
    double dist2 = dp_forceField->distance(at2Idx, at3Idx, pos);

    RDGeom::Point3D p1(pos[3 * at1Idx], pos[3 * at1Idx + 1],
                       pos[3 * at1Idx + 2]);
    RDGeom::Point3D p2(pos[3 * at2Idx], pos[3 * at2Idx + 1],
                       pos[3 * at2Idx + 2]);
    RDGeom::Point3D p3(pos[3 * at3Idx], pos[3 * at3Idx + 1],
                       pos[3 * at3Idx + 2]);
    double *g1 = &(grad[3 * at1Idx]);
    double *g2 = &(grad[3 * at2Idx]);
    double *g3 = &(grad[3 * at3Idx]);

    RDGeom::Point3D p12 = (p1 - p2) / dist1;
    RDGeom::Point3D p32 = (p3 - p2) / dist2;
    double const c5 = MDYNE_A_TO_KCAL_MOL * DEG2RAD;
    double cosTheta = p12.dotProduct(p32);
    clipToOne(cosTheta);
    double sinThetaSq = 1.0 - cosTheta * cosTheta;
    double sinTheta = std::max(sqrt(sinThetaSq), 1.0e-8);
    double angleTerm = RAD2DEG * acos(cosTheta) - theta0;
    double distTerm = RAD2DEG * (forceConstant1 * (dist1 - restLen1) +
                                 forceConstant2 * (dist2 - restLen2));
    double dCos_dS1 = 1.0 / dist1 * (p32.x - cosTheta * p12.x);
    double dCos_dS2 = 1.0 / dist1 * (p32.y - cosTheta * p12.y);
    double dCos_dS3 = 1.0 / dist1 * (p32.z - cosTheta * p12.z);

    double dCos_dS4 = 1.0 / dist2 * (p12.x - cosTheta * p32.x);
    double dCos_dS5 = 1.0 / dist2 * (p12.y - cosTheta * p32.y);
    double dCos_dS6 = 1.0 / dist2 * (p12.z - cosTheta * p32.z);

    g1[0] += c5 * (p12.x * forceConstant1 * angleTerm +
                   dCos_dS1 / (-sinTheta) * distTerm);
    g1[1] += c5 * (p12.y * forceConstant1 * angleTerm +
                   dCos_dS2 / (-sinTheta) * distTerm);
    g1[2] += c5 * (p12.z * forceConstant1 * angleTerm +
                   dCos_dS3 / (-sinTheta) * distTerm);

    g2[0] +=
        c5 * ((-p12.x * forceConstant1 - p32.x * forceConstant2) * angleTerm +
              (-dCos_dS1 - dCos_dS4) / (-sinTheta) * distTerm);
    g2[1] +=
        c5 * ((-p12.y * forceConstant1 - p32.y * forceConstant2) * angleTerm +
              (-dCos_dS2 - dCos_dS5) / (-sinTheta) * distTerm);
    g2[2] +=
        c5 * ((-p12.z * forceConstant1 - p32.z * forceConstant2) * angleTerm +
              (-dCos_dS3 - dCos_dS6) / (-sinTheta) * distTerm);

    g3[0] += c5 * (p32.x * forceConstant2 * angleTerm +
                   dCos_dS4 / (-sinTheta) * distTerm);
    g3[1] += c5 * (p32.y * forceConstant2 * angleTerm +
                   dCos_dS5 / (-sinTheta) * distTerm);
    g3[2] += c5 * (p32.z * forceConstant2 * angleTerm +
                   dCos_dS6 / (-sinTheta) * distTerm);
  }
}
}  // namespace MMFF
}  // namespace ForceFields
