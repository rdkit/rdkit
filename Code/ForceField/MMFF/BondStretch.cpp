//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "BondStretch.h"
#include "Params.h"
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>

namespace ForceFields {
namespace MMFF {
namespace Utils {

double calcBondRestLength(const MMFFBond *mmffBondParams) {
  PRECONDITION(mmffBondParams, "bond parameters not found");

  return mmffBondParams->r0;
}

double calcBondForceConstant(const MMFFBond *mmffBondParams) {
  PRECONDITION(mmffBondParams, "bond parameters not found");

  return mmffBondParams->kb;
}

double calcBondStretchEnergy(const double r0, const double kb,
                             const double distance) {
  double distTerm = distance - r0;
  double distTerm2 = distTerm * distTerm;
  double const c1 = MDYNE_A_TO_KCAL_MOL;
  double const cs = -2.0;
  double const c3 = 7.0 / 12.0;

  return (0.5 * c1 * kb * distTerm2 *
          (1.0 + cs * distTerm + c3 * cs * cs * distTerm2));
}
}  // end of namespace Utils

BondStretchContrib::BondStretchContrib(ForceField *owner) {
  PRECONDITION(owner, "bad owner");


  dp_forceField = owner;

}

void BondStretchContrib::addTerm(const unsigned int idx1,
                                 const unsigned int idx2,
                                 const ForceFields::MMFF::MMFFBond *mmffBondParams) {
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  PRECONDITION(mmffBondParams, "bond parameters not found");

  d_at1Idxs.push_back(idx1);
  d_at2Idxs.push_back(idx2);
  d_r0.push_back(mmffBondParams->r0);
  d_kb.push_back(mmffBondParams->kb);
}

double BondStretchContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  const int numTerms = d_at1Idxs.size();
  double energySum = 0.0;
  for (int i =0; i < numTerms; i++) {
    energySum += Utils::calcBondStretchEnergy(
        d_r0[i], d_kb[i], dp_forceField->distance(d_at1Idxs[i], d_at2Idxs[i], pos));
  }
  return energySum;
}

void BondStretchContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  const int numTerms = d_at1Idxs.size();
  constexpr double cs = -2.0;
  constexpr double c1 = MDYNE_A_TO_KCAL_MOL;
  constexpr double c3 = 7.0 / 12.0;
  for (int termIdx = 0; termIdx < numTerms; termIdx++) {
    const int d_at1Idx = d_at1Idxs[termIdx];
    const int d_at2Idx = d_at2Idxs[termIdx];

    double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);

    double *at1Coords = &(pos[3 * d_at1Idx]);
    double *at2Coords = &(pos[3 * d_at2Idx]);
    double *g1 = &(grad[3 * d_at1Idx]);
    double *g2 = &(grad[3 * d_at2Idx]);

    double distTerm = dist - d_r0[termIdx];
    double dE_dr =
        c1 * d_kb[termIdx] * distTerm *
        (1.0 + 1.5 * cs * distTerm + 2.0 * c3 * cs * cs * distTerm * distTerm);
    double dGrad;
    for (unsigned int i = 0; i < 3; ++i) {
      dGrad = ((dist > 0.0) ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist)
                            : d_kb[termIdx] * 0.01);
      g1[i] += dGrad;
      g2[i] -= dGrad;
    }
  }
}
}  // namespace MMFF
}  // namespace ForceFields
