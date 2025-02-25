//  Copyright (C) 2013-2025 Paolo Tosco and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Nonbonded.h"
#include "Params.h"
#include <cmath>
#include <ForceField/ForceField.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/utils.h>
#include <GraphMol/ForceFieldHelpers/MMFF/AtomTyper.h>

namespace ForceFields {
namespace MMFF {
namespace Utils {
double calcUnscaledVdWMinimum(const MMFFVdWCollection *mmffVdW,
                              const MMFFVdW *mmffVdWParamsIAtom,
                              const MMFFVdW *mmffVdWParamsJAtom) {
  double gamma_ij = (mmffVdWParamsIAtom->R_star - mmffVdWParamsJAtom->R_star) /
                    (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star);

  return (0.5 * (mmffVdWParamsIAtom->R_star + mmffVdWParamsJAtom->R_star) *
          (1.0 +
           (((mmffVdWParamsIAtom->DA == 'D') || (mmffVdWParamsJAtom->DA == 'D'))
                ? 0.0
                : mmffVdW->B *
                      (1.0 - exp(-(mmffVdW->Beta) * gamma_ij * gamma_ij)))));
}

double calcUnscaledVdWWellDepth(double R_star_ij,
                                const MMFFVdW *mmffVdWParamsIAtom,
                                const MMFFVdW *mmffVdWParamsJAtom) {
  double R_star_ij2 = R_star_ij * R_star_ij;
  double const c4 = 181.16;

  return (c4 * mmffVdWParamsIAtom->G_i * mmffVdWParamsJAtom->G_i *
          mmffVdWParamsIAtom->alpha_i * mmffVdWParamsJAtom->alpha_i /
          ((sqrt(mmffVdWParamsIAtom->alpha_i / mmffVdWParamsIAtom->N_i) +
            sqrt(mmffVdWParamsJAtom->alpha_i / mmffVdWParamsJAtom->N_i)) *
           R_star_ij2 * R_star_ij2 * R_star_ij2));
}

double calcVdWEnergy(const double dist, const double R_star_ij,
                     const double wellDepth) {
  double const vdw1 = 1.07;
  double const vdw1m1 = vdw1 - 1.0;
  double const vdw2 = 1.12;
  double const vdw2m1 = vdw2 - 1.0;
  double dist2 = dist * dist;
  double dist7 = dist2 * dist2 * dist2 * dist;
  double aTerm = vdw1 * R_star_ij / (dist + vdw1m1 * R_star_ij);
  double aTerm2 = aTerm * aTerm;
  double aTerm7 = aTerm2 * aTerm2 * aTerm2 * aTerm;
  double R_star_ij2 = R_star_ij * R_star_ij;
  double R_star_ij7 = R_star_ij2 * R_star_ij2 * R_star_ij2 * R_star_ij;
  double bTerm = vdw2 * R_star_ij7 / (dist7 + vdw2m1 * R_star_ij7) - 2.0;
  double res = wellDepth * aTerm7 * bTerm;

  return res;
}

void scaleVdWParams(double &R_star_ij, double &wellDepth,
                    const MMFFVdWCollection *mmffVdW,
                    const MMFFVdW *mmffVdWParamsIAtom,
                    const MMFFVdW *mmffVdWParamsJAtom) {
  if (((mmffVdWParamsIAtom->DA == 'D') && (mmffVdWParamsJAtom->DA == 'A')) ||
      ((mmffVdWParamsIAtom->DA == 'A') && (mmffVdWParamsJAtom->DA == 'D'))) {
    R_star_ij *= mmffVdW->DARAD;
    wellDepth *= mmffVdW->DAEPS;
  }
}

double calcEleEnergy(unsigned int, unsigned int, double dist, double chargeTerm,
                     std::uint8_t dielModel, bool is1_4) {
  double corr_dist = dist + 0.05;
  double const diel = 332.0716;
  double const sc1_4 = 0.75;
  if (dielModel == RDKit::MMFF::DISTANCE) {
    corr_dist *= corr_dist;
  }
  return (diel * chargeTerm / corr_dist * (is1_4 ? sc1_4 : 1.0));
}
}  // namespace Utils

VdWContrib::VdWContrib(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}

void VdWContrib::addTerm(unsigned int idx1, unsigned int idx2, const MMFFVdWRijstarEps *mmffVdWConstants) {
  PRECONDITION(mmffVdWConstants, "bad MMFFVdW parameters");
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  d_at1Idxs.push_back(idx1);
  d_at2Idxs.push_back(idx2);
  d_R_ij_stars.push_back(mmffVdWConstants->R_ij_star);
  d_wellDepths.push_back(mmffVdWConstants->epsilon);
}

double VdWContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  double energySum = 0.0;

  const int numPairs = d_at1Idxs.size();
  for (int i = 0; i < numPairs; ++i) {
    unsigned int d_at1Idx = d_at1Idxs[i];
    unsigned int d_at2Idx = d_at2Idxs[i];
    double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
    double res = Utils::calcVdWEnergy(dist, d_R_ij_stars[i], d_wellDepths[i]);
    energySum += res;
  }
  return energySum;
}

void VdWContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  double const vdw1 = 1.07;
  double const vdw1m1 = vdw1 - 1.0;
  double const vdw2 = 1.12;
  double const vdw2m1 = vdw2 - 1.0;
  double const vdw2t7 = vdw2 * 7.0;


  const int numPairs = d_at1Idxs.size();
  for (int pairIdx = 0; pairIdx < numPairs; ++pairIdx) {
    const int d_at1Idx = d_at1Idxs[pairIdx];
    const int d_at2Idx = d_at2Idxs[pairIdx];
    const double d_R_ij_star = d_R_ij_stars[pairIdx];
    const double d_wellDepth = d_wellDepths[pairIdx];

    double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
    double *at1Coords = &(pos[3 * d_at1Idx]);
    double *at2Coords = &(pos[3 * d_at2Idx]);
    double *g1 = &(grad[3 * d_at1Idx]);
    double *g2 = &(grad[3 * d_at2Idx]);
    double q = dist / d_R_ij_star;
    double q2 = q * q;
    double q6 = q2 * q2 * q2;
    double q7 = q6 * q;
    double q7pvdw2m1 = q7 + vdw2m1;
    double t = vdw1 / (q + vdw1 - 1.0);
    double t2 = t * t;
    double t7 = t2 * t2 * t2 * t;
    double dE_dr = d_wellDepth / d_R_ij_star * t7 *
                   (-vdw2t7 * q6 / (q7pvdw2m1 * q7pvdw2m1) +
                    ((-vdw2t7 / q7pvdw2m1 + 14.0) / (q + vdw1m1)));
    for (unsigned int i = 0; i < 3; ++i) {
      double dGrad;
      dGrad = ((dist > 0.0) ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist)
                            : d_R_ij_star * 0.01);
      g1[i] += dGrad;
      g2[i] -= dGrad;
    }
  }
}

EleContrib::EleContrib(ForceField *owner) {
  PRECONDITION(owner, "bad owner");
  dp_forceField = owner;
}
void EleContrib::addTerm(unsigned int idx1, unsigned int idx2,
             double chargeTerm, std::uint8_t dielModel, bool is1_4) {
  URANGE_CHECK(idx1, dp_forceField->positions().size());
  URANGE_CHECK(idx2, dp_forceField->positions().size());
  d_at1Idxs.push_back(idx1);
  d_at2Idxs.push_back(idx2);
  d_chargeTerms.push_back(chargeTerm);
  d_dielModels.push_back(dielModel);
  d_is_1_4s.push_back(is1_4);
}
double EleContrib::getEnergy(double *pos) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");

  double res = 0.0;
  const int numPairs = d_at1Idxs.size();

  for (int i = 0; i < numPairs; ++i) {
    unsigned int d_at1Idx = d_at1Idxs[i];
    unsigned int d_at2Idx = d_at2Idxs[i];
    double d_chargeTerm = d_chargeTerms[i];
    std::uint8_t d_dielModel = d_dielModels[i];
    bool d_is1_4 = d_is_1_4s[i];

    res += Utils::calcEleEnergy(d_at1Idx,
                                d_at2Idx,
                                dp_forceField->distance(d_at1Idx, d_at2Idx, pos),
                                d_chargeTerm,
                                d_dielModel,
                                d_is1_4);
  }
  return res;
}

void EleContrib::getGrad(double *pos, double *grad) const {
  PRECONDITION(dp_forceField, "no owner");
  PRECONDITION(pos, "bad vector");
  PRECONDITION(grad, "bad vector");

  const int numPairs = d_at1Idxs.size();
  for (int pairIdx = 0; pairIdx < numPairs; ++pairIdx) {
    const int d_at1Idx = d_at1Idxs[pairIdx];
    const int d_at2Idx = d_at2Idxs[pairIdx];
    const double d_chargeTerm = d_chargeTerms[pairIdx];
    const std::uint8_t d_dielModel = d_dielModels[pairIdx];
    const bool d_is1_4 = d_is_1_4s[pairIdx];

    double dist = dp_forceField->distance(d_at1Idx, d_at2Idx, pos);
    double *at1Coords = &(pos[3 * d_at1Idx]);
    double *at2Coords = &(pos[3 * d_at2Idx]);
    double *g1 = &(grad[3 * d_at1Idx]);
    double *g2 = &(grad[3 * d_at2Idx]);
    double corr_dist = dist + 0.05;
    corr_dist *= ((d_dielModel == RDKit::MMFF::DISTANCE) ? corr_dist * corr_dist
                                                         : corr_dist);
    double dE_dr = -332.0716 * (double)(d_dielModel)*d_chargeTerm / corr_dist *
                   (d_is1_4 ? 0.75 : 1.0);
    for (unsigned int i = 0; i < 3; ++i) {
      double dGrad;
      dGrad = ((dist > 0.0) ? (dE_dr * (at1Coords[i] - at2Coords[i]) / dist)
                            : 0.02);
      g1[i] += dGrad;
      g2[i] -= dGrad;
    }
  }
}
}  // namespace MMFF
}  // namespace ForceFields
