//
//  Copyright (C) 2026 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>

constexpr int D = 4;
static const DTYPE PI = 4.0 * atan(1.0);

namespace RDKit {
namespace GaussianShape {
void normalizeQuaternion(DTYPE *q) {
  double magq = sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] * q[3] * q[3]);
  q[0] /= magq;
  q[1] /= magq;
  q[2] /= magq;
  q[3] /= magq;
}

// Extract the quaternion from the optimizer's pos,  normalizing it
// in situ in pos first.
std::array<DTYPE, 4> extractNormalizedQuaternion(double *pos) {
  double magq = sqrt(pos[0] * pos[0] + pos[1] * pos[1] +
                     pos[2] * pos[2] * pos[3] * pos[3]);
  pos[0] /= magq;
  pos[1] /= magq;
  pos[2] /= magq;
  pos[3] /= magq;
  std::array<DTYPE, 4> q;
  q[0] = pos[0];
  q[1] = pos[1];
  q[2] = pos[2];
  q[3] = pos[3];

  return q;
}

std::array<DTYPE, 3> extractTranslation(double *pos, bool reversed) {
  std::array<DTYPE, 3> t;
  if (reversed) {
    t[0] = -pos[4];
    t[1] = -pos[5];
    t[2] = -pos[6];
  } else {
    t[0] = pos[4];
    t[1] = pos[5];
    t[2] = pos[6];
  }
  return t;
}

SingleConformerAlignment::SingleConformerAlignment(
    const DTYPE *molA, DTYPE *molAT, const int *molA_type, int NmolA_shape,
    int NmolA_color, DTYPE shapeVolA, DTYPE colorVolA, const DTYPE *molB,
    DTYPE *molBT, const int *molB_type, int NmolB_shape, int NmolB_color,
    DTYPE shapeVolB, DTYPE colorVolB, OptimMode optimMode, DTYPE mixing_param,
    bool all_carbon_radii, unsigned int maxIts)
    : d_molA(molA),
      d_molAT(molAT),
      d_NmolA_shape(NmolA_shape),
      d_NmolA_color(NmolA_color),
      d_molA_type(molA_type),
      d_shapeVolA(shapeVolA),
      d_colorVolA(colorVolA),
      d_molB(molB),
      d_molBT(molBT),
      d_NmolB_shape(NmolB_shape),
      d_NmolB_color(NmolB_color),
      d_molB_type(molB_type),
      d_shapeVolB(shapeVolB),
      d_colorVolB(colorVolB),
      d_optimMode(optimMode),
      d_mixing_param(mixing_param),
      d_all_carbon_radii(all_carbon_radii),
      d_maxIts(maxIts) {}

std::array<DTYPE, 5> SingleConformerAlignment::calcScores(
    const DTYPE *molA, const DTYPE *molB, bool includeColor) const {
  std::array<DTYPE, 5> scores{0.0, 0.0, 0.0, 0.0, 0.0};
  scores[3] = calcVolAndGrads(molA, d_NmolA_shape, molB, d_NmolB_shape,
                              d_all_carbon_radii, nullptr, nullptr);
  if (d_NmolA_color && d_NmolB_color &&
      (d_optimMode == OptimMode::SHAPE_PLUS_COLOR || includeColor)) {
    scores[4] = calcVolAndGrads(molA + d_NmolA_shape * D, d_NmolA_color,
                                d_molA_type + d_NmolA_shape,
                                molB + d_NmolB_shape * D, d_NmolB_color,
                                d_molB_type + d_NmolB_shape, nullptr, nullptr);
  }
  scores = calcScores(scores[3], scores[4]);
  return scores;
}

std::array<DTYPE, 5> SingleConformerAlignment::calcScores(
    const DTYPE shapeOvVol, const DTYPE colorOvVol) const {
  std::array<DTYPE, 5> scores{0.0, 0.0, 0.0, 0.0, 0.0};
  scores[3] = shapeOvVol;
  scores[4] = colorOvVol;
  scores[1] = scores[3] / (d_shapeVolA + d_shapeVolB - scores[3]);
  if (d_NmolA_color && d_NmolB_color && d_colorVolA > 0.0 &&
      d_colorVolB > 0.0) {
    scores[2] = scores[4] / (d_colorVolA + d_colorVolB - scores[4]);
    scores[0] = scores[1] * (1 - d_mixing_param) + scores[2] * d_mixing_param;
  } else {
    scores[0] = scores[1];
  }
  return scores;
}

void SingleConformerAlignment::calcGradients(const DTYPE *molA,
                                             const DTYPE *molB,
                                             const DTYPE *quat,
                                             DTYPE *gradients) const {
  calcVolAndGrads(molA, d_NmolA_shape, molB, d_NmolB_shape, d_all_carbon_radii,
                  quat, gradients);
}

namespace {
// Set of values to convert the cartesian gradients to quaternion gradients
using Adapt = std::array<DTYPE, 12>;
void makeQuaternionAdapters(const DTYPE *quat, const DTYPE *molB, int numBPts,
                            std::vector<Adapt> &derivAdapters) {
  // for ease of ref
  auto q = quat[0];
  auto r = quat[1];
  auto s = quat[2];
  auto u = quat[3];
  auto coef = 1.0 / (q * q + r * r + s * s + u * u);
  for (int i = 0; i < 4 * numBPts; i += 4) {
    auto x = molB[i];
    auto y = molB[i + 1];
    auto z = molB[i + 2];
    auto dx_dq = coef * 2.0 * (q * x + u * y - s * z);
    auto dx_dr = coef * 2.0 * (r * x + s * y + u * z);
    auto dy_dr = coef * 2.0 * (s * x - r * y + q * z);
    auto dx_du = coef * 2.0 * (-u * x + q * y + r * z);
    auto dz_ds = dx_dq;
    auto dy_du = -dx_dq;
    auto dy_ds = dx_dr;
    auto dz_du = dx_dr;
    auto dx_ds = -dy_dr;
    auto dz_dq = dy_dr;
    auto dy_dq = dx_du;
    auto dz_dr = -dx_du;
    derivAdapters.emplace_back(Adapt{dx_dq, dy_dq, dz_dq, dx_dr, dy_dr, dz_dr,
                                     dx_ds, dy_ds, dz_ds, dx_du, dy_du, dz_du});
  }
}
}  // namespace

// Compute the volume overlap and optionally "quaternion" gradients for the
// overlap volume of A and B, wrt B.  molB is the original coords of B,
// molBT is those subject to any transformation applied by the quaternion
// we're using to optimise the overlap volume.  If gradients is null, they
// won't be calculated.  They are assumed to be initialised correctly.
DTYPE calcVolAndGrads(const DTYPE *molA, int numAPts, const DTYPE *molB,
                      int numBPts, const bool all_carbon_radii,
                      const DTYPE *quat, DTYPE *gradients) {
  std::vector<Adapt> quatAdapters;
  if (gradients) {
    makeQuaternionAdapters(quat, molB, numBPts, quatAdapters);
  }

  static constexpr DTYPE KAPPA = 2.41798793102;
  static const DTYPE CARBON_A = KAPPA / (1.7 * 1.7);
  static const DTYPE CARBON_BIT = pow(PI / (2 * CARBON_A), 1.5);
  DTYPE vol = 0.0;
  for (int i = 0; i < numAPts * 4; i += 4) {
    const auto ai = all_carbon_radii ? CARBON_A : molA[i + 3];
    for (int j = 0, j_idx = 0; j < numBPts * 4; j += 4, j_idx++) {
      auto dx = molA[i] - molB[j];
      auto dy = molA[i + 1] - molB[j + 1];
      auto dz = molA[i + 2] - molB[j + 2];
      auto d2 = dx * dx + dy * dy + dz * dz;
      const auto aj = all_carbon_radii ? CARBON_A : molB[j + 3];
      auto mult = -(ai * aj) / (ai + aj);
      auto kij = exp(mult * d2);

      auto vij = all_carbon_radii ? 8 * kij * CARBON_BIT
                                  : 8 * kij * pow((PI / (ai + aj)), 1.5);
      vol += vij;
      if (gradients) {
        auto r = 2.0 * vij * mult;
        // Use the adapters to calculate the gradients in quaternion space
        gradients[0] +=
            r * (dx * quatAdapters[j_idx][0] + dy * quatAdapters[j_idx][1] +
                 dz * quatAdapters[j_idx][2]);
        gradients[1] +=
            r * (dx * quatAdapters[j_idx][3] + dy * quatAdapters[j_idx][4] +
                 dz * quatAdapters[j_idx][5]);
        gradients[2] +=
            r * (dx * quatAdapters[j_idx][6] + dy * quatAdapters[j_idx][7] +
                 dz * quatAdapters[j_idx][8]);
        gradients[3] +=
            r * (dx * quatAdapters[j_idx][9] + dy * quatAdapters[j_idx][10] +
                 dz * quatAdapters[j_idx][11]);
        gradients[4] += r * dx;
        gradients[5] += r * dy;
        gradients[6] += r * dz;
      }
    }
  }
  return vol;
}

DTYPE calcVolAndGrads(const DTYPE *molA, int numAPts, const int *molATypes,
                      const DTYPE *molB, int numBPts, const int *molBTypes,
                      const DTYPE *quat, DTYPE *gradients) {
  DTYPE vol = 0.0;
  std::vector<Adapt> quatAdapters;
  if (gradients) {
    makeQuaternionAdapters(quat, molB, numBPts, quatAdapters);
  }

  double sum_d2 = 0.0;
  for (int i = 0, i_idx = 0; i < numAPts * 4; i += 4, i_idx++) {
    const auto ai = molA[i + 3];
    const auto aType = molATypes[i_idx];
    for (int j = 0, j_idx = 0; j < numBPts * 4; j += 4, j_idx++) {
      const auto bType = molBTypes[j_idx];
      if (aType != bType) {
        continue;
      }
      auto dx = molA[i] - molB[j];
      auto dy = molA[i + 1] - molB[j + 1];
      auto dz = molA[i + 2] - molB[j + 2];
      auto d2 = dx * dx + dy * dy + dz * dz;
      sum_d2 += d2;
      const auto aj = molB[j + 3];
      auto mult = -(ai * aj) / (ai + aj);
      auto kij = exp(mult * d2);

      auto vij = 8 * kij * pow((PI / (ai + aj)), 1.5);
      vol += vij;
      if (gradients) {
        auto r = 2.0 * vij * mult;
        // Use the adapters to calculate the gradients in quaternion space
        gradients[0] +=
            r * (dx * quatAdapters[j_idx][0] + dy * quatAdapters[j_idx][1] +
                 dz * quatAdapters[j_idx][2]);
        gradients[1] +=
            r * (dx * quatAdapters[j_idx][3] + dy * quatAdapters[j_idx][4] +
                 dz * quatAdapters[j_idx][5]);
        gradients[2] +=
            r * (dx * quatAdapters[j_idx][6] + dy * quatAdapters[j_idx][7] +
                 dz * quatAdapters[j_idx][8]);
        gradients[3] +=
            r * (dx * quatAdapters[j_idx][9] + dy * quatAdapters[j_idx][10] +
                 dz * quatAdapters[j_idx][11]);
        gradients[4] += r * dx;
        gradients[5] += r * dy;
        gradients[6] += r * dz;
      }
    }
  }
  return vol;
}

void SingleConformerAlignment::calcVolumeAndGradients(
    const std::array<DTYPE, 7> &quat, DTYPE &shapeOvlpVol, DTYPE &colorOvlpVol,
    std::array<DTYPE, 7> *gradients) {
  RDGeom::Transform3D transformB;
  // Leave B at the origin, and move A to meet it.
  RDGeom::Point3D translateA{-quat[4], -quat[5], -quat[6]};
  translateShape(d_molA, d_molAT, d_NmolA_shape + d_NmolA_color, translateA);
  double tq[4]{quat[0], quat[1], quat[2], quat[3]};
  transformB.SetRotationFromQuaternion(tq);
  applyTransformToShape(d_molB, d_molBT, d_NmolB_shape + d_NmolB_color,
                        transformB);

  // We assume that d_molAT was once initialised to d_MolA and the same with
  // B so that the types are already there.
  if (gradients) {
    (*gradients)[0] = (*gradients)[1] = (*gradients)[2] = (*gradients)[3] =
        (*gradients)[4] = (*gradients)[5] = (*gradients)[6] = 0.0;
    shapeOvlpVol =
        calcVolAndGrads(d_molAT, d_NmolA_shape, d_molBT, d_NmolB_shape,
                        d_all_carbon_radii, quat.data(), gradients->data());
  } else {
    shapeOvlpVol = calcVolAndGrads(d_molAT, d_NmolA_shape, d_molBT,
                                   d_NmolB_shape, true, quat.data(), nullptr);
  }
  if (d_optimMode == OptimMode::SHAPE_PLUS_COLOR && gradients) {
    std::array<DTYPE, 7> colorGrads{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    colorOvlpVol = calcVolAndGrads(
        d_molAT + 4 * d_NmolA_shape, d_NmolA_color, d_molA_type + d_NmolA_shape,
        d_molBT + 4 * d_NmolB_shape, d_NmolB_color, d_molB_type + d_NmolB_shape,
        quat.data(), colorGrads.data());
    // The color gradients are normally dwarfed by the shape gradients, so
    // normalize them and then mix by the same rule as the final score.
    auto shapeSum = sqrt(std::accumulate(
        gradients->begin(), gradients->end(), 0.0,
        [](const auto init, const auto g) -> DTYPE { return init + g * g; }));
    auto colorSum = sqrt(std::accumulate(
        colorGrads.begin(), colorGrads.end(), 0.0,
        [](const auto init, const auto g) -> DTYPE { return init + g * g; }));
    auto ratio = shapeSum / colorSum;
    std::transform(
        gradients->begin(), gradients->end(), colorGrads.begin(),
        gradients->begin(), [&](const auto g1, const auto g2) -> DTYPE {
          return g1 * (1 - d_mixing_param) + g2 * ratio * d_mixing_param;
        });
  } else {
    colorOvlpVol = 0.0;
  }
}

bool SingleConformerAlignment::doOverlay(std::array<DTYPE, 20> &scores) {
  std::array<double, 7> quatTrans{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  auto res = optimise(quatTrans);

  // Get the final coords for B into d_molBT, and compute the scores
  RDGeom::Transform3D xform;
  xform.SetRotationFromQuaternion(quatTrans.data());
  xform.SetTranslation(
      RDGeom::Point3D{quatTrans[4], quatTrans[5], quatTrans[6]});
  applyTransformToShape(d_molB, d_molBT, d_NmolB_shape + d_NmolB_color, xform);

  auto tscores = calcScores(d_molA, d_molBT, true);
  scores[0] = tscores[0];
  scores[1] = tscores[1];
  scores[2] = tscores[2];
  scores[3] = tscores[3];
  scores[4] = tscores[4];
  scores[5] = d_shapeVolA;
  scores[6] = d_colorVolA;
  scores[7] = d_shapeVolB;
  scores[8] = d_colorVolB;
  scores[9] = quatTrans[0];
  scores[10] = quatTrans[1];
  scores[11] = quatTrans[2];
  scores[12] = quatTrans[3];
  scores[13] = quatTrans[4];
  scores[14] = quatTrans[5];
  scores[15] = quatTrans[6];
  scores[16] = 0.0;
  scores[17] = 0.0;
  scores[18] = 0.0;
  scores[19] = 0.0;
  return res;
}

namespace {
DTYPE oneStep(DTYPE grad, double stepSize, DTYPE quatTrans, DTYPE oldGrad,
              DTYPE oldQuatTrans) {
  DTYPE step = 0.0;
  if (std::signbit(grad) != std::signbit(oldGrad)) {
    step = (((quatTrans * fabs(oldGrad)) + (oldQuatTrans * fabs(grad))) /
            (fabs(oldGrad) + fabs(grad) + fabs(grad))) -
           quatTrans;
    DTYPE newStep = stepSize * grad;
    if (fabs(step) > fabs(newStep)) {
      // This is definitely what the PubChem code says!  I read it as keeping
      // the sign of step, but the value of newStep.
      step *= fabs(newStep / step);
    }
  } else {
    step = stepSize * grad;
  }
  return step;
}

void calcStep(std::array<DTYPE, 7> &grad, DTYPE qStepSize, DTYPE tStepSize,
              std::array<DTYPE, 7> &oldGrad, std::array<DTYPE, 7> &quatTrans,
              std::array<DTYPE, 7> &oldQuatTrans, unsigned int iter,
              std::array<DTYPE, 7> &step) {
  step[0] = 0.0;
  if (iter == 0) {
    // 1st iteration, use default step sizes
    step[1] = qStepSize * grad[1];
    step[2] = qStepSize * grad[2];
    step[3] = qStepSize * grad[3];
    step[4] = tStepSize * grad[4];
    step[5] = tStepSize * grad[5];
    step[6] = tStepSize * grad[6];
  } else {
    step[1] =
        oneStep(grad[1], qStepSize, quatTrans[1], oldGrad[1], oldQuatTrans[1]);
    step[2] =
        oneStep(grad[2], qStepSize, quatTrans[2], oldGrad[2], oldQuatTrans[2]);
    step[3] =
        oneStep(grad[3], qStepSize, quatTrans[3], oldGrad[3], oldQuatTrans[3]);
    step[4] =
        oneStep(grad[4], tStepSize, quatTrans[4], oldGrad[4], oldQuatTrans[4]);
    step[5] =
        oneStep(grad[5], tStepSize, quatTrans[5], oldGrad[5], oldQuatTrans[5]);
    step[6] =
        oneStep(grad[6], tStepSize, quatTrans[6], oldGrad[6], oldQuatTrans[6]);
  }
}

DTYPE constrainStep(DTYPE maxStep, DTYPE *step) {
  DTYPE mStep = std::max({fabs(step[0]), fabs(step[1]), fabs(step[2])});
  if (mStep > maxStep) {
    DTYPE scaleFactor = maxStep / mStep;
    if (fabs(step[0] > maxStep)) {
      step[0] *= scaleFactor;
    }
    if (fabs(step[1] > maxStep)) {
      step[1] *= scaleFactor;
    }
    if (fabs(step[2] > maxStep)) {
      step[2] *= scaleFactor;
    }
  }
  DTYPE quatSquared = step[0] * step[0] + step[1] * step[1] + step[2] * step[2];
  if (quatSquared > 1.0) {
    DTYPE scaleFactor = 1.0 / (2.0 * quatSquared);
    step[0] *= scaleFactor;
    step[1] *= scaleFactor;
    step[2] *= scaleFactor;
  }
  return mStep;
}

std::array<DTYPE, 7> combineQuatTrans(const std::array<DTYPE, 7> &q1,
                                      const std::array<DTYPE, 7> &q2) {
  std::array<DTYPE, 7> res;
  // Multiply the quaternions, which are assumed to be normalised.
  res[0] = q1[0] * q2[0] - q1[1] * q2[1] - q1[2] * q2[2] - q1[3] * q2[3];
  res[1] = q1[0] * q2[1] + q1[1] * q2[0] + q1[2] * q2[3] - q1[3] * q2[2];
  res[2] = q1[0] * q2[2] - q1[1] * q2[3] + q1[2] * q2[0] + q1[3] * q2[1];
  res[3] = q1[0] * q2[3] + q1[1] * q2[2] - q1[2] * q2[1] + q1[3] * q2[0];

  // Add the translations
  res[4] = q1[4] + q2[4];
  res[5] = q1[5] + q2[5];
  res[6] = q1[6] + q2[6];

  return res;
}

DTYPE oneReduceStep(DTYPE grad, DTYPE oldGrad, DTYPE quatTrans,
                    DTYPE oldQuatTrans, DTYPE stepSize, DTYPE step) {
  if (std::signbit(grad) != std::signbit(oldGrad)) {
    step = (((quatTrans * fabs(oldGrad)) + (oldQuatTrans * fabs(grad))) /
            (fabs(oldGrad + fabs(grad)))) -
           quatTrans;
    DTYPE newStep = stepSize * grad;
    if (fabs(step) > fabs(newStep)) {
      step *= fabs(newStep / step);
    }
  } else if (fabs(grad) <= 1.0) {
    step = stepSize * grad;
  } else if (fabs(grad) > fabs(oldGrad)) {
    // Going wrong way relative to other components?
    step += stepSize * grad;
  } else {
    DTYPE delta = grad * (step / (oldGrad - grad));
    if (fabs(delta) > fabs(step * 0.1) && fabs(delta) > 0.001) {
      delta *= 0.0005 / fabs(delta);
    }
    step += delta;
  }
  return step;
}

void reduceStep(std::array<DTYPE, 7> &grad, std::array<DTYPE, 7> &oldGrad,
                std::array<DTYPE, 7> quatTrans,
                std::array<DTYPE, 7> &oldQuatTrans, unsigned int lineIter,
                std::array<DTYPE, 7> &step, DTYPE &qStepSize,
                DTYPE &tStepSize) {
  if (lineIter == 2) {
    qStepSize *= 0.1;
    tStepSize *= 0.1;
    step[1] = qStepSize * oldGrad[1];
    step[2] = qStepSize * oldGrad[2];
    step[3] = qStepSize * oldGrad[3];
    step[4] = tStepSize * oldGrad[4];
    step[5] = tStepSize * oldGrad[5];
    step[6] = tStepSize * oldGrad[6];
    qStepSize *= 5.0;
    tStepSize *= 5.0;
  } else {
    step[1] = oneReduceStep(grad[1], oldGrad[1], quatTrans[1], oldQuatTrans[1],
                            qStepSize, step[1]);
    step[2] = oneReduceStep(grad[2], oldGrad[2], quatTrans[2], oldQuatTrans[2],
                            qStepSize, step[2]);
    step[3] = oneReduceStep(grad[3], oldGrad[3], quatTrans[3], oldQuatTrans[3],
                            qStepSize, step[3]);
    // The original PubChem code used qStepSize for all 6 of these updates, but
    // that would appear to be a cut-and-paste error.
    step[4] = oneReduceStep(grad[4], oldGrad[4], quatTrans[4], oldQuatTrans[4],
                            tStepSize, step[4]);
    step[5] = oneReduceStep(grad[5], oldGrad[5], quatTrans[5], oldQuatTrans[5],
                            tStepSize, step[5]);
    step[6] = oneReduceStep(grad[6], oldGrad[6], quatTrans[6], oldQuatTrans[6],
                            tStepSize, step[6]);
  }
}
}  // namespace

// The optimisation follows closely the procedure used in the PubChem
// code from
// https://github.com/ncbi/pubchem-align3d/blob/main/shape_neighbor.cpp
// Original Authors:  Evan Bolton, Leonid Zaslavsky, Paul Thiessen
bool SingleConformerAlignment::optimise(std::array<DTYPE, 7> &quatTrans) {
  const DTYPE maxQuaternionStep = 0.075;   // Maximum step size for quaternion
  const DTYPE maxTranslationStep = 0.500;  // Maximum step size for translation
  const DTYPE minQuaternionStep =
      0.0002;  // Convergence criteria for quaternion
  const DTYPE minTranslationStep =
      0.0020;                          // Convergence criteria for translation
  const DTYPE minScoreDiff = 0.00002;  // Convergence criteria for comboScore

  std::array<DTYPE, 7> grad;
  std::array<DTYPE, 7> oldGrad;
  std::array<DTYPE, 7> oldQuatTrans;
  std::array<DTYPE, 7> step;
  DTYPE shapeOvlpVol, colorOvlpVol, comboScore = 0.0;
  // Steps for the quaternion and translation.
  DTYPE qStepSize{-0.001}, tStepSize{-0.01};
  bool finished = false;
  for (unsigned iter = 0; iter < d_maxIts; iter++) {
    calcVolumeAndGradients(quatTrans, shapeOvlpVol, colorOvlpVol, &grad);
    auto scores = calcScores(shapeOvlpVol, colorOvlpVol);
    comboScore = scores[0];
    calcStep(grad, qStepSize, tStepSize, oldGrad, quatTrans, oldQuatTrans, iter,
             step);

    // In case we have to backtrack
    DTYPE oldComboScore = comboScore;
    oldQuatTrans = quatTrans;
    oldGrad = grad;

    // What the PubChem code calls "Line search (sort of) loop"
    bool converged = false;
    for (unsigned int lineIter = 0; !converged; lineIter++) {
      // Check that the absolute max step size does not go beyond some
      // reasonable size
      DTYPE mqStep = constrainStep(maxQuaternionStep, step.data() + 1);
      DTYPE mtStep = constrainStep(maxTranslationStep, step.data() + 4);
      if (mqStep <= minQuaternionStep && mtStep <= minTranslationStep) {
        converged = true;
        comboScore = 0.0;  // Make sure we return to the old one.
        break;
      }
      // Calculate the 0th component of the quaternion.  Obviously it
      // relies on the other 3 components being small
      DTYPE quatSquared =
          step[1] * step[1] + step[2] * step[2] + step[3] * step[3];
      step[0] = sqrt(1.0 - quatSquared);
      // Update the quaternion with the step, multiplying them.
      auto newQuatTrans = combineQuatTrans(quatTrans, step);
      calcVolumeAndGradients(newQuatTrans, shapeOvlpVol, colorOvlpVol, &grad);
      auto scores = calcScores(shapeOvlpVol, colorOvlpVol);
      comboScore = scores[0];
      // if we made a good step, keep the quaternion and we're done
      if (comboScore > oldComboScore) {
        quatTrans = newQuatTrans;
        break;
      }
      if (lineIter > 2) {
        converged = true;
        quatTrans = newQuatTrans;
        break;
      }
      // It got worse, so reduce the step.
      reduceStep(grad, oldGrad, newQuatTrans, oldQuatTrans, lineIter, step,
                 qStepSize, tStepSize);
      quatTrans = oldQuatTrans;
    }  // End of line search
    // Did it converge?
    if (converged || minScoreDiff > (comboScore - oldComboScore)) {
      if (oldComboScore > comboScore) {
        // The previous step was better, so keep it.
        comboScore = oldComboScore;
        quatTrans = oldQuatTrans;
      }
      finished = true;
      break;
    }
  }
  return finished;
}
}  // namespace GaussianShape
}  // namespace RDKit
