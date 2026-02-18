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

namespace RDKit {
namespace GaussianShape {

SingleConformerAlignment::SingleConformerAlignment(
    const std::vector<double> &ref, const int *refTypes,
    const std::unique_ptr<boost::dynamic_bitset<>> &refCarbonRadii,
    int nRefShape, int nRefColor, double refShapeVol, double refColorVol,
    const std::vector<double> &fit, const int *fitTypes,
    const std::unique_ptr<boost::dynamic_bitset<>> &fitCarbonRadii,
    int nFitShape, int nFitColor, double fitShapeVol, double fitColorVol,
    const std::array<double, 7> &initQuatTrans, OptimMode optimMode,
    double mixingParam, bool useCutoff, double distCutoff,
    double shapeConvergenceCriterion, unsigned int maxIts)
    : d_ref(ref),
      d_refTypes(refTypes),
      d_refCarbonRadii(refCarbonRadii),
      d_nRefShape(nRefShape),
      d_nRefColor(nRefColor),
      d_refShapeVol(refShapeVol),
      d_refColorVol(refColorVol),
      d_fit(fit),
      d_fitTypes(fitTypes),
      d_fitCarbonRadii(fitCarbonRadii),
      d_nFitShape(nFitShape),
      d_nFitColor(nFitColor),
      d_fitShapeVol(fitShapeVol),
      d_fitColorVol(fitColorVol),
      d_initQuatTrans(initQuatTrans),
      d_optimMode(optimMode),
      d_mixingParam(mixingParam),
      d_useCutoff(useCutoff),
      d_distCutoff2(distCutoff * distCutoff),
      d_shapeConvergenceCriterion(shapeConvergenceCriterion),
      d_maxIts(maxIts) {
  // Move the reference by initialTrans, leaving fit at the origin where
  // the rotations work properly.  Apply the initial rotation to the fit.
  translateShape(d_ref, RDGeom::Point3D{d_initQuatTrans[4], d_initQuatTrans[5],
                                        d_initQuatTrans[6]});
  RDGeom::Transform3D xform;
  xform.SetRotationFromQuaternion(d_initQuatTrans.data());
  applyTransformToShape(d_fit, xform);
  d_refTemp.resize(d_ref.size());
  d_fitTemp.resize(d_fit.size());
  d_gradConverters.resize(d_nFitShape + d_nFitColor);
}

void SingleConformerAlignment::getFinalQuatTrans(
    RDGeom::Transform3D &xform) const {
  RDGeom::Transform3D tmp;
  tmp.SetRotationFromQuaternion(d_quatTrans.data());
  tmp.SetTranslation(
      RDGeom::Point3D{d_quatTrans[4], d_quatTrans[5], d_quatTrans[6]});
  RDGeom::Transform3D reverseInitialTrans;
  reverseInitialTrans.SetTranslation(RDGeom::Point3D{
      -d_initQuatTrans[4], -d_initQuatTrans[5], -d_initQuatTrans[6]});
  RDGeom::Transform3D initialRot;
  initialRot.SetRotationFromQuaternion(d_initQuatTrans.data());
  auto tt = reverseInitialTrans * tmp * initialRot;
  copyTransform(tt, xform);
}

std::array<double, 5> SingleConformerAlignment::calcScores(
    const double *ref, const double *fit, bool includeColor) const {
  std::array<double, 5> scores{0.0, 0.0, 0.0, 0.0, 0.0};
  scores[3] = calcVolAndGrads(ref, d_nRefShape, d_refCarbonRadii, fit,
                              d_nFitShape, d_fitCarbonRadii, d_gradConverters,
                              d_useCutoff, d_distCutoff2, nullptr, nullptr);
  if (d_nRefColor && d_nFitColor &&
      (d_optimMode == OptimMode::SHAPE_PLUS_COLOR || includeColor)) {
    scores[4] = calcVolAndGrads(ref + d_nRefShape * D, d_nRefColor,
                                d_refTypes + d_nRefShape, fit + d_nFitShape * D,
                                d_nFitColor, d_fitTypes + d_nFitShape,
                                d_nFitShape, d_gradConverters, d_useCutoff,
                                d_distCutoff2, nullptr, nullptr);
  }
  scores = calcScores(scores[3], scores[4], includeColor);
  return scores;
}

std::array<double, 5> SingleConformerAlignment::calcScores(bool includeColor) {
  applyQuatTrans(d_quatTrans);
  return calcScores(d_refTemp.data(), d_fitTemp.data(), includeColor);
}

std::array<double, 5> SingleConformerAlignment::calcScores(
    const double shapeOvVol, const double colorOvVol, bool includeColor) const {
  std::array<double, 5> scores{0.0, 0.0, 0.0, 0.0, 0.0};
  scores[3] = shapeOvVol;
  scores[4] = colorOvVol;
  scores[1] = scores[3] / (d_refShapeVol + d_fitShapeVol - scores[3]);
  if (d_nRefColor && d_nFitColor && d_refColorVol > 0.0 &&
      d_fitColorVol > 0.0 && includeColor) {
    scores[2] = scores[4] / (d_refColorVol + d_fitColorVol - scores[4]);
    scores[0] = scores[1] * (1 - d_mixingParam) + scores[2] * d_mixingParam;
  } else {
    scores[0] = scores[1];
  }
  return scores;
}

namespace {
// Set of values to convert the cartesian gradients to quaternion gradients.
// This uses the chain rule: the dV/qQ = (dV/dr) * (dr/dQ) where V is the
// volume overlap and r is the Cartesian space.  Assumes gradConverters
// is already the correct size.
void cartToQuatGrads(const double *quat, const double *mol, int numBPts,
                     std::vector<std::array<double, 12>> &gradConverters,
                     int gradConvOffset) {
  // for ease of ref
  auto q = quat[0];
  auto r = quat[1];
  auto s = quat[2];
  auto u = quat[3];
  auto coef = 1.0 / (q * q + r * r + s * s + u * u);
  for (int i = 0, j = gradConvOffset; i < 4 * numBPts; i += 4, ++j) {
    auto x = mol[i];
    auto y = mol[i + 1];
    auto z = mol[i + 2];
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
    gradConverters[j][0] = dx_dq;
    gradConverters[j][1] = dy_dq;
    gradConverters[j][2] = dz_dq;
    gradConverters[j][3] = dx_dr;
    gradConverters[j][4] = dy_dr;
    gradConverters[j][5] = dz_dr;
    gradConverters[j][6] = dx_ds;
    gradConverters[j][7] = dy_ds;
    gradConverters[j][8] = dz_ds;
    gradConverters[j][9] = dx_du;
    gradConverters[j][10] = dy_du;
    gradConverters[j][11] = dz_du;
  }
}
}  // namespace

// atoms/shape features
double calcVolAndGrads(
    const double *ref, int numRefPts,
    const std::unique_ptr<boost::dynamic_bitset<>> &refCarbonRadii,
    const double *fit, int numFitPts,
    const std::unique_ptr<boost::dynamic_bitset<>> &fitCarbonRadii,
    std::vector<std::array<double, 12>> &gradConverters, bool useCutoff,
    double distCutoff2, const double *quat, double *gradients) {
  if (gradients) {
    cartToQuatGrads(quat, fit, numFitPts, gradConverters, 0);
  }
  static constexpr double KAPPA = 2.41798793102;
  static const double CARBON_A = KAPPA / (1.7 * 1.7);
  static const double CARBON_BIT = 8.0 * pow(PI / (2 * CARBON_A), 1.5);
  double vol = 0.0;
  double vij;
  for (int i = 0, i_idx = 0; i < numRefPts * 4; i += 4, i_idx++) {
    const auto ai = ref[i + 3];
    for (int j = 0, j_idx = 0; j < numFitPts * 4; j += 4, j_idx++) {
      auto dx = ref[i] - fit[j];
      auto dy = ref[i + 1] - fit[j + 1];
      auto dz = ref[i + 2] - fit[j + 2];
      auto d2 = dx * dx + dy * dy + dz * dz;
      if (useCutoff && d2 > distCutoff2) {
        continue;
      }
      const auto aj = fit[j + 3];
      auto mult = -(ai * aj) / (ai + aj);
      auto kij = exp(mult * d2);
      if ((!refCarbonRadii && !fitCarbonRadii) ||
          (refCarbonRadii && fitCarbonRadii && (*refCarbonRadii)[i_idx] &&
           (*fitCarbonRadii)[j_idx])) {
        vij = kij * CARBON_BIT;
      } else {
        vij = 8 * kij * pow((PI / (ai + aj)), 1.5);
      }
      vol += vij;
      if (gradients) {
        auto r = 2.0 * vij * mult;
        // Use the gradient converters to calculate the gradients in quaternion
        // space.
        // The zeroth gradient is never used, so don't waste time calculating
        // it but leave the code here for completeness and possible future use.
        // gradients[0] +=
        //     r * (dx * gradConverters[j_idx][0] + dy *
        //     gradConverters[j_idx][1] +
        //          dz * gradConverters[j_idx][2]);
        gradients[1] +=
            r * (dx * gradConverters[j_idx][3] + dy * gradConverters[j_idx][4] +
                 dz * gradConverters[j_idx][5]);
        gradients[2] +=
            r * (dx * gradConverters[j_idx][6] + dy * gradConverters[j_idx][7] +
                 dz * gradConverters[j_idx][8]);
        gradients[3] += r * (dx * gradConverters[j_idx][9] +
                             dy * gradConverters[j_idx][10] +
                             dz * gradConverters[j_idx][11]);
        gradients[4] += r * dx;
        gradients[5] += r * dy;
        gradients[6] += r * dz;
      }
    }
  }
  return vol;
}

// color features
double calcVolAndGrads(const double *ref, int numRefPts, const int *refTypes,
                       const double *fit, int numFitPts, const int *fitTypes,
                       int numFitShape,
                       std::vector<std::array<double, 12>> &gradConverters,
                       const bool useCutoff, const double distCutoff2,
                       const double *quat, double *gradients) {
  double vol = 0.0;
  if (gradients) {
    cartToQuatGrads(quat, fit, numFitPts, gradConverters, numFitShape);
  }

  for (int i = 0, i_idx = 0; i < numRefPts * 4; i += 4, i_idx++) {
    const auto ai = ref[i + 3];
    const auto aType = refTypes[i_idx];
    for (int j = 0, j_idx = 0; j < numFitPts * 4; j += 4, j_idx++) {
      const auto bType = fitTypes[j_idx];
      if (aType != bType) {
        continue;
      }
      auto dx = ref[i] - fit[j];
      auto dy = ref[i + 1] - fit[j + 1];
      auto dz = ref[i + 2] - fit[j + 2];
      auto d2 = dx * dx + dy * dy + dz * dz;
      if (useCutoff && d2 > distCutoff2) {
        continue;
      }
      const auto aj = fit[j + 3];
      auto mult = -(ai * aj) / (ai + aj);
      auto kij = exp(mult * d2);

      auto pi_ai_aj = PI / (ai + aj);
      auto vij = 8 * kij * pi_ai_aj * std::sqrt(pi_ai_aj);
      vol += vij;
      if (gradients) {
        auto r = 2.0 * vij * mult;
        // Use the converters to calculate the gradients in quaternion space.
        // The zeroth gradient is never used, so don't waste time calculating
        // it but leave the code here for completeness and possible future use.
        // gradients[0] +=
        //     r * (dx * gradConverters[j_idx][0] + dy *
        //     gradConverters[j_idx][1] +
        //          dz * gradConverters[j_idx][2]);
        gradients[1] +=
            r * (dx * gradConverters[j_idx][3] + dy * gradConverters[j_idx][4] +
                 dz * gradConverters[j_idx][5]);
        gradients[2] +=
            r * (dx * gradConverters[j_idx][6] + dy * gradConverters[j_idx][7] +
                 dz * gradConverters[j_idx][8]);
        gradients[3] += r * (dx * gradConverters[j_idx][9] +
                             dy * gradConverters[j_idx][10] +
                             dz * gradConverters[j_idx][11]);
        gradients[4] += r * dx;
        gradients[5] += r * dy;
        gradients[6] += r * dz;
      }
    }
  }
  return vol;
}

void SingleConformerAlignment::applyQuatTrans(
    const std::array<double, 7> &quatTrans) {
  // Leave fit at the origin, and move ref to meet it.
  RDGeom::Point3D translateA{-quatTrans[4], -quatTrans[5], -quatTrans[6]};
  translateShape(d_ref.data(), d_refTemp.data(), d_nRefShape + d_nRefColor,
                 translateA);
  // Rotate fit by quaternion
  // double tq[4]{quatTrans[0], quatTrans[1], quatTrans[2], quatTrans[3]};
  RDGeom::Transform3D transformB;
  transformB.SetRotationFromQuaternion(quatTrans.data());
  applyTransformToShape(d_fit.data(), d_fitTemp.data(),
                        d_nFitShape + d_nFitColor, transformB);
}

void SingleConformerAlignment::calcVolumeAndGradients(
    const std::array<double, 7> &quatTrans, double &shapeOvlpVol,
    double &colorOvlpVol, std::array<double, 7> &gradients) {
  // Set the coords up.
  applyQuatTrans(quatTrans);
  // We assume that d_refTemp was once initialised to d_ref and the same with
  // fit so that the radii are already there.
  gradients[0] = gradients[1] = gradients[2] = gradients[3] = gradients[4] =
      gradients[5] = gradients[6] = 0.0;
  shapeOvlpVol = calcVolAndGrads(
      d_refTemp.data(), d_nRefShape, d_refCarbonRadii, d_fitTemp.data(),
      d_nFitShape, d_fitCarbonRadii, d_gradConverters, d_useCutoff,
      d_distCutoff2, quatTrans.data(), gradients.data());
  if (d_optimMode == OptimMode::SHAPE_PLUS_COLOR) {
    std::array<double, 7> colorGrads{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    colorOvlpVol = calcVolAndGrads(
        d_refTemp.data() + 4 * d_nRefShape, d_nRefColor,
        d_refTypes + d_nRefShape, d_fitTemp.data() + 4 * d_nFitShape,
        d_nFitColor, d_fitTypes + d_nFitShape, d_nFitShape, d_gradConverters,
        d_useCutoff, d_distCutoff2, quatTrans.data(), colorGrads.data());
    // The color gradients are normally dwarfed by the shape gradients, so
    // normalize them and then mix by the same rule as the final score.
    auto shapeSum = sqrt(std::accumulate(
        gradients.begin() + 1, gradients.end(), 0.0,
        [](const auto init, const auto g) -> double { return init + g * g; }));
    auto colorSum = sqrt(std::accumulate(
        colorGrads.begin() + 1, colorGrads.end(), 0.0,
        [](const auto init, const auto g) -> double { return init + g * g; }));
    auto ratio = shapeSum / colorSum;
    std::transform(
        gradients.begin() + 1, gradients.end(), colorGrads.begin(),
        gradients.begin() + 1, [&](const auto g1, const auto g2) -> double {
          return g1 * (1 - d_mixingParam) + g2 * ratio * d_mixingParam;
        });
  } else {
    colorOvlpVol = 0.0;
  }
}

bool SingleConformerAlignment::doOverlay(std::array<double, 20> &scores,
                                         unsigned int cycle) {
  unsigned int maxIters = cycle == 0 ? 10 : d_maxIts - 10;
  auto res = optimise(maxIters);

  // Get the final coords for fit into d_fitTemp, and compute the scores
  RDGeom::Transform3D xform;
  xform.SetRotationFromQuaternion(d_quatTrans.data());
  xform.SetTranslation(
      RDGeom::Point3D{d_quatTrans[4], d_quatTrans[5], d_quatTrans[6]});
  applyTransformToShape(d_fit.data(), d_fitTemp.data(),
                        d_nFitShape + d_nFitColor, xform);

  auto tscores = calcScores(d_ref.data(), d_fitTemp.data(), true);
  scores[0] = tscores[0];
  scores[1] = tscores[1];
  scores[2] = tscores[2];
  scores[3] = tscores[3];
  scores[4] = tscores[4];
  scores[5] = d_refShapeVol;
  scores[6] = d_refColorVol;
  scores[7] = d_fitShapeVol;
  scores[8] = d_fitColorVol;
  scores[9] = d_quatTrans[0];
  scores[10] = d_quatTrans[1];
  scores[11] = d_quatTrans[2];
  scores[12] = d_quatTrans[3];
  scores[13] = d_quatTrans[4];
  scores[14] = d_quatTrans[5];
  scores[15] = d_quatTrans[6];
  scores[16] = 0.0;
  scores[17] = 0.0;
  scores[18] = 0.0;
  scores[19] = 0.0;
  return res;
}

namespace {
double oneStep(double grad, double stepSize, double quatTrans, double oldGrad,
               double oldQuatTrans) {
  double step = 0.0;
  if (std::signbit(grad) != std::signbit(oldGrad)) {
    step = (((quatTrans * fabs(oldGrad)) + (oldQuatTrans * fabs(grad))) /
            (fabs(oldGrad) + fabs(grad) + fabs(grad))) -
           quatTrans;
    double newStep = stepSize * grad;
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

void calcStep(std::array<double, 7> &grad, double qStepSize, double tStepSize,
              std::array<double, 7> &oldGrad, std::array<double, 7> &quatTrans,
              std::array<double, 7> &oldQuatTrans, unsigned int iter,
              std::array<double, 7> &step) {
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

double constrainStep(double maxStep, double *step, bool checkSize) {
  double mStep = std::max({fabs(step[0]), fabs(step[1]), fabs(step[2])});
  if (mStep > maxStep) {
    double scaleFactor = maxStep / mStep;
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
  if (checkSize) {
    double quatSquared =
        step[0] * step[0] + step[1] * step[1] + step[2] * step[2];
    if (quatSquared > 1.0) {
      double scaleFactor = 1.0 / (2.0 * quatSquared);
      step[0] *= scaleFactor;
      step[1] *= scaleFactor;
      step[2] *= scaleFactor;
    }
  }
  return mStep;
}

std::array<double, 7> combineQuatTrans(const std::array<double, 7> &q1,
                                       const std::array<double, 7> &q2) {
  std::array<double, 7> res;
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

double oneReduceStep(double grad, double oldGrad, double quatTrans,
                     double oldQuatTrans, double stepSize, double step) {
  if (std::signbit(grad) != std::signbit(oldGrad)) {
    step = (((quatTrans * fabs(oldGrad)) + (oldQuatTrans * fabs(grad))) /
            (fabs(oldGrad + fabs(grad)))) -
           quatTrans;
    double newStep = stepSize * grad;
    if (fabs(step) > fabs(newStep)) {
      step *= fabs(newStep / step);
    }
  } else if (fabs(grad) <= 1.0) {
    step = stepSize * grad;
  } else if (fabs(grad) > fabs(oldGrad)) {
    // Going wrong way relative to other components?
    step += stepSize * grad;
  } else {
    double delta = grad * (step / (oldGrad - grad));
    if (fabs(delta) > fabs(step * 0.1) && fabs(delta) > 0.001) {
      delta *= 0.0005 / fabs(delta);
    }
    step += delta;
  }
  return step;
}

void reduceStep(std::array<double, 7> &grad, std::array<double, 7> &oldGrad,
                std::array<double, 7> quatTrans,
                std::array<double, 7> &oldQuatTrans, unsigned int lineIter,
                std::array<double, 7> &step, double &qStepSize,
                double &tStepSize) {
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
bool SingleConformerAlignment::optimise(unsigned int maxIters) {
  const double maxQuaternionStep = 0.075;   // Maximum step size for quaternion
  const double maxTranslationStep = 0.500;  // Maximum step size for translation
  const double minQuaternionStep =
      0.0002;  // Convergence criteria for quaternion
  const double minTranslationStep =
      0.0020;  // Convergence criteria for translation

  std::array<double, 7> grad;
  std::array<double, 7> oldGrad{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  std::array<double, 7> oldQuatTrans;
  std::array<double, 7> step;
  double shapeOvlpVol, colorOvlpVol, comboScore = 0.0;
  bool finished = false;
  for (unsigned iter = 0; iter < maxIters; iter++) {
    calcVolumeAndGradients(d_quatTrans, shapeOvlpVol, colorOvlpVol, grad);

    // Note that the combo score will have a zero color score so will be half
    // the shape score unless we're optimising on with color gradients.
    auto scores = calcScores(shapeOvlpVol, colorOvlpVol);
    comboScore = scores[0];
    calcStep(grad, d_qStepSize, d_tStepSize, oldGrad, d_quatTrans, oldQuatTrans,
             iter, step);

    // In case we have to backtrack
    double oldComboScore = comboScore;
    oldQuatTrans = d_quatTrans;
    oldGrad = grad;

    // What the PubChem code calls "Line search (sort of) loop"
    bool converged = false;
    for (unsigned int lineIter = 0; !converged; lineIter++) {
      // Check that the absolute max step size does not go beyond some
      // reasonable size
      double mqStep = constrainStep(maxQuaternionStep, step.data() + 1, true);
      double mtStep = constrainStep(maxTranslationStep, step.data() + 4, false);
      if (mqStep <= minQuaternionStep && mtStep <= minTranslationStep) {
        converged = true;
        comboScore = 0.0;  // Make sure we return to the old one.
        break;
      }
      // Calculate the 0th component of the quaternion.  Obviously it
      // relies on the other 3 components being small
      double quatSquared =
          step[1] * step[1] + step[2] * step[2] + step[3] * step[3];
      step[0] = sqrt(1.0 - quatSquared);
      // Update the quaternion with the step, multiplying them.
      auto newQuatTrans = combineQuatTrans(d_quatTrans, step);
      calcVolumeAndGradients(newQuatTrans, shapeOvlpVol, colorOvlpVol, grad);
      auto scores = calcScores(shapeOvlpVol, colorOvlpVol);
      comboScore = scores[0];
      // if we made a good step, keep the quaternion and we're done
      if (comboScore > oldComboScore) {
        d_quatTrans = newQuatTrans;
        break;
      }
      if (lineIter > 2) {
        converged = true;
        d_quatTrans = newQuatTrans;
        break;
      }
      // It got worse, so reduce the step.
      reduceStep(grad, oldGrad, newQuatTrans, oldQuatTrans, lineIter, step,
                 d_qStepSize, d_tStepSize);
      d_quatTrans = oldQuatTrans;
    }  // End of line search
    // Did it converge?
    if (converged ||
        d_shapeConvergenceCriterion > (comboScore - oldComboScore)) {
      if (oldComboScore > comboScore) {
        // The previous step was better, so keep it.
        comboScore = oldComboScore;
        d_quatTrans = oldQuatTrans;
      }
      finished = true;
      break;
    }
  }
  return finished;
}
}  // namespace GaussianShape
}  // namespace RDKit
