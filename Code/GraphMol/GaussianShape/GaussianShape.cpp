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
// This is an implementation of the Gaussian overlap molecular overlay
// method of Grant, Pickup and Gallardo.
// J. Comp. Chem., 17, 1653-1666 (1996)
// https://doi.org/10.1002/(SICI)1096-987X(19961115)17:14%3C1653::AID-JCC7%3E3.0.CO;2-K
// It uses implementation ideas and some code from the PubChem implementation
// https://github.com/ncbi/pubchem-align3d/blob/main/shape_neighbor.cpp.

#include <cmath>
#include <numbers>

#include <Geometry/Transform3D.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/GaussianShape/GaussianShape.h>

#include "GraphMol/SmilesParse/SmilesWrite.h"

#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>

namespace RDKit {
namespace GaussianShape {

namespace {
// Compute final overlay transform, which applies fitShape's
// initial canonical transformation, followed by the overlay transform and
// finally the inverse of refShape's initial canonical transformation.
RDGeom::Transform3D computeFinalTransform(
    const std::array<double, 3> &inRefTrans,
    const std::array<double, 9> &inRefRot,
    const std::array<double, 3> &inFitTrans,
    const std::array<double, 9> &inFitRot, const RDGeom::Transform3D &ovXform) {
  // Move to fitShape's initial centroid and principal axes
  RDGeom::Transform3D transform0;
  transform0.SetTranslation(
      RDGeom::Point3D{inFitTrans[0], inFitTrans[1], inFitTrans[2]});

  RDGeom::Transform3D transform1;
  transform1.setValUnchecked(0, 0, inFitRot[0]);
  transform1.setValUnchecked(0, 1, inFitRot[1]);
  transform1.setValUnchecked(0, 2, inFitRot[2]);
  transform1.setValUnchecked(1, 0, inFitRot[3]);
  transform1.setValUnchecked(1, 1, inFitRot[4]);
  transform1.setValUnchecked(1, 2, inFitRot[5]);
  transform1.setValUnchecked(2, 0, inFitRot[6]);
  transform1.setValUnchecked(2, 1, inFitRot[7]);
  transform1.setValUnchecked(2, 2, inFitRot[8]);

  RDGeom::Transform3D toRefRefFrame;
  // Rotate by the inverse of the ref shape's canonical rotation and
  // translate by the negative of its canonical translation.
  toRefRefFrame.setValUnchecked(0, 0, inRefRot[0]);
  toRefRefFrame.setValUnchecked(0, 1, inRefRot[3]);
  toRefRefFrame.setValUnchecked(0, 2, inRefRot[6]);
  toRefRefFrame.setValUnchecked(0, 3, -inRefTrans[0]);
  toRefRefFrame.setValUnchecked(1, 0, inRefRot[1]);
  toRefRefFrame.setValUnchecked(1, 1, inRefRot[4]);
  toRefRefFrame.setValUnchecked(1, 2, inRefRot[7]);
  toRefRefFrame.setValUnchecked(1, 3, -inRefTrans[1]);
  toRefRefFrame.setValUnchecked(2, 0, inRefRot[2]);
  toRefRefFrame.setValUnchecked(2, 1, inRefRot[5]);
  toRefRefFrame.setValUnchecked(2, 2, inRefRot[8]);
  toRefRefFrame.setValUnchecked(2, 3, -inRefTrans[2]);

  auto finalTransform = toRefRefFrame * ovXform * transform1 * transform0;
  return finalTransform;
}

// Return the original transformation quaternion for the given index.
// Different optimisation modes have different numbers of starting
// orientations to try. In order these are no transformation, rotate 180
// degrees about each axis and rotate +/- 45 degrees about 2 axes at a time.
std::array<double, 4> getInitialRotationPlain(
    int index, const ShapeInput &refShape, const ShapeInput &fitShape,
    const RDGeom::Point3D &refDisp, const ShapeOverlayOptions &overlayOpts,
    double &score) {
  static const double sinpi_4 = std::sin(std::numbers::pi / 4.0);
  const static std::vector<std::array<double, 4>> quats{
      {1.0, 0.0, 0.0, 0.0},          {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},          {0.0, 0.0, 0.0, 1.0},
      {sinpi_4, -sinpi_4, 0.0, 0.0}, {sinpi_4, sinpi_4, 0.0, 0.0},
      {0.0, 0.0, -sinpi_4, sinpi_4}, {0.0, 0.0, sinpi_4, sinpi_4},
      {sinpi_4, 0.0, 0.0, -sinpi_4}, {0.0, sinpi_4, sinpi_4, 0.0},
      {sinpi_4, 0.0, 0.0, sinpi_4},  {0.0, -sinpi_4, sinpi_4, 0.0},
      {sinpi_4, 0.0, sinpi_4, 0.0},  {0.0, sinpi_4, 0.0, sinpi_4},
      {0.0, -sinpi_4, 0.0, sinpi_4}, {sinpi_4, 0.0, -sinpi_4, 0.0}};
  const bool useColor = overlayOpts.optimMode != OptimMode::SHAPE_ONLY;
  const std::array<double, 7> quatTrans{
      quats[index][0], quats[index][1], quats[index][2], quats[index][3],
      refDisp[0],      refDisp[1],      refDisp[2]};
  SingleConformerAlignment sca(
      refShape.getCoords(), refShape.getAlphas(), refShape.getTypes().data(),
      refShape.getCarbonRadii(), refShape.getNumAtoms(),
      refShape.getNumFeatures(), refShape.getShapeVolume(),
      refShape.getColorVolume(), fitShape.getCoords(), fitShape.getAlphas(),
      fitShape.getTypes().data(), fitShape.getCarbonRadii(),
      fitShape.getNumAtoms(), fitShape.getNumFeatures(),
      fitShape.getShapeVolume(), fitShape.getColorVolume(), quatTrans,
      overlayOpts.optimMode, overlayOpts.simAlpha, overlayOpts.simBeta,
      overlayOpts.optParam, overlayOpts.useDistCutoff, overlayOpts.distCutoff,
      overlayOpts.shapeConvergenceCriterion, overlayOpts.nSteps);
  const auto scores = sca.calcScores(useColor);
  score = scores[0];
  return quats[index];
}

// Return the initial transformation matrix in the manner of the PubChem
// overlay code.  Rotate 180 degrees about each axis, and then
// add +/ ~25 degrees from that.  It is not revealed where that
// angle comes from.
std::array<double, 4> getInitialRotationWiggle(
    int index, const ShapeInput &refShape, const ShapeInput &fitShape,
    const RDGeom::Point3D &refDisp, const ShapeOverlayOptions &overlayOpts,
    double &score) {
  const static double qrot1 = 0.977659114061,
                      qrot = 0.210196709523;  // 0.215 (un-normalized)
  const static std::vector<std::array<double, 4>> quats{
      {1.0, 0.0, 0.0, 0.0},  //  0   X,  Y,  Z
      {qrot1, qrot, 0.0, 0.0}, {qrot1, -qrot, 0.0, 0.0},
      {qrot1, 0.0, qrot, 0.0}, {qrot1, 0.0, -qrot, 0.0},
      {qrot1, 0.0, 0.0, qrot}, {qrot1, 0.0, 0.0, -qrot},
      {0.0, 1.0, 0.0, 0.0},  //  1   X, -Y, -Z
      {qrot, qrot1, 0.0, 0.0}, {qrot, -qrot1, 0.0, 0.0},
      {0.0, qrot1, qrot, 0.0}, {0.0, qrot1, -qrot, 0.0},
      {0.0, qrot1, 0.0, qrot}, {0.0, qrot1, 0.0, -qrot},
      {0.0, 0.0, 0.0, 1.0},  //  2  -X, -Y,  Z
      {qrot, 0.0, 0.0, qrot1}, {qrot, 0.0, 0.0, -qrot1},
      {0.0, qrot, 0.0, qrot1}, {0.0, -qrot, 0.0, qrot1},
      {0.0, 0.0, qrot, qrot1}, {0.0, 0.0, -qrot, qrot1},
      {0.0, 0.0, 1.0, 0.0},  //  3  -X,  Y, -Z
      {qrot, 0.0, qrot1, 0.0}, {qrot, 0.0, -qrot1, 0.0},
      {0.0, qrot, qrot1, 0.0}, {0.0, -qrot, qrot1, 0.0},
      {0.0, 0.0, qrot1, qrot}, {0.0, 0.0, qrot1, -qrot}};
  unsigned int start_quat = index * 7;
  unsigned int bestQuat = 0;
  double bestScore = 0.0;
  bool useColor = overlayOpts.optimMode != OptimMode::SHAPE_ONLY;
  std::array<double, 7> tmpQuatTrans{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  SingleConformerAlignment sca(
      refShape.getCoords(), refShape.getAlphas(), refShape.getTypes().data(),
      refShape.getCarbonRadii(), refShape.getNumAtoms(),
      refShape.getNumFeatures(), refShape.getShapeVolume(),
      refShape.getColorVolume(), fitShape.getCoords(), fitShape.getAlphas(),
      fitShape.getTypes().data(), fitShape.getCarbonRadii(),
      fitShape.getNumAtoms(), fitShape.getNumFeatures(),
      fitShape.getShapeVolume(), fitShape.getColorVolume(), tmpQuatTrans,
      overlayOpts.optimMode, overlayOpts.simAlpha, overlayOpts.simBeta,
      overlayOpts.optParam, overlayOpts.useDistCutoff, overlayOpts.distCutoff,
      overlayOpts.shapeConvergenceCriterion, overlayOpts.nSteps);

  for (unsigned int i = start_quat; i < start_quat + 7; ++i) {
    std::array<double, 7> quatTrans{quats[i][0], quats[i][1], quats[i][2],
                                    quats[i][3], refDisp[0],  refDisp[1],
                                    refDisp[2]};
    sca.setQuatTrans(quatTrans);
    auto scores = sca.calcScores(useColor);
    if (scores[0] > bestScore) {
      bestScore = scores[0];
      bestQuat = i;
    }
  }
  score = bestScore;
  return quats[bestQuat];
}

// Return the translation that puts the extreme of refShape at the
// extreme of the fitShape along the appropriate axis.
RDGeom::Point3D getInitialTranslation(int index, ShapeInput &refShape,
                                      ShapeInput fitShape) {
  auto getDisp = [](ShapeInput &shape, size_t i) -> RDGeom::Point3D {
    const double *coord =
        shape.getCoords().data() + shape.calcExtremes()[i] * 3;
    return RDGeom::Point3D(coord[0], coord[1], coord[2]);
  };
  RDGeom::Point3D disp;
  RDGeom::Point3D refDisp, fitDisp;
  switch (index) {
    case 1:
      refDisp = getDisp(refShape, 0);
      fitDisp = getDisp(fitShape, 0);
      disp = fitDisp - refDisp;
      break;
    case 2:
      refDisp = getDisp(refShape, 1);
      fitDisp = getDisp(fitShape, 1);
      disp = fitDisp - refDisp;
      break;
    case 3:
      refDisp = getDisp(refShape, 2);
      fitDisp = getDisp(fitShape, 2);
      disp = fitDisp - refDisp;
      break;
    case 4:
      refDisp = getDisp(refShape, 3);
      fitDisp = getDisp(fitShape, 3);
      disp = fitDisp - refDisp;
      break;
    case 5:
      refDisp = getDisp(refShape, 4);
      fitDisp = getDisp(fitShape, 4);
      disp = fitDisp - refDisp;
      break;
    case 6:
      refDisp = getDisp(refShape, 5);
      fitDisp = getDisp(fitShape, 5);
      disp = fitDisp - refDisp;
      break;
    default:
      break;
  }
  return disp;
}

// This is how the PubChem code decides between ROTATE_180_WIGGLE and
// ROTATE_45.  I have no clue.
unsigned int calculateQrat(const std::array<double, 3> &eigenValues) {
  double double_ev_oe[3]{eigenValues[1] + eigenValues[2] - eigenValues[0],
                         eigenValues[0] + eigenValues[2] - eigenValues[1],
                         eigenValues[0] + eigenValues[1] - eigenValues[2]};
  std::sort(double_ev_oe, double_ev_oe + 3, std::greater<double>());

  constexpr static double qrat_threshold = 0.7225;  // 0.85*0.85;
  unsigned int qrat = 1000;

  if (double_ev_oe[1] > 0) {
    unsigned int u_rqyx, u_rqzy;
    if (qrat_threshold < double_ev_oe[1] / double_ev_oe[0]) {
      u_rqyx = 1;
    } else {
      u_rqyx = 0;
    }
    if (qrat_threshold < double_ev_oe[2] / double_ev_oe[1]) {
      u_rqzy = 1;
    } else {
      u_rqzy = 0;
    }

    qrat = u_rqyx + u_rqzy;
  }
  return qrat;
}

StartMode decideStartModeFromEigenValues(const ShapeInput &refShape,
                                         const ShapeInput &fitShape) {
  // The PubChem code uses the moments of inertia for this, rather than the
  // canonical transformation.
  const auto rqratwf = calculateQrat(refShape.calcMomentsOfInertia(true));
  const auto fqratwf = calculateQrat(fitShape.calcMomentsOfInertia(true));
  StartMode startModeWF{StartMode::ROTATE_180_WIGGLE};
  if (rqratwf > 0 || fqratwf > 0) {
    startModeWF = StartMode::ROTATE_45;
  }
  return startModeWF;
}

std::array<double, 3> alignShape(ShapeInput &refShape, ShapeInput &fitShape,
                                 RDGeom::Transform3D &bestXform,
                                 const ShapeOverlayOptions &overlayOpts) {
  unsigned int finalRotIndex = 1;
  auto startMode = overlayOpts.startMode;
  if (startMode == StartMode::A_LA_PUBCHEM) {
    startMode = decideStartModeFromEigenValues(refShape, fitShape);
  }

  switch (startMode) {
    case StartMode::ROTATE_0:
    case StartMode::ROTATE_0_FRAGMENT:
      break;
    case StartMode::ROTATE_180:
    case StartMode::ROTATE_180_FRAGMENT:
    case StartMode::ROTATE_180_WIGGLE:
      finalRotIndex = 4;
      break;
    case StartMode::ROTATE_45:
    case StartMode::ROTATE_45_FRAGMENT:
      finalRotIndex = 16;
      break;
    default:
      break;
  }
  unsigned int finalTransIndex = 1;
  if (startMode == StartMode::ROTATE_0_FRAGMENT ||
      startMode == StartMode::ROTATE_45_FRAGMENT ||
      startMode == StartMode::ROTATE_180_FRAGMENT) {
    finalTransIndex = 7;
  }

  std::array<double, 3> bestScore;
  double bestTotal = -1.0;

  // Get together the start transformations.
  std::vector<std::unique_ptr<SingleConformerAlignment>> aligners;
  std::vector<std::pair<double, unsigned int>> bestScoreForStart;
  bestScoreForStart.reserve(finalTransIndex * finalRotIndex);
  unsigned int k = 0;
  for (unsigned int j = 0; j < finalTransIndex; j++) {
    auto refDisp = getInitialTranslation(j, refShape, fitShape);
    std::array<double, 4> quat;
    for (unsigned int i = 0; i < finalRotIndex; i++, k++) {
      double score = 0.0;
      if (startMode == StartMode::ROTATE_180_WIGGLE) {
        quat = getInitialRotationWiggle(i, refShape, fitShape, refDisp,
                                        overlayOpts, score);
      } else {
        quat = getInitialRotationPlain(i, refShape, fitShape, refDisp,
                                       overlayOpts, score);
      }
      std::array<double, 7> initQuat{quat[0],   quat[1],   quat[2],  quat[3],
                                     refDisp.x, refDisp.y, refDisp.z};
      aligners.emplace_back(std::make_unique<SingleConformerAlignment>(
          refShape.getCoords(), refShape.getAlphas(),
          refShape.getTypes().data(), refShape.getCarbonRadii(),
          refShape.getNumAtoms(), refShape.getNumFeatures(),
          refShape.getShapeVolume(), refShape.getColorVolume(),
          fitShape.getCoords(), fitShape.getAlphas(),
          fitShape.getTypes().data(), fitShape.getCarbonRadii(),
          fitShape.getNumAtoms(), fitShape.getNumFeatures(),
          fitShape.getShapeVolume(), fitShape.getColorVolume(), initQuat,
          overlayOpts.optimMode, overlayOpts.simAlpha, overlayOpts.simBeta,
          overlayOpts.optParam, overlayOpts.useDistCutoff,
          overlayOpts.distCutoff, overlayOpts.shapeConvergenceCriterion,
          overlayOpts.nSteps));
      bestScoreForStart.push_back({score, k});
    }
  }

  // Do it in 2 cycles, a quick optimisation first, followed by an additional
  // longer one for those that look like they're going to win.
  for (unsigned int cycle = 0; cycle < 2; cycle++) {
    std::ranges::sort(bestScoreForStart,
                      [](const auto &p1, const auto &p2) -> bool {
                        return p1.first > p2.first;
                      });
    std::vector<std::pair<double, unsigned int>> nextBestScoreForStart;
    nextBestScoreForStart.reserve(finalTransIndex * finalRotIndex);
    for (const auto &[bssf, m] : bestScoreForStart) {
      if (cycle == 1) {
        if (bssf < 0.7 * bestScore[0]) {
          continue;
        }
      }
      std::array<double, 20> outScores;
      aligners[m]->doOverlay(outScores, cycle);
      nextBestScoreForStart.emplace_back(outScores[0], m);
      if (outScores[0] > bestTotal) {
        bestTotal = outScores[0];
        bestScore =
            std::array<double, 3>{outScores[0], outScores[1], outScores[2]};
        aligners[m]->getFinalQuatTrans(bestXform);
      }
    }
    bestScoreForStart = nextBestScoreForStart;
  }
  return bestScore;
}
}  // namespace

std::array<double, 3> AlignShape(const ShapeInput &refShape,
                                 ShapeInput &fitShape,
                                 RDGeom::Transform3D *xform,
                                 const ShapeOverlayOptions &overlayOpts) {
  // The shapes aren't necessarily normalized (it's not done on creation, for
  // example) but they might need to be.
  auto workingRefShape = std::make_unique<ShapeInput>(refShape);
  auto workingFitShape = std::make_unique<ShapeInput>(fitShape);
  const auto inRefTrans = workingRefShape->calcCanonicalTranslation();
  const auto inRefRot = workingRefShape->calcCanonicalRotation();
  const auto inFitTrans = workingFitShape->calcCanonicalTranslation();
  const auto inFitRot = workingFitShape->calcCanonicalRotation();
  // If we're not normalizing, translate both shapes so that the fit
  // is at the origin, so the rotations work.
  RDGeom::Transform3D moveToOrigin;
  RDGeom::Transform3D moveFromOrigin;
  if (overlayOpts.normalize) {
    if (!workingRefShape->getNormalized()) {
      workingRefShape->normalizeCoords();
    }
    if (!workingFitShape->getNormalized()) {
      workingFitShape->normalizeCoords();
    }
  } else {
    const auto &canonTrans = workingFitShape->calcCanonicalTranslation();
    moveToOrigin.SetTranslation(
        RDGeom::Point3D{canonTrans[0], canonTrans[1], canonTrans[2]});
    moveFromOrigin.SetTranslation(
        RDGeom::Point3D{-canonTrans[0], -canonTrans[1], -canonTrans[2]});
    workingFitShape->transformCoords(moveToOrigin);
    workingRefShape->transformCoords(moveToOrigin);
  }

  RDGeom::Transform3D bestXform;
  const auto scores =
      alignShape(*workingRefShape, *workingFitShape, bestXform, overlayOpts);
  if (!overlayOpts.normalize) {
    // Shove it back again.
    bestXform = moveFromOrigin * bestXform * moveToOrigin;
  } else {
    bestXform = computeFinalTransform(inRefTrans, inRefRot, inFitTrans,
                                      inFitRot, bestXform);
  }
  fitShape.transformCoords(bestXform);
  if (xform) {
    *xform = bestXform;
  }

  return scores;
}

std::array<double, 3> AlignMolecule(const ShapeInput &refShape, ROMol &fit,
                                    const ShapeInputOptions &fitOpts,
                                    RDGeom::Transform3D *xform,
                                    const ShapeOverlayOptions &overlayOpts,
                                    int fitConfId) {
  auto fitShape = ShapeInput(fit, fitConfId, fitOpts, overlayOpts);
  RDGeom::Transform3D tmpXform;
  const auto scores = AlignShape(refShape, fitShape, &tmpXform, overlayOpts);
  MolTransforms::transformConformer(fit.getConformer(fitConfId), tmpXform);
  if (xform) {
    *xform = tmpXform;
  }
  return scores;
}

std::array<double, 3> AlignMolecule(const ROMol &ref, ROMol &fit,
                                    const ShapeInputOptions &refOpts,
                                    const ShapeInputOptions &fitOpts,
                                    RDGeom::Transform3D *xform,
                                    const ShapeOverlayOptions &overlayOpts,
                                    int refConfId, int fitConfId) {
  const auto refShape = ShapeInput(ref, refConfId, refOpts, overlayOpts);
  const auto scores =
      AlignMolecule(refShape, fit, fitOpts, xform, overlayOpts, fitConfId);
  return scores;
}

void AlignMoleculesAllConformers(const ROMol &ref, const ROMol &fit,
                                 int &refConfId, int &fitConfId,
                                 std::vector<std::vector<double>> &combScores,
                                 const ShapeInputOptions &refOpts,
                                 const ShapeInputOptions &fitOpts,
                                 const ShapeOverlayOptions &overlayOpts,
                                 RDGeom::Transform3D *xform) {
  // Pruning the shapes wastes time and obviously removes the correspondence
  // between conformers and shapes.
  auto refOptsCp = refOpts;
  refOptsCp.shapePruneThreshold = -1;
  refOptsCp.sortShapes = false;
  auto fitOptsCp = fitOpts;
  fitOptsCp.shapePruneThreshold = -1;
  fitOptsCp.sortShapes = false;
  auto refShape = ShapeInput(ref, -1, refOptsCp, overlayOpts);
  auto fitShape = ShapeInput(fit, -1, fitOptsCp, overlayOpts);
  combScores = std::vector<std::vector<double>>(
      refShape.getNumShapes(), std::vector<double>(fitShape.getNumShapes()));
  double bestScore = -1.0;
  for (unsigned int i = 0; i < refShape.getNumShapes(); i++) {
    refShape.setActiveShape(i);
    for (unsigned int j = 0; j < fitShape.getNumShapes(); j++) {
      fitShape.setActiveShape(j);
      RDGeom::Transform3D thisXform;
      auto scores = AlignShape(refShape, fitShape, &thisXform, overlayOpts);
      combScores[i][j] = scores[0];
      if (scores[0] > bestScore) {
        bestScore = scores[0];
        refConfId = i;
        fitConfId = j;
        if (xform) {
          *xform = thisXform;
        }
      }
    }
  }
}

std::array<double, 3> ScoreShape(const ShapeInput &refShape,
                                 const ShapeInput &fitShape,
                                 const ShapeOverlayOptions &overlayOpts,
                                 std::array<double, 2> *overlapVols) {
  auto refWorking = refShape.getCoords();
  auto fitWorking = fitShape.getCoords();
  std::array<double, 7> quatTrans{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  SingleConformerAlignment sca(
      refShape.getCoords(), refShape.getAlphas(), refShape.getTypes().data(),
      refShape.getCarbonRadii(), refShape.getNumAtoms(),
      refShape.getNumFeatures(), refShape.getShapeVolume(),
      refShape.getColorVolume(), fitShape.getCoords(), fitShape.getAlphas(),
      fitShape.getTypes().data(), fitShape.getCarbonRadii(),
      fitShape.getNumAtoms(), fitShape.getNumFeatures(),
      fitShape.getShapeVolume(), fitShape.getColorVolume(), quatTrans,
      overlayOpts.optimMode, overlayOpts.simAlpha, overlayOpts.simBeta,
      overlayOpts.optParam, overlayOpts.useDistCutoff, overlayOpts.distCutoff,
      overlayOpts.shapeConvergenceCriterion, overlayOpts.nSteps);
  const bool includeColor = overlayOpts.optimMode != OptimMode::SHAPE_ONLY;
  const auto scores = sca.calcScores(refShape.getCoords().data(),
                                     fitShape.getCoords().data(), includeColor);
  if (overlapVols) {
    (*overlapVols)[0] = scores[3];
    (*overlapVols)[1] = scores[4];
  }
  return std::array{scores[0], scores[1], scores[2]};
}

std::array<double, 3> ScoreMolecule(const ShapeInput &refShape,
                                    const ROMol &fit,
                                    const ShapeInputOptions &fitOpts,
                                    const ShapeOverlayOptions &overlayOpts,
                                    int fitConfId,
                                    std::array<double, 2> *overlapVols) {
  const auto fitShape = ShapeInput(fit, fitConfId, fitOpts, overlayOpts);
  return ScoreShape(refShape, fitShape, overlayOpts, overlapVols);
}

std::array<double, 3> ScoreMolecule(const ROMol &ref, const ROMol &fit,
                                    const ShapeInputOptions &refOpts,
                                    const ShapeInputOptions &fitOpts,
                                    const ShapeOverlayOptions &overlayOpts,
                                    int refConfId, int fitConfId,
                                    std::array<double, 2> *overlapVols) {
  ShapeOverlayOptions tmpOpts = overlayOpts;
  tmpOpts.normalize = false;
  tmpOpts.startMode = StartMode::ROTATE_0;
  ShapeInputOptions tmpRefOpts = refOpts;
  auto refShape = ShapeInput(ref, refConfId, refOpts, tmpOpts);

  ShapeInputOptions tmpFitOpts = fitOpts;
  const auto fitShape = ShapeInput(fit, fitConfId, fitOpts, tmpOpts);

  return ScoreShape(refShape, fitShape, tmpOpts, overlapVols);
}
}  // namespace GaussianShape
}  // namespace RDKit