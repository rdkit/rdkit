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

#include <Geometry/Transform3D.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

namespace RDKit {
namespace GaussianShape {

namespace {
// Compute final overlay transform, which applies fitShape's
// canonical transformation, followed by the overlay transform and
// finally the inverse of refShape's canonical transformation.
RDGeom::Transform3D computeFinalTransform(const ShapeInput &refShape,
                                          const ShapeInput &fitShape,
                                          RDGeom::Transform3D &ovXform) {
  // Move fitConf to fitShape's initial centroid and principal axes
  RDGeom::Transform3D transform0;
  transform0.SetTranslation(
      RDGeom::Point3D{fitShape.getCanonicalTranslation()[0],
                      fitShape.getCanonicalTranslation()[1],
                      fitShape.getCanonicalTranslation()[2]});

  RDGeom::Transform3D transform1;
  transform1.setValUnchecked(0, 0, fitShape.getCanonicalRotation()[0]);
  transform1.setValUnchecked(0, 1, fitShape.getCanonicalRotation()[1]);
  transform1.setValUnchecked(0, 2, fitShape.getCanonicalRotation()[2]);
  transform1.setValUnchecked(1, 0, fitShape.getCanonicalRotation()[3]);
  transform1.setValUnchecked(1, 1, fitShape.getCanonicalRotation()[4]);
  transform1.setValUnchecked(1, 2, fitShape.getCanonicalRotation()[5]);
  transform1.setValUnchecked(2, 0, fitShape.getCanonicalRotation()[6]);
  transform1.setValUnchecked(2, 1, fitShape.getCanonicalRotation()[7]);
  transform1.setValUnchecked(2, 2, fitShape.getCanonicalRotation()[8]);

  RDGeom::Transform3D toRefRefFrame;
  // Rotate by the inverse of the ref shape's canonical rotation and
  // translate by the negative of its canonical translation.
  toRefRefFrame.setValUnchecked(0, 0, refShape.getCanonicalRotation()[0]);
  toRefRefFrame.setValUnchecked(0, 1, refShape.getCanonicalRotation()[3]);
  toRefRefFrame.setValUnchecked(0, 2, refShape.getCanonicalRotation()[6]);
  toRefRefFrame.setValUnchecked(0, 3, -refShape.getCanonicalTranslation()[0]);
  toRefRefFrame.setValUnchecked(1, 0, refShape.getCanonicalRotation()[1]);
  toRefRefFrame.setValUnchecked(1, 1, refShape.getCanonicalRotation()[4]);
  toRefRefFrame.setValUnchecked(1, 2, refShape.getCanonicalRotation()[7]);
  toRefRefFrame.setValUnchecked(1, 3, -refShape.getCanonicalTranslation()[1]);
  toRefRefFrame.setValUnchecked(2, 0, refShape.getCanonicalRotation()[2]);
  toRefRefFrame.setValUnchecked(2, 1, refShape.getCanonicalRotation()[5]);
  toRefRefFrame.setValUnchecked(2, 2, refShape.getCanonicalRotation()[8]);
  toRefRefFrame.setValUnchecked(2, 3, -refShape.getCanonicalTranslation()[2]);

  auto finalTransform = toRefRefFrame * ovXform * transform1 * transform0;
  return finalTransform;
}

// Return the original transformation matrix for the given index.
// Different optimisation modes have different numbers of starting
// orientations to try. In order these are no transformation, rotate 180
// degrees about each axis and rotate +/- 45 degrees about 2 axes at a time.
std::array<double, 4> getInitialRotationPlain(
    int index, const ShapeInput &refShape, const ShapeInput &fitShape,
    const RDGeom::Point3D &refDisp, const ShapeOverlayOptions &overlayOpts,
    std::vector<double> &workingRef, std::vector<double> &workingFit,
    RDGeom::Transform3D &initXform, double &score) {
  static const double sinpi_4 = std::sin(std::atan(1.0));
  const static std::vector<std::array<double, 4>> quats{
      {1.0, 0.0, 0.0, 0.0},          {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},          {0.0, 0.0, 0.0, 1.0},
      {sinpi_4, -sinpi_4, 0.0, 0.0}, {sinpi_4, sinpi_4, 0.0, 0.0},
      {0.0, 0.0, -sinpi_4, sinpi_4}, {0.0, 0.0, sinpi_4, sinpi_4},
      {sinpi_4, 0.0, 0.0, -sinpi_4}, {0.0, sinpi_4, sinpi_4, 0.0},
      {sinpi_4, 0.0, 0.0, sinpi_4},  {0.0, -sinpi_4, sinpi_4, 0.0},
      {sinpi_4, 0.0, sinpi_4, 0.0},  {0.0, sinpi_4, 0.0, sinpi_4},
      {0.0, -sinpi_4, 0.0, sinpi_4}, {sinpi_4, 0.0, -sinpi_4, 0.0}};
  auto refCoords = refShape.getCoords();
  translateShape(refCoords, refDisp);
  bool useColor = overlayOpts.optimMode != OptimMode::SHAPE_ONLY;
  auto fitCoords = fitShape.getCoords();
  initXform.setToIdentity();
  initXform.SetRotationFromQuaternion(quats[index].data());
  applyTransformToShape(fitCoords, initXform);
  SingleConformerAlignment sca(
      refCoords.data(), workingRef.data(), refShape.getTypes().data(),
      &refShape.getCarbonRadii(), refShape.getNumAtoms(),
      refShape.getNumFeatures(), refShape.getShapeVolume(),
      refShape.getColorVolume(), fitCoords.data(), workingFit.data(),
      fitShape.getTypes().data(), &fitShape.getCarbonRadii(),
      fitShape.getNumAtoms(), fitShape.getNumFeatures(),
      fitShape.getShapeVolume(), fitShape.getColorVolume(),
      overlayOpts.optimMode, overlayOpts.optParam, overlayOpts.useDistCutoff,
      overlayOpts.distCutoff, overlayOpts.nSteps);
  auto scores = sca.calcScores(refCoords.data(), fitCoords.data(), useColor);
  std::cout << "     index " << index << " : " << scores[0] << ", " << scores[1]
            << ", " << scores[2] << ", " << scores[3] << ", " << scores[4]
            << std::endl;
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
    std::vector<double> &workingRef, std::vector<double> &workingFit,
    RDGeom::Transform3D &initXform, double &score) {
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
  auto refCoords = refShape.getCoords();
  translateShape(refCoords, refDisp);
  unsigned int start_quat = index * 7;
  unsigned int bestQuat = 0;
  double bestScore = 0.0;
  bool useColor = overlayOpts.optimMode != OptimMode::SHAPE_ONLY;

  for (unsigned int i = start_quat; i < start_quat + 7; ++i) {
    auto fitCoords = fitShape.getCoords();
    RDGeom::Transform3D xForm;
    xForm.SetRotationFromQuaternion(quats[i].data());
    applyTransformToShape(fitCoords, xForm);
    SingleConformerAlignment sca(
        refCoords.data(), workingRef.data(), refShape.getTypes().data(),
        &refShape.getCarbonRadii(), refShape.getNumAtoms(),
        refShape.getNumFeatures(), refShape.getShapeVolume(),
        refShape.getColorVolume(), fitCoords.data(), workingFit.data(),
        fitShape.getTypes().data(), &fitShape.getCarbonRadii(),
        fitShape.getNumAtoms(), fitShape.getNumFeatures(),
        fitShape.getShapeVolume(), fitShape.getColorVolume(),
        overlayOpts.optimMode, overlayOpts.optParam, overlayOpts.useDistCutoff,
        overlayOpts.distCutoff, overlayOpts.nSteps);
    auto scores = sca.calcScores(refCoords.data(), fitCoords.data(), useColor);
    if (scores[0] > bestScore) {
      bestScore = scores[0];
      bestQuat = i;
    }
  }
  initXform.setToIdentity();
  initXform.SetRotationFromQuaternion(quats[bestQuat].data());
  score = bestScore;
  return quats[bestQuat];
}

// Return the translation that puts the extreme of refShape at the
// extreme of the fitShape along the appropriate axis.
RDGeom::Point3D getInitialTranslation(int index, const ShapeInput &refShape,
                                      const ShapeInput fitShape) {
  auto getDisp = [](const ShapeInput &shape, size_t i) -> RDGeom::Point3D {
    const double *coord = shape.getCoords().data() + shape.getExtremes()[i] * 4;
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
  const static double qrat_threshold = 0.7225;  // 0.85*0.85;
  unsigned int qrat = 1000;
  unsigned int u_rqyx, u_rqzy;
  double double_ev_oe[3]{eigenValues[1] + eigenValues[2] - eigenValues[0],
                         eigenValues[0] + eigenValues[2] - eigenValues[1],
                         eigenValues[0] + eigenValues[1] - eigenValues[2]};
  if (double_ev_oe[1] > 0) {
    if (qrat_threshold < (double_ev_oe[1] / double_ev_oe[0])) {
      u_rqyx = 1;
    } else {
      u_rqyx = 0;
    }
    if (qrat_threshold < (double_ev_oe[2] / double_ev_oe[1])) {
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
  auto rqrat = calculateQrat(refShape.getEigenValues());
  auto fqrat = calculateQrat(fitShape.getEigenValues());
  if (rqrat > 0 || fqrat > 0) {
    return StartMode::ROTATE_45;
  }
  return StartMode::ROTATE_180_WIGGLE;
}

std::array<double, 3> alignShape(const ShapeInput &refShape,
                                 const ShapeInput &fitShape,
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
  RDGeom::Transform3D initialRot, initialTrans, initialXform;
  // For working values of the coordinates.
  std::vector<double> workingRef(refShape.getCoords());
  std::vector<double> workingFit(fitShape.getCoords());

  // Get together the start transformations.
  std::vector<RDGeom::Point3D> initialTranslations;
  std::vector<RDGeom::Transform3D> initialRotations;
  std::vector<std::array<double, 7>> initialQuats;
  std::vector<std::pair<double, unsigned int>> bestScoreForStart;
  initialTranslations.reserve(finalTransIndex * finalRotIndex);
  initialRotations.reserve(finalTransIndex * finalRotIndex);
  initialQuats.reserve(finalTransIndex * finalRotIndex);
  bestScoreForStart.reserve(finalTransIndex * finalRotIndex);
  unsigned int k = 0;
  for (unsigned int j = 0; j < finalTransIndex; j++) {
    auto refDisp = getInitialTranslation(j, refShape, fitShape);
    std::array<double, 4> quat;
    for (unsigned int i = 0; i < finalRotIndex; i++, k++) {
      double score = 0.0;
      if (startMode == StartMode::ROTATE_180_WIGGLE) {
        quat = getInitialRotationWiggle(i, refShape, fitShape, refDisp,
                                        overlayOpts, workingRef, workingFit,
                                        initialRot, score);
      } else {
        quat =
            getInitialRotationPlain(i, refShape, fitShape, refDisp, overlayOpts,
                                    workingRef, workingFit, initialRot, score);
      }
      initialTranslations.push_back(refDisp);
      initialRotations.push_back(initialRot);
      initialQuats.push_back({quat[0], quat[1], quat[2], quat[3], refDisp.x,
                              refDisp.y, refDisp.z});
      bestScoreForStart.push_back({score, k});
    }
  }

  // Do it in 2 cycles, a quick optimisation first, followed by an additional
  // longer one for those that look like they're going to win.
  for (unsigned int cycle = 0; cycle < 2; cycle++) {
    std::cout << "__--------------------------  " << cycle << std::endl;
    std::ranges::sort(bestScoreForStart,
                      [](const auto &p1, const auto &p2) -> bool {
                        return p1.first > p2.first;
                      });
    std::vector<std::pair<double, unsigned int>> nextBestScoreForStart;
    nextBestScoreForStart.reserve(finalTransIndex * finalRotIndex);
    for (const auto &[bssf, k] : bestScoreForStart) {
      std::cout << "case " << k << " with bssf: " << bssf << std::endl;
      RDGeom::Point3D refDisp{initialQuats[k][4], initialQuats[k][5],
                              initialQuats[k][6]};
      if (cycle == 1) {
        if (bssf < 0.7 * bestScore[0]) {
          std::cout << "skipping because " << bssf << " for " << k << " vs "
                    << 0.7 * bestScore[0] << " from " << bestScore[0]
                    << std::endl;
          continue;
        }
      }
      std::vector<double> startFit(fitShape.getCoords());
      // initialRot.setToIdentity();
      // initialRot.SetRotationFromQuaternion(initialQuats[k].data());
      // if (cycle == 1) {
      //   // Add in where the optimisation started from in cycle 0.
      //   copyTransform(initialRotations[k] * initialRot, initialRot);
      //   refDisp += initialTranslations[k];
      // }
      // applyTransformToShape(startFit, initialRot);
      std::vector<double> startRef(refShape.getCoords());
      // // Move the reference by initialTrans, leaving fit at the origin where
      // // the rotations work properly.
      // translateShape(startRef, refDisp);
      std::array<double, 20> outScores;
      SingleConformerAlignment sca(
          startRef.data(), workingRef.data(), refShape.getTypes().data(),
          &refShape.getCarbonRadii(), refShape.getNumAtoms(),
          refShape.getNumFeatures(), refShape.getShapeVolume(),
          refShape.getColorVolume(), startFit.data(), workingFit.data(),
          fitShape.getTypes().data(), &fitShape.getCarbonRadii(),
          fitShape.getNumAtoms(), fitShape.getNumFeatures(),
          fitShape.getShapeVolume(), fitShape.getColorVolume(),
          overlayOpts.optimMode, overlayOpts.optParam,
          overlayOpts.useDistCutoff, overlayOpts.distCutoff,
          overlayOpts.nSteps);
      sca.doOverlay(outScores, initialQuats[k], cycle);
      std::cout << cycle << " : " << k << " : " << outScores[0] << std::endl;
      nextBestScoreForStart.emplace_back(outScores[0], k);
      if (outScores[0] > bestTotal) {
        bestTotal = outScores[0];
        bestScore =
            std::array<double, 3>{outScores[0], outScores[1], outScores[2]};
        RDGeom::Transform3D tmp;
        tmp.SetRotationFromQuaternion(outScores.data() + 9);
        tmp.SetTranslation(
            RDGeom::Point3D{outScores[13], outScores[14], outScores[15]});
        RDGeom::Transform3D reverseInitialTrans;
        reverseInitialTrans.SetTranslation(-initialTranslations[k]);
        auto tt = reverseInitialTrans * tmp * initialRotations[k];
        copyTransform(tt, bestXform);
        copyTransform(tmp, bestXform);
      }
      initialQuats[k] = std::array<double, 7>{
          outScores[9],  outScores[10], outScores[11], outScores[12],
          outScores[13], outScores[14], outScores[15]};
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
    moveToOrigin.SetTranslation(
        RDGeom::Point3D{workingFitShape->getCanonicalTranslation()[0],
                        workingFitShape->getCanonicalTranslation()[1],
                        workingFitShape->getCanonicalTranslation()[2]});
    moveFromOrigin.SetTranslation(
        RDGeom::Point3D{-fitShape.getCanonicalTranslation()[0],
                        -fitShape.getCanonicalTranslation()[1],
                        -fitShape.getCanonicalTranslation()[2]});
    workingFitShape->transformCoords(moveToOrigin);
    workingRefShape->transformCoords(moveToOrigin);
  }

  RDGeom::Transform3D bestXform;
  auto scores =
      alignShape(*workingRefShape, *workingFitShape, bestXform, overlayOpts);
  if (!overlayOpts.normalize) {
    // Shove it back again.
    auto finalXform = moveFromOrigin * bestXform * moveToOrigin;
    copyTransform(finalXform, bestXform);
  } else {
    auto finalXform =
        computeFinalTransform(*workingRefShape, *workingFitShape, bestXform);
    copyTransform(finalXform, bestXform);
    fitShape.transformCoords(bestXform);
  }
  if (xform) {
    copyTransform(bestXform, *xform);
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
  auto scores = AlignShape(refShape, fitShape, &tmpXform, overlayOpts);
  MolTransforms::transformConformer(fit.getConformer(fitConfId), tmpXform);
  if (xform) {
    copyTransform(tmpXform, *xform);
  }
  return scores;
}

std::array<double, 3> AlignMolecule(const ROMol &ref, ROMol &fit,
                                    const ShapeInputOptions &refOpts,
                                    const ShapeInputOptions &fitOpts,
                                    RDGeom::Transform3D *xform,
                                    const ShapeOverlayOptions &overlayOpts,
                                    int refConfId, int fitConfId) {
  auto refShape = ShapeInput(ref, refConfId, refOpts, overlayOpts);
  RDGeom::Transform3D tmpXform;
  auto scores =
      AlignMolecule(refShape, fit, fitOpts, xform, overlayOpts, fitConfId);
  return scores;
}

std::array<double, 3> ScoreShape(const ShapeInput &refShape,
                                 const ShapeInput &fitShape,
                                 const ShapeOverlayOptions &overlayOpts) {
  auto refWorking = refShape.getCoords();
  auto fitWorking = fitShape.getCoords();
  SingleConformerAlignment sca(
      refShape.getCoords().data(), refWorking.data(),
      refShape.getTypes().data(), &refShape.getCarbonRadii(),
      refShape.getNumAtoms(), refShape.getNumFeatures(),
      refShape.getShapeVolume(), refShape.getColorVolume(),
      fitShape.getCoords().data(), fitWorking.data(),
      fitShape.getTypes().data(), &fitShape.getCarbonRadii(),
      fitShape.getNumAtoms(), fitShape.getNumFeatures(),
      fitShape.getShapeVolume(), fitShape.getColorVolume(),
      overlayOpts.optimMode, overlayOpts.optParam, overlayOpts.useDistCutoff,
      overlayOpts.distCutoff, overlayOpts.nSteps);
  bool includeColor = overlayOpts.optimMode != OptimMode::SHAPE_ONLY;
  auto scores = sca.calcScores(refShape.getCoords().data(),
                               fitShape.getCoords().data(), includeColor);
  return std::array{scores[0], scores[1], scores[2]};
}

std::array<double, 3> ScoreMolecule(const ShapeInput &refShape,
                                    const ROMol &fit,
                                    const ShapeInputOptions &fitOpts,
                                    const ShapeOverlayOptions &overlayOpts,
                                    int fitConfId) {
  auto fitShape = ShapeInput(fit, fitConfId, fitOpts, overlayOpts);
  return ScoreShape(refShape, fitShape, overlayOpts);
}

std::array<double, 3> ScoreMolecule(const ROMol &ref, const ROMol &fit,
                                    const ShapeInputOptions &refOpts,
                                    const ShapeInputOptions &fitOpts,
                                    const ShapeOverlayOptions &overlayOpts,
                                    int refConfId, int fitConfId) {
  ShapeOverlayOptions tmpOpts = overlayOpts;
  tmpOpts.normalize = false;
  tmpOpts.startMode = StartMode::ROTATE_0;
  ShapeInputOptions tmpRefOpts = refOpts;
  auto refShape = ShapeInput(ref, refConfId, refOpts, tmpOpts);

  ShapeInputOptions tmpFitOpts = fitOpts;
  auto fitShape = ShapeInput(fit, fitConfId, fitOpts, tmpOpts);

  return ScoreShape(refShape, fitShape, tmpOpts);
}
}  // namespace GaussianShape
}  // namespace RDKit