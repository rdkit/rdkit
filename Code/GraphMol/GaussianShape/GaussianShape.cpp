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
// degrees about each axis and rotate +/- 90 degrees about each axis.
void getInitialRotation(int index, RDGeom::Transform3D &initXform) {
  static const DTYPE sinpi_2 = std::sin(2.0 * std::atan(1.0));
  const static std::vector<std::array<DTYPE, 4>> quats{
      {1.0, 0.0, 0.0, 0.0},          {0.0, 1.0, 0.0, 0.0},
      {0.0, 0.0, 1.0, 0.0},          {0.0, 0.0, 0.0, 1.0},
      {sinpi_2, -sinpi_2, 0.0, 0.0}, {sinpi_2, sinpi_2, 0.0, 0.0},
      {0.0, 0.0, -sinpi_2, sinpi_2}, {0.0, 0.0, sinpi_2, sinpi_2},
      {sinpi_2, 0.0, 0.0, -sinpi_2}, {0.0, sinpi_2, sinpi_2, 0.0},
      {sinpi_2, 0.0, 0.0, sinpi_2},  {0.0, -sinpi_2, sinpi_2, 0.0},
      {sinpi_2, 0.0, sinpi_2, 0.0},  {0.0, sinpi_2, 0.0, sinpi_2},
      {0.0, -sinpi_2, 0.0, sinpi_2}, {sinpi_2, 0.0, -sinpi_2, 0.0}};
  initXform.setToIdentity();
  initXform.SetRotationFromQuaternion(quats[index].data());
}

// Return the initial transformation matrix in the manner of the PubChem
// overlay code.  Rotate 180 degrees about each axis, and then
// add +/ ~25 degrees from that.  It is not revealed where that
// angle comes from.
void getInitialRotation(int index, const ShapeInput &refShape,
                        const ShapeInput &fitShape,
                        const RDGeom::Point3D &refDisp,
                        const ShapeOverlayOptions &overlayOpts,
                        std::vector<DTYPE> &workingRef,
                        std::vector<DTYPE> &workingFit,
                        RDGeom::Transform3D &initXForm) {
  const static double qrot1 = 0.977659114061,
                      qrot = 0.210196709523;  // 0.215 (un-normalized)
  const static std::vector<std::array<DTYPE, 4>> quats{
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
        refShape.getNumFeatures(), refShape.getSelfOverlapVol(),
        refShape.getSelfOverlapColor(), fitCoords.data(), workingFit.data(),
        fitShape.getTypes().data(), &fitShape.getCarbonRadii(),
        fitShape.getNumAtoms(), fitShape.getNumFeatures(),
        fitShape.getSelfOverlapVol(), fitShape.getSelfOverlapColor(),
        overlayOpts.optimMode, overlayOpts.optParam, overlayOpts.useDistCutoff,
        overlayOpts.distCutoff, overlayOpts.nSteps);
    auto scores = sca.calcScores(refCoords.data(), fitCoords.data(), useColor);
    if (scores[0] > bestScore) {
      bestScore = scores[0];
      bestQuat = i;
    }
  }
  initXForm.setToIdentity();
  initXForm.SetRotationFromQuaternion(quats[bestQuat].data());
}

// Return the translation that puts the extreme of refShape at the
// extreme of the fitShape along the appropriate axis.
RDGeom::Point3D getInitialTranslation(int index, const ShapeInput &refShape,
                                      const ShapeInput fitShape) {
  auto getDisp = [](const ShapeInput &shape, size_t i) -> RDGeom::Point3D {
    const DTYPE *coord = shape.getCoords().data() + shape.getExtremes()[i] * 4;
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

std::pair<double, double> alignShape(const ShapeInput &refShape,
                                     const ShapeInput &fitShape,
                                     RDGeom::Transform3D &bestXform,
                                     const ShapeOverlayOptions &overlayOpts) {
  unsigned int finalRotIndex = 1;
  switch (overlayOpts.startMode) {
    case StartMode::ROTATE_0:
    case StartMode::ROTATE_0_FRAGMENT:
      break;
    case StartMode::ROTATE_180:
    case StartMode::ROTATE_180_FRAGMENT:
    case StartMode::ROTATE_180_WIGGLE:
      finalRotIndex = 4;
      break;
    case StartMode::ROTATE_90:
    case StartMode::ROTATE_90_FRAGMENT:
      finalRotIndex = 16;
      break;
  }

  unsigned int finalTransIndex = 1;
  if (overlayOpts.startMode == StartMode::ROTATE_0_FRAGMENT ||
      overlayOpts.startMode == StartMode::ROTATE_90_FRAGMENT ||
      overlayOpts.startMode == StartMode::ROTATE_180_FRAGMENT) {
    finalTransIndex = 7;
  }

  std::pair<double, double> bestScore;
  double bestTotal = -1.0;
  RDGeom::Transform3D initialRot, initialTrans, initialXform;

  // For working values of the coordinates.
  std::vector<DTYPE> workingRef(refShape.getCoords());
  std::vector<DTYPE> workingFit(fitShape.getCoords());
  for (unsigned int j = 0; j < finalTransIndex; j++) {
    auto refDisp = getInitialTranslation(j, refShape, fitShape);
    for (unsigned int i = 0; i < finalRotIndex; i++) {
      if (overlayOpts.startMode == StartMode::ROTATE_180_WIGGLE) {
        getInitialRotation(i, refShape, fitShape, refDisp, overlayOpts,
                           workingRef, workingFit, initialRot);
      } else {
        getInitialRotation(i, initialRot);
      }
      std::vector<DTYPE> startFit(fitShape.getCoords());
      applyTransformToShape(startFit, initialRot);
      std::vector<DTYPE> startRef(refShape.getCoords());
      // Move the reference by initialTrans, leaving fit at the origin where
      // the rotations work properly.
      translateShape(startRef, refDisp);
      std::array<DTYPE, 20> outScores;
      SingleConformerAlignment sca(
          startRef.data(), workingRef.data(), refShape.getTypes().data(),
          &refShape.getCarbonRadii(), refShape.getNumAtoms(),
          refShape.getNumFeatures(), refShape.getSelfOverlapVol(),
          refShape.getSelfOverlapColor(), startFit.data(), workingFit.data(),
          fitShape.getTypes().data(), &fitShape.getCarbonRadii(),
          fitShape.getNumAtoms(), fitShape.getNumFeatures(),
          fitShape.getSelfOverlapVol(), fitShape.getSelfOverlapColor(),
          overlayOpts.optimMode, overlayOpts.optParam,
          overlayOpts.useDistCutoff, overlayOpts.distCutoff,
          overlayOpts.nSteps);
      sca.doOverlay(outScores);
      if (outScores[0] > bestTotal) {
        bestTotal = outScores[0];
        bestScore = std::make_pair(outScores[1], outScores[2]);
        RDGeom::Transform3D tmp;
        tmp.SetRotationFromQuaternion(outScores.data() + 9);
        tmp.SetTranslation(
            RDGeom::Point3D{outScores[13], outScores[14], outScores[15]});
        RDGeom::Transform3D reverseInitialTrans;
        reverseInitialTrans.SetTranslation(-refDisp);
        auto tt = reverseInitialTrans * tmp * initialRot;
        copyTransform(tt, bestXform);
      }
    }
  }
  return bestScore;
}

}  // namespace

std::pair<double, double> AlignShape(const ShapeInput &refShape,
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

std::pair<double, double> AlignMolecule(const ShapeInput &refShape, ROMol &fit,
                                        const ShapeInputOptions &fitOpts,
                                        RDGeom::Transform3D *xform,
                                        const ShapeOverlayOptions &overlayOpts,
                                        int fitConfId) {
  auto fitShape = ShapeInput(fit, fitConfId, fitOpts, overlayOpts);
  RDGeom::Transform3D tmpXform;
  auto tcs = AlignShape(refShape, fitShape, &tmpXform, overlayOpts);
  MolTransforms::transformConformer(fit.getConformer(fitConfId), tmpXform);
  if (xform) {
    copyTransform(tmpXform, *xform);
  }
  return tcs;
}

std::pair<double, double> AlignMolecule(const ROMol &ref, ROMol &fit,
                                        const ShapeInputOptions &refOpts,
                                        const ShapeInputOptions &fitOpts,
                                        RDGeom::Transform3D *xform,
                                        const ShapeOverlayOptions &overlayOpts,
                                        int refConfId, int fitConfId) {
  auto refShape = ShapeInput(ref, refConfId, refOpts, overlayOpts);
  RDGeom::Transform3D tmpXform;
  auto tcs =
      AlignMolecule(refShape, fit, fitOpts, xform, overlayOpts, fitConfId);
  return tcs;
}

std::pair<double, double> ScoreShape(const ShapeInput &refShape,
                                     const ShapeInput &fitShape,
                                     const ShapeOverlayOptions &overlayOpts) {
  auto refWorking = refShape.getCoords();
  auto fitWorking = fitShape.getCoords();
  SingleConformerAlignment sca(
      refShape.getCoords().data(), refWorking.data(),
      refShape.getTypes().data(), &refShape.getCarbonRadii(),
      refShape.getNumAtoms(), refShape.getNumFeatures(),
      refShape.getSelfOverlapVol(), refShape.getSelfOverlapColor(),
      fitShape.getCoords().data(), fitWorking.data(),
      fitShape.getTypes().data(), &fitShape.getCarbonRadii(),
      fitShape.getNumAtoms(), fitShape.getNumFeatures(),
      fitShape.getSelfOverlapVol(), fitShape.getSelfOverlapColor(),
      overlayOpts.optimMode, overlayOpts.optParam, overlayOpts.useDistCutoff,
      overlayOpts.distCutoff, overlayOpts.nSteps);
  bool includeColor = overlayOpts.optimMode != OptimMode::SHAPE_ONLY;
  auto scores = sca.calcScores(refShape.getCoords().data(),
                               fitShape.getCoords().data(), includeColor);
  return std::make_pair(scores[1], scores[2]);
}

std::pair<double, double> ScoreMolecule(const ShapeInput &refShape,
                                        const ROMol &fit,
                                        const ShapeInputOptions &fitOpts,
                                        const ShapeOverlayOptions &overlayOpts,
                                        int fitConfId) {
  auto fitShape = ShapeInput(fit, fitConfId, fitOpts, overlayOpts);
  return ScoreShape(refShape, fitShape, overlayOpts);
}

std::pair<double, double> ScoreMolecule(const ROMol &ref, const ROMol &fit,
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