//
//  Copyright (C) 2021-2025 David Cosgrove and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// Original author: David Cosgrove (CozChemIx Limited)
//
// This is the implementation of the functions for using the Roshambo2 backend
// to perform shape-based molecule alignments and scoring.

#include <cmath>

#include <Geometry/Transform3D.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include "Roshambo2Shape.hpp"
#include "ShapeInput.h"
#include "roshambo2/roshambo2/backends/cpp_src/cpp_helper_functions.h"

#include <ranges>

namespace RDKit {
namespace ShapeAlign {

// Apply the final overlay transform to the conformer.
void TransformConformer(Conformer &fitConf, const ShapeInput &refShape,
                        const ShapeInput &fitShape,
                        RDGeom::Transform3D &ovXform) {
  // Move fitConf to fitShape's initial centroid and principal axes
  RDGeom::Transform3D transform0;
  transform0.SetTranslation(
      RDGeom::Point3D{fitShape.getCanonicalTranslation()[0],
                      fitShape.getCanonicalTranslation()[1],
                      fitShape.getCanonicalTranslation()[2]});

  RDGeom::Transform3D transform1;
  transform1.setToIdentity();
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
  toRefRefFrame.setToIdentity();
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
  MolTransforms::transformConformer(fitConf, finalTransform);
}

namespace {
std::array<double, 4> axisAngleToQuaternion(std::array<double, 4> axisAngle) {
  std::array<double, 4> quat;
  double hangle = axisAngle[3] * 0.5;
  quat[0] = cos(hangle);
  double sinHangle = sin(hangle);
  quat[1] = axisAngle[0] * sinHangle;
  quat[2] = axisAngle[1] * sinHangle;
  quat[3] = axisAngle[2] * sinHangle;
  return quat;
}

// Return the transformation matrix for the given method and index.  Different
// methods have different numbers of starting orientations to try.  The setup
// is very much modelled on the code in Roshambo2's (in cpp_functions.cpp).
void getInitialTransformation(int index, RDGeom::Transform3D &initXform) {
  static const std::vector<std::array<double, 4>> axisAngle{
      {1, 0, 0, 0},      {1, 0, 0, M_PI},    {0, 1, 0, M_PI},
      {0, 0, 1, M_PI},   {1, 0, 0, M_PI_2},  {0, 1, 0, M_PI_2},
      {0, 0, 1, M_PI_2}, {1, 0, 0, -M_PI_2}, {0, 1, 0, -M_PI_2},
      {0, 0, 1, -M_PI_2}};
  auto quat = axisAngleToQuaternion(axisAngle[index]);
  initXform.setToIdentity();
  initXform.SetRotationFromQuaternion(quat.data());
}

}  // namespace

std::pair<double, double> AlignShape(const ShapeInput &refShape,
                                     ShapeInput &fitShape,
                                     const ShapeOverlayOptions &overlayOpts,
                                     RDGeom::Transform3D *xform) {
  int finalIndex = 1;
  switch (overlayOpts.d_mode) {
    case StartMode::AS_IS:
      break;
    case StartMode::ROTATE_180:
      finalIndex = 4;
      break;
    case StartMode::ROTATE_90:
      finalIndex = 10;
      break;
  }

  std::pair<double, double> bestScore;
  double bestTotal = 0.0;
  RDGeom::Transform3D bestXform;
  RDGeom::Transform3D initialXform;

  for (int i = 0; i < finalIndex; i++) {
    getInitialTransformation(i, initialXform);
    std::vector<DTYPE> startFit(fitShape.getCoords());
    applyTransformToShape(startFit, initialXform);
    std::vector<DTYPE> workingFit(startFit);
    std::vector<DTYPE> workingRef(refShape.getCoords());
    std::array<DTYPE, 20> outScores;
    single_conformer_optimiser(
        refShape.getCoords().data(), workingRef.data(),
        refShape.getTypes().data(),
        refShape.getNumAtoms() + refShape.getNumFeatures(),
        refShape.getNumAtoms(), refShape.getNumFeatures(),
        refShape.getSelfOverlapVol(), refShape.getSelfOverlapColor(),
        startFit.data(), workingFit.data(), fitShape.getTypes().data(),
        fitShape.getNumAtoms() + fitShape.getNumFeatures(),
        fitShape.getNumAtoms(), fitShape.getNumFeatures(),
        overlayOpts.d_rMat.data(), overlayOpts.d_pMat.data(),
        overlayOpts.d_nTypes, overlayOpts.d_useColors, overlayOpts.d_optParam,
        overlayOpts.lr_q, overlayOpts.lr_t, overlayOpts.d_nSteps,
        outScores.data());
    auto tmp =
        quatTransToTransform(outScores.data() + 9, outScores.data() + 13);
    applyTransformToShape(startFit, tmp);
    if (outScores[0] > bestTotal) {
      bestTotal = outScores[0];
      bestScore = std::make_pair(outScores[1], outScores[2]);
      auto tmp =
          quatTransToTransform(outScores.data() + 9, outScores.data() + 13);
      auto tt = tmp * initialXform;
      copyTransform(tt, bestXform);
    }
  }
  if (xform) {
    copyTransform(bestXform, *xform);
  }
  return bestScore;
}

std::pair<double, double> AlignMolecule(const ROMol &ref, ROMol &fit,
                                        RDGeom::Transform3D *xform,
                                        const ShapeOverlayOptions &overlayOpts,
                                        int refConfId, int fitConfId) {
  auto refShape = ShapeInput(ref, refConfId, overlayOpts);
  auto fitShape = ShapeInput(fit, fitConfId, overlayOpts);
  RDGeom::Transform3D tmpXform;
  auto tcs = AlignShape(refShape, fitShape, overlayOpts, &tmpXform);
  TransformConformer(fit.getConformer(fitConfId), refShape, fitShape, tmpXform);
  if (xform) {
    copyTransform(tmpXform, *xform);
  }
  return tcs;
}

std::pair<double, double> ScoreShape(const ShapeInput &refShape,
                                     const ShapeInput &fitShape,
                                     const ShapeOverlayOptions &overlayOpts) {
  DTYPE ovVol = volume(refShape.getCoords().data(), refShape.getNumAtoms(),
                       fitShape.getCoords().data(), fitShape.getNumAtoms());
  auto ts = ovVol / (refShape.getSelfOverlapVol() +
                     fitShape.getSelfOverlapVol() - ovVol);
  DTYPE tc = 0.0;
  if (overlayOpts.d_useColors) {
    DTYPE ovCol =
        volume_color(refShape.getCoords().data() + 4 * refShape.getNumAtoms(),
                     refShape.getNumFeatures(),
                     refShape.getTypes().data() + refShape.getNumAtoms(),
                     fitShape.getCoords().data() + 4 * fitShape.getNumAtoms(),
                     fitShape.getNumFeatures(),
                     fitShape.getTypes().data() + fitShape.getNumAtoms(),
                     overlayOpts.d_rMat.data(), overlayOpts.d_pMat.data(),
                     overlayOpts.d_nTypes);
    tc = ovCol / (refShape.getSelfOverlapColor() +
                  fitShape.getSelfOverlapColor() - ovCol);
  }

  return std::make_pair(ts, tc);
}

std::pair<double, double> ScoreMolecule(const ShapeInput &refShape,
                                        const ROMol &fit,
                                        const ShapeOverlayOptions &overlayOpts,
                                        int fitConfId) {
  ShapeOverlayOptions tmpOpts = overlayOpts;
  tmpOpts.d_normalize = false;
  auto fitShape = ShapeInput(fit, fitConfId, tmpOpts);
  return ScoreShape(refShape, fitShape, tmpOpts);
}

std::pair<double, double> ScoreMolecule(const ROMol &ref, const ROMol &fit,
                                        const ShapeOverlayOptions &overlayOpts,
                                        int refConfId, int fitConfId) {
  ShapeOverlayOptions tmpOpts = overlayOpts;
  tmpOpts.d_normalize = false;
  auto refShape = ShapeInput(ref, refConfId, tmpOpts);
  auto fitShape = ShapeInput(fit, fitConfId, tmpOpts);
  return ScoreShape(refShape, fitShape, tmpOpts);
}
}  // namespace ShapeAlign
}  // namespace RDKit