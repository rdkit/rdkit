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
// This is the class that does optimises a moving molecule (fit)
// to maximise its Gaussian overlap with the reference molecule (ref).
// The optimiser is a modified BFGS taken in large part, but re-arranged
// for readability, from the PubChem shape overlay code
// https://github.com/ncbi/pubchem-align3d/blob/main/shape_neighbor.cpp

#ifndef RDKIT_SINGLECONFORMERALIGNMENT_GUARD
#define RDKIT_SINGLECONFORMERALIGNMENT_GUARD

#include <array>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/export.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>

namespace RDKit {
namespace GaussianShape {
struct RDKIT_GAUSSIANSHAPE_EXPORT SingleConformerAlignment {
  SingleConformerAlignment() = delete;
  /// @brief Do the overlay for a single conformer of fit against a single
  /// conformer of ref.  The output in scores is the rotation and translation
  /// that moves fit to optimise its score with ref.
  /// @param ref - the query molecule as 1D array of 4 * N entries. Each
  /// block of 4 is the coords and atom radius
  /// @param refTypes - the feature types for molecule ref
  /// @param refCarbonRadii - whether each atom has a carbon radius
  /// @param nRefShape - the number of atoms in ref
  /// @param nRefColor - the number of features in ref
  /// @param refShapeVol - overlap volume of ref with itself
  /// @param refColorVol - color overlap of ref with itself
  /// @param fit - the fit molecule as 1D array of 4 * N entries. Each
  /// block of 4 is the coords and atom radius.
  /// @param fitTypes - the feature types for fit molecule
  /// @param fitCarbonRadii - whether each atom has a carbon radius
  /// @param nFitShape - the number of atoms in fit
  /// @param nFitColor - the number of features in fit
  /// @param fitShapeVol - overlap volume of fit with itself
  /// @param fitColorVol - color overlap of fit with itself
  /// @param optimMode - optimisation mode
  /// @param mixingParam - how to mix the 2 tanimoto values
  /// @param useCutoff - whether to use a distance cutoff in the volume
  /// calculation
  /// @param distCutoff - the cutoff to use if we're doing it.
  /// carbon.  This makes it faster but less correct.
  /// @param maxIts - maximum number of iterations for optimiser
  /// of optimiser
  SingleConformerAlignment(
      const std::vector<double> &ref, const int *refTypes,
      const std::unique_ptr<boost::dynamic_bitset<>> &refCarbonRadii,
      int nRefShape, int nRefColor, double refShapeVol, double refColorVol,
      const std::vector<double> &fit, const int *fitTypes,
      const std::unique_ptr<boost::dynamic_bitset<>> &fitCarbonRadii,
      int nFitShape, int nFitColor, double fitShapeVol, double fitColorVol,
      const std::array<double, 7> &initQuatTrans, OptimMode optimMode,
      double mixingParam, bool useCutoff, double distCutoff,
      double shapeConvergenceCriterion, unsigned int maxIts);

  SingleConformerAlignment(const SingleConformerAlignment &other) = delete;
  SingleConformerAlignment(SingleConformerAlignment &&other) = delete;
  SingleConformerAlignment &operator=(const SingleConformerAlignment &other) =
      delete;
  SingleConformerAlignment &operator=(SingleConformerAlignment &&other) =
      delete;
  ~SingleConformerAlignment() = default;

  void setQuatTrans(const std::array<double, 7> &quatTrans) {
    d_quatTrans = quatTrans;
  }

  // Get the final transformation by adding the initial transformation
  // and the optimised final answer.
  void getFinalQuatTrans(RDGeom::Transform3D &xform) const;

  // Calculate the combined, shape, and color tanimotos as appropriate,
  // plus the volume of the shape and color overlaps, in that order.
  // Assumes that ref and fit are already in the correct configurations.
  // If includeColor is passed in true, it will compute the color score
  // irrespective of the value in d_optimMode.  We still want the color
  // score even if doing SHAPE_ONLY optimisation, for example.
  std::array<double, 5> calcScores(const double *ref, const double *fit,
                                   bool includeColor = false) const;
  // This one applies the current quatTrans to the coords and then calculates
  // the score.
  std::array<double, 5> calcScores(bool includeColor = false);
  // This one computes the scores from the given overlap volumes.  Color score
  // only calculated if the color volumes are non-zero.
  std::array<double, 5> calcScores(const double shapeOvVol,
                                   const double colorOvVol,
                                   bool includeColor = true) const;

  // Apply the quatTrans to the ref and fit shapes and put the results
  // into their tmp equivalents.  Ref is translated by the -ve of the
  // translation, fit is rotated by the rotation bit.
  void applyQuatTrans(const std::array<double, 7> &quatTrans);

  // Calculate the overlap volume between A and B after the given "quaternion"
  // has been applied.  The "quaternion" is 7 elements, the first 4 the
  // quaternion the last 3 the translation that currently form the
  // transformation that overlays B onto A.
  void calcVolumeAndGradients(const std::array<double, 7> &quatTrans,
                              double &shapeOvlpVol, double &colorOvlpVol,
                              std::array<double, 7> &gradients);

  /// @brief Do the overlay, feeding the results into scores.
  /// @return scores - the output scores and transformation to reproduce the
  /// overlay - an array of size 20. Only the first 16 are used here. They are:
  /// 0 - the combo score
  /// 1 - the shape tanimoto
  /// 2 - the color tanimoto
  /// 3 - the shape overlap volume
  /// 4 - the color overlap volume
  /// 5 - the shape volume of fit
  /// 6 - the shape volume of ref
  /// 7 - the color volume of fit
  /// 8 - the color volume of ref
  /// 9-12 - the quaternion to rotate fit onto ref. Applied first.
  /// 13-15 - the translation to move fit onto ref. Applied second.
  /// 16-19 - not used at present, returned as zeros.
  /// Returns false if it didn't finish with the allowed maximum number of
  /// iterations.
  bool doOverlay(std::array<double, 20> &scores, unsigned int cycle);

  // Find the quaternion and translation that maximises the volume
  // overlap appropriate to d_optimMode.  Returns false if it didn't finish with
  // the allowed maximum number of iterations.
  bool optimise(unsigned int maxIters);

  std::vector<double> d_ref;
  std::vector<double> d_refTemp;
  const int *d_refTypes;
  const std::unique_ptr<boost::dynamic_bitset<>> &d_refCarbonRadii;
  const int d_nRefShape;
  const int d_nRefColor;
  const double d_refShapeVol;
  const double d_refColorVol;
  std::vector<double> d_fit;
  std::vector<double> d_fitTemp;
  const int *d_fitTypes;
  const std::unique_ptr<boost::dynamic_bitset<>> &d_fitCarbonRadii;
  const int d_nFitShape;
  const int d_nFitColor;
  double d_fitShapeVol;
  double d_fitColorVol;
  std::array<double, 7> d_initQuatTrans;
  const OptimMode d_optimMode;
  const double d_mixingParam;
  const bool d_useCutoff;
  const double d_distCutoff2;
  const double d_shapeConvergenceCriterion;
  const unsigned int d_maxIts;
  // The quaternion/translation as the optimisation proceeds
  std::array<double, 7> d_quatTrans{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // The step sizes of the quaternion and translation during the
  // optimisition.  Taken from the PubChem code.
  double d_qStepSize{-0.001};
  double d_tStepSize{-0.01};
  // Scratch space for the gradients dr/dQ of the fit molecule.
  mutable std::vector<std::array<double, 12>> d_gradConverters;
};

// Compute the volume overlap and optionally "quaternion" gradients for the
// overlap volume of ref and fit, wrt fit.  fit is the original coords of
// the fit molecule, fitTemp is those subject to any transformation applied
// by the quaternion we're using to optimise the overlap volume.  If
// gradients is null, they won't be calculated.  They are assumed to be
// initialised correctly.
// This is for the atoms/shape features.
double calcVolAndGrads(
    const double *ref, int numRefPts,
    const std::unique_ptr<boost::dynamic_bitset<>> &refCarbonRadii,
    const double *fit, int numFitPts,
    const std::unique_ptr<boost::dynamic_bitset<>> &fitCarbonRadii,
    std::vector<std::array<double, 12>> &gradConverters, const bool useCutoff,
    const double distCutoff2, const double *quat = nullptr,
    double *gradients = nullptr);
// This one is for the features, and only calculates values if the types
// of 2 features match.
double calcVolAndGrads(const double *ref, int numRefPts, const int *refTypes,
                       const double *fit, int numFitPts, const int *fitTypes,
                       int numFitShape,
                       std::vector<std::array<double, 12>> &gradConverters,
                       const bool useCutoff, const double distCutoff2,
                       const double *quat, double *gradients);

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SINGLECONFORMERALIGNMENT_GUARD
