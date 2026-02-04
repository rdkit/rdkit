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

#ifndef RDKIT_SINGLECONFORMERALIGNMENT_GUARD
#define RDKIT_SINGLECONFORMERALIGNMENT_GUARD

#include <array>

#include <RDGeneral/export.h>
#include "ShapeOverlayOptions.h"

namespace RDKit {
namespace GaussianShape {
struct RDKIT_GAUSSIANSHAPE_EXPORT SingleConformerAlignment {
  SingleConformerAlignment() = delete;
  /// @brief Do the overlay for a single conformer of A against a single
  /// conformer of B.  The output in scores is the rotation and translation that
  /// moves B to optimise its score with A.
  /// @param molA - the query molecule as 1D array of 4 * NmolA entries. Each
  /// block of 4 is the coords and w parameter
  /// @param molAT - the working copy of A
  /// @param molA_type - the features types for molecule A
  /// @param NmolA - the number of atoms and features in A
  /// @param NmolA_real - the number of atoms in A
  /// @param NmolA_color - the number of features in A
  /// @param self_overlap_A - overlap volume of A with itself
  /// @param self_overlap_A_color - color overlap of A with itself
  /// @param molB - the target molecule as 1D array of 4 * NmolB entries. Each
  /// block of 4 is the coords and w parameter
  /// @param molBT - the working copy of B - the final overlaid coordinates will
  /// be left here
  /// @param molB_type - the features types for molecule B
  /// @param NmolB - the number of atoms and features in B
  /// @param NmolB_real - the number of atoms in B
  /// @param NmolB_color - the number of features in B
  /// @param optimMode - optimisation mode
  /// @param mixing_param - how to mix the 2 tanimoto values
  /// @param maxIts - maximum number of iterations for optimiser
  SingleConformerAlignment(const DTYPE *molA, DTYPE *molAT,
                           const int *molA_type, int NmolA_shape,
                           int NmolA_color, DTYPE shapeVolA, DTYPE colorVolA,
                           const DTYPE *molB, DTYPE *molBT,
                           const int *molB_type, int NmolB_shape,
                           int NmolB_color, DTYPE shapeVolB, DTYPE colorVolB,
                           OptimMode optimMode, DTYPE mixing_param,
                           bool all_carbon_radii, unsigned int maxIts);

  SingleConformerAlignment(const SingleConformerAlignment &other) = delete;
  SingleConformerAlignment(SingleConformerAlignment &&other) = delete;
  SingleConformerAlignment &operator=(const SingleConformerAlignment &other) =
      delete;
  SingleConformerAlignment &operator=(SingleConformerAlignment &&other) =
      delete;
  ~SingleConformerAlignment() = default;

  // Calculate the combined, shape, and color tanimotos as appropriate,
  // plus the volume of the shape and color overlaps, in that order.
  // Assumes that molA and molB are already in the correct configurations.
  // If includeColor is passed in true, it will compute the color score
  // irrespective of the value in d_optimMode.  We still want the color
  // score even if doing SHAPE_ONLY optimisation, for example.
  std::array<DTYPE, 5> calcScores(const DTYPE *molA, const DTYPE *molB,
                                  bool includeColor = false) const;
  // This one computes the scores from the given overlap volumes.  Color score
  // only calculated if the color volumes are non-zero.
  std::array<DTYPE, 5> calcScores(const DTYPE shapeOvVol,
                                  const DTYPE colorOvVol) const;

  // Calculate the gradients as appropriate.  Assumes that molA and molB are
  // already in the correct configurations.
  void calcGradients(const DTYPE *molA, const DTYPE *molB, const DTYPE *quat,
                     DTYPE *gradients) const;

  // Calculate the overlap volume between A and B after the given "quaternion"
  // has been applied.  The "quaternion" is 7 elements, the first 4 the
  // quaternion the last 3 the translation that currently form the
  // transformation that overlays B onto A.  If gradients is non-null,
  // they will be computed as well.
  void calcVolumeAndGradients(const std::array<DTYPE, 7> &quat,
                              DTYPE &shapeOvlpVol, DTYPE &colorOvlpVol,
                              std::array<DTYPE, 7> *gradients = nullptr);

  /// @brief Do the overlay, feeding the results into scores.
  /// @return scores - the output scores and transformation to reproduce the
  /// overlay - an array of size 20. Only the first 16 are used here. They are:
  /// 0 - the combo score
  /// 1 - the shape tanimoto
  /// 2 - the color tanimoto
  /// 3 - the overlap volume
  /// 4 - the color overlap volume
  /// 5 - the volume of A
  /// 6 - the volume of B
  /// 7 - the color volume of A
  /// 8 - the color volume of B
  /// 9-12 - the quaternion to rotate B onto A. Applied first.
  /// 13-15 - the translation to move B onto A. Applied second.
  /// 16-19 - not used at present, returned as zeros.
  /// Returns false if it didn't finish with the allowed maximum number of
  /// iterations.
  bool doOverlay(std::array<DTYPE, 20> &scores);

  // Find the quaternion and translation that maximises the volume
  // overlap.  Assume it is set to 1,0,0,0,0,0,0 on input.  Returns
  // false if it didn't finish with the allowed maximum number of
  // iterations.
  bool optimise(std::array<DTYPE, 7> &quatTrans);

  const DTYPE *d_molA;
  DTYPE *d_molAT;
  const int d_NmolA_shape;
  const int d_NmolA_color;
  const int *d_molA_type;
  const DTYPE d_shapeVolA;
  const DTYPE d_colorVolA;
  const DTYPE *d_molB;
  DTYPE *d_molBT;
  const int d_NmolB_shape;
  const int d_NmolB_color;
  const int *d_molB_type;
  DTYPE d_shapeVolB;
  DTYPE d_colorVolB;
  const OptimMode d_optimMode;
  const DTYPE d_mixing_param;
  const bool d_all_carbon_radii;
  const unsigned int d_maxIts;
};

DTYPE calcVolAndGrads(const DTYPE *molA, int numAPts, const DTYPE *molB,
                      int numBPts, const bool all_carbon_radii = true,
                      const DTYPE *quat = nullptr, DTYPE *gradients = nullptr);
// This one if for the features, and only calculates values if the types
// of 2 features match.
DTYPE calcVolAndGrads(const DTYPE *molA, int numAPts, const int *molATypes,
                      const DTYPE *molB, int numBPts, const int *molBTypes,
                      const DTYPE *quat, DTYPE *gradients);
}  // namespace ShapeAlign
}  // namespace RDKit

#endif  // RDKIT_SINGLECONFORMERALIGNMENT_GUARD
