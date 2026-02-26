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
// Options for the Roshambo2-based shape overlay.

#ifndef RDKIT_SHAPEOVERLAYOPTIONS_GUARD
#define RDKIT_SHAPEOVERLAYOPTIONS_GUARD

#include <RDGeneral/export.h>

using DTYPE = float;

namespace RDKit {
class ROMol;
namespace ShapeAlign {

enum class RDKIT_ROSHAMBO2SHAPE_EXPORT StartMode {
  AS_IS,
  ROTATE_180,
  ROTATE_90
};

struct RDKIT_ROSHAMBO2SHAPE_EXPORT ShapeOverlayOptions {
  ShapeOverlayOptions();
  ShapeOverlayOptions(const ShapeOverlayOptions &) = default;
  ShapeOverlayOptions(ShapeOverlayOptions &&) = default;
  ShapeOverlayOptions &operator=(const ShapeOverlayOptions &) = default;
  ShapeOverlayOptions &operator=(ShapeOverlayOptions &&) = default;

  ~ShapeOverlayOptions() = default;

  // Whether to normalise the shape by putting into
  // its canonical conformation (centred at the origin,
  // aligned along its principal axes).
  bool d_normalize{true};
  // Different modes for starting the optimisation.  Default is as originally
  // defined by Grant and Pickup - 180 rotations about the x, y and z axes,
  // although here the molecules are normalized so the principal axes are along
  // the cartesian axes rather than the shape quadrupole axes.
  StartMode d_startMode{StartMode::ROTATE_180};
  bool d_useColors{true};  // Whether to use the features/colors as part of the
  // optimisation.
  DTYPE d_optParam{0.5};  // If using colors, the relative weights of shape and
  // color scores.
  DTYPE lr_q{0.1};    // Learning rate for optimising quaternion
  DTYPE lr_t{0.1};    // Learning rate for optimising translation
  int d_nSteps{100};  // Number of steps for optimiser to take.

  // Patterns for assigning features. Uses the normal RDKit ones by default.
  std::vector<std::vector<std::shared_ptr<ROMol>>> d_ph4Patterns;
  unsigned int d_nTypes{0};

  // The r and p matrices for the features in the shapes.  Stored
  // as linearised square matrices, of size nFeats * nFeats, where nFeats
  // is d_ph4Patterns->size(). Currently both parameters are set to 1.0
  // for all types.
  std::vector<DTYPE> d_rMat;
  std::vector<DTYPE> d_pMat;

  void buildPh4Patterns();
  void buildInteractionMatrices();
};
}  // namespace ShapeAlign
}  // namespace RDKit

#endif  // RDKIT_SHAPEOVERLAYOPTIONS_GUARD
