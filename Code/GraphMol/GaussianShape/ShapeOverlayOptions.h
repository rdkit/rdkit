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

#include <memory>
#include <vector>

#include <RDGeneral/export.h>

using DTYPE = double;

namespace RDKit {
class ROMol;
namespace GaussianShape {

enum class RDKIT_GAUSSIANSHAPE_EXPORT StartMode {
  ROTATE_0,    //! No rotation, just normalization if requested
  ROTATE_180,  //! ROTATE_0 plus rotate by 180 degrees about each of x, y, z
  ROTATE_180_WIGGLE,  //! ROTATE_180 plus, as the PubChem code does, rotate by a
                      //! small amount (~25 degrees) about each axis and use the
                      //! highest scoring of those as the start point for that
                      //! rotation angle
  ROTATE_90,  //! ROTATE_180 plus rotate by 90 degrees about each of x, y, z
  ROTATE_0_FRAGMENT,    //! No rotation, translate probe to each end of ref
  ROTATE_180_FRAGMENT,  //! Translate probe to each end of ref and then
                        //! ROTATE_180
  ROTATE_90_FRAGMENT,  //! Translate probe to each end of ref and then ROTATE_90
};

enum class RDKIT_GAUSSIANSHAPE_EXPORT OptimMode {
  SHAPE_ONLY,              //! Drive the optimisation by shape overlap only.
  SHAPE_PLUS_COLOR_SCORE,  //! Drive the optimisation by shape, but include
                           //! color in the score to determine the best
                           //! solution.  Color never used in the optimisation
                           //! stage.
  SHAPE_PLUS_COLOR,  //! Drive the optimisation by overlap of shape and color
                     //! features.
};

struct RDKIT_GAUSSIANSHAPE_EXPORT ShapeOverlayOptions {
  ShapeOverlayOptions();
  ShapeOverlayOptions(const ShapeOverlayOptions &) = default;
  ShapeOverlayOptions(ShapeOverlayOptions &&) = default;
  ShapeOverlayOptions &operator=(const ShapeOverlayOptions &) = default;
  ShapeOverlayOptions &operator=(ShapeOverlayOptions &&) = default;

  ~ShapeOverlayOptions() = default;

  // Whether to normalise the shape by putting into
  // its canonical conformation (centred at the origin,
  // aligned along its principal axes).
  bool normalize{true};
  // Different modes for starting the optimisation.  Default is as originally
  // defined by Grant and Pickup - 180 rotations about the x, y and z axes,
  // although here the molecules are normalized so the principal axes are along
  // the cartesian axes rather than the shape quadrupole axes.
  StartMode startMode{StartMode::ROTATE_180};
  OptimMode optimMode{OptimMode::SHAPE_PLUS_COLOR_SCORE};  // Optimisation mode.
  // Whether to use carbon radii for all atoms (which is quicker) or
  // vdw radii appropriate for the elements.
  bool all_carbon_radii{true};

  DTYPE optParam{0.5};  // If using colors, the relative weights of shape and
  // color scores.
  int nSteps{100};  // Number of steps for optimiser to take.

  // Patterns for assigning features. Uses the normal RDKit ones by default.
  std::vector<std::vector<std::shared_ptr<ROMol>>> d_ph4Patterns;
  unsigned int d_nTypes{0};

  void buildPh4Patterns();
};
}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SHAPEOVERLAYOPTIONS_GUARD
