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
  ROTATE_45,  //! ROTATE_180 plus rotate by 45 degrees about pairs of each of x,
              //! y, z
  ROTATE_0_FRAGMENT,    //! No rotation, translate probe to each end of ref
  ROTATE_180_FRAGMENT,  //! Translate probe to each end of ref and then
                        //! ROTATE_180
  ROTATE_45_FRAGMENT,  //! Translate probe to each end of ref and then ROTATE_90
  A_LA_PUBCHEM,  //! Uses the eigenvalues of the principal vectors do decide
                 //! whether to do ROTATE_180_WIGGLE or ROTATE_45
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
  //! Different modes for starting the optimisation.  Default is as used by the
  //! PubChem code.  The molecules are normalized so the principal axes are
  //! along the cartesian axes rather than the shape quadrupole axes as Grant et
  //! al. did.
  StartMode startMode{StartMode::A_LA_PUBCHEM};
  OptimMode optimMode{
      OptimMode::SHAPE_PLUS_COLOR_SCORE};  //! Optimisation mode.
  double optParam{0.5};  //! If using colors, the relative weights of shape and
                         //! color scores.
  int nSteps{100};       //! Maximum number of steps for optimiser to take.

  bool normalize{
      true};  //! Whether to normalise the shapes by putting them into their
              //! canonical conformations (centred at the origin, aligned along
              //! its principal axes) before starting.
  bool useDistCutoff{
      true};  //! Whether to use a distance cutoff fo the volume calculation.
  double distCutoff{4.5};  //! The distance cutoff.  If 2 atoms are more than
                           //! this distance apart, they are not included in the
                           //! volume calculation. A smaller value is faster but
                           //! less precise.
  double shapeConvergenceCriterion{
      0.001};  //! Optimisation stops when the shape tanimoto changes by less
               //! than this amount.  A larger number is faster but less
               //! precise.
};
}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SHAPEOVERLAYOPTIONS_GUARD
