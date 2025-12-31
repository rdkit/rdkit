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
// This is the interface for the functions for using the Roshambo2 backend
// to perform shape-based molecule alignments and scoring.

#ifndef RDKIT_ROSHAMBO2SHAPE_GUARD
#define RDKIT_ROSHAMBO2SHAPE_GUARD

#include <Geometry/Transform3D.h>

#include "ShapeOverlayOptions.h"

namespace RDKit {
class ROMol;
namespace ShapeAlign {

class ShapeInput;

//! Apply the final overlay transform to the conformer so that it is overlaid
//  onto refShape.  This comprises 3 transformations -
//  the conformation (assumed to be the one that fitShape is derived from)
//  is moved to its centroid and principal axes, the ovXform is applied
//  and it is then moved with the inverse of the refShape normalization
//  transformation so it is finally in refShape's coordinate frame.
RDKIT_ROSHAMBO2SHAPE_EXPORT void TransformConformer(
    Conformer &fitConf, const ShapeInput &refShape, const ShapeInput &fitShape,
    RDGeom::Transform3D &ovXform);

//! Align a shape onto a reference shape.  Assumes the shapes are both
//! centred on the origin and aligned along their principal axes.
/*!
  \param refShape      the reference shape
  \param fitShape      the shape to align
  \param overlayOpts   options for the overlay
  \param xform         if passed in as non-null, will be populated with the
                       transformation matrix that aligns fit onto ref, with the
                       latter on its input position.

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
  if useColors is false)
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> AlignShape(
    const ShapeInput &refShape, ShapeInput &fitShape,
    const ShapeOverlayOptions &overlayOpts,
    RDGeom::Transform3D *xform = nullptr);

//! Align a molecule to a reference molecule
/*!
  \param ref           the reference molecule
  \param fit           the molecule to align
  \param xform         if passed in as non-null, will be populated with the
                       transformation matrix that aligns fit onto ref, with the
                       latter on its input position.
  \param overlayOpts   options for setting up and running the overlay
  \param refConfId     (optional) the conformer to use for the reference
                       molecule
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> AlignMolecule(
    const ROMol &ref, ROMol &fit, RDGeom::Transform3D *xform = nullptr,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int refConfId = -1, int fitConfId = -1);

//! Score the overlap of a shape to a reference shape without moving
//  either.  Note that if you take the output from one of the Align...
//  functions and feed it into a Score... function you won't get
//  exactly the same answer.  This is because the formula for the
//  tanimoto uses the fit volume and that is calculated once at the
//  start before the rotations and translations that form the
//  optimisation.  Floating point cruft moves that atoms by small
//  amounts relative to each other which means that the final
//  calculated volume differs slightly from the one calculated at
//  the start.  The fixed scoring obviously doesn't have this effect.
/*!
  \param refShape      the reference shape
  \param fitShape      the shape to score
  \param overlayOpts   options for setting up the shapes

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> ScoreShape(
    const ShapeInput &refShape, const ShapeInput &fitShape,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());

//! Score the overlap of a molecule to a reference shape without moving
//  either.
/*!
  \param ref           the reference shape
  \param fit           the molecule to score
  \param overlayOpts   options for setting up the shapes
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> ScoreMolecule(
    const ShapeInput &refShape, const ROMol &fit,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int fitConfId = -1);

//! Score the overlap of a molecule to a reference molecule without moving
//  either.
/*!
  \param ref           the reference molecule
  \param fit           the molecule to score
  \param overlayOpts   options for setting up the shapes
  \param refConfId     (optional) the conformer to use for the reference
                       molecule
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_PUBCHEMSHAPE_EXPORT std::pair<double, double> ScoreMolecule(
    const ROMol &ref, const ROMol &fit,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int refConfId = -1, int fitConfId = -1);

}  // namespace ShapeAlign
}  // namespace RDKit

#endif  // RDKIT_ROSHAMBO2SHAPE_GUARD
