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

#ifndef RDKIT_GAUSSIANSHAPE_GUARD
#define RDKIT_GAUSSIANSHAPE_GUARD

#include <RDGeneral/export.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>

namespace RDKit {
class ROMol;
class Conformer;

namespace GaussianShape {

//! Align a shape onto a reference shape.  Assumes the shapes are both
//! centred on the origin and aligned along their principal axes.
/*!
  \param refShape      the reference shape
  \param fitShape      the shape to align
  \param overlayOpts   options for the overlay
  \param xform         if passed in as non-null, will be populated with the
                       transformation matrix that aligns fit onto ref.

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
  if useColors is false)
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::pair<double, double> AlignShape(
    const ShapeInput &refShape, ShapeInput &fitShape,
    RDGeom::Transform3D *xform = nullptr,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());

//! Align a molecule to a reference shape
/*!
  \param refShape      the reference shape
  \param fit           the molecule to align
  \param fitOpts       the options for creating the fit shape
  \param xform         if passed in as non-null, will be populated with the
                       transformation matrix that aligns fit onto ref.
  \param overlayOpts   options for setting up and running the overlay
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::pair<double, double> AlignMolecule(
    const ShapeInput &refShape, ROMol &fit,
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    RDGeom::Transform3D *xform = nullptr,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int fitConfId = -1);

//! Align a molecule to a reference molecule
/*!
  \param ref           the reference molecule
  \param fit           the molecule to align
  \param refOpts       the options for creating the ref shape
  \param fitOpts       the options for creating the fit shape
  \param xform         if passed in as non-null, will be populated with the
                       transformation matrix that aligns fit onto ref.
  \param overlayOpts   options for setting up and running the overlay
  \param refConfId     (optional) the conformer to use for the reference
                       molecule
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::pair<double, double> AlignMolecule(
    const ROMol &ref, ROMol &fit,
    const ShapeInputOptions &refOpts = ShapeInputOptions(),
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    RDGeom::Transform3D *xform = nullptr,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int refConfId = -1, int fitConfId = -1);

//! Score the overlap of a shape to a reference shape without moving
//  either.  Note that if you take the output from one of the Align...
//  functions and feed it into a Score... function you won't get
//  exactly the same answer.  This is because the formula for the
//  tanimoto uses the fit volume and that is calculated once at the
//  start before the rotations and translations that form the
//  optimisation.  Floating point cruft moves the atoms by small
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
RDKIT_GAUSSIANSHAPE_EXPORT std::pair<double, double> ScoreShape(
    const ShapeInput &refShape, const ShapeInput &fitShape,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());

//! Score the overlap of a molecule to a reference shape without moving
//  either.
/*!
  \param ref           the reference shape
  \param fit           the molecule to score
  \param fitOpts       the options for creating the fit shape
  \param overlayOpts   options for setting up the shapes
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::pair<double, double> ScoreMolecule(
    const ShapeInput &refShape, const ROMol &fit,
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int fitConfId = -1);

//! Score the overlap of a molecule to a reference molecule without moving
//  either.
/*!
  \param ref           the reference molecule
  \param fit           the molecule to score
  \param refOpts       the options for creating the ref shape
  \param fitOpts       the options for creating the fit shape
  \param overlayOpts   options for setting up the shapes
  \param refConfId     (optional) the conformer to use for the reference
                       molecule
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule

  \return a pair of the shape Tanimoto value and the color Tanimoto value (zero
          if opt_param is 1.0.)  If useColors is True, uses RDKit pphore
          types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::pair<double, double> ScoreMolecule(
    const ROMol &ref, const ROMol &fit,
    const ShapeInputOptions &refOpts = ShapeInputOptions(),
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int refConfId = -1, int fitConfId = -1);

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_GAUSSIANSHAPE_GUARD
