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
// This is the interface for the functions to perform shape-based molecule
// alignments and scoring.  It is experimental code and the API and/or
// results may change in future releases.

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
RDKIT_GAUSSIANSHAPE_EXPORT

//! Align a shape onto a reference shape.
/*!
  \param refShape      the reference shape
  \param fitShape      the shape to align
  \param xform         if passed in as non-null, will be populated with the
                       transformation matrix that aligns fit onto ref.
  \param overlayOpts   options for the overlay

\return an array of the combination score of the shape Tversky value and the
  color Tversky value (zero if colors not used) and the individual values.  If
  using color features, defaults to RDKit pharmacophore types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::array<double, 3> AlignShape(
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

\return an array of the combination score of the shape Tversky value and the
  color Tversky value (zero if colors not used) and the individual values.  If
  using color features, defaults to RDKit pharmacophore types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::array<double, 3> AlignMolecule(
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

  \return an array of the combination score of the shape Tversky value and the
  color Tversky value (zero if colors not used) and the individual values.  If
  using color features, defaults to RDKit pharmacophore types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::array<double, 3> AlignMolecule(
    const ROMol &ref, ROMol &fit,
    const ShapeInputOptions &refOpts = ShapeInputOptions(),
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    RDGeom::Transform3D *xform = nullptr,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int refConfId = -1, int fitConfId = -1);

//! Calculate scores for the alignment of all conformers of one molecule
//! onto another.  Returns a matrix of the combination scores, the conformer
//! numbers of the two molecules that gave the best overlay and the
//! transformation matrix for that overlay if requested.  The molecules
//! themselves are not altered.  scores[0][1] is the score of aligning
//! fit conformation 1 onto ref conformation 0
/*!
  \param ref           the reference molecule
  \param fit           the molecule to align
  \param refConfId     returns the reference conformer for the best scoring
                       overlay
  \param fitConfId     returns the fit conformer for the best scoring overlay
  \param combScores    the scores for all the overlays.  Will be returned sized
                       by the number of conformers of the ref and fit molecules.
                       combScores[i][j] will be the score for the jth fit
                       conformer onto the ith ref conformer.
  \param refOpts       the options for creating the ref shape
  \param fitOpts       the options for creating the fit shape
  \param overlayOpts   options for setting up and running the overlay
  \param xform         if passed in as non-null, will be populated with the
                       transformation matrix that gives the best-scoring
                       overlay.
 */
RDKIT_GAUSSIANSHAPE_EXPORT void ScoreMoleculeAllConformers(
    const ROMol &ref, const ROMol &fit, int &refConfId, int &fitConfId,
    std::vector<std::vector<double>> &combScores,
    const ShapeInputOptions &refOpts = ShapeInputOptions(),
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    RDGeom::Transform3D *xform = nullptr);

//! Score the overlap of a shape to a reference shape without moving
//  either.
/*!
  \param refShape      the reference shape
  \param fitShape      the shape to score
  \param overlayOpts   options for controlling the volume calculation
  \param overlapVols  if not-null, is filled with the raw overlap volumes

\return an array of the combination score of the shape Tversky value and the
  color Tversky value (zero if colors not used) and the individual values.  If
  using color features, defaults to RDKit pharmacophore types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::array<double, 3> ScoreShape(
    const ShapeInput &refShape, const ShapeInput &fitShape,
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    std::array<double, 2> *overlapVols = nullptr);

//! Score the overlap of a molecule to a reference shape without moving
//  either.
/*!
  \param refShape      the reference shape
  \param fit           the molecule to score
  \param fitOpts       the options for creating the fit shape
  \param overlayOpts   options for controlling the volume calculation
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule
  \param overlapVols  if not-null, is filled with the raw overlap volumes

\return an array of the combination score of the shape Tversky value and the
  color Tversky value (zero if colors not used) and the individual values.  If
  using color features, defaults to RDKit pharmacophore types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::array<double, 3> ScoreMolecule(
    const ShapeInput &refShape, const ROMol &fit,
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int fitConfId = -1, std::array<double, 2> *overlapVols = nullptr);

//! Score the overlap of a molecule to a reference molecule without moving
//  either.
/*!
  \param ref           the reference molecule
  \param fit           the molecule to score
  \param refOpts       the options for creating the ref shape
  \param fitOpts       the options for creating the fit shape
  \param overlayOpts   options for controlling the volume calculation
  \param refConfId     (optional) the conformer to use for the reference
                       molecule
  \param fitConfId     (optional) the conformer to use for the fit
                       molecule
  \param overlapVols  if not-null, is filled with the raw overlap volumes

\return an array of the combination score of the shape Tverksy value and the
  color Tversky value (zero if colors not used) and the individual values.  If
  using color features, defaults to RDKit pharmacophore types for the features.
*/
RDKIT_GAUSSIANSHAPE_EXPORT std::array<double, 3> ScoreMolecule(
    const ROMol &ref, const ROMol &fit,
    const ShapeInputOptions &refOpts = ShapeInputOptions(),
    const ShapeInputOptions &fitOpts = ShapeInputOptions(),
    const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions(),
    int refConfId = -1, int fitConfId = -1,
    std::array<double, 2> *overlapVols = nullptr);

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_GAUSSIANSHAPE_GUARD
