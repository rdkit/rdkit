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

#ifndef RDKIT_SHAPEINPUT_GUARD
#define RDKIT_SHAPEINPUT_GUARD

#include <vector>

#include <RDGeneral/export.h>
#include <Geometry/point.h>

#include "ShapeOverlayOptions.h"

using DTYPE = float;

namespace RDKit {
class ROMol;
namespace ShapeAlign {

// Data for Roshambo2 shape alignment code
class RDKIT_ROSHAMBO2SHAPE_EXPORT ShapeInput {
 public:
  ShapeInput(const ROMol &mol, int confId = -1,
             const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());
  ShapeInput(const ShapeInput &other);
  ShapeInput(ShapeInput &&other) = default;
  ShapeInput &operator=(const ShapeInput &other);
  ShapeInput &operator=(ShapeInput &&other) = default;
  ~ShapeInput() = default;

  const std::vector<DTYPE> &getCoords() const { return d_coords; }
  bool getNormalized() const { return d_normalized; }
  const std::vector<int> &getTypes() const { return d_types; }
  int getNumAtoms() const { return d_numAtoms; }
  int getNumFeatures() const { return d_numFeats; }
  DTYPE getSelfOverlapVol() const { return d_selfOverlapVol; }
  DTYPE getSelfOverlapColor() const { return d_selfOverlapColor; }
  const std::array<double, 9> &getCanonicalRotation() const {
    return *d_canonRot;
  }
  const std::array<double, 3> &getCanonicalTranslation() const {
    return *d_centroid;
  }
  // Align the principal axes to the cartesian axes and centre on the origin.
  // Doesn't require that the shape was created from a molecule.  Creates
  // the necessary transformation if not already done.
  void normalizeCoords();

  void transformCoords(RDGeom::Transform3D &xform);

 private:
  void extractAtoms(const ROMol &mol, int confId);
  // Extract the features for the color scores, using RDKit pphore features
  // for now.  Other options to be added later.
  void extractFeatures(const ROMol &mol, int confId,
                       const ShapeOverlayOptions &shapeOpts);
  // Calculate the rotation and translation that will align the principal axes
  // to the cartesian axes and centre on the origin, using the conformer.
  void calcNormalization(const ROMol &mol, int confId);

  void calcExtremes();

  std::vector<DTYPE> d_coords;  // The coordinates and w parameter for the
  // atoms and features, packed as 4 floats per
  // item - x, y, z and w.  At present, w is
  // always 1.0 for atoms and features.  In future
  // it might be used to give different Gaussian
  // widths to different atoms (per the Roshambo2
  // code.)
  std::vector<int> d_types;  // The feature types.  The size is the same
  // as the number of coordinates, padded with 0
  // for the atoms.
  int d_numAtoms;                 // The number of atoms
  int d_numFeats;                 // The number of features
  DTYPE d_selfOverlapVol{0.0};    // Shape volume
  DTYPE d_selfOverlapColor{0.0};  // Color volume
  // These are the points at the extremes of the x, y and z axes.
  // they are min_x, min_y, min_z and max_x, max_y, max_z.
  std::array<size_t, 6> d_extreme_points;

  // This is the rotation and translation to align the principal axes of the
  // shape with cartesian axes.  If d_normalized is true, it has been applied
  // to the coordinates.
  bool d_normalized{false};
  std::unique_ptr<std::array<double, 9>> d_canonRot;
  std::unique_ptr<std::array<double, 3>> d_centroid;
};

// Calculate the mean position of the given atoms.
RDKIT_ROSHAMBO2SHAPE_EXPORT RDGeom::Point3D computeFeaturePos(
    const ROMol &mol, int confId, const std::vector<unsigned int> &ats);

RDKIT_ROSHAMBO2SHAPE_EXPORT void writeCoords(const std::vector<DTYPE> &shape,
                                             const std::string &label,
                                             char lineEnd = '\n');

RDKIT_ROSHAMBO2SHAPE_EXPORT RDGeom::Transform3D quatTransToTransform(
    const DTYPE *quat, const DTYPE *trans);

RDKIT_ROSHAMBO2SHAPE_EXPORT void copyTransform(const RDGeom::Transform3D &src,
                                               RDGeom::Transform3D &dest);

// Apply the transformation to the coordinates assumed to be in
// ShapeInput.d_coords form.
RDKIT_ROSHAMBO2SHAPE_EXPORT void applyTransformToShape(
    std::vector<DTYPE> &shape, RDGeom::Transform3D &xform);

}  // namespace ShapeAlign
}  // namespace RDKit

#endif  // RDKIT_SHAPEINPUT_GUARD
