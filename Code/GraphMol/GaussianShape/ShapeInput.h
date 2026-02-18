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

#include <array>
#include <vector>

#include <RDGeneral/export.h>
#include <Geometry/Transform3D.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>

namespace RDKit {
class ROMol;
namespace GaussianShape {

const double PI = 4 * std::atan(1.0);
// From Grant et al.
const double P = 2.7;
const double KAPPA = 2.41798793102;
using CustomFeatures =
    std::vector<std::tuple<unsigned int, RDGeom::Point3D, double>>;

struct ShapeInputOptions {
  ShapeInputOptions() = default;
  ShapeInputOptions(const ShapeInputOptions &) = default;
  ShapeInputOptions(ShapeInputOptions &&) = default;
  ShapeInputOptions &operator=(const ShapeInputOptions &) = default;
  ShapeInputOptions &operator=(ShapeInputOptions &&) = default;

  ~ShapeInputOptions() = default;

  // By default, it will create features using the RDKit pharmacophore
  // definitions.
  bool useColors{true};
  // Custom color features used verbatim.  A vector of
  // tuples of integer type, Point3D coords, double radius.
  CustomFeatures customFeatures;
  // Whether to use carbon radii for all atoms (which is quicker) or
  // vdw radii appropriate for the elements.
  bool allCarbonRadii{true};
};

// Data for shape alignment code
class RDKIT_GAUSSIANSHAPE_EXPORT ShapeInput {
 public:
  //! Create the ShapeInput object.
  //! @param mol: The molecule of interest
  //! @param confId: The conformer to use
  //! @param opts: Options for setting up the shape
  ShapeInput(const ROMol &mol, int confId = -1,
             const ShapeInputOptions &opts = ShapeInputOptions(),
             const ShapeOverlayOptions &shapeOpts = ShapeOverlayOptions());
  ShapeInput(const ShapeInput &other);
  ShapeInput(ShapeInput &&other) = default;
  ShapeInput &operator=(const ShapeInput &other);
  ShapeInput &operator=(ShapeInput &&other) = default;
  ~ShapeInput() = default;

  const std::vector<double> &getCoords() const { return d_coords; }
  bool getNormalized() const { return d_normalized; }
  const std::vector<int> &getTypes() const { return d_types; }
  unsigned int getNumAtoms() const { return d_numAtoms; }
  unsigned int getNumFeatures() const { return d_numFeats; }
  double getShapeVolume() const { return d_selfOverlapVol; }
  double getColorVolume() const { return d_selfOverlapColor; }
  const std::unique_ptr<boost::dynamic_bitset<>> &getCarbonRadii() const {
    return d_carbonRadii;
  }
  const std::array<double, 9> &getCanonicalRotation() const {
    return *d_canonRot;
  }
  const std::array<double, 3> &getCanonicalTranslation() const {
    return *d_centroid;
  }
  const std::unique_ptr<std::array<double, 3>> &getEigenValues() const {
    return d_eigenValues;
  }
  const std::array<size_t, 6> &getExtremes() const { return d_extremePoints; }

  // Align the principal axes to the cartesian axes and centre on the origin.
  // Doesn't require that the shape was created from a molecule.  Creates
  // the necessary transformation if not already done.
  void normalizeCoords();

  void transformCoords(RDGeom::Transform3D &xform);

 private:
  void extractAtoms(const ROMol &mol, int confId,
                    const ShapeInputOptions &opts);
  // Extract the features for the color scores, using RDKit pphore features
  // for now.  Other options to be added later.
  void extractFeatures(const ROMol &mol, int confId,
                       const ShapeInputOptions &shapeOpts);
  // Calculate the rotation and translation that will align the principal axes
  // to the cartesian axes and centre on the origin, using the conformer.
  void calcNormalization(const ROMol &mol, int confId);

  void calcExtremes();

  std::vector<double> d_coords;  // The coordinates and alpha parameter for the
  // atoms and features, packed as 4 floats per
  // item - x, y, z and alpha. alpha is KAPPA / (r * r) where r is the radius
  // of the atom.  This is not used if using all_atoms_carbon mode.
  std::vector<int> d_types;  // The feature types.  The size is the same
  // as the number of coordinates, padded with 0
  // for the atoms.
  int d_numAtoms;                  // The number of atoms
  int d_numFeats;                  // The number of features
  double d_selfOverlapVol{0.0};    // Shape volume
  double d_selfOverlapColor{0.0};  // Color volume
  // These are the points at the extremes of the x, y and z axes.
  // they are min_x, min_y, min_z and max_x, max_y, max_z.
  std::array<size_t, 6> d_extremePoints;
  std::unique_ptr<boost::dynamic_bitset<>>
      d_carbonRadii;  // Flags those atoms with a carbon radius, for faster
  // calculation later.

  // This is the rotation and translation to align the principal axes of the
  // shape with cartesian axes.  If d_normalized is true, it has been applied
  // to the coordinates.
  bool d_normalized{false};
  std::unique_ptr<std::array<double, 9>> d_canonRot;
  std::unique_ptr<std::array<double, 3>> d_centroid;
  // The sorted eigenvalues of the principal axes.
  std::unique_ptr<std::array<double, 3>> d_eigenValues;
};

// Calculate the mean position of the given atoms.
RDKIT_GAUSSIANSHAPE_EXPORT RDGeom::Point3D computeFeaturePos(
    const ROMol &mol, int confId, const std::vector<unsigned int> &ats);

RDKIT_GAUSSIANSHAPE_EXPORT void writeCoords(const std::vector<double> &shape,
                                            const std::string &label,
                                            char lineEnd = '\n');
RDKIT_GAUSSIANSHAPE_EXPORT void writeCoords(const double *shape,
                                            unsigned int numPts,
                                            const std::string &label,
                                            char lineEnd = '\n');

RDKIT_GAUSSIANSHAPE_EXPORT RDGeom::Transform3D quatTransToTransform(
    const double *quat, const double *trans);

RDKIT_GAUSSIANSHAPE_EXPORT void copyTransform(const RDGeom::Transform3D &src,
                                              RDGeom::Transform3D &dest);

// Apply the transformation to the coordinates assumed to be in
// ShapeInput.d_coords form.
RDKIT_GAUSSIANSHAPE_EXPORT void applyTransformToShape(
    std::vector<double> &shape, RDGeom::Transform3D &xform);
RDKIT_GAUSSIANSHAPE_EXPORT void applyTransformToShape(
    const double *inShape, double *outShape, size_t numPoints,
    RDGeom::Transform3D &xform);
RDKIT_GAUSSIANSHAPE_EXPORT void translateShape(
    std::vector<double> &shape, const RDGeom::Point3D &translation);
RDKIT_GAUSSIANSHAPE_EXPORT void translateShape(
    const double *inShape, double *outShape, size_t numPoints,
    const RDGeom::Point3D &translation);

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SHAPEINPUT_GUARD
