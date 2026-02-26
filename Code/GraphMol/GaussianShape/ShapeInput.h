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
#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/unique_ptr.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>

// The code below was provided by Claude (Sonnet 4.6).
// If first tried to get me to use boost/serialization/dynamic_bitset.hpp
// and then admitted that it had made that up.
namespace boost {
namespace serialization {

template <class Archive, typename Block, typename Allocator>
void serialize(Archive &ar, boost::dynamic_bitset<Block, Allocator> &bs,
               const unsigned int /*version*/) {
  size_t num_bits = bs.size();
  ar & num_bits;

  std::vector<Block> blocks;

  if (Archive::is_saving::value) {
    to_block_range(bs, std::back_inserter(blocks));
  }

  ar & blocks;

  if (Archive::is_loading::value) {
    bs.resize(num_bits);
    from_block_range(blocks.begin(), blocks.end(), bs);
    bs.resize(num_bits);  // trim any excess bits
  }
}

}  // namespace serialization
}  // namespace boost

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
  bool useColors{
      true};  //! Whether to build the color features.  By default, it will
              //! create features using the RDKit pharmacophore definitions.

  CustomFeatures customFeatures;  //! Custom color features used verbatim.  A
                                  //! vector of tuples of integer type, Point3D
                                  //! coords, double radius.
  std::vector<unsigned int>
      atomSubset;  //! If not empty, use just these atoms in the molecule to
                   //! form the ShapeInput object.
  std::vector<std::pair<unsigned int, double>>
      atomRadii;  //! Use these non-standard radii for these atoms. The int is
                  //! for the atom index in the molecule, not the atomic number.
                  //! Not all atoms need be specified, just some radii can be
                  //! over-ridden, with the rest left as standard.
  bool allCarbonRadii{
      true};  //! Whether to use carbon radii for all atoms (which is quicker
              //! but less accurate) or vdw radii appropriate for the elements.
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
             const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());
  ShapeInput(const std::string &str) {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    std::stringstream ss(str);
    boost::archive::text_iarchive ia(ss);
    ia &*this;
#endif
  }
  ShapeInput(const ShapeInput &other);
  ShapeInput(ShapeInput &&other) = default;
  ShapeInput &operator=(const ShapeInput &other);
  ShapeInput &operator=(ShapeInput &&other) = default;
  virtual ~ShapeInput() = default;

  std::string toString() const {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    std::stringstream ss;
    boost::archive::text_oarchive oa(ss);
    oa &*this;
    return ss.str();
#endif
  }

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
    return d_canonRot;
  }
  const std::array<double, 3> &getCanonicalTranslation() const {
    return d_canonTrans;
  }
  const std::array<double, 3> &getEigenValues() const { return d_eigenValues; }
  const std::array<size_t, 6> &getExtremes() const { return d_extremePoints; }

  // Align the principal axes to the cartesian axes and centre on the origin.
  // Doesn't require that the shape was created from a molecule.  Creates
  // the necessary transformation if not already done.
  void normalizeCoords();

  void transformCoords(RDGeom::Transform3D &xform);

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int) {
    ar & d_coords;
    ar & d_types;
    ar & d_numAtoms;
    ar & d_numFeats;
    ar & d_selfOverlapVol;
    ar & d_selfOverlapColor;
    ar & d_extremePoints;
    ar & d_carbonRadii;
    ar & d_normalized;
    ar & d_canonRot;
    ar & d_canonTrans;
    ar & d_eigenValues;
  }
#endif

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
  std::array<double, 9> d_canonRot;
  std::array<double, 3> d_canonTrans;
  // The sorted eigenvalues of the principal axes.
  std::array<double, 3> d_eigenValues;
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
