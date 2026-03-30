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
class RWMol;
namespace GaussianShape {

constexpr double CARBON_RAD = 1.70;
constexpr double DUMMY_RAD = 2.16;  // same as Xe
// From Grant et al.
constexpr double P = 2.7;
constexpr double KAPPA = 2.41798793102;

struct CustomFeature {
  CustomFeature(
      unsigned int t, const RDGeom::Point3D &p, double r,
      const std::vector<unsigned int> &a = std::vector<unsigned int>())
      : type(t), pos(p), rad(r), atoms(a) {}
  unsigned int type;
  RDGeom::Point3D pos;
  double rad;
  std::vector<unsigned int>
      atoms;  // That the feature was derived from.  May be left empty.
};

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

  std::vector<std::vector<CustomFeature>>
      customFeatures;  //! Custom color features used verbatim.  A
                       //! vector of vectors of tuples of integer type, Point3D
                       //! coords, double radius for each conformation in the
                       //! input molecule.  One outer vector for each
                       //! conformation in the molecule.
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
  double shapePruneThreshold{-1.0};  //! If there is more than 1 conformer for
                                     //! the input molecule, prune the shapes so
                                     //! that none of them are more similar to
                                     //! each other than the threshold.  Default
                                     //! -1.0 means no pruning.
  bool includeDummies{true};         //! Whether to include dummy atoms or not.
};

// Data for shape alignment code
class RDKIT_GAUSSIANSHAPE_EXPORT ShapeInput {
 public:
  //! Create the ShapeInput object.
  //! @param mol: The molecule of interest
  //! @param confId: The conformer to use.  If -1, uses all conformers.
  //! @param opts: Options for setting up the shape
  ShapeInput(const ROMol &mol, int confId = -1,
             const ShapeInputOptions &opts = ShapeInputOptions(),
             const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());
  //! Create a ShapeInput object with a single shape copied from
  //! other.
  //! @param other: the ShapeInput that supplies the shape
  //! @param shapeNum: the number of the shape of interest.
  ShapeInput(const ShapeInput &other, unsigned int shapeNum);
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
  ~ShapeInput() = default;

  // Merge the other ShapeInput, assuming it has the correct number
  // of atoms etc.  Empties other, unless they can't be merged in which case
  // it returns unscathed.
  void merge(ShapeInput &other);

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

  std::string getSmiles() const { return d_smiles; }
  unsigned int getActiveShape() const { return d_activeShape; }
  //! Set the currently active conformation to the new value.
  //! @param newShape: the number of the conformation to be used
  //!                 for future calculations.  Counts from 0,
  //!                 obviously.  If invalid, throws a runtime
  //!                 error.
  void setActiveShape(unsigned int newShape);
  // Note that the coords returned is a vector size 3*getNumAtoms()
  const std::vector<double> &getCoords() const {
    return d_coords[d_activeShape];
  }
  const std::vector<double> &getAlphas() const { return d_alphas; }
  //! Fetch the coordinates of the atoms and optionally features.
  std::vector<RDGeom::Point3D> getAtomPoints(bool includeColors = false) const;
  bool getNormalized() const { return d_normalizeds[d_activeShape]; }
  const std::vector<int> &getTypes() const { return d_types; }
  unsigned int getNumAtoms() const { return d_numAtoms; }
  unsigned int getNumFeatures() const { return d_numFeats; }
  unsigned int getNumShapes() const { return d_coords.size(); }
  double getShapeVolume() const {
    return d_selfOverlapShapeVols[d_activeShape];
  }
  double getShapeVolume(unsigned int shapeNum) const;
  double getColorVolume() const {
    return d_selfOverlapColorVols[d_activeShape];
  }
  double getColorVolume(unsigned int shapeNum) const;
  const boost::dynamic_bitset<> *getCarbonRadii() const {
    return d_carbonRadii.get();
  }
  // These functions use cached values if available.
  const std::array<double, 9> &calcCanonicalRotation();
  const std::array<double, 3> &calcCanonicalTranslation();
  const std::array<double, 3> &calcEigenValues();
  const std::array<size_t, 6> &calcExtremes();
  // Return the principal moments of inertia, if Eigen3 is available, and the
  // eigenvalues of the canonical transformation if not.
  std::array<double, 3> calcMomentsOfInertia(bool includeColors = false) const;

  // Align the principal axes to the cartesian axes and centre on the origin.
  // Doesn't require that the shape was created from a molecule.  Creates
  // the necessary transformation if not already done.
  void normalizeCoords();

  void transformCoords(RDGeom::Transform3D &xform);

  // Make a molecule from the shape.  If required, features are added
  // as xenon atoms.  If withBonds is false, just makes a molecule from
  // the atoms, otherwise builds a full molecule.
  std::unique_ptr<RWMol> shapeToMol(bool includeColors = false,
                                    bool withBonds = true) const;

  // Find the best similarity score between all shapes in this shape and the
  // other one. Stops as soon as it gets something above the threshold.
  // The score runs between 0.0 and 1.0, so the default threshold of -1.0
  // means no threshold. Fills in the shape numbers of the two that were
  // responsible if there is something above the threshold, and the
  // transformation that did it. Returns -1.0 for the similarity if there was
  // nothing above the threshold.  Note that the shape numbers are not
  // necessarily the same as the original molecule conformation numbers.
  double bestSimilarity(
      ShapeInput &fitShape, unsigned int &bestThisShape,
      unsigned int &bestFitShape, RDGeom::Transform3D &bestXform,
      double threshold = -1.0,
      const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());

  // Return the maximum similarity achievable between the 2 shapes.  The
  // maximum similarity is when one shape is entirely inside the other.  This
  // returns the similarity in that case, which is the upper bound on what
  // is achievable between these 2 shapes.
  double maxPossibleSimilarity(
      const ShapeInput &fitShape,
      const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions()) const;

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int);
#endif

 private:
  void extractAtoms(const ROMol &mol, int confId, const ShapeInputOptions &opts,
                    bool fillAlphas);
  // Extract the features for the color scores, using RDKit pphore features
  // for now.  Other options to be added later.
  void extractFeatures(const ROMol &mol, unsigned int confId,
                       const ShapeInputOptions &shapeOpts, bool fillAlphas);
  // Calculate the rotation and translation that will align the principal axes
  // to the cartesian axes and centre on the origin.
  void calcNormalization();

  void calculateExtremes();

  unsigned int d_activeShape;

  std::vector<std::vector<double>>
      d_coords;  // The coordinates for the atoms and features,
  // packed as 3 floats per item - x, y, z
  std::vector<double> d_alphas;  // The alpha values for the atoms and features.
  // alpha is KAPPA / (r * r) where r is the radius
  // of the atom.  This is not used if using all_atoms_carbon mode.
  std::vector<int> d_types;  // The feature types.  The size is the same
  // as the number of coordinates, padded with 0
  // for the atoms.
  unsigned int d_numAtoms;                     // The number of atoms
  unsigned int d_numFeats;                     // The number of features
  std::vector<double> d_selfOverlapShapeVols;  // Shape volume
  std::vector<double> d_selfOverlapColorVols;  // Color volume
  // These are the points at the extremes of the x, y and z axes.
  // they are min_x, min_y, min_z and max_x, max_y, max_z.
  std::vector<std::array<size_t, 6>> d_extremePointss;
  std::unique_ptr<boost::dynamic_bitset<>>
      d_carbonRadii;  // Flags those atoms with a carbon radius, for faster
  // calculation later.
  std::string d_smiles;  // The SMILES string of the input molecule

  // These are the rotation and translation matrices to align the principal
  // axes of the shape with cartesian axes.  If d_normalized is true, it has
  // been applied to the coordinates.
  boost::dynamic_bitset<> d_normalizeds;
  // If the shape is moved, the normalization matrices are no longer valid.
  // This flags that so it is re-computed as required.
  boost::dynamic_bitset<> d_normalizationOKs;

  std::vector<std::array<double, 9>> d_canonRots;
  std::vector<std::array<double, 3>> d_canonTranss;
  // The sorted eigenvalues of the principal axes.
  std::vector<std::array<double, 3>> d_eigenValuess;

  // Prune the shapes so none a more similar to each other than
  // the threshold.
  void pruneShapes(double simThreshold);
  void selectConformations(const std::vector<int> &picks);
  void calculateSelfOverlaps(const ShapeOverlayOptions &overlayOpts);
  // Sort the shapes in descending order of the sume of the shape
  // and color volumes.
  void sortShapesByVolumes();
};

#ifdef RDK_USE_BOOST_SERIALIZATION
template <class Archive>
void ShapeInput::serialize(Archive &ar, const unsigned int) {
  ar & d_activeShape;
  ar & d_coords;
  ar & d_alphas;
  ar & d_types;
  ar & d_numAtoms;
  ar & d_numFeats;
  ar & d_selfOverlapShapeVols;
  ar & d_selfOverlapColorVols;
  ar & d_extremePointss;
  ar & d_carbonRadii;
  ar & d_smiles;
  ar & d_normalizeds;
  ar & d_normalizationOKs;
  ar & d_canonRots;
  ar & d_canonTranss;
  ar & d_eigenValuess;
}
#endif

// Extract the features from the molecule, optionally just for the subset
// of atoms.
RDKIT_GAUSSIANSHAPE_EXPORT void findFeatures(
    const ROMol &mol, int confId, std::vector<CustomFeature> &features,
    const std::vector<unsigned int> &atomSubset = std::vector<unsigned int>());

// Calculate the mean position of the given atoms.
RDKIT_GAUSSIANSHAPE_EXPORT RDGeom::Point3D computeFeaturePos(
    const ROMol &mol, int confId, const std::vector<unsigned int> &ats);

RDKIT_GAUSSIANSHAPE_EXPORT RDGeom::Transform3D quatTransToTransform(
    const double *quat, const double *trans);

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

// Maximum possible score of the 2 shape (v[12]) and color (c[12]) volumes
double maxScore(double v1, double v2, double c1, double c2,
                const ShapeOverlayOptions &overlayOpts);

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // RDKIT_SHAPEINPUT_GUARD
