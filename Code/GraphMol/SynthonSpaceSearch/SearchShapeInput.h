//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// A class derived from the GaussianShape::ShapeInput object with more
// stuff.

#ifndef SEARCHSHAPEINPUT_H
#define SEARCHSHAPEINPUT_H

#include <GraphMol/GaussianShape/ShapeInput.h>

namespace RDKit {
namespace GaussianShape {
// Make a subclass of ShapeInput with some extra info, including allowing
// for multiple conformations of the same atoms.  ShapeInput is in the
// GaussianShape namespace, so it's easier if this is too.
class RDKIT_SYNTHONSPACESEARCH_EXPORT SearchShapeInput : public ShapeInput {
 public:
  SearchShapeInput() = default;
  // Creates the shapes for all conformations in the molecule.  Prunes them
  // so that they are all at least pruneThreshold comboScore apart and
  // sorts them into descending order of total volume.
  SearchShapeInput(
      const ROMol &mol, double pruneThreshold,
      const ShapeInputOptions &opts = ShapeInputOptions(),
      const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());

  SearchShapeInput(const std::string &str);

  SearchShapeInput(const ShapeInput &si);
  SearchShapeInput(const SearchShapeInput &other) = default;
  SearchShapeInput(SearchShapeInput &&other) = default;
  SearchShapeInput &operator=(const SearchShapeInput &other) = default;
  SearchShapeInput &operator=(SearchShapeInput &&other) = default;
  ~SearchShapeInput() override = default;

  using ShapeInput::getColorVolume;
  using ShapeInput::getShapeVolume;
  double getShapeVolume(unsigned int shapeNum) const;
  double getColorVolume(unsigned int shapenum) const;
  double getDummyVolume(unsigned int shapeNum) const;

  // Merge the other ShapeInputSet, assuming it has the correct number
  // of atoms etc.  Empties the multiple conformation parts of other,
  // unless they can't be merged in which case it returns unscathed.
  // Always leaves the base elements of other unchanged.
  void merge(SearchShapeInput &other);

  size_t getNumShapes() const { return d_confCoords.size(); }
  const std::vector<std::vector<double>> &getConfCoords() const {
    return d_confCoords;
  }

  // Make a single ShapeInput from the given shape number.
  // If shapeNum is out of range, use the first shape.
  ShapeInput makeSingleShape(unsigned int shapeNum) const;

  void setActiveShape(unsigned int shapeNum);

  std::string toString() const;

  // Cut the conformations down so that none of them have shapes that
  // are more than the simThreshold in similarity with each other.
  void pruneShapes(double simThreshold);

  // Sort the shapes in descending order of shape and color volume.
  void sortShapesByScore();

  // Find the best similarity score between all shapes this shape and the
  // other one. Stops as soon as it gets something above the threshold.
  // The threshold is applied to the combo score. The similarity is between 0.0
  // and 1.0 so the default threshold of 2.0 effectively means no threshold.
  // Fills in the conformations of the two that were responsible if there is
  // something above the threshold.  Returns -1.0 for
  // the similarity if there was nothing above the threshold.
  double bestSimilarity(
      SearchShapeInput &fitShape, unsigned int &bestThisConf,
      unsigned int &bestFitConf, double threshold = 2.0,
      const ShapeOverlayOptions &overlayOpts = ShapeOverlayOptions());

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int);
#endif

 private:
  unsigned int d_numDummies{0};
  double d_dummyVol{0.0};
  unsigned int d_actConf{0};
  // The base class has 4 floats per atom for the coords, the 4th
  // being the radius.  That's a bit of a waste of space since the
  // radii are constant, so these are just 3 floats per atom.  It'll
  // make the copying a bit slower.
  std::vector<std::vector<double>> d_confCoords;
  std::vector<unsigned int>
      d_molConfs;  // the conformer from the input molecule
  // that this shape refers to.  The shapes
  // are pruned and sorted so they may end
  // up not in the original order.
  // All the other information that varies from conformation to conformation.
  std::vector<double> d_dummyVolumes;
  std::vector<double> d_shapeVolumes;
  std::vector<double> d_colorVolumes;
  std::vector<std::array<size_t, 6>> d_extremePointss;
  std::vector<std::array<double, 9>> d_canonRots;
  std::vector<std::array<double, 3>> d_centroids;
  std::vector<std::array<double, 3>> d_eigenValuess;

  void initializeFromBase();
  // Copy the conf coords for shapeNum into coords, in 4-float format
  // for use in the base ShapeInput.
  void confCoordsToShapeCoords(unsigned int shapeNum,
                               std::vector<double> &coords) const;

  void selectConformations(const std::vector<int> &picks);
};

// Make a SearchShapeInput from all conformations of a molecule and then
// prune them at the given threshold. so that all selected shapes have
// a similarity score with each other that's less than the threshold.
// Assumes that shapeOpts.atomRadii only includes the dummy atoms, and
// doesn't use non-standard radii for other atoms.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<SearchShapeInput>
PrepareConformers(const RDKit::ROMol &mol,
                  const ShapeInputOptions &shapeOpts = ShapeInputOptions(),
                  double pruneThreshold = 1.9);

// Mock a molecule up from the shape for visual inspection.  No bonds.
// Atoms are C, features are N.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::unique_ptr<RDKit::RWMol> shapeToMol(
    const ShapeInput &shape);

#ifdef RDK_USE_BOOST_SERIALIZATION
template <class Archive>
void SearchShapeInput::serialize(Archive &ar, const unsigned int) {
  ar &boost::serialization::base_object<ShapeInput>(*this);
  ar & d_numDummies;
  ar & d_dummyVol;
  ar & d_actConf;
  ar & d_confCoords;
  ar & d_molConfs;
  ar & d_dummyVolumes;
  ar & d_shapeVolumes;
  ar & d_colorVolumes;
  ar & d_extremePointss;
  ar & d_canonRots;
  ar & d_centroids;
  ar & d_eigenValuess;
}
#endif

}  // namespace GaussianShape
}  // namespace RDKit

#endif  // SEARCHSHAPEINPUT_H
