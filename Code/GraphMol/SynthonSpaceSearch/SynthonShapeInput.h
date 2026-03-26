//
// Copyright (C) David Cosgrove 2026.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// A class building on ShapeInput for Synthon shape searching.

#ifndef RDKIT_SYNTHONSHAPEINPUT_H
#define RDKIT_SYNTHONSHAPEINPUT_H

#include <GraphMol/RWMol.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/ShapeOverlayOptions.h>

#include <RDGeneral/BoostStartInclude.h>
#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/utility.hpp>  // for std::pair
#include <boost/serialization/vector.hpp>
#include <boost/serialization/split_member.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace SynthonSpaceSearch {

class SynthonShapeInput {
 public:
  SynthonShapeInput(const ROMol &mol, int confId = -1,
                    const GaussianShape::ShapeInputOptions &opts =
                        GaussianShape::ShapeInputOptions(),
                    const GaussianShape::ShapeOverlayOptions &overlayOpts =
                        GaussianShape::ShapeOverlayOptions());
  SynthonShapeInput(const std::string &str) {
#ifndef RDK_USE_BOOST_SERIALIZATION
    PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
    std::stringstream ss(str);
    boost::archive::text_iarchive ia(ss);
    ia &*this;
#endif
  }
  SynthonShapeInput(const SynthonShapeInput &) = default;
  SynthonShapeInput(SynthonShapeInput &&) = default;
  SynthonShapeInput &operator=(const SynthonShapeInput &) = default;
  SynthonShapeInput &operator=(SynthonShapeInput &&) = default;
  ~SynthonShapeInput() = default;

  // Merge the other SynthonShapeInput, assuming it has the correct number
  // of atoms etc.  Empties other, unless they can't be merged in which case
  // it returns unscathed.
  void merge(SynthonShapeInput &other);

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

  GaussianShape::ShapeInput &getShapes() const;

  double getDummyVolume(unsigned int shapeNum) const;
  unsigned int getNumDummyAtoms() const;
  unsigned int getNumDummyAtomNbors() const {
    return d_dummyAtomsAndNbrs.size();
  }
  const std::vector<std::pair<unsigned int, unsigned int>> &
  getDummyAtomsAndNbrs() const {
    return d_dummyAtomsAndNbrs;
  }

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int) {
    if (Archive::is_saving::value) {
      std::string shape = d_shapes->toString();
      ar & shape;
    } else {
      std::string shape;
      ar & shape;
      d_shapes.reset(new GaussianShape::ShapeInput(shape));
    }
    ar & d_dummyVolumes;
    ar & d_dummyAtomsAndNbrs;
  }
#endif

 private:
  std::unique_ptr<GaussianShape::ShapeInput> d_shapes;
  std::vector<double>
      d_dummyVolumes;  // The volumes of the dummy atoms in the shapes
  // These pairs are the dummy atoms and neighbour atoms.  On rare occasions
  // a dummy may have 2 neighbours, so the first numbers in the pairs
  // may not be unique.
  std::vector<std::pair<unsigned int, unsigned int>> d_dummyAtomsAndNbrs;

  void buildDummyAtomsAndNbrs();
  void calculateDummyVols(const GaussianShape::ShapeOverlayOptions &opts);
};

}  // namespace SynthonSpaceSearch
}  // namespace RDKit

#endif  // RDKIT_SYNTHONSHAPEINPUT_H
