//
// Copyright (C) David Cosgrove 2026.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/SynthonSpaceSearch/SynthonShapeInput.h>

namespace RDKit {
namespace SynthonSpaceSearch {
SynthonShapeInput::SynthonShapeInput(
    const ROMol &mol, int confId, const GaussianShape::ShapeInputOptions &opts,
    const GaussianShape::ShapeOverlayOptions &overlayOpts) {
  d_shapes = std::make_unique<GaussianShape::ShapeInput>(mol, confId, opts,
                                                         overlayOpts);
  buildDummyAtomsAndNbrs();
  calculateDummyVols(overlayOpts);
}

void SynthonShapeInput::merge(SynthonShapeInput &other) {
  d_shapes->merge(*other.d_shapes);
  d_dummyVolumes.reserve(d_dummyVolumes.size() + other.d_dummyVolumes.size());
  d_dummyVolumes.insert(d_dummyVolumes.end(),
                        std::make_move_iterator(other.d_dummyVolumes.begin()),
                        std::make_move_iterator(other.d_dummyVolumes.end()));
  other.d_dummyVolumes.clear();
  d_dummyAtomsAndNbrs.reserve(d_dummyAtomsAndNbrs.size() +
                              other.d_dummyAtomsAndNbrs.size());
  d_dummyAtomsAndNbrs.insert(
      d_dummyAtomsAndNbrs.end(),
      std::make_move_iterator(other.d_dummyAtomsAndNbrs.begin()),
      std::make_move_iterator(other.d_dummyAtomsAndNbrs.end()));
  other.d_dummyAtomsAndNbrs.clear();
}

GaussianShape::ShapeInput &SynthonShapeInput::getShapes() const {
  return *d_shapes;
}

double SynthonShapeInput::getDummyVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < getShapes().getNumShapes(),
               "Invalid shape number (" + std::to_string(shapeNum) + " vs " +
                   std::to_string(getShapes().getNumShapes()) +
                   ")"
                   ").");
  return d_dummyVolumes[shapeNum];
}

unsigned int SynthonShapeInput::getNumDummyAtoms() const {
  std::vector<unsigned int> dummyAtoms;
  std::ranges::transform(d_dummyAtomsAndNbrs, std::back_inserter(dummyAtoms),
                         [](const auto &p) -> unsigned int { return p.first; });
  std::ranges::sort(dummyAtoms);
  dummyAtoms.erase(std::unique(dummyAtoms.begin(), dummyAtoms.end()),
                   dummyAtoms.end());
  return dummyAtoms.size();
}

void SynthonShapeInput::buildDummyAtomsAndNbrs() {
  auto tmpMol = d_shapes->shapeToMol(false, true);
  for (auto atom : tmpMol->atoms()) {
    if (!atom->getAtomicNum()) {
      for (auto nbr : tmpMol->atomNeighbors(atom)) {
        d_dummyAtomsAndNbrs.emplace_back(atom->getIdx(), nbr->getIdx());
      }
    }
  }
}

void SynthonShapeInput::calculateDummyVols(
    const GaussianShape::ShapeOverlayOptions &opts) {
  d_dummyVolumes.reserve(d_shapes->getNumShapes());
  // Set the alphas for the dummy atoms to -1.0 so they are skipped in
  // the volume calculation.
  auto alphas = d_shapes->getAlphas();
  for (const auto &p : d_dummyAtomsAndNbrs) {
    alphas[p.first] = -1.0;
  }

  std::vector<double> gradConverters(
      12 * (d_shapes->getNumAtoms() + d_shapes->getNumFeatures()));
  for (unsigned int i = 0; i < d_shapes->getNumShapes(); i++) {
    d_shapes->setActiveShape(i);
    auto fullVol = d_shapes->getShapeVolume();
    auto partVol = GaussianShape::calcVolAndGrads(
        d_shapes->getCoords().data(), alphas.data(), d_shapes->getNumAtoms(),
        d_shapes->getCarbonRadii(), d_shapes->getCoords().data(), alphas.data(),
        d_shapes->getNumAtoms(), d_shapes->getCarbonRadii(), gradConverters,
        opts.useDistCutoff, opts.distCutoff * opts.distCutoff);
    d_dummyVolumes.emplace_back(fullVol - partVol);
  }
}
}  // namespace SynthonSpaceSearch
}  // namespace RDKit
