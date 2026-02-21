//
// Copyright (C) David Cosgrove 2025.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/Abbreviations/Abbreviations.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SearchShapeInput.h>
#include <SimDivPickers/LeaderPicker.h>

namespace RDKit {
namespace GaussianShape {
SearchShapeInput::SearchShapeInput(const ROMol &mol, double pruneThreshold,
                                   const ShapeInputOptions &opts,
                                   const ShapeOverlayOptions &overlayOpts)
    : ShapeInput(mol, 0, opts, overlayOpts) {
  // The ShapeInput c'tor puts a PRECONDITION on conformers being available.
  initializeFromBase();

  // build the shapes for conformations 1 onwards.
  d_confCoords.reserve(mol.getNumConformers());
  d_molConfs.reserve(mol.getNumConformers());
  d_dummyVolumes.reserve(mol.getNumConformers());
  d_shapeVolumes.reserve(mol.getNumConformers());
  d_colorVolumes.reserve(mol.getNumConformers());
  d_extremePointss.reserve(mol.getNumConformers());
  d_canonRots.reserve(mol.getNumConformers());
  d_centroids.reserve(mol.getNumConformers());
  d_eigenValuess.reserve(mol.getNumConformers());
  for (unsigned int cn = 0; cn < mol.getNumConformers(); ++cn) {
    auto shape = ShapeInput(mol, cn, opts, overlayOpts);
    const auto &cdsAndRads = shape.getCoords();
    std::vector<double> cdsCoords(3 * (getNumAtoms() + getNumFeatures()));
    for (unsigned int i = 0; i < getNumAtoms() + getNumFeatures(); ++i) {
      cdsCoords[3 * i] = cdsAndRads[4 * i];
      cdsCoords[3 * i + 1] = cdsAndRads[4 * i + 1];
      cdsCoords[3 * i + 1] = cdsAndRads[4 * i + 1];
    }
    d_confCoords.emplace_back(std::move(cdsCoords));
    d_molConfs.push_back(cn);
    d_dummyVolumes.push_back(-1);
    d_shapeVolumes.push_back(shape.getShapeVolume());
    d_colorVolumes.push_back(shape.getColorVolume());
    d_extremePointss.push_back(shape.getExtremes());
    d_canonRots.push_back(shape.getCanonicalRotation());
    d_centroids.push_back(shape.getCanonicalTranslation());
    d_eigenValuess.push_back(shape.getEigenValues());
  }
  pruneShapes(pruneThreshold);
  sortShapesByScore();
}

SearchShapeInput::SearchShapeInput(const std::string &str) {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  std::stringstream ss(str);
  boost::archive::text_iarchive ia(ss);
  ia &*this;
#endif
}

SearchShapeInput::SearchShapeInput(const ShapeInput &si) : ShapeInput(si) {
  initializeFromBase();
}

double SearchShapeInput::getShapeVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_confCoords.size(), "Bad shapeNum");
  return d_shapeVolumes[shapeNum];
}

double SearchShapeInput::getColorVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_confCoords.size(), "Bad shapeNum");
  return d_colorVolumes[shapeNum];
}

double SearchShapeInput::getDummyVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_confCoords.size(), "Bad shapeNum");
  return d_dummyVolumes[shapeNum];
}

void SearchShapeInput::merge(SearchShapeInput &other) {
  if (!d_confCoords.empty() &&
      d_confCoords.front().size() != other.d_confCoords.front().size()) {
    BOOST_LOG(rdWarningLog) << "Can't merge shapes as different sizes.\n";
    return;
  }
  if (other.d_confCoords.empty()) {
    return;
  }
  d_confCoords.reserve(d_confCoords.size() + other.d_confCoords.size());
  d_confCoords.insert(d_confCoords.end(),
                      std::make_move_iterator(other.d_confCoords.begin()),
                      std::make_move_iterator(other.d_confCoords.end()));
  d_molConfs.reserve(d_molConfs.size() + other.d_molConfs.size());
  d_molConfs.insert(d_molConfs.end(),
                    std::make_move_iterator(other.d_molConfs.begin()),
                    std::make_move_iterator(other.d_molConfs.end()));
  d_dummyVolumes.reserve(d_dummyVolumes.size() + other.d_dummyVolumes.size());
  d_dummyVolumes.insert(d_dummyVolumes.end(),
                        std::make_move_iterator(other.d_dummyVolumes.begin()),
                        std::make_move_iterator(other.d_dummyVolumes.end()));
  d_shapeVolumes.reserve(d_shapeVolumes.size() + other.d_shapeVolumes.size());
  d_shapeVolumes.insert(d_shapeVolumes.end(),
                        std::make_move_iterator(other.d_shapeVolumes.begin()),
                        std::make_move_iterator(other.d_shapeVolumes.end()));
  d_colorVolumes.reserve(d_colorVolumes.size() + other.d_colorVolumes.size());
  d_colorVolumes.insert(d_colorVolumes.end(),
                        std::make_move_iterator(other.d_colorVolumes.begin()),
                        std::make_move_iterator(other.d_colorVolumes.end()));
  d_extremePointss.reserve(d_extremePointss.size() +
                           other.d_extremePointss.size());
  d_extremePointss.insert(
      d_extremePointss.end(),
      std::make_move_iterator(other.d_extremePointss.begin()),
      std::make_move_iterator(other.d_extremePointss.end()));
  d_canonRots.reserve(d_canonRots.size() + other.d_canonRots.size());
  d_canonRots.insert(d_canonRots.end(),
                     std::make_move_iterator(other.d_canonRots.begin()),
                     std::make_move_iterator(other.d_canonRots.end()));
  d_centroids.reserve(d_centroids.size() + other.d_centroids.size());
  d_centroids.insert(d_centroids.end(),
                     std::make_move_iterator(other.d_centroids.begin()),
                     std::make_move_iterator(other.d_centroids.end()));
  d_eigenValuess.reserve(d_eigenValuess.size() + other.d_eigenValuess.size());
  d_eigenValuess.insert(d_eigenValuess.end(),
                        std::make_move_iterator(other.d_eigenValuess.begin()),
                        std::make_move_iterator(other.d_eigenValuess.end()));

  other.d_confCoords.clear();
  other.d_molConfs.clear();
  other.d_dummyVolumes.clear();
  other.d_shapeVolumes.clear();
  other.d_colorVolumes.clear();
  other.d_extremePointss.clear();
  other.d_canonRots.clear();
  other.d_centroids.clear();
  other.d_eigenValuess.clear();
}

ShapeInput SearchShapeInput::makeSingleShape(unsigned int shapeNum) const {
  if (shapeNum >= d_confCoords.size()) {
    shapeNum = 0;
  }
  SearchShapeInput cp(*this);
  cp.setActiveShape(shapeNum);
  auto baseString = cp.ShapeInput::toString();
  ShapeInput shape(baseString);
  return shape;
}

void SearchShapeInput::setActiveShape(unsigned int shapeNum) {
  PRECONDITION(shapeNum < d_confCoords.size(), "confNum is out of bounds");
  d_actConf = shapeNum;
  std::vector<double> coords(getCoords().size());
  confCoordsToShapeCoords(d_actConf, coords);
  setCoords(coords);
  d_dummyVol = d_dummyVolumes[shapeNum];
  setShapeVolume(d_shapeVolumes[shapeNum]);
  setColorVolume(d_colorVolumes[shapeNum]);
  setExtremes(d_extremePointss[shapeNum]);
  setCanonicalRotation(d_canonRots[shapeNum]);
  setCanonicalTranslation(d_centroids[shapeNum]);
  setEigenValues(d_eigenValuess[shapeNum]);
}

std::string SearchShapeInput::toString() const {
#ifndef RDK_USE_BOOST_SERIALIZATION
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  std::stringstream ss;
  boost::archive::text_oarchive oa(ss);
  oa &*this;
  return ss.str();
#endif
}

void SearchShapeInput::pruneShapes(double simThreshold) {
  if (d_confCoords.size() < 2) {
    return;
  }
  class DistFunctor {
   public:
    DistFunctor(const SearchShapeInput &shapes) : d_shapes(shapes) {
      d_shapei = std::make_unique<ShapeInput>(d_shapes.makeSingleShape(0));
      d_shapej = std::make_unique<ShapeInput>(d_shapes.makeSingleShape(1));
    }
    ~DistFunctor() = default;
    double operator()(unsigned int i, unsigned int j) {
      d_shapei->setCoords(d_shapes.d_confCoords[i]);
      d_shapej->setCoords(d_shapes.d_confCoords[j]);
      auto scores = AlignShape(*d_shapei, *d_shapej);
      double res = scores[0];
      return res;
    }
    const SearchShapeInput &d_shapes;
    std::unique_ptr<ShapeInput> d_shapei, d_shapej;
  };
  RDPickers::LeaderPicker leaderPicker;
  DistFunctor distFunctor(*this);
  auto picks = leaderPicker.lazyPick(distFunctor, d_confCoords.size(), 0,
                                     2.0 - simThreshold);
  // Allow for mysterious LeaderPicker behaviour where it returns a full vector
  // of 0s when it should only pick 1 shape.
  if (picks.size() == d_confCoords.size()) {
    std::ranges::sort(picks);
    auto [first, last] = std::ranges::unique(picks);
    picks.erase(first, last);
  }
  selectConformations(picks);
}

void SearchShapeInput::sortShapesByScore() {
  std::vector<std::pair<double, size_t>> vals;
  vals.reserve(d_confCoords.size());
  for (size_t i = 0; i < d_confCoords.size(); i++) {
    vals.push_back(std::make_pair(d_shapeVolumes[i] + d_colorVolumes[i], i));
  }
  std::ranges::sort(vals,
                    [](const std::pair<double, size_t> &a,
                       const std::pair<double, size_t> &b) -> bool {
                      return a.first > b.first;
                    });

  std::vector<int> picks(vals.size());
  std::ranges::transform(vals, begin(picks),
                         [](const std::pair<double, size_t> &a) -> int {
                           return static_cast<int>(a.second);
                         });
  selectConformations(picks);
}

double SearchShapeInput::bestSimilarity(
    SearchShapeInput &fitShape, unsigned int &bestRefConf,
    unsigned int &bestFitConf, double threshold,
    const ShapeOverlayOptions &overlayOpts) {
  // The best score achievable is when the smaller volume is entirely inside
  // the larger volume.  The Shape tanimoto is the fraction of volume in
  // common.  The scores for the different conformations are sorted in
  // descending order.  So try the smallest refShape score against the
  // largest fitShape score and vice versa.
  auto maxScore = [](double v1, double v2, double f1, double f2) -> double {
    double maxSt = std::min(v1, v2) / std::max(v1, v2);
    double maxCt = std::min(f1, f2) / std::max(f1, f2);
    return maxSt + maxCt;
  };
  if (maxScore(fitShape.d_shapeVolumes.front(), d_shapeVolumes.back(),
               fitShape.d_colorVolumes.front(),
               d_colorVolumes.back()) < threshold &&
      maxScore(fitShape.d_shapeVolumes.back(), d_shapeVolumes.front(),
               fitShape.d_colorVolumes.back(),
               d_colorVolumes.front()) < threshold) {
    return -1.0;
  }

  double best_combo_t = -1.0;
  for (size_t i = 0; i < getNumAtoms() + getNumFeatures(); i++) {
    setActiveShape(i);
    for (size_t j = 0; j < fitShape.getNumAtoms() + fitShape.getNumFeatures();
         j++) {
      auto maxSim = maxScore(d_shapeVolumes[i], fitShape.d_shapeVolumes[j],
                             d_colorVolumes[i], fitShape.d_colorVolumes[j]);
      if (maxSim > threshold) {
        fitShape.setActiveShape(j);
        auto scores = AlignShape(*this, fitShape, nullptr, overlayOpts);
        if (scores[0] > best_combo_t) {
          best_combo_t = scores[0];
          bestRefConf = i;
          bestFitConf = j;
        }
      }
    }
  }
  return best_combo_t;
}

std::unique_ptr<RWMol> shapeToMol(const ShapeInput &shape) {
  auto mol = std::make_unique<RWMol>();
  for (unsigned int i = 0; i < shape.getNumAtoms(); i++) {
    RDKit::Atom *atom = new RDKit::Atom(6);
    mol->addAtom(atom, true, true);
  }
  for (unsigned int i = 0; i < shape.getNumFeatures(); i++) {
    RDKit::Atom *atom = new RDKit::Atom(7);
    mol->addAtom(atom, true, true);
  }
  RDKit::Conformer *conf =
      new RDKit::Conformer(shape.getNumAtoms() + shape.getNumFeatures());
  const auto &shapeCds = shape.getCoords();
  for (unsigned int i = 0; i < shape.getNumAtoms() + shape.getNumFeatures();
       i++) {
    auto &pos = conf->getAtomPos(i);
    pos.x = shapeCds[4 * i];
    pos.y = shapeCds[4 * i + 1];
    pos.z = shapeCds[4 * i + 2];
  }
  mol->addConformer(conf, true);
  return mol;
}

void SearchShapeInput::initializeFromBase() {
  d_confCoords = std::vector<std::vector<double>>(
      1, std::vector<double>(3 * (getNumAtoms() + getNumFeatures())));
  const auto &cdsAndRads = getCoords();
  for (unsigned int i = 0; i < getNumAtoms() + getNumFeatures(); ++i) {
    d_confCoords[0][3 * i] = cdsAndRads[4 * i];
    d_confCoords[0][3 * i + 1] = cdsAndRads[4 * i + 1];
    d_confCoords[0][3 * i + 1] = cdsAndRads[4 * i + 1];
  }
  // We don't actually know what conformer it came from
  d_molConfs = std::vector<unsigned int>(1, 0);
  // Keep the dummyVols in synch with the other data, but
  // clearly flag it as not calculated.
  d_dummyVolumes = std::vector<double>(1, -1);
  d_shapeVolumes = std::vector<double>(1, getShapeVolume());
  d_colorVolumes = std::vector<double>(1, getColorVolume());
  d_extremePointss = std::vector<std::array<size_t, 6>>(1, getExtremes());
  d_canonRots = std::vector<std::array<double, 9>>(1, getCanonicalRotation());
  d_centroids =
      std::vector<std::array<double, 3>>(1, getCanonicalTranslation());
  d_eigenValuess = std::vector<std::array<double, 3>>(1, getEigenValues());
}

void SearchShapeInput::confCoordsToShapeCoords(
    unsigned int shapeNum, std::vector<double> &coords) const {
  coords.resize(getCoords().size());
  for (unsigned int i = 0; i < getNumAtoms() + getNumFeatures(); i++) {
    coords[4 * i] = d_confCoords[shapeNum][3 * i];
    coords[4 * i + 1] = d_confCoords[shapeNum][3 * i + 1];
    coords[4 * i + 2] = d_confCoords[shapeNum][3 * i + 2];
    coords[4 * i + 3] = getCoords()[4 * i + 3];
  }
}

void SearchShapeInput::selectConformations(const std::vector<int> &picks) {
  std::vector<std::vector<double>> newCoords;
  newCoords.reserve(picks.size());
  std::vector<unsigned int> newMolConfs;
  newMolConfs.reserve(picks.size());
  std::vector<double> newDummyVols;
  newDummyVols.reserve(picks.size());
  std::vector<double> newShapeVolumes;
  newShapeVolumes.reserve(picks.size());
  std::vector<double> newColorVolumes;
  newColorVolumes.reserve(picks.size());
  std::vector<std::array<size_t, 6>> newExtremePointss;
  newExtremePointss.reserve(picks.size());
  std::vector<std::array<double, 9>> newCanRots;
  newCanRots.reserve(picks.size());
  std::vector<std::array<double, 3>> newCentroids;
  newCentroids.reserve(picks.size());
  std::vector<std::array<double, 3>> newEigenValuess;
  newEigenValuess.reserve(picks.size());
  for (auto p : picks) {
    newCoords.push_back(std::move(d_confCoords[p]));
    newMolConfs.push_back(d_molConfs[p]);
    newDummyVols.push_back(d_dummyVolumes[p]);
    newShapeVolumes.push_back(d_shapeVolumes[p]);
    newColorVolumes.push_back(d_colorVolumes[p]);
    newExtremePointss.push_back(std::move(d_extremePointss[p]));
    newCanRots.push_back(std::move(d_canonRots[p]));
    newCentroids.push_back(std::move(d_centroids[p]));
    newEigenValuess.push_back(std::move(d_eigenValuess[p]));
  }
  d_confCoords = std::move(newCoords);
  d_molConfs = std::move(newMolConfs);
  d_dummyVolumes = std::move(newDummyVols);
  d_shapeVolumes = std::move(newShapeVolumes);
  d_colorVolumes = std::move(newColorVolumes);
  d_extremePointss = std::move(newExtremePointss);
  d_eigenValuess = std::move(newEigenValuess);
}
}  // namespace GaussianShape
}  // namespace RDKit
