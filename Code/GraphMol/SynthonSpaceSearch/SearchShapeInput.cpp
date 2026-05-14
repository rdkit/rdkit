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
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SynthonSpaceSearch/SearchShapeInput.h>

#include <GraphMol/PeriodicTable.h>
#include <GraphMol/MolTransforms/MolTransforms.h>

#include <SimDivPickers/LeaderPicker.h>
#include <boost/locale/date_time.hpp>

namespace RDKit {
namespace GaussianShape {

SearchShapeInput::SearchShapeInput(const ROMol &mol, double pruneThreshold,
                                   const ShapeInputOptions &opts,
                                   const ShapeOverlayOptions &overlayOpts)
    : ShapeInput(mol, 0, opts, overlayOpts) {
  PRECONDITION(pruneThreshold <= 1.0,
               "Shape prune threshold can't be above 1.0");
  // The ShapeInput c'tor puts a PRECONDITION on conformers being available.
  initializeFromBase();

  // Flag the dummy atoms.  The coords are in smilesAtomOutputOrder of
  // the molecule built from getSmiles() so the indices should be correct.
  v2::SmilesParse::SmilesParserParams smilesParams;
  smilesParams.sanitize = false;
  auto bitMol = v2::SmilesParse::MolFromSmiles(getSmiles(), smilesParams);
  d_dummyAtoms = boost::dynamic_bitset<>(getNumAtoms() + getNumFeatures());
  for (auto atom : bitMol->atoms()) {
    if (atom->getAtomicNum() == 0) {
      d_dummyAtoms[atom->getIdx()] = true;
    }
  }
  // build the shapes for conformations 1 onwards.
  d_confCoords.reserve(mol.getNumConformers());
  d_molConfs.reserve(mol.getNumConformers());
  d_dummyVolumes.reserve(mol.getNumConformers());
  d_shapeVolumes.reserve(mol.getNumConformers());
  d_colorVolumes.reserve(mol.getNumConformers());
  d_extremePointss.reserve(mol.getNumConformers());
  d_canonRots.reserve(mol.getNumConformers());
  d_canonTranss.reserve(mol.getNumConformers());
  d_eigenValuess.reserve(mol.getNumConformers());
  for (unsigned int cn = 1; cn < mol.getNumConformers(); ++cn) {
    auto shape = ShapeInput(mol, cn, opts, overlayOpts);
    const auto &cdsAndRads = shape.getCoords();
    std::vector<double> cdsCoords(3 * (getNumAtoms() + getNumFeatures()));
    for (unsigned int i = 0; i < getNumAtoms() + getNumFeatures(); ++i) {
      cdsCoords[3 * i] = cdsAndRads[4 * i];
      cdsCoords[3 * i + 1] = cdsAndRads[4 * i + 1];
      cdsCoords[3 * i + 2] = cdsAndRads[4 * i + 2];
    }
    d_confCoords.emplace_back(std::move(cdsCoords));
    d_molConfs.push_back(cn);
    d_dummyVolumes.push_back(-1);
    d_shapeVolumes.push_back(shape.getShapeVolume());
    d_colorVolumes.push_back(shape.getColorVolume());
    d_extremePointss.push_back(shape.calcExtremes());
    d_canonRots.push_back(shape.calcCanonicalRotation());
    d_canonTranss.push_back(shape.calcCanonicalTranslation());
    d_eigenValuess.push_back(shape.calcEigenValues());
  }
  if (pruneThreshold > 0.0 && mol.getNumConformers() > 1) {
    pruneShapes(pruneThreshold);
  }
  sortShapesByScore();
  calculateDummyVolumes(overlayOpts);
  for (unsigned int i = 0; i < getNumShapes(); ++i) {
    setActiveShape(i);
  }
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

unsigned int SearchShapeInput::getMolConf(unsigned int shapeNum) {
  PRECONDITION(shapeNum < d_confCoords.size(),
               "Bad shapeNum in getMolConf : " + std::to_string(shapeNum) +
                   " vs " + std::to_string(d_confCoords.size()));
  return d_molConfs[shapeNum];
}

double SearchShapeInput::getShapeVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_confCoords.size(),
               "Bad shapeNum in getShapeVolume : " + std::to_string(shapeNum) +
                   " vs " + std::to_string(d_confCoords.size()));
  return d_shapeVolumes[shapeNum];
}

double SearchShapeInput::getColorVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_confCoords.size(),
               "Bad shapeNum in getColorVolume : " + std::to_string(shapeNum) +
                   " vs " + std::to_string(d_confCoords.size()));
  return d_colorVolumes[shapeNum];
}

double SearchShapeInput::getDummyVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_confCoords.size(),
               "Bad shapeNum in getDummyVolume : " + std::to_string(shapeNum) +
                   " vs " + std::to_string(d_confCoords.size()));
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
  d_canonTranss.reserve(d_canonTranss.size() + other.d_canonTranss.size());
  d_canonTranss.insert(d_canonTranss.end(),
                       std::make_move_iterator(other.d_canonTranss.begin()),
                       std::make_move_iterator(other.d_canonTranss.end()));
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
  other.d_canonTranss.clear();
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
  PRECONDITION(shapeNum < d_confCoords.size(),
               "confNum " + std::to_string(shapeNum) + " is out of bounds vs " +
                   std::to_string(d_confCoords.size()) + ".");
  d_actConf = shapeNum;
  std::vector<double> coords(getCoords().size());
  confCoordsToShapeCoords(d_actConf, coords);
  setCoords(coords);
  d_dummyVol = d_dummyVolumes[shapeNum];
  setShapeVolume(d_shapeVolumes[shapeNum]);
  setColorVolume(d_colorVolumes[shapeNum]);
  setExtremes(d_extremePointss[shapeNum]);
  setCanonicalRotation(d_canonRots[shapeNum]);
  setCanonicalTranslation(d_canonTranss[shapeNum]);
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
      d_shapei->setShapeVolume(d_shapes.d_shapeVolumes[i]);
      d_shapei->setColorVolume(d_shapes.d_colorVolumes[i]);
      d_shapej->setCoords(d_shapes.d_confCoords[j]);
      d_shapej->setShapeVolume(d_shapes.d_shapeVolumes[j]);
      d_shapej->setColorVolume(d_shapes.d_colorVolumes[j]);
      auto scores = AlignShape(*d_shapei, *d_shapej);
      // The picker works on distances.
      return 1.0 - scores[0];
    }
    const SearchShapeInput &d_shapes;
    std::unique_ptr<ShapeInput> d_shapei, d_shapej;
  };
  RDPickers::LeaderPicker leaderPicker;
  DistFunctor distFunctor(*this);
  // The picker works on distances.
  auto picks = leaderPicker.lazyPick(distFunctor, d_confCoords.size(), 0,
                                     1.0 - simThreshold);
  // Allow for mysterious LeaderPicker behaviour where it returns a
  // vector of 0s when it should only pick 1 shape.
  std::ranges::sort(picks);
  auto [first, last] = std::ranges::unique(picks);
  picks.erase(first, last);
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

namespace {
// Maximum possible score of the 2 shape (v[12]) and color (c[12]) volumes
double maxScore(double v1, double v2, double c1, double c2) {
  double maxSt = std::min(v1, v2) / std::max(v1, v2);
  double maxCt = std::min(c1, c2) / std::max(c1, c2);
  return maxSt + maxCt;
};

}  // namespace
double SearchShapeInput::bestSimilarity(
    SearchShapeInput &fitShape, unsigned int &bestThisShape,
    unsigned int &bestFitShape, RDGeom::Transform3D &bestXform,
    double threshold, const ShapeOverlayOptions &overlayOpts) {
  if (maxSimilarity(fitShape) < threshold) {
    return -1.0;
  }

  double best_combo_t = -1.0;
  RDGeom::Transform3D xform;
  for (size_t i = 0; i < getNumShapes(); i++) {
    setActiveShape(i);
    for (size_t j = 0; j < fitShape.getNumShapes(); j++) {
      auto maxSim = maxScore(d_shapeVolumes[i], fitShape.d_shapeVolumes[j],
                             d_colorVolumes[i], fitShape.d_colorVolumes[j]);
      if (maxSim > threshold) {
        fitShape.setActiveShape(j);
        auto scores = AlignShape(*this, fitShape, &xform, overlayOpts);
        if (scores[0] > best_combo_t) {
          best_combo_t = scores[0];
          bestThisShape = i;
          bestFitShape = j;
          bestXform = xform;
        }
      }
    }
  }
  return best_combo_t;
}

double SearchShapeInput::maxSimilarity(const SearchShapeInput &fitShape) const {
  return std::max(
      maxScore(fitShape.d_shapeVolumes.front(), d_shapeVolumes.back(),
               fitShape.d_colorVolumes.front(), d_colorVolumes.back()),
      maxScore(fitShape.d_shapeVolumes.back(), d_shapeVolumes.front(),
               fitShape.d_colorVolumes.back(), d_colorVolumes.front()));
}

const double *SearchShapeInput::getDummyCoords(unsigned int dummyNumber) const {
  PRECONDITION(dummyNumber < d_dummyAtoms.count(), "Bad dummy atom number.");
  for (unsigned int i = 0, j = 0; i < d_dummyAtoms.size(); i++) {
    if (d_dummyAtoms[i]) {
      if (j == dummyNumber) {
        return getCoords().data() + 4 * i;
      }
      ++j;
    }
  }
  // To keep the compiler happy.  It should never get here because dummyNumber
  // should always be valid.
  return nullptr;
}

void SearchShapeInput::getDummyAndNbr(unsigned int dummyNumber,
                                      unsigned int &dummyIdx,
                                      unsigned int &nbrIdx) const {
  PRECONDITION(dummyNumber < d_dummyAtoms.count(), "Bad dummy atom number.");
  // The coords are in smilesAtomOutputOrder of the molecule built from
  // getSmiles() so the indices should be correct.
  dummyIdx = -1;
  nbrIdx = -1;
  for (unsigned int i = 0, j = 0; i < d_dummyAtoms.size(); i++) {
    if (d_dummyAtoms[i]) {
      if (j == dummyNumber) {
        dummyIdx = i;
        v2::SmilesParse::SmilesParserParams smilesParams;
        smilesParams.sanitize = false;
        auto bitMol = v2::SmilesParse::MolFromSmiles(getSmiles(), smilesParams);
        auto dummy = bitMol->getAtomWithIdx(dummyIdx);
        PRECONDITION(!dummy->getAtomicNum(), "Bad dummy atom number.");
        for (auto nbr : bitMol->atomNeighbors(dummy)) {
          // There should only be 1 neighbour.
          nbrIdx = nbr->getIdx();
          return;
        }
      }
      ++j;
    }
  }
}

void SearchShapeInput::initializeFromBase() {
  d_confCoords = std::vector<std::vector<double>>(
      1, std::vector<double>(3 * (getNumAtoms() + getNumFeatures())));
  const auto &cdsAndRads = getCoords();
  for (unsigned int i = 0; i < getNumAtoms() + getNumFeatures(); ++i) {
    d_confCoords[0][3 * i] = cdsAndRads[4 * i];
    d_confCoords[0][3 * i + 1] = cdsAndRads[4 * i + 1];
    d_confCoords[0][3 * i + 2] = cdsAndRads[4 * i + 2];
  }
  // We don't actually know what conformer it came from
  d_molConfs = std::vector<unsigned int>(1, 0);
  // Keep the dummyVols in synch with the other data, but
  // clearly flag it as not calculated.
  d_dummyVol = -1.0;
  d_dummyVolumes = std::vector<double>(1, -1);
  d_shapeVolumes = std::vector<double>(1, getShapeVolume());
  d_colorVolumes = std::vector<double>(1, getColorVolume());
  d_extremePointss = std::vector<std::array<size_t, 6>>(1, calcExtremes());
  d_canonRots = std::vector<std::array<double, 9>>(1, calcCanonicalRotation());
  d_canonTranss =
      std::vector<std::array<double, 3>>(1, calcCanonicalTranslation());
  d_eigenValuess = std::vector<std::array<double, 3>>(1, calcEigenValues());
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
    newCentroids.push_back(std::move(d_canonTranss[p]));
    newEigenValuess.push_back(std::move(d_eigenValuess[p]));
  }
  d_confCoords = std::move(newCoords);
  d_molConfs = std::move(newMolConfs);
  d_dummyVolumes = std::move(newDummyVols);
  d_shapeVolumes = std::move(newShapeVolumes);
  d_colorVolumes = std::move(newColorVolumes);
  d_extremePointss = std::move(newExtremePointss);
  d_canonRots = std::move(newCanRots);
  d_canonTranss = std::move(newCentroids);
  d_eigenValuess = std::move(newEigenValuess);
}

void SearchShapeInput::calculateDummyVolumes(
    const ShapeOverlayOptions &overlayOpts) {
  // Set the rad value of the coords to -ve for the dummy atoms, to
  // flag them as to be skipped.
  auto tmpCoords = getCoords();
  for (unsigned int i = 0; i < getNumShapes(); ++i) {
    confCoordsToShapeCoords(i, tmpCoords);
    // We need to do this every time, as the radius part is taken from
    // the base shape each time.
    for (unsigned int i = 0; i < getNumAtoms(); i++) {
      if (d_dummyAtoms[i]) {
        tmpCoords[4 * i + 3] *= -1;
      }
    }
    std::vector<double> gradConverters(12 * (getNumAtoms() + getNumFeatures()));

    auto vol =
        calcVolAndGrads(tmpCoords.data(), getNumAtoms(), getCarbonRadii(),
                        tmpCoords.data(), getNumAtoms(), getCarbonRadii(),
                        gradConverters, overlayOpts.useDistCutoff,
                        overlayOpts.distCutoff * overlayOpts.distCutoff);
    d_dummyVolumes[i] = d_shapeVolumes[i] - vol;
  }
}

}  // namespace GaussianShape
}  // namespace RDKit
