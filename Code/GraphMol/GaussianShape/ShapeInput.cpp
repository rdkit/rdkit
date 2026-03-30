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

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/GaussianShape/GaussianShape.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <SimDivPickers/LeaderPicker.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

#ifdef RDK_HAS_EIGEN3
#include <Eigen/Dense>
#endif

namespace RDKit {
namespace GaussianShape {

// Bondi radii
// You can find more of these in Table 12 of this publication:
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
const std::map<unsigned int, double> vdw_radii = {
    {0, DUMMY_RAD},   // Dummy, same as Xe.
    {1, 1.10},        // H
    {2, 1.40},        // He
    {3, 1.81},        // Li
    {4, 1.53},        // Be
    {5, 1.92},        // B
    {6, CARBON_RAD},  // C
    {7, 1.55},        // N
    {8, 1.52},        // O
    {9, 1.47},        // F
    {10, 1.54},       // Ne
    {11, 2.27},       // Na
    {12, 1.73},       // Mg
    {13, 1.84},       // Al
    {14, 2.10},       // Si
    {15, 1.80},       // P
    {16, 1.80},       // S
    {17, 1.75},       // Cl
    {18, 1.88},       // Ar
    {19, 2.75},       // K
    {20, 2.31},       // Ca
    {31, 1.87},       // Ga
    {32, 2.11},       // Ge
    {33, 1.85},       // As
    {34, 1.90},       // Se
    {35, 1.83},       // Br
    {36, 2.02},       // Kr
    {37, 3.03},       // Rb
    {38, 2.49},       // Sr
    {49, 1.93},       // In
    {50, 2.17},       // Sn
    {51, 2.06},       // Sb
    {52, 2.06},       // Te
    {53, 1.98},       // I
    {54, 2.16},       // Xe
    {55, 3.43},       // Cs
    {56, 2.68},       // Ba
    {81, 1.96},       // Tl
    {82, 2.02},       // Pb
    {83, 2.07},       // Bi
    {84, 1.97},       // Po
    {85, 2.02},       // At
    {86, 2.20},       // Rn
    {87, 3.48},       // Fr
    {88, 2.83},       // Ra
};
constexpr double radius_color =
    1.08265;  // same radius for all feature/color "atoms", as used by the
// PubChem code.

namespace {
// Throw out any atoms we don't want.  copySubsetMol does more work than
// is necessary for this and seemed to leave the molecule in an odd state
// that didn't play well with extractAtoms.
std::unique_ptr<RWMol> extractSubset(
    const RWMol &mol, const std::vector<unsigned int> &atomsToKeep) {
  auto retMol = std::make_unique<RWMol>(mol);
  boost::dynamic_bitset<> keepAtoms(mol.getNumAtoms());
  for (const auto atk : atomsToKeep) {
    keepAtoms[atk] = true;
  }
  retMol->beginBatchEdit();
  for (auto atom : retMol->atoms()) {
    if (!keepAtoms[atom->getIdx()]) {
      retMol->removeAtom(atom->getIdx());
    }
  }
  retMol->commitBatchEdit();
  return retMol;
}
}  // namespace

ShapeInput::ShapeInput(const ROMol &mol, int confId,
                       const ShapeInputOptions &opts,
                       const ShapeOverlayOptions &overlayOpts) {
  PRECONDITION(mol.getNumConformers() > 0,
               "ShapeInput object needs the molecule to have conformers.  " +
                   mol.getProp<std::string>("_Name") + "  " + MolToSmiles(mol));
  std::unique_ptr<RWMol> tmpMol;
  // Subsetting the molecule makes any bespoke atom radii, identified
  // by atom index, incorrect so stash them as atom properties.
  for (const auto &[idx, radius] : opts.atomRadii) {
    auto atom = mol.getAtomWithIdx(idx);
    atom->setProp<double>("BespokeRadius", radius);
  }

  if (!opts.atomSubset.empty()) {
    tmpMol = extractSubset(mol, opts.atomSubset);
  } else {
    tmpMol.reset(new RWMol(mol));
  }
  // The input molecule may have been in Kekule form, so fix that so
  // aromatic features are found.
  MolOps::setAromaticity(*tmpMol);
  d_smiles = MolToSmiles(*tmpMol);
  std::vector<unsigned int> atOrder;
  tmpMol->getProp(common_properties::_smilesAtomOutputOrder, atOrder);
  tmpMol.reset(dynamic_cast<RWMol *>(MolOps::renumberAtoms(*tmpMol, atOrder)));

  if (!tmpMol->getRingInfo()->isInitialized()) {
    // Query molecules don't seem to have the ring info generated on creation.
    MolOps::findSSSR(*tmpMol);
  }

  d_normalizeds = boost::dynamic_bitset<>(tmpMol->getNumConformers());
  d_normalizationOKs = boost::dynamic_bitset<>(tmpMol->getNumConformers());

  auto processConf = [&](ROMol &m, unsigned int ci, bool fa) {
    extractAtoms(m, ci, opts, fa);
    if (opts.useColors) {
      extractFeatures(m, ci, opts, fa);
    }
  };

  if (confId >= 0) {
    processConf(*tmpMol, confId, true);
  } else {
    for (unsigned int i = 0; i < tmpMol->getNumConformers(); ++i) {
      processConf(*tmpMol, i, i == 0);
    }
  }
  d_activeShape = 0;
  calcNormalization();
  calcExtremes();
  calculateSelfOverlaps(overlayOpts);
  if (opts.shapePruneThreshold > 0.0 && tmpMol->getNumConformers() > 1) {
    pruneShapes(opts.shapePruneThreshold);
  }
  sortShapesByVolumes();
}

ShapeInput::ShapeInput(const ShapeInput &other, unsigned int shapeNum) {
  PRECONDITION(
      shapeNum < other.getNumShapes(),
      "Invalid shape number in makeSingleShape : " + std::to_string(shapeNum) +
          " vs " + std::to_string(other.getNumShapes()) + ".");
  d_activeShape = 0;
  d_coords.emplace_back(other.d_coords[shapeNum]);
  d_alphas = other.d_alphas;
  d_types = other.d_types;
  d_numAtoms = other.d_numAtoms;
  d_numFeats = other.d_numFeats;
  d_selfOverlapShapeVols.emplace_back(other.d_selfOverlapShapeVols[shapeNum]);
  d_selfOverlapColorVols.emplace_back(other.d_selfOverlapColorVols[shapeNum]);
  d_extremePointss.emplace_back(other.d_extremePointss[shapeNum]);
  if (other.d_carbonRadii) {
    d_carbonRadii.reset(new boost::dynamic_bitset<>(*other.d_carbonRadii));
  }
  d_smiles = other.d_smiles;
  d_normalizeds = boost::dynamic_bitset<>(1);
  d_normalizeds[0] = other.d_normalizeds[shapeNum];
  d_normalizationOKs = boost::dynamic_bitset<>(1);
  d_normalizationOKs[0] = other.d_normalizationOKs[shapeNum];
  d_canonRots.emplace_back(other.d_canonRots[shapeNum]);
  d_canonTranss.emplace_back(other.d_canonTranss[shapeNum]);
  d_eigenValuess.emplace_back(other.d_eigenValuess[shapeNum]);
}

ShapeInput::ShapeInput(const ShapeInput &other)
    : d_activeShape(other.d_activeShape),
      d_coords(other.d_coords),
      d_alphas(other.d_alphas),
      d_types(other.d_types),
      d_numAtoms(other.d_numAtoms),
      d_numFeats(other.d_numFeats),
      d_selfOverlapShapeVols(other.d_selfOverlapShapeVols),
      d_selfOverlapColorVols(other.d_selfOverlapColorVols),
      d_extremePointss(other.d_extremePointss),
      d_smiles(other.d_smiles),
      d_normalizeds(other.d_normalizeds),
      d_normalizationOKs(other.d_normalizationOKs),
      d_canonRots(other.d_canonRots),
      d_canonTranss(other.d_canonTranss),
      d_eigenValuess(other.d_eigenValuess) {
  if (other.d_carbonRadii) {
    d_carbonRadii.reset(new boost::dynamic_bitset<>(*other.d_carbonRadii));
  }
}

ShapeInput &ShapeInput::operator=(const ShapeInput &other) {
  if (this == &other) {
    return *this;
  }
  d_activeShape = other.d_activeShape;
  d_coords = other.d_coords;
  d_alphas = other.d_alphas;
  d_types = other.d_types;
  d_numAtoms = other.d_numAtoms;
  d_numFeats = other.d_numFeats;
  d_selfOverlapShapeVols = other.d_selfOverlapShapeVols;
  d_selfOverlapColorVols = other.d_selfOverlapColorVols;
  d_extremePointss = other.d_extremePointss;
  d_smiles = other.d_smiles;
  d_normalizeds = other.d_normalizeds;
  d_normalizationOKs = other.d_normalizationOKs;
  d_canonRots = other.d_canonRots;
  d_canonTranss = other.d_canonTranss;
  d_eigenValuess = other.d_eigenValuess;
  if (other.d_carbonRadii) {
    d_carbonRadii.reset(new boost::dynamic_bitset<>(*other.d_carbonRadii));
  } else {
    d_carbonRadii.reset();
  }
  return *this;
}

void ShapeInput::merge(ShapeInput &other) {
  if (!d_coords.empty() &&
      d_coords.front().size() != other.d_coords.front().size()) {
    BOOST_LOG(rdWarningLog) << "Can't merge shapes as different sizes.\n";
    return;
  }
  if (other.d_coords.empty()) {
    return;
  }
  d_coords.reserve(d_coords.size() + other.d_coords.size());
  d_coords.insert(d_coords.end(),
                  std::make_move_iterator(other.d_coords.begin()),
                  std::make_move_iterator(other.d_coords.end()));
  d_selfOverlapShapeVols.reserve(d_selfOverlapShapeVols.size() +
                                 other.d_selfOverlapShapeVols.size());
  d_selfOverlapShapeVols.insert(
      d_selfOverlapShapeVols.end(),
      std::make_move_iterator(other.d_selfOverlapShapeVols.begin()),
      std::make_move_iterator(other.d_selfOverlapShapeVols.end()));
  d_selfOverlapColorVols.reserve(d_selfOverlapColorVols.size() +
                                 other.d_selfOverlapColorVols.size());
  d_selfOverlapColorVols.insert(
      d_selfOverlapColorVols.end(),
      std::make_move_iterator(other.d_selfOverlapColorVols.begin()),
      std::make_move_iterator(other.d_selfOverlapColorVols.end()));
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

  other.d_coords.clear();
  other.d_selfOverlapShapeVols.clear();
  other.d_selfOverlapColorVols.clear();
  other.d_extremePointss.clear();
  other.d_canonRots.clear();
  other.d_canonTranss.clear();
  other.d_eigenValuess.clear();
}

void ShapeInput::setActiveShape(unsigned int newShape) {
  PRECONDITION(newShape < d_coords.size(),
               "Invalid shape number (" + std::to_string(newShape) + " vs " +
                   std::to_string(d_coords.size()) + ").");
  if (d_activeShape != newShape) {
    d_activeShape = newShape;
  }
}

std::vector<RDGeom::Point3D> ShapeInput::getAtomPoints(
    bool includeColors) const {
  std::vector<RDGeom::Point3D> atomPoints;
  unsigned int numPoints = getNumAtoms();
  if (includeColors) {
    numPoints += getNumFeatures();
  }
  atomPoints.reserve(numPoints);
  for (unsigned int i = 0; i < 3 * numPoints; i += 3) {
    atomPoints.emplace_back(RDGeom::Point3D(d_coords[d_activeShape][i],
                                            d_coords[d_activeShape][i + 1],
                                            d_coords[d_activeShape][i + 2]));
  }
  return atomPoints;
}

double ShapeInput::getShapeVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_coords.size(),
               "Invalid shape number (" + std::to_string(shapeNum) + " vs " +
                   std::to_string(d_coords.size()) + ").");
  return d_selfOverlapShapeVols[shapeNum];
}

double ShapeInput::getColorVolume(unsigned int shapeNum) const {
  PRECONDITION(shapeNum < d_coords.size(),
               "Invalid shape number (" + std::to_string(shapeNum) + " vs " +
                   std::to_string(d_coords.size()) + ").");
  return d_selfOverlapColorVols[shapeNum];
}

const std::array<double, 9> &ShapeInput::calcCanonicalRotation() {
  if (!d_normalizationOKs[d_activeShape]) {
    calcNormalization();
  }
  return d_canonRots[d_activeShape];
}

const std::array<double, 3> &ShapeInput::calcCanonicalTranslation() {
  if (!d_normalizationOKs[d_activeShape]) {
    calcNormalization();
  }
  return d_canonTranss[d_activeShape];
}

const std::array<double, 3> &ShapeInput::calcEigenValues() {
  if (!d_normalizationOKs[d_activeShape]) {
    calcNormalization();
  }
  return d_eigenValuess[d_activeShape];
}

const std::array<size_t, 6> &ShapeInput::calcExtremes() {
  if (!d_normalizationOKs[d_activeShape]) {
    calculateExtremes();
  }
  return d_extremePointss[d_activeShape];
}

std::array<double, 3> ShapeInput::calcMomentsOfInertia(
    bool includeColors) const {
  auto tmpMol = shapeToMol(includeColors);
  std::array<double, 3> eVals;
#if RDK_HAS_EIGEN3
  Eigen::Matrix3d axes;
  Eigen::Vector3d moments;
  MolTransforms::computePrincipalAxesAndMoments(tmpMol->getConformer(), axes,
                                                moments);
  eVals[0] = moments[0];
  eVals[1] = moments[1];
  eVals[2] = moments[2];
#else
  std::unique_ptr<RDGeom::Transform3D> canonXform(
      MolTransforms::computeCanonicalTransform(tmpMol->getConformer(), nullptr,
                                               false, true, eVals.data()));
#endif
  return eVals;
}

void ShapeInput::normalizeCoords() {
  if (d_normalizeds[d_activeShape]) {
    return;
  }
  if (!d_normalizationOKs[d_activeShape]) {
    calcNormalization();
  }
  RDGeom::Transform3D canonRot;
  for (unsigned int i = 0, k = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j, ++k) {
      canonRot.setValUnchecked(i, j, d_canonRots[d_activeShape][k]);
    }
  }
  RDGeom::Point3D trans{d_canonTranss[d_activeShape][0],
                        d_canonTranss[d_activeShape][1],
                        d_canonTranss[d_activeShape][2]};
  canonRot.TransformPoint(trans);
  canonRot.SetTranslation(trans);

  transformCoords(canonRot);
  d_normalizeds[d_activeShape] = true;
  // Recalculate the extremes now we've changed the coordinates.
  calcExtremes();
}

void ShapeInput::transformCoords(RDGeom::Transform3D &xform) {
  applyTransformToShape(d_coords[d_activeShape], xform);
  d_normalizeds[d_activeShape] = false;
  d_normalizationOKs[d_activeShape] = false;
}

std::unique_ptr<RWMol> ShapeInput::shapeToMol(bool includeColors,
                                              bool withBonds) const {
  std::unique_ptr<RWMol> mol;
  if (withBonds) {
    // The SMILES string and the atom coordinates should be in the same
    // order.
    v2::SmilesParse::SmilesParserParams params;
    params.sanitize = false;
    mol = v2::SmilesParse::MolFromSmiles(d_smiles, params);
  } else {
    mol.reset(new RWMol());
    for (unsigned int i = 0; i < getNumAtoms(); i++) {
      Atom *atom = new Atom(6);
      mol->addAtom(atom, true, true);
    }
  }
  if (includeColors) {
    for (unsigned int i = 0; i < getNumFeatures(); i++) {
      Atom *atom = new Atom(54);
      mol->addAtom(atom, true, true);
    }
  }
  unsigned int num = getNumAtoms();
  if (includeColors) {
    num += getNumFeatures();
  }
  Conformer *conf = new Conformer(num);
  const auto &shapeCds = getCoords();
  for (unsigned int i = 0; i < num; i++) {
    auto &pos = conf->getAtomPos(i);
    pos.x = shapeCds[3 * i];
    pos.y = shapeCds[3 * i + 1];
    pos.z = shapeCds[3 * i + 2];
  }
  mol->addConformer(conf, true);
  return mol;
}

double ShapeInput::bestSimilarity(ShapeInput &fitShape,
                                  unsigned int &bestThisShape,
                                  unsigned int &bestFitShape,
                                  RDGeom::Transform3D &bestXform,
                                  double threshold,
                                  const ShapeOverlayOptions &overlayOpts) {
  bestThisShape = -1;
  bestFitShape = -1;
  if (maxPossibleSimilarity(fitShape) < threshold) {
    return -1.0;
  }

  double bestSim = -1.0;
  RDGeom::Transform3D xform;
  for (size_t i = 0; i < getNumShapes(); i++) {
    setActiveShape(i);
    for (size_t j = 0; j < fitShape.getNumShapes(); j++) {
      fitShape.setActiveShape(j);
      auto maxSim =
          maxScore(getShapeVolume(), fitShape.getShapeVolume(),
                   getColorVolume(), fitShape.getColorVolume(), overlayOpts);

      if (maxSim > threshold) {
        auto scores = AlignShape(*this, fitShape, &xform, overlayOpts);
        if (scores[0] > bestSim) {
          bestSim = scores[0];
          bestThisShape = i;
          bestFitShape = j;
          bestXform = xform;
        }
      }
      // Floating point cruft means we sometimes get a similarity slightly
      // above 1.0.  1.0 is the maximum possible, so stop if we hit it.
      if (bestSim > 1.0 || fabs(bestSim - 1.0) < 1.0e-6) {
        return bestSim;
      }
    }
  }
  return bestSim;
}

double ShapeInput::maxPossibleSimilarity(
    const ShapeInput &fitShape, const ShapeOverlayOptions &overlayOpts) const {
  // The best score achievable is when the smaller volume is entirely inside
  // the larger volume.
  double maxSim = 0.0;
  for (unsigned int i = 0; i < getNumShapes(); i++) {
    for (unsigned int j = 0; j < fitShape.getNumShapes(); j++) {
      auto sim = maxScore(fitShape.d_selfOverlapShapeVols[j],
                          d_selfOverlapShapeVols[i],
                          fitShape.d_selfOverlapColorVols[j],
                          d_selfOverlapColorVols[i], overlayOpts);
      if (sim > maxSim) {
        maxSim = sim;
      }
    }
  }
  return maxSim;
}

namespace {
double getStandardAtomRadius(unsigned int atomicNum) {
  // Mostly they will be carbons, so just return that without lookup.
  if (atomicNum == 6) {
    return CARBON_RAD;
  }
  if (auto rad = vdw_radii.find(static_cast<unsigned int>(atomicNum));
      rad != vdw_radii.end()) {
    return rad->second;
  }
  throw ValueErrorException("No VdW radius for atom with Z=" +
                            std::to_string(atomicNum));
}

}  // namespace
void ShapeInput::extractAtoms(const ROMol &mol, int confId,
                              const ShapeInputOptions &opts, bool fillAlphas) {
  if (!opts.allCarbonRadii) {
    d_carbonRadii.reset(new boost::dynamic_bitset<>(
        !opts.atomSubset.empty() ? opts.atomSubset.size() : mol.getNumAtoms()));
  }
  auto conf = mol.getConformer(confId);
  std::vector<double> theseCoords;
  theseCoords.reserve(mol.getNumAtoms() * 3);
  for (const auto atom : mol.atoms()) {
    // Ignore H atoms except deuterium and tritium which are treated elsewhere
    // as explicit atoms, and do use dummies if requested.  The latter two
    // are set to be ignored in the volume calculation, however.  This is
    // just to keep the atom numberings correct.  Dummy atoms can also have
    // isotope number 1.
    if (atom->getAtomicNum() != 1 ||
        (atom->getAtomicNum() == 1 && atom->getIsotope() > 1)) {
      if (!opts.includeDummies && !atom->getAtomicNum()) {
        continue;
      }
      const auto atIdx = atom->getIdx();
      auto &pos = conf.getAtomPos(atIdx);
      theseCoords.push_back(pos.x);
      theseCoords.push_back(pos.y);
      theseCoords.push_back(pos.z);
      if (fillAlphas) {
        if (opts.allCarbonRadii && !d_carbonRadii) {
          d_alphas.push_back(KAPPA / (CARBON_RAD * CARBON_RAD));
        } else {
          double rad = 0.0;
          if (opts.atomRadii.empty()) {
            if ((opts.allCarbonRadii && atom->getAtomicNum()) ||
                atom->getAtomicNum() == 6) {
              (*d_carbonRadii)[atIdx] = true;
            }
            rad = getStandardAtomRadius(atom->getAtomicNum());
          } else {
            if (atom->hasProp("BespokeRadius")) {
              rad = atom->getProp<double>("BespokeRadius");
            } else {
              rad = getStandardAtomRadius(atom->getAtomicNum());
              if ((opts.allCarbonRadii && atom->getAtomicNum()) ||
                  atom->getAtomicNum() == 6) {
                (*d_carbonRadii)[atIdx] = true;
              }
            }
          }
          // TODO: For now, check there hasn't been a logic error leaving a rad
          // of 0.0.  Take it out later when it's been tried enough we can be
          // confident it's ok.
          if (rad == 0.0) {
            throw ValueErrorException("Atom has radius 0.0.");
          }
          d_alphas.push_back(KAPPA / (rad * rad));
        }
        if (atom->getAtomicNum() == 1) {
          // So it's not included in the volume calcs.
          d_alphas.back() = -1.0;
        }
      }
    }
  }
  d_activeShape = d_coords.size();
  d_coords.push_back(theseCoords);
  if (fillAlphas) {
    d_numAtoms = d_coords.front().size() / 3;
    d_types.resize(d_numAtoms);
    d_numFeats = 0;
  }
}

namespace {
class ss_matcher {
 public:
  ss_matcher(const std::string &pattern) : m_pattern(pattern) {
    m_needCopies = (pattern.find_first_of("$") != std::string::npos);
    RWMol *p = SmartsToMol(pattern);
    m_matcher = p;
    POSTCONDITION(m_matcher, "no matcher");
  };
  const ROMol *getMatcher() const { return m_matcher; };
  ~ss_matcher() { delete m_matcher; };

 private:
  ss_matcher() : m_pattern("") {};
  std::string m_pattern;
  bool m_needCopies{false};
  const ROMol *m_matcher{nullptr};
};
}  // namespace

// This came from the original PubChemShape.cpp
typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
                         boost::flyweights::no_tracking>
    pattern_flyweight;
// Definitions for feature points adapted from:
// Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998)
const std::vector<std::vector<std::string>> smartsPatterns = {
    {"[$([N;!H0;v3,v4&+1]),\
$([O,S;H1;+0]),\
n&H1&+0]"},                                // Donor
    {"[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),\
$([O,S;H0;v2]),\
$([O,S;-]),\
$([N;v3;!$(N-*=[O,N,P,S])]),\
n&H0&+0,\
$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]"},  // Acceptor
    {
        "[r]1[r][r]1",
        "[r]1[r][r][r]1",
        "[r]1[r][r][r][r]1",
        "[r]1[r][r][r][r][r]1",
        "[r]1[r][r][r][r][r][r]1",
    },  // rings
        //    "[a]",                                                  //
        //    Aromatic
        //    "[F,Cl,Br,I]",                                          // Halogen
    {"[#7;+,\
$([N;H2&+0][$([C,a]);!$([C,a](=O))]),\
$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),\
$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]"},  // Basic
    {"[$([C,S](=[O,S,P])-[O;H1,-1])]"}                       // Acidic
};
std::vector<std::vector<const ROMol *>> *getPh4Patterns() {
  static std::unique_ptr<std::vector<std::vector<const ROMol *>>> patterns;
  if (!patterns) {
    patterns.reset(new std::vector<std::vector<const ROMol *>>());
    for (const auto &smartsV : smartsPatterns) {
      std::vector<const ROMol *> v;
      for (const auto &smarts : smartsV) {
        const ROMol *matcher = pattern_flyweight(smarts).get().getMatcher();
        CHECK_INVARIANT(matcher, "bad smarts");
        v.push_back(matcher);
      }
      patterns->push_back(std::move(v));
    }
  }

  return patterns.get();
}

// Extract the features for the color scores, using RDKit pphore features
// for now.  Other options to be added later.  This should be working on
// the parent molecule, but only extracting features that are due to
// atoms entirely within any atom subset.
void ShapeInput::extractFeatures(const ROMol &mol, unsigned int confId,
                                 const ShapeInputOptions &opts,
                                 bool fillAlphas) {
  std::vector<CustomFeature> feats;
  if (opts.customFeatures.empty()) {
    findFeatures(mol, confId, feats);
  } else {
    PRECONDITION(confId < opts.customFeatures.size(),
                 "Conformer number " + std::to_string(confId) +
                     " too large for custom features size " +
                     std::to_string(opts.customFeatures.size()));
    feats = opts.customFeatures[confId];
  }
  // Just copy them directly except that the last is a radius, so convert
  // to alpha.
  for (const auto &f : feats) {
    const auto &pos = f.pos;
    d_coords[d_activeShape].push_back(pos.x);
    d_coords[d_activeShape].push_back(pos.y);
    d_coords[d_activeShape].push_back(pos.z);
    if (fillAlphas) {
      auto rad = f.rad;
      d_alphas.push_back(KAPPA / (rad * rad));
      d_types.push_back(f.type);
      d_numFeats++;
    }
  }
}

void ShapeInput::calcNormalization() {
  if (d_eigenValuess.size() < d_coords.size()) {
    d_eigenValuess.resize(d_coords.size());
    d_canonRots.resize(d_coords.size());
    d_canonTranss.resize(d_coords.size());
    d_eigenValuess.resize(d_coords.size());
    d_extremePointss.resize(d_coords.size());
  }
  d_eigenValuess[d_activeShape] = std::array<double, 3>{0.0, 0.0, 0.0};
  // Build a "molecule" just of the shape points, not the color features
  // with which to calculate the canonical transformation.  Doesn't ever
  // use the input molecule in case the shape was built from a subset of
  // atoms in that molecule.
  auto tmpMol = shapeToMol(false, false);
  std::unique_ptr<RDGeom::Transform3D> canonXform(
      MolTransforms::computeCanonicalTransform(
          tmpMol->getConformer(), nullptr, false, true,
          d_eigenValuess[d_activeShape].data()));
  d_canonRots[d_activeShape] =
      std::array<double, 9>{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  for (unsigned int i = 0, k = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j, ++k) {
      d_canonRots[d_activeShape][k] = canonXform->getValUnchecked(i, j);
    }
  }
  d_canonTranss[d_activeShape] = std::array<double, 3>{0.0, 0.0, 0.0};
  for (unsigned int i = 0; i < 3 * d_numAtoms; i += 3) {
    d_canonTranss[d_activeShape][0] -= d_coords[d_activeShape][i];
    d_canonTranss[d_activeShape][1] -= d_coords[d_activeShape][i + 1];
    d_canonTranss[d_activeShape][2] -= d_coords[d_activeShape][i + 2];
  }
  d_canonTranss[d_activeShape][0] /= d_numAtoms;
  d_canonTranss[d_activeShape][1] /= d_numAtoms;
  d_canonTranss[d_activeShape][2] /= d_numAtoms;
  d_normalizationOKs[d_activeShape] = true;
}

void ShapeInput::calculateExtremes() {
  if (d_extremePointss.size() < d_coords.size()) {
    d_extremePointss.resize(d_coords.size());
  }
  d_extremePointss[d_activeShape] = std::array<size_t, 6>{0, 0, 0, 0, 0, 0};
  const std::vector<double> &theseCoords = d_coords[d_activeShape];
  for (size_t i = 0, j = 0; i < theseCoords.size(); i += 3, ++j) {
    if (theseCoords[i] < theseCoords[3 * d_extremePointss[d_activeShape][0]]) {
      d_extremePointss[d_activeShape][0] = j;
    }
    if (theseCoords[i] > theseCoords[3 * d_extremePointss[d_activeShape][3]]) {
      d_extremePointss[d_activeShape][3] = j;
    }

    if (theseCoords[i + 1] <
        theseCoords[3 * d_extremePointss[d_activeShape][1] + 1]) {
      d_extremePointss[d_activeShape][1] = j;
    }
    if (theseCoords[i + 1] >
        theseCoords[3 * d_extremePointss[d_activeShape][4] + 1]) {
      d_extremePointss[d_activeShape][4] = j;
    }

    if (theseCoords[i + 2] <
        theseCoords[3 * d_extremePointss[d_activeShape][2] + 2]) {
      d_extremePointss[d_activeShape][2] = j;
    }
    if (theseCoords[i + 2] >
        theseCoords[3 * d_extremePointss[d_activeShape][5] + 2]) {
      d_extremePointss[d_activeShape][5] = j;
    }
  }
}

void ShapeInput::pruneShapes(double simThreshold) {
  if (d_coords.size() < 2) {
    return;
  }
  class DistFunctor {
   public:
    DistFunctor(ShapeInput &shapes) : d_shapes(shapes), d_shapesCp(shapes) {}
    ~DistFunctor() = default;
    double operator()(unsigned int i, unsigned int j) {
      d_shapes.setActiveShape(i);
      d_shapesCp.setActiveShape(j);
      auto scores = AlignShape(d_shapes, d_shapesCp);
      // The picker works on distances.
      return 1.0 - scores[0];
    }
    ShapeInput d_shapes;
    ShapeInput d_shapesCp;
  };
  RDPickers::LeaderPicker leaderPicker;
  DistFunctor distFunctor(*this);
  // The picker works on distances.
  auto picks = leaderPicker.lazyPick(distFunctor, d_coords.size(), 0,
                                     1.0 - simThreshold);
  // Allow for mysterious LeaderPicker behaviour where it returns a
  // vector of 0s when it should only pick 1 shape.
  std::ranges::sort(picks);
  auto [first, last] = std::ranges::unique(picks);
  picks.erase(first, last);
  selectConformations(picks);
}

void ShapeInput::selectConformations(const std::vector<int> &picks) {
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
  boost::dynamic_bitset<> newNormalizeds(picks.size());
  boost::dynamic_bitset<> newNormalizationOKs(picks.size());
  std::vector<std::array<double, 9>> newCanRots;
  newCanRots.reserve(picks.size());
  std::vector<std::array<double, 3>> newCentroids;
  newCentroids.reserve(picks.size());
  std::vector<std::array<double, 3>> newEigenValuess;
  newEigenValuess.reserve(picks.size());
  unsigned int i = 0;
  for (auto p : picks) {
    newCoords.push_back(std::move(d_coords[p]));
    newShapeVolumes.push_back(d_selfOverlapShapeVols[p]);
    newColorVolumes.push_back(d_selfOverlapColorVols[p]);
    newExtremePointss.push_back(std::move(d_extremePointss[p]));
    newNormalizeds[i] = d_normalizeds[p];
    newNormalizationOKs[i] = d_normalizationOKs[p];
    newCanRots.push_back(std::move(d_canonRots[p]));
    newCentroids.push_back(std::move(d_canonTranss[p]));
    newEigenValuess.push_back(std::move(d_eigenValuess[p]));
    ++i;
  }
  d_coords = std::move(newCoords);
  d_selfOverlapShapeVols = std::move(newShapeVolumes);
  d_selfOverlapColorVols = std::move(newColorVolumes);
  d_extremePointss = std::move(newExtremePointss);
  d_normalizeds = newNormalizeds;
  d_normalizationOKs = newNormalizationOKs;
  d_canonRots = std::move(newCanRots);
  d_canonTranss = std::move(newCentroids);
  d_eigenValuess = std::move(newEigenValuess);
}

void ShapeInput::calculateSelfOverlaps(const ShapeOverlayOptions &overlayOpts) {
  d_selfOverlapShapeVols.reserve(getNumShapes());
  d_selfOverlapColorVols.reserve(getNumShapes());
  std::vector<double> gradConverters(12 * (d_numAtoms + d_numFeats));
  for (unsigned int i = 0; i < getNumShapes(); i++) {
    d_selfOverlapShapeVols.push_back(calcVolAndGrads(
        d_coords[i].data(), d_alphas.data(), d_numAtoms, d_carbonRadii.get(),
        d_coords[i].data(), d_alphas.data(), d_numAtoms, d_carbonRadii.get(),
        gradConverters, overlayOpts.useDistCutoff,
        overlayOpts.distCutoff * overlayOpts.distCutoff));
    d_selfOverlapColorVols.push_back(calcVolAndGrads(
        d_coords[i].data() + 3 * d_numAtoms, d_alphas.data() + d_numAtoms,
        d_numFeats, d_types.data() + d_numAtoms,
        d_coords[i].data() + 3 * d_numAtoms, d_alphas.data() + d_numAtoms,
        d_numFeats, d_types.data() + d_numAtoms, d_numAtoms, gradConverters,
        overlayOpts.useDistCutoff,
        overlayOpts.distCutoff * overlayOpts.distCutoff, nullptr, nullptr));
  }
}

void ShapeInput::sortShapesByVolumes() {
  std::vector<std::pair<double, size_t>> vals;
  vals.reserve(d_coords.size());
  for (size_t i = 0; i < d_coords.size(); i++) {
    vals.push_back(std::make_pair(
        d_selfOverlapShapeVols[i] + d_selfOverlapColorVols[i], i));
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

void findFeatures(const ROMol &mol, int confId,
                  std::vector<CustomFeature> &features,
                  const std::vector<unsigned int> &atomSubset) {
  unsigned pattIdx = 1;
  const auto pattVects = getPh4Patterns();
  for (const auto &patts : *pattVects) {
    for (const auto &patt : patts) {
      std::vector<MatchVectType> matches;
      // Once, recursive queries weren't thread safe but Greg assures
      // me that they now are.
      matches = SubstructMatch(mol, *patt);
      for (const auto &match : matches) {
        std::vector<unsigned int> ats;
        bool featOk = true;
        for (const auto &pr : match) {
          // make sure all the atoms are in the subset, if there is one
          if (!atomSubset.empty()) {
            if (std::ranges::find_if(atomSubset, [pr](const auto &p) -> bool {
                  return p == static_cast<unsigned int>(pr.second);
                }) == atomSubset.end()) {
              featOk = false;
              break;
            }
          }
          ats.push_back(pr.second);
        }
        if (!featOk) {
          continue;
        }
        auto featPos = computeFeaturePos(mol, confId, ats);
        features.emplace_back(
            CustomFeature{pattIdx, featPos, radius_color, ats});
      }
    }
    ++pattIdx;
  }
}

RDGeom::Point3D computeFeaturePos(const ROMol &mol, int confId,
                                  const std::vector<unsigned int> &ats) {
  RDGeom::Point3D featPos;
  auto &conf = mol.getConformer(confId);
  for (const auto at : ats) {
    featPos += conf.getAtomPos(at);
  }
  featPos /= ats.size();
  return featPos;
}

void applyTransformToShape(std::vector<double> &shape,
                           RDGeom::Transform3D &xform) {
  for (size_t i = 0; i < shape.size(); i += 3) {
    RDGeom::Point3D pos{shape[i], shape[i + 1], shape[i + 2]};
    xform.TransformPoint(pos);
    shape[i] = pos.x;
    shape[i + 1] = pos.y;
    shape[i + 2] = pos.z;
  }
}

void applyTransformToShape(const double *inShape, double *outShape,
                           size_t numPoints, RDGeom::Transform3D &xform) {
  for (size_t i = 0; i < 3 * numPoints; i += 3) {
    RDGeom::Point3D pos{inShape[i], inShape[i + 1], inShape[i + 2]};
    xform.TransformPoint(pos);
    outShape[i] = pos.x;
    outShape[i + 1] = pos.y;
    outShape[i + 2] = pos.z;
  }
}

void translateShape(std::vector<double> &shape,
                    const RDGeom::Point3D &translation) {
  for (size_t i = 0; i < shape.size(); i += 3) {
    shape[i] += translation.x;
    shape[i + 1] += translation.y;
    shape[i + 2] += translation.z;
  }
}

void translateShape(const double *inShape, double *outShape, size_t numPoints,
                    const RDGeom::Point3D &translation) {
  for (size_t i = 0; i < 3 * numPoints; i += 3) {
    outShape[i] = inShape[i] + translation.x;
    outShape[i + 1] = inShape[i + 1] + translation.y;
    outShape[i + 2] = inShape[i + 2] + translation.z;
  }
}

// Maximum possible score of the 2 shape (v[12]) and color (c[12]) volumes
double maxScore(double v1, double v2, double c1, double c2,
                const ShapeOverlayOptions &overlayOpts) {
  // We're dealing with a Tversky score
  // s = O / (A * (R - O) + B * (F - O) + O)
  // There are 2 cases to handle, where v1 < v2 in which case the max overlap
  // is v1, and the opposite.
  auto maxPart = [](double p1, double p2,
                    const ShapeOverlayOptions &overlayOpts) -> double {
    if (p1 < p2) {
      return p1 / (overlayOpts.simBeta * (p2 - p1) + p1);
    }
    return p2 / (overlayOpts.simAlpha * (p1 - p2) + p2);
  };
  auto maxSt = maxPart(v1, v2, overlayOpts);
  auto maxCt = maxPart(c1, c2, overlayOpts);
  return maxSt * (1.0 - overlayOpts.optParam) + maxCt * overlayOpts.optParam;
};

}  // namespace GaussianShape
}  // namespace RDKit
