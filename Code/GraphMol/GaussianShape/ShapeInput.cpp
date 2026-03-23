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
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
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

std::mutex mtx;

namespace RDKit {
namespace GaussianShape {

// Bondi radii
// You can find more of these in Table 12 of this publication:
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
// The dummy atom radius (atomic number 0) is set to
// 2.16.
const std::map<unsigned int, double> vdw_radii = {
    {0, 2.16},   // Dummy, same as Xe.
    {1, 1.10},   // H
    {2, 1.40},   // He
    {3, 1.81},   // Li
    {4, 1.53},   // Be
    {5, 1.92},   // B
    {6, 1.70},   // C
    {7, 1.55},   // N
    {8, 1.52},   // O
    {9, 1.47},   // F
    {10, 1.54},  // Ne
    {11, 2.27},  // Na
    {12, 1.73},  // Mg
    {13, 1.84},  // Al
    {14, 2.10},  // Si
    {15, 1.80},  // P
    {16, 1.80},  // S
    {17, 1.75},  // Cl
    {18, 1.88},  // Ar
    {19, 2.75},  // K
    {20, 2.31},  // Ca
    {31, 1.87},  // Ga
    {32, 2.11},  // Ge
    {33, 1.85},  // As
    {34, 1.90},  // Se
    {35, 1.83},  // Br
    {36, 2.02},  // Kr
    {37, 3.03},  // Rb
    {38, 2.49},  // Sr
    {49, 1.93},  // In
    {50, 2.17},  // Sn
    {51, 2.06},  // Sb
    {52, 2.06},  // Te
    {53, 1.98},  // I
    {54, 2.16},  // Xe
    {55, 3.43},  // Cs
    {56, 2.68},  // Ba
    {81, 1.96},  // Tl
    {82, 2.02},  // Pb
    {83, 2.07},  // Bi
    {84, 1.97},  // Po
    {85, 2.02},  // At
    {86, 2.20},  // Rn
    {87, 3.48},  // Fr
    {88, 2.83},  // Ra
};
constexpr double radius_color =
    1.08265;  // same radius for all feature/color "atoms", as used by the
              // PubChem code.

ShapeInput::ShapeInput(const ROMol &mol, int confId,
                       const ShapeInputOptions &opts,
                       const ShapeOverlayOptions &overlayOpts) {
  PRECONDITION(mol.getNumConformers() > 0,
               "ShapeInput object needs the molecule to have conformers.  " +
                   mol.getProp<std::string>("_Name") + "  " + MolToSmiles(mol));

  if (opts.allCarbonRadii && !opts.atomRadii.empty()) {
    BOOST_LOG(rdWarningLog)
        << "Specifying allCarbonRadii and providing custom atom radii doesn't"
           " make sense.  Ignoring the radii."
        << std::endl;
  }
  if (confId >= 0) {
    extractAtoms(mol, confId, opts, true);
    if (opts.useColors) {
      extractFeatures(mol, confId, opts, true);
    }
  } else {
    for (unsigned int i = 0; i < mol.getNumConformers(); ++i) {
      extractAtoms(mol, i, opts, i == 0);
      if (opts.useColors) {
        extractFeatures(mol, i, opts, i == 0);
      }
    }
  }
  calcNormalization();
  calcExtremes();
  calculateSelfOverlaps(overlayOpts);
  if (opts.shapePruneThreshold > 0.0 && mol.getNumConformers() > 1) {
    pruneShapes(opts.shapePruneThreshold);
  }
  sortShapesByVolumes();
  d_activeShape = 0;
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
      d_normalized(other.d_normalized),
      d_normalizationOK(other.d_normalizationOK),
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
  d_normalized = other.d_normalized;
  d_normalizationOK = other.d_normalizationOK;
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

void ShapeInput::setActiveShape(unsigned int newShape) {
  PRECONDITION(newShape < d_coords.size(),
               "Invalid conformation number (" + std::to_string(newShape) +
                   " vs " + std::to_string(d_coords.size()) + ").");
  if (d_activeShape != newShape) {
    d_activeShape = newShape;
    d_normalizationOK = false;
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

const std::array<double, 9> &ShapeInput::calcCanonicalRotation() {
  if (!d_normalizationOK) {
    calcNormalization();
  }
  return d_canonRots[d_activeShape];
}

const std::array<double, 3> &ShapeInput::calcCanonicalTranslation() {
  if (!d_normalizationOK) {
    calcNormalization();
  }
  return d_canonTranss[d_activeShape];
}

const std::array<double, 3> &ShapeInput::calcEigenValues() {
  if (!d_normalizationOK) {
    calcNormalization();
  }
  return d_eigenValuess[d_activeShape];
}

const std::array<size_t, 6> &ShapeInput::calcExtremes() {
  if (!d_normalizationOK) {
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
  if (d_normalized) {
    return;
  }
  if (!d_normalizationOK) {
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
  d_normalized = true;
  // Recalculate the extremes now we've changed the coordinates.
  calcExtremes();
}

void ShapeInput::transformCoords(RDGeom::Transform3D &xform) {
  applyTransformToShape(d_coords[d_activeShape], xform);
  d_normalized = false;
  d_normalizationOK = false;
}

std::unique_ptr<RWMol> ShapeInput::shapeToMol(bool includeColors) const {
  auto mol = std::make_unique<RWMol>();
  for (unsigned int i = 0; i < getNumAtoms(); i++) {
    Atom *atom = new Atom(6);
    mol->addAtom(atom, true, true);
  }
  if (includeColors) {
    for (unsigned int i = 0; i < getNumFeatures(); i++) {
      Atom *atom = new Atom(7);
      mol->addAtom(atom, true, true);
    }
  }
  unsigned int num = getNumAtoms();
  if (includeColors) {
    num += getNumFeatures();
  }
  Conformer *conf = new Conformer(num);
  const auto &shapeCds = getCoords();
  for (unsigned int i = 0, j = 0; i < 3 * num; i += 3, ++j) {
    auto &pos = conf->getAtomPos(j);
    pos.x = shapeCds[i];
    pos.y = shapeCds[i + 1];
    pos.z = shapeCds[i + 2];
  }
  mol->addConformer(conf, true);
  return mol;
}

namespace {
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
}  // namespace

double ShapeInput::bestSimilarity(ShapeInput &fitShape,
                                  unsigned int &bestThisShape,
                                  unsigned int &bestFitShape,
                                  RDGeom::Transform3D &bestXform,
                                  double threshold,
                                  const ShapeOverlayOptions &overlayOpts) {
  bestThisShape = -1;
  bestFitShape = -1;
  if (maxSimilarity(fitShape) < threshold) {
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
        std::cout << i << " -> " << j << " = " << scores[0] << ", " << scores[1]
                  << ", " << scores[2] << std::endl;
        if (scores[0] > bestSim) {
          bestSim = scores[0];
          bestThisShape = i;
          bestFitShape = j;
          bestXform = xform;
          std::cout << fabs(bestSim - 1.0) << std::endl;
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

double ShapeInput::maxSimilarity(const ShapeInput &fitShape,
                                 const ShapeOverlayOptions &overlayOpts) const {
  // The best score achievable is when the smaller volume is entirely inside
  // the larger volume.
  double maxSim = 0.0;
  for (unsigned int i = 0; i < getNumShapes(); i++) {
    for (unsigned int j = 0; j < fitShape.getNumShapes(); j++) {
      auto sim = maxScore(fitShape.d_selfOverlapShapeVols[i],
                          d_selfOverlapShapeVols[j],
                          fitShape.d_selfOverlapColorVols[i],
                          d_selfOverlapColorVols[j], overlayOpts);
      if (sim > maxSim) {
        maxSim = sim;
      }
    }
  }
  return maxSim;
}

namespace {
double getStandardAtomRadius(unsigned int atomicNum) {
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
    const auto atomIdx = atom->getIdx();
    if (!opts.atomSubset.empty()) {
      if (auto it = std::ranges::find_if(
              opts.atomSubset,
              [atomIdx](const auto &p) -> bool { return p == atomIdx; });
          it == opts.atomSubset.end()) {
        continue;
      }
    }
    if (atom->getAtomicNum() > 1) {
      auto atIdx = atom->getIdx();
      auto &pos = conf.getAtomPos(atIdx);
      theseCoords.push_back(pos.x);
      theseCoords.push_back(pos.y);
      theseCoords.push_back(pos.z);
      if (fillAlphas) {
        if (opts.allCarbonRadii) {
          d_alphas.push_back(KAPPA / (1.7 * 1.7));
        } else {
          double rad = 0.0;
          if (opts.atomRadii.empty()) {
            if (atom->getAtomicNum() == 6) {
              rad = 1.7;
              (*d_carbonRadii)[atomIdx] = true;
            } else {
              rad = getStandardAtomRadius(atom->getAtomicNum());
            }
          } else {
            auto it = std::ranges::find_if(
                opts.atomRadii,
                [atIdx](const auto &p) -> bool { return p.first == atIdx; });
            if (it == opts.atomRadii.end()) {
              rad = getStandardAtomRadius(atom->getAtomicNum());
            } else {
              rad = it->second;
            }
          }
          d_alphas.push_back(KAPPA / (rad * rad));
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
// for now.  Other options to be added later.
void ShapeInput::extractFeatures(const ROMol &mol, int confId,
                                 const ShapeInputOptions &opts,
                                 bool fillAlphas) {
  if (opts.customFeatures.empty()) {
    unsigned pattIdx = 1;
    const auto pattVects = getPh4Patterns();
    for (const auto &patts : *pattVects) {
      for (const auto &patt : patts) {
        std::vector<MatchVectType> matches;
        {
          // Once, recursive queries weren't thread safe but Greg assures
          // me that they now are.
          matches = SubstructMatch(mol, *patt);
        }
        for (const auto &match : matches) {
          std::vector<unsigned int> ats;
          bool featOk = true;
          for (const auto &pr : match) {
            // make sure all the atoms are in the subset, if there is one
            if (!opts.atomSubset.empty()) {
              if (std::ranges::find_if(
                      opts.atomSubset, [pr](const auto &p) -> bool {
                        return p == static_cast<unsigned int>(pr.second);
                      }) == opts.atomSubset.end()) {
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
          d_types.push_back(pattIdx);
          d_coords[d_activeShape].push_back(featPos.x);
          d_coords[d_activeShape].push_back(featPos.y);
          d_coords[d_activeShape].push_back(featPos.z);
          if (fillAlphas) {
            d_alphas.push_back(KAPPA / (radius_color * radius_color));
            d_numFeats++;
          }
        }
      }
      ++pattIdx;
    }
  } else {
    // Just copy them directly except that the last is a radius, so convert
    // to alpha.
    for (const auto &f : opts.customFeatures) {
      d_types.push_back(std::get<0>(f));
      d_numFeats++;
      const auto &pos = std::get<1>(f);
      d_coords[d_activeShape].push_back(pos.x);
      d_coords[d_activeShape].push_back(pos.y);
      d_coords[d_activeShape].push_back(pos.z);
      if (fillAlphas) {
        auto rad = std::get<2>(f);
        d_alphas.push_back(KAPPA / (rad * rad));
      }
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
  auto tmpMol = shapeToMol(false);
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
  d_normalizationOK = true;
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

void writeCoords(const std::vector<double> &shape, const std::string &label,
                 char lineEnd) {
  std::cout << label << " :: ";
  if (lineEnd == '\n') {
    std::cout << lineEnd;
  }
  for (size_t i = 0; i < shape.size(); i += 3) {
    std::cout << shape[i] << "," << shape[i + 1] << "," << shape[i + 2]
              << lineEnd;
  }
}

void writeCoords(const double *shape, unsigned int numPts,
                 const std::string &label, char lineEnd) {
  std::cout << label << " :: ";
  if (lineEnd == '\n') {
    std::cout << lineEnd;
  }
  for (unsigned int i = 0; i < numPts * 3; i += 3) {
    std::cout << shape[i] << "," << shape[i + 1] << "," << shape[i + 2]
              << lineEnd;
  }
}

void copyTransform(const RDGeom::Transform3D &src, RDGeom::Transform3D &dest) {
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      dest.setValUnchecked(i, j, src.getValUnchecked(i, j));
    }
  }
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
  std::vector<std::array<double, 9>> newCanRots;
  newCanRots.reserve(picks.size());
  std::vector<std::array<double, 3>> newCentroids;
  newCentroids.reserve(picks.size());
  std::vector<std::array<double, 3>> newEigenValuess;
  newEigenValuess.reserve(picks.size());
  for (auto p : picks) {
    newCoords.push_back(std::move(d_coords[p]));
    newShapeVolumes.push_back(d_selfOverlapShapeVols[p]);
    newColorVolumes.push_back(d_selfOverlapColorVols[p]);
    newExtremePointss.push_back(std::move(d_extremePointss[p]));
    newCanRots.push_back(std::move(d_canonRots[p]));
    newCentroids.push_back(std::move(d_canonTranss[p]));
    newEigenValuess.push_back(std::move(d_eigenValuess[p]));
  }
  d_coords = std::move(newCoords);
  d_selfOverlapShapeVols = std::move(newShapeVolumes);
  d_selfOverlapColorVols = std::move(newColorVolumes);
  d_extremePointss = std::move(newExtremePointss);
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

}  // namespace GaussianShape
}  // namespace RDKit
