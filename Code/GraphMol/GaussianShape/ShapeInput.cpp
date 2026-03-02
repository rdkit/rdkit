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

#include <mutex>
#include <cmath>

#include <Geometry/point.h>
#include <Geometry/Transform3D.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/RWMol.h>
#include <GraphMol/GaussianShape/ShapeInput.h>
#include <GraphMol/GaussianShape/SingleConformerAlignment.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>

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
  extractAtoms(mol, confId, opts);
  if (opts.useColors) {
    extractFeatures(mol, confId, opts);
  }
  calcNormalization();
  calcExtremes();
  std::vector<double> gradConverters(12 * (d_numAtoms + d_numFeats));
  d_selfOverlapVol = calcVolAndGrads(
      d_coords.data(), d_numAtoms, d_carbonRadii, d_coords.data(), d_numAtoms,
      d_carbonRadii, gradConverters, overlayOpts.useDistCutoff,
      overlayOpts.distCutoff * overlayOpts.distCutoff);
  d_selfOverlapColor = calcVolAndGrads(
      d_coords.data() + 4 * d_numAtoms, d_numFeats, d_types.data() + d_numAtoms,
      d_coords.data() + 4 * d_numAtoms, d_numFeats, d_types.data() + d_numAtoms,
      d_numAtoms, gradConverters, overlayOpts.useDistCutoff,
      overlayOpts.distCutoff * overlayOpts.distCutoff, nullptr, nullptr);
}

ShapeInput::ShapeInput(const ShapeInput &other)
    : d_coords(other.d_coords),
      d_types(other.d_types),
      d_numAtoms(other.d_numAtoms),
      d_numFeats(other.d_numFeats),
      d_selfOverlapVol(other.d_selfOverlapVol),
      d_selfOverlapColor(other.d_selfOverlapColor),
      d_extremePoints(other.d_extremePoints),
      d_normalized(other.d_normalized),
      d_normalizationOK(other.d_normalizationOK),
      d_canonRot(other.d_canonRot),
      d_canonTrans(other.d_canonTrans),
      d_eigenValues(other.d_eigenValues) {
  if (other.d_carbonRadii) {
    d_carbonRadii.reset(new boost::dynamic_bitset<>(*other.d_carbonRadii));
  }
}

ShapeInput &ShapeInput::operator=(const ShapeInput &other) {
  if (this == &other) {
    return *this;
  }
  d_coords = other.d_coords;
  d_types = other.d_types;
  d_numAtoms = other.d_numAtoms;
  d_numFeats = other.d_numFeats;
  d_selfOverlapVol = other.d_selfOverlapVol;
  d_selfOverlapColor = other.d_selfOverlapColor;
  d_extremePoints = other.d_extremePoints;
  d_normalized = other.d_normalized;
  d_normalizationOK = other.d_normalizationOK;
  d_canonRot = other.d_canonRot;
  d_canonTrans = other.d_canonTrans;
  d_eigenValues = other.d_eigenValues;
  if (other.d_carbonRadii) {
    d_carbonRadii.reset(new boost::dynamic_bitset<>(*other.d_carbonRadii));
  } else {
    d_carbonRadii.reset();
  }
  d_canonRot = other.d_canonRot;
  d_canonTrans = other.d_canonTrans;
  return *this;
}

std::vector<RDGeom::Point3D> ShapeInput::getAtomPoints(
    bool includeColors) const {
  std::vector<RDGeom::Point3D> atomPoints;
  unsigned int numPoints = getNumAtoms();
  if (includeColors) {
    numPoints += getNumFeatures();
  }
  atomPoints.reserve(numPoints);
  for (unsigned int i = 0; i < 4 * numPoints; i += 4) {
    atomPoints.emplace_back(
        RDGeom::Point3D(d_coords[i], d_coords[i + 1], d_coords[i + 2]));
  }
  return atomPoints;
}

const std::array<double, 9> &ShapeInput::getCanonicalRotation() {
  if (!d_normalizationOK) {
    calcNormalization();
  }
  return d_canonRot;
}

const std::array<double, 3> &ShapeInput::getCanonicalTranslation() {
  if (!d_normalizationOK) {
    calcNormalization();
  }
  return d_canonTrans;
}

const std::array<double, 3> &ShapeInput::getEigenValues() {
  if (!d_normalizationOK) {
    calcNormalization();
  }
  return d_eigenValues;
}

const std::array<size_t, 6> &ShapeInput::getExtremes() {
  if (!d_normalizationOK) {
    calcNormalization();
  }
  return d_extremePoints;
}

std::array<double, 3> ShapeInput::getMomentsOfInertia(
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
  if (!d_normalizationOK) {
    calcNormalization();
  }
  RDGeom::Transform3D canonRot;
  for (unsigned int i = 0, k = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j, ++k) {
      canonRot.setValUnchecked(i, j, d_canonRot[k]);
    }
  }
  RDGeom::Point3D trans{d_canonTrans[0], d_canonTrans[1], d_canonTrans[2]};
  canonRot.TransformPoint(trans);
  canonRot.SetTranslation(trans);

  transformCoords(canonRot);
  d_normalized = true;
  // Recalculate the extremes now we've changed the coordinates.
  calcExtremes();
}

void ShapeInput::transformCoords(RDGeom::Transform3D &xform) {
  applyTransformToShape(d_coords, xform);
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
  for (unsigned int i = 0; i < num; i++) {
    auto &pos = conf->getAtomPos(i);
    pos.x = shapeCds[4 * i];
    pos.y = shapeCds[4 * i + 1];
    pos.z = shapeCds[4 * i + 2];
  }
  mol->addConformer(conf, true);
  return mol;
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
                              const ShapeInputOptions &opts) {
  d_coords.reserve(mol.getNumAtoms() * 4);
  if (!opts.allCarbonRadii) {
    d_carbonRadii.reset(new boost::dynamic_bitset<>(
        !opts.atomSubset.empty() ? opts.atomSubset.size() : mol.getNumAtoms()));
  }
  auto conf = mol.getConformer(confId);
  // Index of atoms that have been added to the shape.
  unsigned int idx = 0;
  for (const auto atom : mol.atoms()) {
    if (!opts.atomSubset.empty()) {
      const auto atomIdx = atom->getIdx();
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
      d_coords.push_back(pos.x);
      d_coords.push_back(pos.y);
      d_coords.push_back(pos.z);
      if (opts.allCarbonRadii) {
        d_coords.push_back(KAPPA / (1.7 * 1.7));
      } else {
        double rad = 0.0;
        if (opts.atomRadii.empty()) {
          if (atom->getAtomicNum() == 6) {
            rad = 1.7;
            (*d_carbonRadii)[idx] = true;
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
        d_coords.push_back(KAPPA / (rad * rad));
      }
    }
    ++idx;
  }
  d_numAtoms = d_coords.size() / 4;
  d_types.resize(d_numAtoms);
  d_numFeats = 0;
}

namespace {
class ss_matcher {
 public:
  ss_matcher(const std::string &pattern) : m_pattern(pattern) {
    m_needCopies = (pattern.find_first_of("$") != std::string::npos);
    RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    m_matcher = p;
    POSTCONDITION(m_matcher, "no matcher");
  };
  const RDKit::ROMol *getMatcher() const { return m_matcher; };
  unsigned int countMatches(const RDKit::ROMol &mol) const {
    PRECONDITION(m_matcher, "no matcher");
    std::vector<RDKit::MatchVectType> matches;
    // This is an ugly one. Recursive queries aren't thread safe.
    // Unfortunately we have to take a performance hit here in order
    // to guarantee thread safety
    if (m_needCopies) {
      const RDKit::ROMol nm(*(m_matcher), true);
      RDKit::SubstructMatch(mol, nm, matches);
    } else {
      const RDKit::ROMol &nm = *m_matcher;
      RDKit::SubstructMatch(mol, nm, matches);
    }
    return matches.size();
  }
  ~ss_matcher() { delete m_matcher; };

 private:
  ss_matcher() : m_pattern("") {};
  std::string m_pattern;
  bool m_needCopies{false};
  const RDKit::ROMol *m_matcher{nullptr};
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
                                 const ShapeInputOptions &opts) {
  if (opts.customFeatures.empty()) {
    unsigned pattIdx = 1;
    const auto pattVects = getPh4Patterns();
    for (const auto &patts : *pattVects) {
      for (const auto &patt : patts) {
        std::vector<MatchVectType> matches;
        {
          // recursive queries aren't thread safe.
          std::unique_lock<std::mutex> lock(mtx);
          matches = SubstructMatch(mol, *patt);
        }
        for (auto match : matches) {
          std::vector<unsigned int> ats;
          bool featOk = true;
          for (const auto &pr : match) {
            // make sure all the atoms are in the subset, if there is one
            if (!opts.atomSubset.empty()) {
              if (auto it = std::ranges::find_if(
                      opts.atomSubset,
                      [pr](const auto &p) -> bool {
                        return p == static_cast<unsigned int>(pr.second);
                      });
                  it == opts.atomSubset.end()) {
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
          d_coords.push_back(featPos.x);
          d_coords.push_back(featPos.y);
          d_coords.push_back(featPos.z);
          d_coords.push_back(KAPPA / (radius_color * radius_color));
          d_numFeats++;
        }
      }
      ++pattIdx;
    }
  } else {
    // Just copy them directly
    for (const auto &f : opts.customFeatures) {
      d_types.push_back(std::get<0>(f));
      d_numFeats++;
      const auto &pos = std::get<1>(f);
      d_coords.push_back(pos.x);
      d_coords.push_back(pos.y);
      d_coords.push_back(pos.z);
      d_coords.push_back(std::get<2>(f));
    }
  }
}

void ShapeInput::calcNormalization() {
  d_eigenValues = std::array<double, 3>{0.0, 0.0, 0.0};
  // Build a "molecule" just of the shape points, not the color features
  // with which to calculate the canonical transformation.  Doesn't ever
  // use the input molecule in case the shape was built from a subset of
  // atoms in that molecule.
  auto tmpMol = shapeToMol(false);
  std::unique_ptr<RDGeom::Transform3D> canonXform(
      MolTransforms::computeCanonicalTransform(
          tmpMol->getConformer(), nullptr, false, true, d_eigenValues.data()));
  d_canonRot =
      std::array<double, 9>{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  for (unsigned int i = 0, k = 0; i < 3; ++i) {
    for (unsigned int j = 0; j < 3; ++j, ++k) {
      d_canonRot[k] = canonXform->getValUnchecked(i, j);
    }
  }
  d_canonTrans = std::array<double, 3>{0.0, 0.0, 0.0};
  for (unsigned int i = 0; i < 4 * d_numAtoms; i += 4) {
    d_canonTrans[0] -= d_coords[i];
    d_canonTrans[1] -= d_coords[i + 1];
    d_canonTrans[2] -= d_coords[i + 2];
  }
  d_canonTrans[0] /= d_numAtoms;
  d_canonTrans[1] /= d_numAtoms;
  d_canonTrans[2] /= d_numAtoms;
  d_normalizationOK = true;
}

void ShapeInput::calcExtremes() {
  d_extremePoints = std::array<size_t, 6>{0, 0, 0, 0, 0, 0};
  for (size_t i = 0, j = 0; i < d_coords.size(); i += 4, ++j) {
    if (d_coords[i] < d_coords[4 * d_extremePoints[0]]) {
      d_extremePoints[0] = j;
    }
    if (d_coords[i] > d_coords[4 * d_extremePoints[3]]) {
      d_extremePoints[3] = j;
    }

    if (d_coords[i + 1] < d_coords[4 * d_extremePoints[1] + 1]) {
      d_extremePoints[1] = j;
    }
    if (d_coords[i + 1] > d_coords[4 * d_extremePoints[4] + 1]) {
      d_extremePoints[4] = j;
    }

    if (d_coords[i + 2] < d_coords[4 * d_extremePoints[2] + 2]) {
      d_extremePoints[2] = j;
    }
    if (d_coords[i + 2] > d_coords[4 * d_extremePoints[5] + 2]) {
      d_extremePoints[5] = j;
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
  for (size_t i = 0; i < shape.size(); i += 4) {
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
  for (unsigned int i = 0; i < numPts * 4; i += 4) {
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
  for (size_t i = 0; i < shape.size(); i += 4) {
    RDGeom::Point3D pos{shape[i], shape[i + 1], shape[i + 2]};
    xform.TransformPoint(pos);
    shape[i] = pos.x;
    shape[i + 1] = pos.y;
    shape[i + 2] = pos.z;
  }
}

void applyTransformToShape(const double *inShape, double *outShape,
                           size_t numPoints, RDGeom::Transform3D &xform) {
  for (size_t i = 0; i < 4 * numPoints; i += 4) {
    RDGeom::Point3D pos{inShape[i], inShape[i + 1], inShape[i + 2]};
    xform.TransformPoint(pos);
    outShape[i] = pos.x;
    outShape[i + 1] = pos.y;
    outShape[i + 2] = pos.z;
    outShape[i + 3] = inShape[i + 3];
  }
}

void translateShape(std::vector<double> &shape,
                    const RDGeom::Point3D &translation) {
  for (size_t i = 0; i < shape.size(); i += 4) {
    shape[i] += translation.x;
    shape[i + 1] += translation.y;
    shape[i + 2] += translation.z;
  }
}

void translateShape(const double *inShape, double *outShape, size_t numPoints,
                    const RDGeom::Point3D &translation) {
  for (size_t i = 0; i < 4 * numPoints; i += 4) {
    outShape[i] = inShape[i] + translation.x;
    outShape[i + 1] = inShape[i + 1] + translation.y;
    outShape[i + 2] = inShape[i + 2] + translation.z;
    outShape[i + 3] = inShape[i + 3];
  }
}

}  // namespace GaussianShape
}  // namespace RDKit
