#include <ranges>
#include <sstream>
#include <stdexcept>

#include <GraphMol/RWMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include "pubchem-align3d/shape_functions.hpp"
#include "PubChemShape.hpp"

constexpr auto pubchemFeatureName = "PUBCHEM_PHARMACOPHORE_FEATURES";

// #define DEBUG_MSG(msg_stream) cout << msg_stream << '\n'
#define DEBUG_MSG(msg_stream)

using namespace RDKit;

// Bondi radii
//  can find more of these in Table 12 of this publication:
//   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3658832/
// The dummy atom radius (atomic number 0) is set to
// 2.16 in ShapeInputOptions and may be varied there, as
// may all the other radii if required, including the
// addition of atoms not covered here.
const std::map<unsigned int, double> vdw_radii = {
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
    1.08265;  // same radius for all feature/color "atoms"

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

namespace {
std::vector<std::pair<std::vector<unsigned int>, unsigned int>> extractFeatures(
    const ROMol &mol, const ShapeInputOptions &shapeOpts) {
  // unpack features (PubChem-specific property from SDF)
  // NOTE: this unpacking assumes that RWMol-atom-index = SDF-atom-number - 1
  //   e.g. RWMol uses [0..N-1] and SDF uses [1..N], with atoms in the same
  //   order
  // If there are no PubChem features, falls back on RDKit pphore types.

  std::vector<std::pair<std::vector<unsigned int>, unsigned int>>
      feature_idx_type;
  if (shapeOpts.useColors) {
    std::string features;
    if (mol.getPropIfPresent(pubchemFeatureName, features)) {
      // regular atoms have type 0; feature "atoms" (features represented by a
      // single point+radius) must have type > 0
      static const std::map<std::string, unsigned int> atomTypes = {
          {"acceptor", 1}, {"anion", 2},      {"cation", 3},
          {"donor", 4},    {"hydrophobe", 5}, {"rings", 6},
      };

      std::istringstream iss(features);
      std::string line;
      unsigned int n = 0;
      while (std::getline(iss, line)) {
        if (n == 0) {
          feature_idx_type.resize(stoul(line));
        }

        else {
          unsigned int f = n - 1;
          if (f >= feature_idx_type.size()) {
            throw ValueErrorException("Too many features");
          }

          std::istringstream iss2(line);
          std::string token;
          unsigned int t = 0;
          while (std::getline(iss2, token, ' ')) {
            if (t == 0) {
              feature_idx_type[f].first.resize(stoul(token));
            } else if (t <= feature_idx_type[f].first.size()) {
              feature_idx_type[f].first[t - 1] = stoul(token) - 1;
            } else {
              auto type = atomTypes.find(token);
              if (type == atomTypes.end()) {
                throw ValueErrorException("Invalid feature type: " + token);
              }
              feature_idx_type[f].second = type->second;
            }
            ++t;
          }
          if (t != (feature_idx_type[f].first.size() + 2)) {
            throw ValueErrorException("Wrong number of tokens in feature");
          }
        }

        ++n;
      }
      if (n != (feature_idx_type.size() + 1)) {
        throw ValueErrorException("Wrong number of features");
      }

      DEBUG_MSG("# features: " << feature_idx_type.size());
    } else {
      const auto pattVects = getPh4Patterns();
      feature_idx_type.clear();

      unsigned pattIdx = 1;
      for (const auto &patts : *pattVects) {
        for (const auto patt : patts) {
          auto matches = SubstructMatch(mol, *patt);
          for (auto match : matches) {
            std::vector<unsigned int> ats;
            for (const auto &pr : match) {
              ats.push_back(pr.second);
            }
            feature_idx_type.emplace_back(ats, pattIdx);
          }
        }
        ++pattIdx;
      }
    }
  }
  return feature_idx_type;
}

bool atomInSubset(unsigned int atomIdx, const ShapeInputOptions &shapeOpts) {
  if (shapeOpts.atomSubset.empty()) {
    return true;
  }
  return std::ranges::find(shapeOpts.atomSubset, atomIdx) !=
         shapeOpts.atomSubset.end();
}
bool atomAllowedInColor(unsigned int atomIdx,
                        const ShapeInputOptions &shapeOpts) {
  return std::ranges::find(shapeOpts.notColorAtoms, atomIdx) ==
         shapeOpts.notColorAtoms.end();
}

double getAtomRadius(unsigned int atomIdx, const ShapeInputOptions &shapeOpts) {
  auto it = std::ranges::find_if(
      shapeOpts.atomRadii,
      [atomIdx](const auto &p) -> bool { return p.first == atomIdx; });
  return it == shapeOpts.atomRadii.end() ? -1.0 : it->second;
}

// Get the atom radii.  rad_vector is expected to be big enough to hold them
// all.  Also computes the average coordinates of the selected atoms.
void extractAtomRadii(const Conformer &conformer, unsigned int nAtoms,
                      const ShapeInputOptions &shapeOpts, RDGeom::Point3D &ave,
                      unsigned int &nSelectedAtoms,
                      std::vector<double> &rad_vector) {
  nSelectedAtoms = 0;
  for (unsigned int i = 0u; i < nAtoms; ++i) {
    if (!atomInSubset(i, shapeOpts)) {
      continue;
    }
    double rad = getAtomRadius(i, shapeOpts);
    if (rad > 0.0) {
      rad_vector[nSelectedAtoms++] = rad;
      ave += conformer.getAtomPos(i);
    } else {
      unsigned int Z =
          conformer.getOwningMol().getAtomWithIdx(i)->getAtomicNum();
      if (Z > 1) {
        ave += conformer.getAtomPos(i);
        if (auto rad = vdw_radii.find(Z); rad != vdw_radii.end()) {
          rad_vector[nSelectedAtoms++] = rad->second;
        } else {
          throw ValueErrorException("No VdW radius for atom with Z=" +
                                    std::to_string(Z));
        }
      } else if (shapeOpts.includeDummies && Z == 0) {
        ave += conformer.getAtomPos(i);
        rad_vector[nSelectedAtoms++] = shapeOpts.dummyRadius;
      }
    }
  }
  ave /= nSelectedAtoms;
}

void extractAtomCoords(const Conformer &conformer, const unsigned int nAtoms,
                       const ShapeInputOptions &shapeOpts,
                       const RDGeom::Point3D &ave, std::vector<float> &coords) {
  for (unsigned i = 0, j = 0; i < nAtoms; ++i) {
    if (!atomInSubset(i, shapeOpts)) {
      continue;
    }
    // use only non-H for alignment, optionally with dummy atoms.
    unsigned int Z = conformer.getOwningMol().getAtomWithIdx(i)->getAtomicNum();
    if (Z > 1 || (shapeOpts.includeDummies && Z == 0)) {
      RDGeom::Point3D pos = conformer.getAtomPos(i);
      if (shapeOpts.normalize) {
        pos -= ave;
      }
      coords[j * 3] = pos.x;
      coords[(j * 3) + 1] = pos.y;
      coords[(j * 3) + 2] = pos.z;
      ++j;
    }
  }
}

void extractFeatureCoords(
    const Conformer &conformer, const unsigned int nAtoms,
    const unsigned int nSelectedAtoms,
    const std::vector<std::pair<std::vector<unsigned int>, unsigned int>>
        &feature_idx_type,
    const ShapeInputOptions &shapeOpts, const RDGeom::Point3D &ave,
    unsigned int &numFeatures, ShapeInput &res,
    std::vector<double> &rad_vector) {
  // get feature coordinates - simply the average of coords of all atoms in the
  // feature
  for (unsigned i = 0; i < feature_idx_type.size(); ++i) {
    RDGeom::Point3D floc;
    unsigned int nSel = 0;
    for (unsigned int j = 0; j < feature_idx_type[i].first.size(); ++j) {
      unsigned int idx = feature_idx_type[i].first[j];
      if (!atomInSubset(idx, shapeOpts) ||
          !atomAllowedInColor(idx, shapeOpts)) {
        continue;
      }
      if (idx >= nAtoms ||
          conformer.getOwningMol().getAtomWithIdx(idx)->getAtomicNum() <= 1) {
        throw ValueErrorException("Invalid feature atom index");
      }
      floc += conformer.getAtomPos(idx);
      ++nSel;
    }
    if (nSel == feature_idx_type[i].first.size()) {
      floc /= nSel;
      if (shapeOpts.normalize) {
        floc -= ave;
      }
      DEBUG_MSG("feature type " << feature_idx_type[i].second << " (" << floc
                                << ")");

      auto array_idx = nSelectedAtoms + numFeatures;
      res.coord[array_idx * 3] = floc.x;
      res.coord[(array_idx * 3) + 1] = floc.y;
      res.coord[(array_idx * 3) + 2] = floc.z;
      rad_vector[array_idx] = radius_color;
      res.atom_type_vector[array_idx] = feature_idx_type[i].second;
      ++numFeatures;
    }
  }
}

void extractCustomFeatureCoords(const unsigned int nSelectedAtoms,
                                const ShapeInputOptions &shapeOpts,
                                const RDGeom::Point3D &ave,
                                unsigned int &numFeatures, ShapeInput &res,
                                std::vector<double> &rad_vector) {
  for (const auto &feature : shapeOpts.customFeatures) {
    unsigned int feature_type = std::get<0>(feature);
    RDGeom::Point3D floc = std::get<1>(feature);
    double radius = std::get<2>(feature);
    floc -= ave;
    DEBUG_MSG("custom feature type " << feature_type << " (" << floc << ")");

    auto array_idx = nSelectedAtoms + numFeatures;
    res.coord[array_idx * 3] = floc.x;
    res.coord[(array_idx * 3) + 1] = floc.y;
    res.coord[(array_idx * 3) + 2] = floc.z;
    rad_vector[array_idx] = radius;
    res.atom_type_vector[array_idx] = feature_type;
    ++numFeatures;
  }
}

}  // namespace
// The conformer is left where it is, the shape is translated to the origin.
ShapeInput PrepareConformer(const ROMol &mol, int confId,
                            const ShapeInputOptions &shapeOpts) {
  Align3D::setUseCutOff(true);
  ShapeInput res;

  auto feature_idx_type = extractFeatures(mol, shapeOpts);

  auto &conformer = mol.getConformer(confId);
  if (!conformer.is3D()) {
    throw ValueErrorException("Conformer must be 3D");
  }
  unsigned int nAtoms = mol.getNumAtoms();
  // DEBUG_MSG("num atoms: " << nAtoms);

  // Start with the arrays as large as they will possibly have to be.
  // They will be re-sized later.
  unsigned int nAlignmentAtoms = nAtoms;
  if (shapeOpts.customFeatures.empty()) {
    nAlignmentAtoms += feature_idx_type.size();
  } else {
    nAlignmentAtoms += shapeOpts.customFeatures.size();
  }
  std::vector<double> rad_vector(nAlignmentAtoms);
  res.atom_type_vector.resize(nAlignmentAtoms, 0);

  RDGeom::Point3D ave;
  unsigned int nSelectedAtoms = 0;
  extractAtomRadii(conformer, nAtoms, shapeOpts, ave, nSelectedAtoms,
                   rad_vector);

  // translate steric center to origin
  DEBUG_MSG("steric center: (" << ave << ")");
  res.shift = {-ave.x, -ave.y, -ave.z};
  res.coord.resize(nAlignmentAtoms * 3);

  extractAtomCoords(conformer, nAtoms, shapeOpts, ave, res.coord);
  unsigned int numFeatures = 0;
  if (shapeOpts.customFeatures.empty()) {
    extractFeatureCoords(conformer, nAtoms, nSelectedAtoms, feature_idx_type,
                         shapeOpts, ave, numFeatures, res, rad_vector);
  } else if (shapeOpts.useColors) {
    extractCustomFeatureCoords(nSelectedAtoms, shapeOpts, ave, numFeatures, res,
                               rad_vector);
  }

  // Now cut the final vectors down to the actual number of atoms and
  // features used.
  nAlignmentAtoms = nSelectedAtoms + numFeatures;
  res.coord.resize(nAlignmentAtoms * 3);
  rad_vector.resize(nAlignmentAtoms);
  res.atom_type_vector.resize(nAlignmentAtoms);
  Align3D::setAlpha(rad_vector.data(), rad_vector.size(), res.alpha_vector);

  // regular atom self overlap
  Align3D::getVolumeAtomIndexVector(res.atom_type_vector.data(),
                                    res.atom_type_vector.size(),
                                    res.volumeAtomIndexVector);
  res.sov = Align3D::ComputeShapeOverlap(
      res.coord.data(), res.alpha_vector, res.volumeAtomIndexVector,
      res.coord.data(), res.alpha_vector, res.volumeAtomIndexVector);
  DEBUG_MSG("sov: " << res.sov);

  // feature self overlap
  if (feature_idx_type.size() > 0) {
    Align3D::getColorAtomType2IndexVectorMap(res.atom_type_vector.data(),
                                             res.atom_type_vector.size(),
                                             res.colorAtomType2IndexVectorMap);
    res.sof = Align3D::ComputeFeatureOverlap(
        res.coord.data(), res.alpha_vector, res.colorAtomType2IndexVectorMap,
        res.coord.data(), res.alpha_vector, res.colorAtomType2IndexVectorMap);
    DEBUG_MSG("sof: " << res.sof);
  }
  return res;
}

std::pair<double, double> AlignShape(const ShapeInput &refShape,
                                     ShapeInput &fitShape,
                                     std::vector<float> &matrix,
                                     double opt_param,
                                     unsigned int max_preiters,
                                     unsigned int max_postiters) {
  std::set<unsigned int> jointColorAtomTypeSet;
  Align3D::getJointColorTypeSet(
      refShape.atom_type_vector.data(), refShape.atom_type_vector.size(),
      fitShape.atom_type_vector.data(), fitShape.atom_type_vector.size(),
      jointColorAtomTypeSet);
  auto mapCp = refShape.colorAtomType2IndexVectorMap;
  Align3D::restrictColorAtomType2IndexVectorMap(mapCp, jointColorAtomTypeSet);
  // Take copy of the color atom mappings so as not to alter the input shape
  // which might be re-used.
  auto fitMapCp = fitShape.colorAtomType2IndexVectorMap;
  Align3D::restrictColorAtomType2IndexVectorMap(fitMapCp,
                                                jointColorAtomTypeSet);

  DEBUG_MSG("Running alignment...");
  double nbr_st = 0.0;
  double nbr_ct = 0.0;
  Align3D::Neighbor_Conformers(
      refShape.coord.data(), refShape.alpha_vector,
      refShape.volumeAtomIndexVector, mapCp, refShape.sov, refShape.sof,
      fitShape.coord.data(), fitShape.alpha_vector,
      fitShape.volumeAtomIndexVector, fitMapCp, fitShape.sov, fitShape.sof,
      !jointColorAtomTypeSet.empty(), max_preiters, max_postiters, opt_param,
      matrix.data(), nbr_st, nbr_ct);

  DEBUG_MSG("Done!");
  DEBUG_MSG("nbr_st: " << nbr_st);
  DEBUG_MSG("nbr_ct: " << nbr_ct);

  std::vector<float> transformed(fitShape.coord.size());
  Align3D::VApplyRotTransMatrix(transformed.data(), fitShape.coord.data(),
                                fitShape.coord.size() / 3, matrix.data());
  fitShape.coord = transformed;
  return std::make_pair(nbr_st, nbr_ct);
}

void TransformConformer(const std::vector<double> &finalTrans,
                        const std::vector<float> &matrix, ShapeInput &fitShape,
                        Conformer &fitConf) {
  // we reuse/modify the coord member of fitShape in order to avoid memory
  // allocations
  const unsigned int nAtoms = fitConf.getOwningMol().getNumAtoms();
  if (nAtoms > fitShape.volumeAtomIndexVector.size()) {
    // Hs are missing... make sure we have space for them
    fitShape.coord.resize(3 * nAtoms);
  }
  // initialize the fitShape coords with the starting atomic positions from
  // the conformer shifted to the center of "mass" coordinates.
  for (unsigned int i = 0; i < nAtoms; ++i) {
    const auto &pos = fitConf.getAtomPos(i);
    fitShape.coord[i * 3] = pos.x + fitShape.shift[0];
    fitShape.coord[i * 3 + 1] = pos.y + fitShape.shift[1];
    fitShape.coord[i * 3 + 2] = pos.z + fitShape.shift[2];
  }

  std::vector<float> transformed(nAtoms * 3);
  Align3D::VApplyRotTransMatrix(transformed.data(), fitShape.coord.data(),
                                nAtoms, matrix.data());

  // now set the coordinates in the conformer; undo the shift of the reference
  // shape to center of "mass" coordinates
  for (unsigned i = 0; i < nAtoms; ++i) {
    RDGeom::Point3D &pos = fitConf.getAtomPos(i);
    pos.x = transformed[i * 3] - finalTrans[0];
    pos.y = transformed[(i * 3) + 1] - finalTrans[1];
    pos.z = transformed[(i * 3) + 2] - finalTrans[2];
  }
}

std::pair<double, double> AlignMolecule(
    const ShapeInput &refShape, ROMol &fit, std::vector<float> &matrix,
    const ShapeInputOptions &shapeOpts, int fitConfId, double opt_param,
    unsigned int max_preiters, unsigned int max_postiters, bool applyRefShift) {
  PRECONDITION(matrix.size() == 12, "bad matrix size");
  Align3D::setUseCutOff(true);

  DEBUG_MSG("Fit details:");
  auto fitShape = PrepareConformer(fit, fitConfId, shapeOpts);
  auto tanis = AlignShape(refShape, fitShape, matrix, opt_param, max_preiters,
                          max_postiters);

  // transform fit coords
  Conformer &fit_conformer = fit.getConformer(fitConfId);
  std::vector<double> finalTrans{0.0, 0.0, 0.0};
  if (applyRefShift) {
    finalTrans = refShape.shift;
  }
  TransformConformer(finalTrans, matrix, fitShape, fit_conformer);
  fit.setProp("shape_align_shape_tanimoto", tanis.first);
  fit.setProp("shape_align_color_tanimoto", tanis.second);

  return tanis;
}

std::pair<double, double> AlignMolecule(
    const ShapeInput &refShape, ROMol &fit, std::vector<float> &matrix,
    int fitConfId, bool useColors, double opt_param, unsigned int max_preiters,
    unsigned int max_postiters, bool applyRefShift) {
  ShapeInputOptions shapeOpts;
  shapeOpts.useColors = useColors;
  return AlignMolecule(refShape, fit, matrix, shapeOpts, fitConfId, opt_param,
                       max_preiters, max_postiters, applyRefShift);
}

std::pair<double, double> AlignMolecule(
    const ROMol &ref, ROMol &fit, std::vector<float> &matrix,
    const ShapeInputOptions &refShapeOpts,
    const ShapeInputOptions &probeShapeOpts, int refConfId, int fitConfId,
    double opt_param, unsigned int max_preiters, unsigned int max_postiters) {
  Align3D::setUseCutOff(true);

  if (refShapeOpts.useColors != probeShapeOpts.useColors) {
    BOOST_LOG(rdWarningLog)
        << "useColor values inconsistent between the reference and fit molecules. Colors will not be used in the alignment."
        << std::endl;
  }

  DEBUG_MSG("Reference details:");
  auto refShape = PrepareConformer(ref, refConfId, refShapeOpts);
  bool applyRefShift = true;
  auto scores =
      AlignMolecule(refShape, fit, matrix, probeShapeOpts, fitConfId, opt_param,
                    max_preiters, max_postiters, applyRefShift);
  return scores;
}

std::pair<double, double> AlignMolecule(const ROMol &ref, ROMol &fit,
                                        std::vector<float> &matrix,
                                        int refConfId, int fitConfId,
                                        bool useColors, double opt_param,
                                        unsigned int max_preiters,
                                        unsigned int max_postiters) {
  Align3D::setUseCutOff(true);

  DEBUG_MSG("Reference details:");
  ShapeInputOptions shapeOpts;
  shapeOpts.useColors = useColors;
  return AlignMolecule(ref, fit, matrix, shapeOpts, shapeOpts, refConfId,
                       fitConfId, opt_param, max_preiters, max_postiters);
}

std::pair<double, double> ScoreShape(const ShapeInput &shape1,
                                     ShapeInput &shape2, bool useColors) {
  double shapeOverlap = Align3D::ComputeShapeOverlap(
      shape1.coord.data(), shape1.alpha_vector, shape1.volumeAtomIndexVector,
      shape2.coord.data(), shape2.alpha_vector, shape2.volumeAtomIndexVector);

  double shape = shapeOverlap / (shape1.sov + shape2.sov - shapeOverlap);

  std::set<unsigned int> jointColorAtomTypeSet;
  Align3D::getJointColorTypeSet(
      shape1.atom_type_vector.data(), shape1.atom_type_vector.size(),
      shape2.atom_type_vector.data(), shape2.atom_type_vector.size(),
      jointColorAtomTypeSet);
  // Take copy of the color atom mappings so as not to alter the input shape
  // which might be re-used.
  double color = 0.0;
  if (useColors) {
    auto shape1MapCp = shape1.colorAtomType2IndexVectorMap;
    Align3D::restrictColorAtomType2IndexVectorMap(shape1MapCp,
                                                  jointColorAtomTypeSet);
    auto shape2MapCp = shape2.colorAtomType2IndexVectorMap;
    Align3D::restrictColorAtomType2IndexVectorMap(shape2MapCp,
                                                  jointColorAtomTypeSet);
    double featureOverlap = Align3D::ComputeFeatureOverlap(
        shape1.coord.data(), shape1.alpha_vector, shape1MapCp,
        shape2.coord.data(), shape2.alpha_vector, shape2MapCp, nullptr);
    color = featureOverlap / (shape1.sof + shape2.sof - featureOverlap);
  }
  return std::make_pair(shape, color);
}

std::pair<double, double> ScoreMolecule(const ShapeInput &shape,
                                        RDKit::ROMol &mol,
                                        const ShapeInputOptions &molShapeOpts,
                                        int molConfId) {
  ShapeInputOptions shapeOpts = molShapeOpts;
  shapeOpts.normalize = false;
  auto shape2 = PrepareConformer(mol, molConfId, shapeOpts);
  return ScoreShape(shape, shape2, shapeOpts.useColors);
}

std::pair<double, double> ScoreMolecule(const RDKit::ROMol &mol1,
                                        RDKit::ROMol &mol2,
                                        const ShapeInputOptions &mol1ShapeOpts,
                                        const ShapeInputOptions &mol2ShapeOpts,
                                        int mol1ConfId, int mol2ConfId) {
  ShapeInputOptions shapeOpts = mol1ShapeOpts;
  shapeOpts.normalize = false;
  auto shape1 = PrepareConformer(mol1, mol1ConfId, shapeOpts);
  shapeOpts = mol2ShapeOpts;
  shapeOpts.normalize = false;
  auto shape2 = PrepareConformer(mol2, mol2ConfId, shapeOpts);
  // If either shapeOpts has useColors of false, we have to go with that.
  bool useColors = mol1ShapeOpts.useColors && mol2ShapeOpts.useColors;
  return ScoreShape(shape1, shape2, useColors);
}
