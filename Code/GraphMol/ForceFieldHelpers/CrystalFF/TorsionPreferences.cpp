//
//  Copyright (C) 2017-2026 Sereina Riniker and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "TorsionPreferences.h"
#include <GraphMol/RDKitBase.h>
#include <Geometry/Utils.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/DistGeomHelpers/BoundsMatrixBuilderDetails.h>
#include <GraphMol/ForceFieldHelpers/CrystalFF/GaussianTorsionAngleContribs.h>
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/StreamOps.h>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

#include <sstream>
#include <numbers>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ranges>

namespace ForceFields {
namespace CrystalFF {
using namespace RDKit;

// the "macrocycle" patterns for ETKDGv3 use a minimum ring size of 9
constexpr unsigned int MIN_MACROCYCLE_SIZE = 9;

// parameters for ETKDGv4
const std::vector<double> AROMATIC_HEIGHTS = {2.305};
const std::vector<double> AROMATIC_POSITIONS = {0.0};
const std::vector<double> TRANS_POSITION = {std::numbers::pi_v<double>};
const std::vector<double> CIS_POSITION = {0.0};
const std::vector<double> AROMATIC_WIDTHS = {0.251};
constexpr std::size_t AMIDE_TRANS_IDX = 2000;
constexpr std::size_t AMIDE_CIS_IDX = 2001;
constexpr std::size_t AROM_TORSION_IDX = 2002;
constexpr std::size_t LOOKUP_GRID_SIZE = 180;

/* SMARTS patterns for experimental torsion angle preferences
 * Version 1 taken from J. Med. Chem. 56, 1026-2028 (2013)
 * Version 2 taken from J. Chem. Inf. Model. 56, 1 (2016)
 *
 * torsion-angle potential form:
 * V = V1*(1 + s1*cos(1x)) + V2*(1 + s2*cos(2x)) + V3*(1 + s3*cos(3x))
 *     + V4*(1 + s4*cos(4x)) + V5*(1 + s5*cos(5x)) + V6*(1 + s6*cos(6x))
 *
 * format: [SMARTS, s1, V1, s2, V2, s3, V3, s4, V4, s5, V5, s6, V6]
 */
#include "torsionData/torsionPreferences_v1.in"
#include "torsionData/torsionPreferences_v2.in"
#include "torsionData/torsionPreferences_smallrings.in"
#include "torsionData/torsionPreferences_macrocycles.in"

#include "torsionData/gaussTorsionPreferences_v4.h"
#include "torsionData/gaussTorsionPreferences_v4_macrocycles.h"
#include "torsionData/gaussTorsionPreferences_v4_smallrings.h"

// class to store the experimental torsion angles
template <TorsionAngleType T>
class ExpTorsionAngleCollection {
 public:
  typedef std::vector<T> ParamsVect;
  static const ExpTorsionAngleCollection<T> *getParams(
      unsigned int version, bool useSmallRingTorsions,
      bool useMacrocycleTorsions, const std::string &paramData = "") requires
      std::is_same_v<T, ExpTorsionAngle>;
  static const ExpTorsionAngleCollection<T> *getParams(
      unsigned int version, bool useSmallRingTorsions,
      bool useMacrocycleTorsions,
      const GaussianTorsionAngleRaw &paramData = {}) requires
      std::is_same_v<T, GaussianExpTorsionAngle>;
  ParamsVect::const_iterator begin() const { return d_params.begin(); };
  ParamsVect::const_iterator end() const { return d_params.end(); };
  ExpTorsionAngleCollection(const std::string &paramData);
  ExpTorsionAngleCollection(
      const std::string &paramData) requires std::is_same_v<T, ExpTorsionAngle>;
  ExpTorsionAngleCollection(const GaussianTorsionAngleRaw &paramData) requires
      std::is_same_v<T, GaussianExpTorsionAngle>;

 private:
  ParamsVect d_params;  //!< the parameters
};

using param_flyweight = boost::flyweight<
    boost::flyweights::key_value<std::string,
                                 ExpTorsionAngleCollection<ExpTorsionAngle>>,
    boost::flyweights::no_tracking>;

using GaussianParamFlyweight =
    boost::flyweight<boost::flyweights::key_value<
                         GaussianTorsionAngleRaw,
                         ExpTorsionAngleCollection<GaussianExpTorsionAngle>>,
                     boost::flyweights::no_tracking>;

template <TorsionAngleType T>
const ExpTorsionAngleCollection<T> *ExpTorsionAngleCollection<T>::getParams(
    unsigned int version, bool useSmallRingTorsions, bool useMacrocycleTorsions,
    const std::string &paramData) requires std::is_same_v<T, ExpTorsionAngle> {
  std::string params;
  if (paramData.empty()) {
    switch (version) {
      case 1:
        params = torsionPreferencesV1;
        break;
      case 2:
        params = torsionPreferencesV2;
        break;
      default:
        throw ValueErrorException("ETversion must be 1 or 2.");
    }
  } else {
    params = paramData;
  }

  if (useSmallRingTorsions) {
    params += torsionPreferencesSmallRings;
  }

  if (useMacrocycleTorsions) {
    params += torsionPreferencesMacrocycles;
  }

  return &(param_flyweight(params).get());
}
template <TorsionAngleType T>
const ExpTorsionAngleCollection<T> *ExpTorsionAngleCollection<T>::getParams(
    unsigned int version, bool useSmallRingTorsions, bool useMacrocycleTorsions,
    const GaussianTorsionAngleRaw &paramData) requires
    std::is_same_v<T, GaussianExpTorsionAngle> {
  if (version != 4) {
    throw ValueErrorException(
        "Currently, only ETVersions 4 is supported together with the Gaussian fit functional form.");
  }
  if (!paramData.empty()) {
    return &(GaussianParamFlyweight(paramData).get());
  }
  namespace GT = GaussTorsions::V4;
  auto data = GT::Regular;
  if (useSmallRingTorsions) {
    std::ranges::copy(GaussTorsions::SmallRings, std::back_inserter(data));
  }
  if (useMacrocycleTorsions) {
    std::ranges::copy(GT::MacroCycles, std::back_inserter(data));
  }
  return &(GaussianParamFlyweight(data).get());
}

template <TorsionAngleType T>
void addPattern(T &angle) {
  angle.dp_pattern.reset(SmartsToMol(angle.smarts));
  // get the atom indices for atom 1, 2, 3, 4 in the pattern
  for (unsigned int i = 0; i < (angle.dp_pattern.get())->getNumAtoms(); ++i) {
    Atom const *atom = (angle.dp_pattern.get())->getAtomWithIdx(i);
    int num;
    if (atom->getPropIfPresent("molAtomMapNumber", num)) {
      if (num > 0 && num < 5) {
        angle.idx[num - 1] = i;
      }
    }
  }
}

template <TorsionAngleType T>
ExpTorsionAngleCollection<T>::ExpTorsionAngleCollection(
    const std::string &paramData) requires std::is_same_v<T, ExpTorsionAngle> {
  boost::char_separator<char> tabSep(" ", "", boost::drop_empty_tokens);
  std::istringstream inStream(paramData);

  std::string inLine = RDKit::getLine(inStream);
  std::size_t torsionIdx = 0;
  while (!inStream.eof()) {
    if (inLine[0] != '#') {
      ExpTorsionAngle angle;
      tokenizer tokens(inLine, tabSep);
      tokenizer::iterator token = tokens.begin();
      angle.smarts = *token;
      angle.torsionIdx = torsionIdx++;
      ++token;
      for (std::size_t i = 0; i < 12; i += 2) {
        angle.signs.push_back(boost::lexical_cast<int>(*token));
        ++token;
        angle.V.push_back(boost::lexical_cast<double>(*token));
        ++token;
      }
      addPattern(angle);
      d_params.push_back(std::move(angle));
    }
    inLine = RDKit::getLine(inStream);
  }  // while loop
  // std::cerr << "Exp. torsion angles = " << d_params.size() << " "
  //    << d_params[d_params.size()-1].smarts << std::endl;
}

template <TorsionAngleType T>
ExpTorsionAngleCollection<T>::ExpTorsionAngleCollection(
    const GaussianTorsionAngleRaw &paramData) requires
    std::is_same_v<T, GaussianExpTorsionAngle> {
  for (std::size_t i = 0; i < paramData.size(); ++i) {
    auto &param = paramData[i];
    GaussianExpTorsionAngle angle;
    angle.torsionIdx = i;
    angle.smarts = std::get<0>(param);
    angle.heights = std::get<1>(param);
    angle.positions = std::get<2>(param);
    angle.widths = std::get<3>(param);
    addPattern(angle);
    d_params.push_back(std::move(angle));
  }
}

template <TorsionParamType T>
void getExperimentalTorsions(
    const RDKit::ROMol &mol, CrystalFFDetails<T> &details,
    std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                           const MappedAngle_T<T> *>> &torsionBonds,
    bool useExpTorsions, bool useSmallRingTorsions, bool useMacrocycleTorsions,
    bool useBasicKnowledge, unsigned int version, bool verbose) {
  using Angle_T = MappedAngle_T<T>;
  torsionBonds.clear();
  unsigned int nb = mol.getNumBonds();
  unsigned int na = mol.getNumAtoms();
  if (!na) {
    throw ValueErrorException("molecule has no atoms");
  }

  RDLog::LogStateSetter logs(
      RDLog::RDLoggerList({rdInfoLog, rdErrorLog, rdWarningLog}));
  // check that vectors are empty
  details.expTorsionAtoms.clear();
  details.expTorsionAngles.clear();
  details.improperAtoms.clear();

  unsigned int aid1, aid2, aid3, aid4;
  unsigned int bid2;

  // exclude bonds in bridged ring systems
  boost::dynamic_bitset<> excludedBonds(nb);
  const RingInfo *rinfo = mol.getRingInfo();
  const VECT_INT_VECT &bondRings = rinfo->bondRings();
  for (auto rii = bondRings.begin(); rii != bondRings.end(); ++rii) {
    boost::dynamic_bitset<> rs1(nb);  // bitset for ring 1
    for (auto riiv : *rii) {
      rs1[riiv] = 1;
    }
    for (auto rjj = rii + 1; rjj != bondRings.end(); ++rjj) {
      // we don't worry about the overlap if both rings are macrocycles:
      if (rii->size() >= MIN_MACROCYCLE_SIZE &&
          rjj->size() >= MIN_MACROCYCLE_SIZE) {
        continue;
      }
      unsigned int nInCommon = 0;
      for (auto rjj_i : *rjj) {
        if (rs1[rjj_i]) {
          if (++nInCommon > 1) {
            break;
          }
        }
      }
      if (nInCommon > 1) {  // more than one bond in common
        // exclude bonds from non-macrocycles:
        if (rii->size() < MIN_MACROCYCLE_SIZE) {
          for (unsigned int i = 0; i < rii->size(); i++) {
            excludedBonds[(*rii)[i]] = 1;  // exclude all bonds of ring 1
          }
        }
        if (rjj->size() < MIN_MACROCYCLE_SIZE) {
          for (unsigned int i = 0; i < rjj->size(); i++) {
            excludedBonds[(*rjj)[i]] = 1;  // exclude all bonds of ring 2
          }
        }
      }
    }
  }

  boost::dynamic_bitset<> doneBonds(nb);
  // apply basic knowledge such as flat aromatic rings, other sp2-centers,
  // straight triple bonds, etc.
  if (useBasicKnowledge) {
    // torsions for forced trans amides / esters
    auto is_forced_cis_or_trans = [](const auto &config) {
      if (!config.type.isForced) {
        return false;
      }
      return config.type.type == DGeomHelpers::TorsionType::TRANS ||
             config.type.type == DGeomHelpers::TorsionType::CIS;
    };
    for (const auto &config :
         details.path14Configs | std::views::filter(is_forced_cis_or_trans)) {
      const auto i = config.aid1;
      const auto j = config.aid2;
      const auto k = config.aid3;
      const auto l = config.aid4;
      const auto bndIdx = mol.getBondBetweenAtoms(j, k)->getIdx();

      if (excludedBonds[bndIdx] ||
          mol.getRingInfo()->numBondRings(bndIdx) > 3) {
        doneBonds[bndIdx] = 1;
      }
      if (doneBonds[bndIdx]) {
        continue;
      }

      if (!details.constrainedAtoms.empty() && details.constrainedAtoms[i] &&
          details.constrainedAtoms[j] && details.constrainedAtoms[k] &&
          details.constrainedAtoms[l]) {
        continue;
      }
      const bool isCis = config.type.type == DGeomHelpers::TorsionType::CIS;
      details.expTorsionAtoms.push_back(
          {static_cast<int>(i), static_cast<int>(j), static_cast<int>(k),
           static_cast<int>(l)});
      const bool isAIO =
          std::fabs(details.forceConsts.etTermScaling - 1.0) > 1e-3;
      details.torsionIdx.push_back(isCis ? AMIDE_CIS_IDX : AMIDE_TRANS_IDX);
      if constexpr (std::is_same_v<T, CosineExp_T>) {
        std::vector<double> V(6, 0.0);
        std::vector<int> signs(6, 1);

        V[0] = 75.0;
        if (isAIO) {
          V[0] = 4.0;
        }
        if (isCis) {
          signs[0] = -1;
        }
        details.expTorsionAngles.emplace_back(signs, V);
      } else if constexpr (std::is_same_v<T, GaussianExp_T>) {
        details.expTorsionAngles.emplace_back(
            AROMATIC_HEIGHTS, isCis ? CIS_POSITION : TRANS_POSITION,
            AROMATIC_WIDTHS, isAIO ? 0.05 : 1.0);
      }
    }

  }  // if useBasicKnowledge
  if (useExpTorsions) {
    // we set the torsion angles with experimental data
    const auto *params = ExpTorsionAngleCollection<Angle_T>::getParams(
        version, useSmallRingTorsions, useMacrocycleTorsions);
    CHECK_INVARIANT(params, "no parameters available");
    // loop over patterns
    for (const auto &param : *params) {
      std::vector<MatchVectType> matches;
      SubstructMatch(mol, *(param.dp_pattern.get()), matches, false, true);
      // loop over matches
      for (const auto &match : matches) {
        // get bond indices
        aid1 = match[param.idx[0]].second;
        aid2 = match[param.idx[1]].second;
        aid3 = match[param.idx[2]].second;
        aid4 = match[param.idx[3]].second;
        const auto bnd = mol.getBondBetweenAtoms(aid2, aid3);
        CHECK_INVARIANT(bnd, "bond between central atoms not found")
        bid2 = bnd->getIdx();

        // check that a bond is part of maximum one ring
        if (excludedBonds[bid2] || mol.getRingInfo()->numBondRings(bid2) > 3) {
          doneBonds[bid2] = 1;
        }
        if (doneBonds[bid2]) {
          continue;
        }
        // do not add ET terms between constrained atoms
        // REVIEW: do we really need to check all 4 atoms?
        if (!details.constrainedAtoms.empty() &&
            details.constrainedAtoms[aid1] && details.constrainedAtoms[aid2] &&
            details.constrainedAtoms[aid3] && details.constrainedAtoms[aid4]) {
          continue;
        }
        std::vector<unsigned int> aids{aid1, aid2, aid3, aid4};
        torsionBonds.emplace_back(bid2, aids, &param);
        doneBonds[bid2] = 1;
        std::vector<int> atoms(4);
        atoms[0] = aid1;
        atoms[1] = aid2;
        atoms[2] = aid3;
        atoms[3] = aid4;
        details.expTorsionAtoms.push_back(atoms);
        details.torsionIdx.push_back(param.torsionIdx);
        if constexpr (std::is_same_v<T, CosineExp_T>) {
          std::vector<double> vals(param.V);
          for (auto &val : vals) {
            val *= details.forceConsts.etTermScaling;
          }
          details.expTorsionAngles.emplace_back(param.signs, vals);
          if (verbose) {
            // using the stringstream seems redundant, but we don't want the
            // extra formatting provided by the logger after every entry;
            std::stringstream sstr;
            sstr << param.smarts << ": " << aid1 << " " << aid2 << " " << aid3
                 << " " << aid4 << ", [";
            for (unsigned int i = 0; i < param.V.size() - 1; ++i) {
              sstr << "(" << param.signs[i] << " " << param.V[i] << "), ";
            }
            sstr << "(" << param.signs.back() << " " << param.V.back() << ")] ";
            BOOST_LOG(rdInfoLog) << sstr.str() << std::endl;
          }
        } else if constexpr (std::is_same_v<T, GaussianExp_T>) {
          details.expTorsionAngles.emplace_back(
              param.heights, param.positions, param.widths,
              details.forceConsts.etTermScaling);
          if (verbose) {
            // using the stringstream seems redundant, but we don't want the
            // extra formatting provided by the logger after every entry;
            std::stringstream sstr;
            sstr << param.smarts << ": " << aid1 << " " << aid2 << " " << aid3
                 << " " << aid4 << ", [";
            BOOST_LOG(rdInfoLog) << sstr.str() << std::endl;
          }
        }
      }  // if not donePaths
    }    // end loop over matches
  }      // end loop over patterns
  if (useBasicKnowledge) {
    boost::dynamic_bitset<> doneAtoms(na);

    // inversion terms (improper torsions / out-of-plane bends / inversion)
    // loop over atoms
    for (aid2 = 0; aid2 < na; ++aid2) {
      if (!(doneAtoms[aid2])) {
        std::vector<int> atoms(4, -1);
        atoms[1] = aid2;
        const Atom *atom2 = mol.getAtomWithIdx(atoms[1]);
        int at2AtomicNum = atom2->getAtomicNum();

        // if atom is a N,O or C, SP2-hybridized, and has three neighbors
        if (((at2AtomicNum == 6) || (at2AtomicNum == 7) ||
             (at2AtomicNum == 8)) &&
            (atom2->getHybridization() == Atom::SP2) &&
            mol.getAtomDegree(atom2) == 3) {
          unsigned int i = 0;
          unsigned int isBoundToSP2O = 0;  // false
          for (const auto atomX : mol.atomNeighbors(atom2)) {
            atoms[i] = atomX->getIdx();
            // if the central atom is sp2 carbon and is bound to sp2 oxygen,
            // set a flag
            if (!isBoundToSP2O) {
              isBoundToSP2O =
                  ((at2AtomicNum == 6) && (atomX->getAtomicNum() == 8) &&
                   (atomX->getHybridization() == Atom::SP2));
            }
            if (!i) {
              ++i;
            }
            ++i;
          }
          atoms.push_back(at2AtomicNum);
          atoms.push_back(isBoundToSP2O);
          details.improperAtoms.push_back(atoms);
          /*if (verbose) {
            std::cout << "out-of-plane bend: " << atoms[0] << " " << atoms[1]
          << " "
                << atoms[2] << " " << atoms[3] << std::endl;
          }*/
        }
      }  // if atom is a N,O or C and SP2-hybridized
    }

    // torsions for flat rings
    const RingInfo *rinfo = mol.getRingInfo();
    CHECK_INVARIANT(rinfo, "no ring info");
    CHECK_INVARIANT(rinfo->isInitialized(), "ring info not initialized");
    for (const auto &atomRing : rinfo->atomRings()) {
      std::size_t rSize = atomRing.size();
      // we don't need to deal with 3 membered rings
      // and we do not treat rings greater than 6
      if (rSize < 4 || rSize > 6) {
        continue;
      }
      // loop over ring atoms
      for (std::size_t i = 0; i < rSize; ++i) {
        // proper torsions
        aid1 = atomRing[i];
        aid2 = atomRing[(i + 1) % rSize];
        aid3 = atomRing[(i + 2) % rSize];
        aid4 = atomRing[(i + 3) % rSize];
        bid2 = mol.getBondBetweenAtoms(aid2, aid3)->getIdx();
        // if all 4 atoms are SP2, add torsion
        if (!(doneBonds[bid2]) &&
            (mol.getAtomWithIdx(aid1)->getHybridization() == Atom::SP2) &&
            (mol.getAtomWithIdx(aid2)->getHybridization() == Atom::SP2) &&
            (mol.getAtomWithIdx(aid3)->getHybridization() == Atom::SP2) &&
            (mol.getAtomWithIdx(aid4)->getHybridization() == Atom::SP2)) {
          doneBonds[bid2] = 1;
          std::vector<int> atoms(4);
          atoms[0] = aid1;
          atoms[1] = aid2;
          atoms[2] = aid3;
          atoms[3] = aid4;
          details.expTorsionAtoms.push_back(atoms);
          details.torsionIdx.push_back(AROM_TORSION_IDX);

          if constexpr (std::is_same_v<T, CosineExp_T>) {
            std::vector<int> signs(6, 1);
            signs[1] = -1;  // MMFF sign for m = 2
            std::vector<double> fconsts(6, 0.0);
            fconsts[1] = details.forceConsts
                             .kTermTorsion;  // 7.0 is MMFF force constants
                                             // for aromatic rings
            details.expTorsionAngles.emplace_back(signs, fconsts);
          } else if constexpr (std::is_same_v<T, GaussianExp_T>) {
            details.expTorsionAngles.emplace_back(
                AROMATIC_HEIGHTS, AROMATIC_POSITIONS, AROMATIC_WIDTHS,
                details.forceConsts.kTermTorsion);
          }
          /*if (verbose) {
            std::cout << "SP2 ring: " << aid1 << " " << aid2 << " " << aid3 <<
          " " << aid4 << std::endl;
          }*/
        }

      }  // loop over atoms in ring
    }    // loop over rings
  }      // end function
}

template <TorsionParamType T>
void getExperimentalTorsions(const RDKit::ROMol &mol,
                             CrystalFFDetails<T> &details, bool useExpTorsions,
                             bool useSmallRingTorsions,
                             bool useMacrocycleTorsions, bool useBasicKnowledge,
                             unsigned int version, bool verbose) {
  std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                         const MappedAngle_T<T> *>>
      torsionBonds;
  getExperimentalTorsions(mol, details, torsionBonds, useExpTorsions,
                          useSmallRingTorsions, useMacrocycleTorsions,
                          useBasicKnowledge, version, verbose);
}

void populateRefTable(CrystalFFDetails<GaussianExp_T> &details) {
  std::vector<std::size_t> sorted = details.torsionIdx;
  std::ranges::sort(sorted);
  auto [first, last] = std::ranges::unique(sorted);
  sorted.erase(first, last);
  std::ranges::transform(
      details.torsionIdx, details.torsionIdx.begin(), [&sorted](int x) {
        return std::ranges::distance(sorted.begin(),
                                     std::ranges::lower_bound(sorted, x));
      });
  details.cosPhiToEnergy.resize(sorted.size(),
                                std::vector<double>(LOOKUP_GRID_SIZE));
  details.cosPhiToGrad.resize(sorted.size(),
                              std::vector<double>(LOOKUP_GRID_SIZE));
  for (std::size_t torsionIdx = 0; torsionIdx < sorted.size(); ++torsionIdx) {
    auto it = std::ranges::find(details.torsionIdx, torsionIdx);
    std::size_t termIdx = std::distance(details.torsionIdx.begin(), it);
    auto &heights = std::get<0>(details.expTorsionAngles[termIdx]);
    auto &positions = std::get<1>(details.expTorsionAngles[termIdx]);
    auto &widths = std::get<2>(details.expTorsionAngles[termIdx]);
    for (std::size_t gridPoint = 0; gridPoint < LOOKUP_GRID_SIZE; ++gridPoint) {
      const double phi = gridPoint * std::numbers::pi / (LOOKUP_GRID_SIZE - 1);
      details.cosPhiToEnergy[torsionIdx][gridPoint] =
          getEnergy(heights, positions, widths, phi);
      details.cosPhiToGrad[torsionIdx][gridPoint] =
          getdEdPhi(heights, positions, widths, phi);
    }
  }
}

// Some explicit instantiations
template struct CrystalFFDetails<CosineExp_T>;
template struct CrystalFFDetails<GaussianExp_T>;
template class ExpTorsionAngleCollection<ExpTorsionAngle>;
template class ExpTorsionAngleCollection<GaussianExpTorsionAngle>;

template RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &, CrystalFFDetails<CosineExp_T> &, bool, bool, bool,
    bool, unsigned int, bool);
template RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &, CrystalFFDetails<GaussianExp_T> &, bool, bool, bool,
    bool, unsigned int, bool);
template RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &, CrystalFFDetails<CosineExp_T> &,
    std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                           const ExpTorsionAngle *>> &,
    bool, bool, bool, bool, unsigned int, bool);
template RDKIT_FORCEFIELDHELPERS_EXPORT void getExperimentalTorsions(
    const RDKit::ROMol &, CrystalFFDetails<GaussianExp_T> &,
    std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                           const GaussianExpTorsionAngle *>> &,
    bool, bool, bool, bool, unsigned int, bool);

}  // namespace CrystalFF
}  // namespace ForceFields
