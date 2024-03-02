//
//  Copyright (C) 2017-2023 Sereina Riniker and other RDKit contributors
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
#include <RDGeneral/utils.h>
#include <RDGeneral/RDLog.h>
#include <RDGeneral/Exceptions.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <RDGeneral/StreamOps.h>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
typedef boost::tokenizer<boost::char_separator<char>> tokenizer;
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>

namespace ForceFields {
namespace CrystalFF {
using namespace RDKit;

// the "macrocycle" patterns for ETKDGv3 use a minimum ring size of 9
const unsigned int MIN_MACROCYCLE_SIZE = 9;

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
#include "torsionPreferences_v1.in"
#include "torsionPreferences_v2.in"
#include "torsionPreferences_smallrings.in"
#include "torsionPreferences_macrocycles.in"

// class to store the experimental torsion angles
class ExpTorsionAngleCollection {
 public:
  typedef std::vector<ExpTorsionAngle> ParamsVect;
  static const ExpTorsionAngleCollection *getParams(
      unsigned int version, bool useSmallRingTorsions,
      bool useMacrocycleTorsions, const std::string &paramData = "");
  ParamsVect::const_iterator begin() const { return d_params.begin(); };
  ParamsVect::const_iterator end() const { return d_params.end(); };
  ExpTorsionAngleCollection(const std::string &paramData);

 private:
  ParamsVect d_params;  //!< the parameters
};

typedef boost::flyweight<
    boost::flyweights::key_value<std::string, ExpTorsionAngleCollection>,
    boost::flyweights::no_tracking>
    param_flyweight;

const ExpTorsionAngleCollection *ExpTorsionAngleCollection::getParams(
    unsigned int version, bool useSmallRingTorsions, bool useMacrocycleTorsions,
    const std::string &paramData) {
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

ExpTorsionAngleCollection::ExpTorsionAngleCollection(
    const std::string &paramData) {
  boost::char_separator<char> tabSep(" ", "", boost::drop_empty_tokens);
  std::istringstream inStream(paramData);

  std::string inLine = RDKit::getLine(inStream);
  unsigned int torsionIdx = 0;
  while (!inStream.eof()) {
    if (inLine[0] != '#') {
      ExpTorsionAngle angle;
      tokenizer tokens(inLine, tabSep);
      tokenizer::iterator token = tokens.begin();
      angle.smarts = *token;
      angle.torsionIdx = torsionIdx++;
      ++token;
      for (unsigned int i = 0; i < 12; i += 2) {
        angle.signs.push_back(boost::lexical_cast<int>(*token));
        ++token;
        angle.V.push_back(boost::lexical_cast<double>(*token));
        ++token;
      }
      angle.dp_pattern.reset(SmartsToMol(angle.smarts));
      // get the atom indices for atom 1, 2, 3, 4 in the pattern
      for (unsigned int i = 0; i < (angle.dp_pattern.get())->getNumAtoms();
           ++i) {
        Atom const *atom = (angle.dp_pattern.get())->getAtomWithIdx(i);
        int num;
        if (atom->getPropIfPresent("molAtomMapNumber", num)) {
          if (num > 0 && num < 5) {
            angle.idx[num - 1] = i;
          }
        }
      }
      d_params.push_back(std::move(angle));
    }
    inLine = RDKit::getLine(inStream);
  }  // while loop
  // std::cerr << "Exp. torsion angles = " << d_params.size() << " "
  //    << d_params[d_params.size()-1].smarts << std::endl;
}

void getExperimentalTorsions(
    const RDKit::ROMol &mol, CrystalFFDetails &details,
    std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                           const ExpTorsionAngle *>> &torsionBonds,
    bool useExpTorsions, bool useSmallRingTorsions, bool useMacrocycleTorsions,
    bool useBasicKnowledge, unsigned int version, bool verbose) {
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

  if (useExpTorsions) {
    // we set the torsion angles with experimental data
    const auto *params = ExpTorsionAngleCollection::getParams(
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
        if (!doneBonds[bid2]) {
          // do not add ET terms between constrained atoms
          // REVIEW: do we really need to check all 4 atoms?
          if (!details.constrainedAtoms.empty() &&
              details.constrainedAtoms[aid1] &&
              details.constrainedAtoms[aid2] &&
              details.constrainedAtoms[aid3] &&
              details.constrainedAtoms[aid4]) {
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
          details.expTorsionAngles.emplace_back(param.signs, param.V);
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
        }  // if not donePaths
      }    // end loop over matches

    }  // end loop over patterns
  }

  // apply basic knowledge such as flat aromatic rings, other sp2-centers,
  // straight triple bonds, etc.
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
      unsigned int rSize = atomRing.size();
      // we don't need to deal with 3 membered rings
      // and we do not treat rings greater than 6
      if (rSize < 4 || rSize > 6) {
        continue;
      }
      // loop over ring atoms
      for (unsigned int i = 0; i < rSize; ++i) {
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

          std::vector<int> signs(6, 1);
          signs[1] = -1;  // MMFF sign for m = 2
          std::vector<double> fconsts(6, 0.0);
          fconsts[1] = 100.0;  // 7.0 is MMFF force constants for aromatic rings
          details.expTorsionAngles.emplace_back(signs, fconsts);
          /*if (verbose) {
            std::cout << "SP2 ring: " << aid1 << " " << aid2 << " " << aid3 <<
          " " << aid4 << std::endl;
          }*/
        }

      }  // loop over atoms in ring
    }    // loop over rings
  }      // if useBasicKnowledge

}  // end function

void getExperimentalTorsions(const RDKit::ROMol &mol, CrystalFFDetails &details,
                             bool useExpTorsions, bool useSmallRingTorsions,
                             bool useMacrocycleTorsions, bool useBasicKnowledge,
                             unsigned int version, bool verbose) {
  std::vector<std::tuple<unsigned int, std::vector<unsigned int>,
                         const ExpTorsionAngle *>>
      torsionBonds;
  getExperimentalTorsions(mol, details, torsionBonds, useExpTorsions,
                          useSmallRingTorsions, useMacrocycleTorsions,
                          useBasicKnowledge, version, verbose);
}

}  // namespace CrystalFF
}  // namespace ForceFields
