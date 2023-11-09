//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <list>
#include <algorithm>
#include <cmath>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <iostream>
#include <sstream>
#include "SubstructMatchCustom.h"
#include "MaximumCommonSubgraph.h"
#include <GraphMol/QueryOps.h>

namespace RDKit {

namespace {
struct cmpCStr {
  bool operator()(const char* a, const char* b) const {
    return std::strcmp(a, b) < 0;
  }
};
}  // namespace

void MCSParameters::setMCSAtomTyperFromEnum(AtomComparator atomComp) {
  switch (atomComp) {
    case AtomCompareAny:
      AtomTyper = MCSAtomCompareAny;
      break;
    case AtomCompareElements:
      AtomTyper = MCSAtomCompareElements;
      break;
    case AtomCompareIsotopes:
      AtomTyper = MCSAtomCompareIsotopes;
      break;
    case AtomCompareAnyHeavyAtom:
      AtomTyper = MCSAtomCompareAnyHeavyAtom;
      break;
    default:
      throw ValueErrorException("Unknown AtomComparator");
  }
}

void MCSParameters::setMCSAtomTyperFromConstChar(const char* atomComp) {
  PRECONDITION(atomComp, "atomComp must not be NULL");
  static const std::map<const char*, AtomComparator, cmpCStr>
      atomCompStringToEnum = {{"Any", AtomCompareAny},
                              {"Elements", AtomCompareElements},
                              {"Isotopes", AtomCompareIsotopes},
                              {"AnyHeavy", AtomCompareAnyHeavyAtom}};
  const auto it = atomCompStringToEnum.find(atomComp);
  // we accept "def" as a no-op
  if (it != atomCompStringToEnum.end()) {
    setMCSAtomTyperFromEnum(it->second);
  }
}

void MCSParameters::setMCSBondTyperFromEnum(BondComparator bondComp) {
  switch (bondComp) {
    case BondCompareAny:
      BondTyper = MCSBondCompareAny;
      break;
    case BondCompareOrder:
      BondTyper = MCSBondCompareOrder;
      break;
    case BondCompareOrderExact:
      BondTyper = MCSBondCompareOrderExact;
      break;
    default:
      throw ValueErrorException("Unknown BondComparator");
  }
}

void MCSParameters::setMCSBondTyperFromConstChar(const char* bondComp) {
  PRECONDITION(bondComp, "bondComp must not be NULL");
  static const std::map<const char*, BondComparator, cmpCStr>
      bondCompStringToEnum = {{"Any", BondCompareAny},
                              {"Order", BondCompareOrder},
                              {"OrderExact", BondCompareOrderExact}};
  const auto it = bondCompStringToEnum.find(bondComp);
  // we accept "def" as a no-op
  if (it != bondCompStringToEnum.end()) {
    setMCSBondTyperFromEnum(it->second);
  }
}

void parseMCSParametersJSON(const char* json, MCSParameters* params) {
  if (!params || !json || !strlen(json)) {
    return;
  }
  std::istringstream ss;
  ss.str(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);

  auto& p = *params;
  p.MaximizeBonds = pt.get<bool>("MaximizeBonds", p.MaximizeBonds);
  p.Threshold = pt.get<double>("Threshold", p.Threshold);
  p.Timeout = pt.get<unsigned int>("Timeout", p.Timeout);
  p.AtomCompareParameters.MatchValences =
      pt.get<bool>("MatchValences", p.AtomCompareParameters.MatchValences);
  p.AtomCompareParameters.MatchChiralTag =
      pt.get<bool>("MatchChiralTag", p.AtomCompareParameters.MatchChiralTag);
  p.AtomCompareParameters.MatchFormalCharge = pt.get<bool>(
      "MatchFormalCharge", p.AtomCompareParameters.MatchFormalCharge);
  p.AtomCompareParameters.RingMatchesRingOnly = pt.get<bool>(
      "RingMatchesRingOnly", p.AtomCompareParameters.RingMatchesRingOnly);
  p.AtomCompareParameters.MaxDistance =
      pt.get<double>("MaxDistance", p.AtomCompareParameters.MaxDistance);
  p.BondCompareParameters.RingMatchesRingOnly = pt.get<bool>(
      "RingMatchesRingOnly", p.BondCompareParameters.RingMatchesRingOnly);
  p.AtomCompareParameters.RingMatchesRingOnly = pt.get<bool>(
      "AtomRingMatchesRingOnly", p.AtomCompareParameters.RingMatchesRingOnly);
  p.BondCompareParameters.RingMatchesRingOnly = pt.get<bool>(
      "BondRingMatchesRingOnly", p.BondCompareParameters.RingMatchesRingOnly);
  p.BondCompareParameters.CompleteRingsOnly = pt.get<bool>(
      "CompleteRingsOnly", p.BondCompareParameters.CompleteRingsOnly);
  p.AtomCompareParameters.CompleteRingsOnly = pt.get<bool>(
      "AtomCompleteRingsOnly", p.AtomCompareParameters.CompleteRingsOnly);
  p.BondCompareParameters.CompleteRingsOnly = pt.get<bool>(
      "BondCompleteRingsOnly", p.BondCompareParameters.CompleteRingsOnly);
  p.BondCompareParameters.MatchFusedRings =
      pt.get<bool>("MatchFusedRings", p.BondCompareParameters.MatchFusedRings);
  p.BondCompareParameters.MatchFusedRingsStrict = pt.get<bool>(
      "MatchFusedRingsStrict", p.BondCompareParameters.MatchFusedRingsStrict);
  p.BondCompareParameters.MatchStereo =
      pt.get<bool>("MatchStereo", p.BondCompareParameters.MatchStereo);
  p.StoreAll = pt.get<bool>("StoreAll", p.StoreAll);

  p.setMCSAtomTyperFromConstChar(
      pt.get<std::string>("AtomCompare", "def").c_str());
  p.setMCSBondTyperFromConstChar(
      pt.get<std::string>("BondCompare", "def").c_str());

  p.InitialSeed = pt.get<std::string>("InitialSeed", "");
}

MCSResult findMCS(const std::vector<ROMOL_SPTR>& mols,
                  const MCSParameters* params) {
  MCSParameters p;
  if (nullptr == params) {
    params = &p;
  }
  RDKit::FMCS::MaximumCommonSubgraph fmcs(params);
  return fmcs.find(mols);
}

MCSResult findMCS_P(const std::vector<ROMOL_SPTR>& mols,
                    const char* params_json) {
  MCSParameters p;
  parseMCSParametersJSON(params_json, &p);
  return findMCS(mols, &p);
}

MCSResult findMCS(const std::vector<ROMOL_SPTR>& mols, bool maximizeBonds,
                  double threshold, unsigned int timeout, bool verbose,
                  bool matchValences, bool ringMatchesRingOnly,
                  bool completeRingsOnly, bool matchChiralTag,
                  AtomComparator atomComp, BondComparator bondComp) {
  return findMCS(mols, maximizeBonds, threshold, timeout, verbose,
                 matchValences, ringMatchesRingOnly, completeRingsOnly,
                 matchChiralTag, atomComp, bondComp, IgnoreRingFusion);
}

MCSResult findMCS(const std::vector<ROMOL_SPTR>& mols, bool maximizeBonds,
                  double threshold, unsigned int timeout, bool verbose,
                  bool matchValences, bool ringMatchesRingOnly,
                  bool completeRingsOnly, bool matchChiralTag,
                  AtomComparator atomComp, BondComparator bondComp,
                  RingComparator ringComp) {
  MCSParameters ps;
  ps.MaximizeBonds = maximizeBonds;
  ps.Threshold = threshold;
  ps.Timeout = timeout;
  ps.Verbose = verbose;
  ps.setMCSAtomTyperFromEnum(atomComp);
  ps.AtomCompareParameters.MatchValences = matchValences;
  ps.AtomCompareParameters.MatchChiralTag = matchChiralTag;
  ps.AtomCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  ps.setMCSBondTyperFromEnum(bondComp);
  ps.BondCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  ps.BondCompareParameters.CompleteRingsOnly = completeRingsOnly;
  ps.BondCompareParameters.MatchFusedRings = (ringComp != IgnoreRingFusion);
  ps.BondCompareParameters.MatchFusedRingsStrict =
      (ringComp == StrictRingFusion);
  return findMCS(mols, &ps);
}

bool MCSProgressCallbackTimeout(const MCSProgressData&,
                                const MCSParameters& params, void* userData) {
  PRECONDITION(userData, "userData must not be NULL");
  auto t0 = static_cast<unsigned long long*>(userData);
  unsigned long long t = nanoClock();
  return !params.Timeout || (t - *t0 <= params.Timeout * 1000000ULL);
}

// PREDEFINED FUNCTORS:

//=== ATOM COMPARE ========================================================
bool checkAtomRingMatch(const MCSAtomCompareParameters& p, const ROMol& mol1,
                        unsigned int atom1, const ROMol& mol2,
                        unsigned int atom2) {
  if (p.RingMatchesRingOnly) {
    const auto ri1 = mol1.getRingInfo();
    const auto ri2 = mol2.getRingInfo();
    bool atom1inRing = (ri1->numAtomRings(atom1) > 0);
    bool atom2inRing = (ri2->numAtomRings(atom2) > 0);
    return atom1inRing == atom2inRing;
  } else {
    return true;
  }
}

bool checkAtomCharge(const MCSAtomCompareParameters&, const ROMol& mol1,
                     unsigned int atom1, const ROMol& mol2,
                     unsigned int atom2) {
  const auto a1 = mol1.getAtomWithIdx(atom1);
  const auto a2 = mol2.getAtomWithIdx(atom2);
  return a1->getFormalCharge() == a2->getFormalCharge();
}

bool checkAtomChirality(const MCSAtomCompareParameters&, const ROMol& mol1,
                        unsigned int atom1, const ROMol& mol2,
                        unsigned int atom2) {
  const auto a1 = mol1.getAtomWithIdx(atom1);
  const auto a2 = mol2.getAtomWithIdx(atom2);
  const auto ac1 = a1->getChiralTag();
  const auto ac2 = a2->getChiralTag();
  if (ac1 == Atom::CHI_TETRAHEDRAL_CW || ac1 == Atom::CHI_TETRAHEDRAL_CCW) {
    return (ac2 == Atom::CHI_TETRAHEDRAL_CW ||
            ac2 == Atom::CHI_TETRAHEDRAL_CCW);
  }
  return true;
}

bool checkAtomDistance(const MCSAtomCompareParameters& p, const ROMol& mol1,
                       unsigned int atom1, const ROMol& mol2,
                       unsigned int atom2) {
  const auto& ci1 = mol1.getConformer();
  const auto& ci2 = mol2.getConformer();
  const auto& pos1 = ci1.getAtomPos(atom1);
  const auto& pos2 = ci2.getAtomPos(atom2);
  bool withinRange = (pos1 - pos2).length() <= p.MaxDistance;
  return withinRange;
}

bool MCSAtomCompareAny(const MCSAtomCompareParameters& p, const ROMol& mol1,
                       unsigned int atom1, const ROMol& mol2,
                       unsigned int atom2, void*) {
  if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MaxDistance > 0 && !checkAtomDistance(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.RingMatchesRingOnly) {
    return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
  }

  return true;
}

bool MCSAtomCompareElements(const MCSAtomCompareParameters& p,
                            const ROMol& mol1, unsigned int atom1,
                            const ROMol& mol2, unsigned int atom2, void*) {
  const auto a1 = mol1.getAtomWithIdx(atom1);
  const auto a2 = mol2.getAtomWithIdx(atom2);
  if (a1->getAtomicNum() != a2->getAtomicNum()) {
    return false;
  }
  if (p.MatchValences && a1->getTotalValence() != a2->getTotalValence()) {
    return false;
  }
  if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MaxDistance > 0 && !checkAtomDistance(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.RingMatchesRingOnly) {
    return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
  }
  return true;
}

bool MCSAtomCompareIsotopes(const MCSAtomCompareParameters& p,
                            const ROMol& mol1, unsigned int atom1,
                            const ROMol& mol2, unsigned int atom2, void*) {
  // ignore everything except isotope information:
  const auto a1 = mol1.getAtomWithIdx(atom1);
  const auto a2 = mol2.getAtomWithIdx(atom2);
  if (a1->getIsotope() != a2->getIsotope()) {
    return false;
  }
  if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MaxDistance > 0 && !checkAtomDistance(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.RingMatchesRingOnly) {
    return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
  }
  return true;
}

bool MCSAtomCompareAnyHeavyAtom(const MCSAtomCompareParameters& p,
                                const ROMol& mol1, unsigned int atom1,
                                const ROMol& mol2, unsigned int atom2, void*) {
  const auto a1 = mol1.getAtomWithIdx(atom1);
  const auto a2 = mol2.getAtomWithIdx(atom2);
  // Any atom, including H, matches another atom of the same type,  according to
  // the other flags
  if (a1->getAtomicNum() == a2->getAtomicNum() ||
      (a1->getAtomicNum() > 1 && a2->getAtomicNum() > 1)) {
    return MCSAtomCompareAny(p, mol1, atom1, mol2, atom2, nullptr);
  }
  return false;
}

//=== BOND COMPARE ========================================================

class BondMatchOrderMatrix {
  bool MatchMatrix[Bond::ZERO + 1][Bond::ZERO + 1];

 public:
  BondMatchOrderMatrix(bool ignoreAromatization) {
    memset(MatchMatrix, 0, sizeof(MatchMatrix));
    // fill cells of the same and unspecified type
    for (size_t i = 0; i <= Bond::ZERO; ++i) {
      MatchMatrix[i][i] = true;
      MatchMatrix[Bond::UNSPECIFIED][i] = MatchMatrix[i][Bond::UNSPECIFIED] =
          true;
      MatchMatrix[Bond::ZERO][i] = MatchMatrix[i][Bond::ZERO] = true;
    }
    if (ignoreAromatization) {
      MatchMatrix[Bond::SINGLE][Bond::AROMATIC] =
          MatchMatrix[Bond::AROMATIC][Bond::SINGLE] = true;
      MatchMatrix[Bond::SINGLE][Bond::ONEANDAHALF] =
          MatchMatrix[Bond::ONEANDAHALF][Bond::SINGLE] = true;
      MatchMatrix[Bond::DOUBLE][Bond::TWOANDAHALF] =
          MatchMatrix[Bond::TWOANDAHALF][Bond::DOUBLE] = true;
      MatchMatrix[Bond::TRIPLE][Bond::THREEANDAHALF] =
          MatchMatrix[Bond::THREEANDAHALF][Bond::TRIPLE] = true;
      MatchMatrix[Bond::QUADRUPLE][Bond::FOURANDAHALF] =
          MatchMatrix[Bond::FOURANDAHALF][Bond::QUADRUPLE] = true;
      MatchMatrix[Bond::QUINTUPLE][Bond::FIVEANDAHALF] =
          MatchMatrix[Bond::FIVEANDAHALF][Bond::QUINTUPLE] = true;
    }
  }
  inline bool isEqual(unsigned int i, unsigned int j) const {
    return MatchMatrix[i][j];
  }
};

bool checkBondStereo(const MCSBondCompareParameters&, const ROMol& mol1,
                     unsigned int bond1, const ROMol& mol2,
                     unsigned int bond2) {
  const auto b1 = mol1.getBondWithIdx(bond1);
  const auto b2 = mol2.getBondWithIdx(bond2);
  auto bs1 = b1->getStereo();
  auto bs2 = b2->getStereo();
  if (b1->getBondType() == Bond::DOUBLE && b2->getBondType() == Bond::DOUBLE) {
    if (bs1 > Bond::STEREOANY && !(bs2 > Bond::STEREOANY)) {
      return false;
    } else {
      return bs1 == bs2;
    }
  }
  return true;
}

bool havePairOfCompatibleRings(const MCSBondCompareParameters&,
                               const ROMol& mol1, unsigned int bond1,
                               const ROMol& mol2, unsigned int bond2) {
  const auto ri1 = mol1.getRingInfo();
  const auto ri2 = mol2.getRingInfo();
  const auto& bondRings1 = ri1->bondRings();
  const auto& bondRings2 = ri2->bondRings();
  for (unsigned int ringIdx1 : ri1->bondMembers(bond1)) {
    const auto& ring1 = bondRings1.at(ringIdx1);
    bool isRing1Fused = ri1->isRingFused(ringIdx1);
    for (unsigned int ringIdx2 : ri2->bondMembers(bond2)) {
      const auto& ring2 = bondRings2.at(ringIdx2);
      if (ring1.size() == ring2.size()) {
        return true;
      }
      if (isRing1Fused && ring2.size() > ring1.size()) {
        return true;
      }
      bool isRing2Fused = ri2->isRingFused(ringIdx2);
      if (isRing2Fused && ring1.size() > ring2.size()) {
        return true;
      }
    }
  }
  return false;
}

bool checkBondRingMatch(const MCSBondCompareParameters& p, const ROMol& mol1,
                        unsigned int bond1, const ROMol& mol2,
                        unsigned int bond2) {
  const auto ri1 = mol1.getRingInfo();
  const auto ri2 = mol2.getRingInfo();
  // indices of rings in the query molecule
  const auto& ringIndices1 = ri1->bondMembers(bond1);
  // indices of rings in the target molecule
  const auto& ringIndices2 = ri2->bondMembers(bond2);
  bool bond1inRing = !ringIndices1.empty();
  bool bond2inRing = !ringIndices2.empty();
  bool res = (bond1inRing == bond2inRing);
  // if rings should be complete, we need to check upfront that there
  // is at least one pair of compatible rings; if there isn't, there
  // will never be a chance of complete match, so we should fail early
  if (p.CompleteRingsOnly && bond1inRing && bond2inRing) {
    res = havePairOfCompatibleRings(p, mol1, bond1, mol2, bond2);
  }

  // bond are both either in a ring or not
  return res;
}

bool MCSBondCompareAny(const MCSBondCompareParameters& p, const ROMol& mol1,
                       unsigned int bond1, const ROMol& mol2,
                       unsigned int bond2, void*) {
  if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
    return false;
  }
  if (p.RingMatchesRingOnly) {
    return checkBondRingMatch(p, mol1, bond1, mol2, bond2);
  }
  return true;
}

bool MCSBondCompareOrder(const MCSBondCompareParameters& p, const ROMol& mol1,
                         unsigned int bond1, const ROMol& mol2,
                         unsigned int bond2, void*) {
  static const BondMatchOrderMatrix match(true);  // ignore Aromatization
  const auto b1 = mol1.getBondWithIdx(bond1);
  const auto b2 = mol2.getBondWithIdx(bond2);
  auto t1 = b1->getBondType();
  auto t2 = b2->getBondType();
  if (match.isEqual(t1, t2)) {
    if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
      return false;
    }
    if (p.RingMatchesRingOnly) {
      return checkBondRingMatch(p, mol1, bond1, mol2, bond2);
    }
    return true;
  }
  return false;
}

bool MCSBondCompareOrderExact(const MCSBondCompareParameters& p,
                              const ROMol& mol1, unsigned int bond1,
                              const ROMol& mol2, unsigned int bond2, void*) {
  static const BondMatchOrderMatrix match(false);  // AROMATIC != SINGLE
  const auto b1 = mol1.getBondWithIdx(bond1);
  const auto b2 = mol2.getBondWithIdx(bond2);
  auto t1 = b1->getBondType();
  auto t2 = b2->getBondType();
  if (match.isEqual(t1, t2)) {
    if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
      return false;
    }
    if (p.RingMatchesRingOnly) {
      return checkBondRingMatch(p, mol1, bond1, mol2, bond2);
    }
    return true;
  }
  return false;
}

//=== RING COMPARE ========================================================
namespace {
// there are 3 bitsets for each ring
// isMCSRingBond: bits are set in correspondence of ring bond indices which are
// part of MCS.
// isMCSRingBondFused: bits are set in correspondence of ring bond
// indices which are part of MCS and are fused.
// isMCSRingBondNonFused: bits are set in correspondence of ring bond indices
// which are part of MCS and are not fused.
class RingBondCountVect {
 public:
  RingBondCountVect(const ROMol& mol) :
    d_mol(mol),
    d_ringInfo(mol.getRingInfo()) {
    d_ringBondCountVect.resize(d_ringInfo->numRings());
    d_isMCSBond.resize(mol.getNumBonds());
    for (auto& ringBondCount : d_ringBondCountVect) {
      ringBondCount.isMCSRingBond.resize(mol.getNumBonds());
      ringBondCount.isMCSRingBondNonFused.resize(mol.getNumBonds());
      ringBondCount.isMCSRingBondFused.resize(mol.getNumBonds());
    }
  }
  // In the 1st pass, we set fused/non-fused bits simply based on
  // numBondRings.
  void setMCSBondBitsPass1(unsigned int beginAtomIdx, unsigned int endAtomIdx,
                          const std::uint32_t c[],
                          const FMCS::Graph& graph) {
    const auto bond =
        d_mol.getBondBetweenAtoms(graph[c[beginAtomIdx]], graph[c[endAtomIdx]]);
    CHECK_INVARIANT(bond, "");
    const auto bi = bond->getIdx();
    d_isMCSBond.set(bi);
    if (!d_ringInfo->numBondRings(bi)) {
      return;
    }
    if (d_ringInfo->numBondRings(bi) == 1) {
      const auto ringIdx = d_ringInfo->bondMembers(bi).front();
      d_ringBondCountVect[ringIdx].isMCSRingBond.set(bi);
      d_ringBondCountVect[ringIdx].isMCSRingBondNonFused.set(bi);
    } else {
      for (const auto& ringIdx : d_ringInfo->bondMembers(bi)) {
        d_ringBondCountVect[ringIdx].isMCSRingBond.set(bi);
        d_ringBondCountVect[ringIdx].isMCSRingBondFused.set(bi);
      }
    }
  }
  // In the 2nd pass, we refine the above as certain bonds originally
  // marked as fused can be relabelled as non-fused
  void setMCSBondBitsPass2() {
    for (auto& ringBondCount : d_ringBondCountVect) {
      ringBondCount.nonFusedCountPass1 = ringBondCount.isMCSRingBondNonFused.count();
      ringBondCount.fusedCountPass1 = ringBondCount.isMCSRingBondFused.count();
    }
    for (unsigned int bi = 0; bi < d_mol.getNumBonds(); ++bi) {
      if (!d_isMCSBond.test(bi)) {
        continue;
      }
      int fusedBondRingIdx = -1;
      unsigned int fusedBondCount = 0;
      for (auto& ringBondCount : d_ringBondCountVect) {
        if (!ringBondCount.isMCSRingBondFused.test(bi)) {
          continue;
        }
        auto ringIdx = &ringBondCount - &d_ringBondCountVect.front();
        if (ringBondCount.nonFusedCountPass1 == 0 &&
            ringBondCount.fusedCountPass1 <
                d_ringInfo->bondRings().at(ringIdx).size()) {
          ringBondCount.isMCSRingBondFused.set(bi, false);
        } else {
          ++fusedBondCount;
          fusedBondRingIdx = ringIdx;
        }
      }
      if (fusedBondCount == 1) {
        d_ringBondCountVect[fusedBondRingIdx].isMCSRingBondNonFused.set(bi);
        d_ringBondCountVect[fusedBondRingIdx].isMCSRingBondFused.set(bi, false);
      }
    }
  }
  bool isRingFusionHonored() {
    for (const auto& ringBondCount : d_ringBondCountVect) {
      unsigned int ringIdx = &ringBondCount - &d_ringBondCountVect.front();
      const auto& bondRings = d_ringInfo->bondRings().at(ringIdx);
      // if all or no bonds of this ring are part of MCS, no need to do further
      // checks
      const auto numRingBondsInMCS = ringBondCount.isMCSRingBond.count();
      if (!numRingBondsInMCS || numRingBondsInMCS == bondRings.size()) {
        continue;
      }
      // check ring bonds:
      // if they are non-fused, it's OK
      // if they are fused but classified as non-fused, it's OK
      // otherwise count a missing fused bond
      // if the sum of missing fused bonds + fused bonds in MCS + non-fused bonds
      // in MCS equals to the ring size, then we have failed the check
      const auto numNonFusedRingBondsInMCS =
          ringBondCount.isMCSRingBondNonFused.count();
      const auto numFusedRingBondsInMCS =
          ringBondCount.isMCSRingBondFused.count();
      const auto numMissingFusedBonds =
          std::count_if(bondRings.begin(), bondRings.end(),
                        [this, &ringBondCount](const auto bi) {
                          return (d_ringInfo->numBondRings(bi) > 1 &&
                                  !ringBondCount.isMCSRingBondNonFused.test(bi) &&
                                  !ringBondCount.isMCSRingBondFused.test(bi));
                        });
      if (numMissingFusedBonds + numFusedRingBondsInMCS +
              numNonFusedRingBondsInMCS ==
          bondRings.size()) {
        return false;
      }
    }
    return true;
  }
 private:
  struct BondCount {
    boost::dynamic_bitset<> isMCSRingBond;
    boost::dynamic_bitset<> isMCSRingBondNonFused;
    unsigned int nonFusedCountPass1 = 0;
    boost::dynamic_bitset<> isMCSRingBondFused;
    unsigned int fusedCountPass1 = 0;
  };
  const ROMol& d_mol;
  const RingInfo *d_ringInfo;
  std::vector<BondCount> d_ringBondCountVect;
  boost::dynamic_bitset<> d_isMCSBond;
};
}  // end of anonymous namespace

inline bool ringFusionCheck(const std::uint32_t c1[], const std::uint32_t c2[],
                            const ROMol& mol1, const FMCS::Graph& query,
                            const ROMol& mol2, const FMCS::Graph& target,
                            const MCSParameters& p) {
  /*
  Consider this case: MCS between
  2-methylbicyclo[4.3.0]nonane and 1-methylbicyclo[3.1.0]hexane

                \__
                /  \_                          \___
                \__/ \                         / \/
                   \ /                         \ /
                    C                           C
                    H2                          H2

  In permissive mode, we are happy for methylcyclohexane to be the MCS.
  in strict mode, we don't want methylcyclohexane to be the MCS.

                                      \__
                                      /  \
                                      \__/

  When methylcyclohexane is checked against 2-methylbicyclo[4.3.0]nonane
  there is no missing fused bond. This is OK for permissive mode.
  In strict mode, we also need to check against 1-methylbicyclo[3.1.0]hexane,
  where there is indeed a missing fused bond.
  Basically, in permissive mode one of two molecules is allowed to fail the
  match, but not both. In strict mode, none is.
  */
  bool res = true;
  if (boost::num_edges(target) < boost::num_edges(query)) {
    return true;
  }
  RingBondCountVect mol1RingBondCountVect(mol1);
  RingBondCountVect mol2RingBondCountVect(mol2);
  auto queryEdges = boost::edges(query);
  std::for_each(
      queryEdges.first, queryEdges.second,
      [&c1, &c2, &query, &target, &mol1RingBondCountVect,
       &mol2RingBondCountVect](const auto& edge) {
        const auto beginAtomIdx = boost::source(edge, query);
        const auto endAtomIdx = boost::target(edge, query);
        mol1RingBondCountVect.setMCSBondBitsPass1(beginAtomIdx, endAtomIdx, c1, query);
        mol2RingBondCountVect.setMCSBondBitsPass1(beginAtomIdx, endAtomIdx, c2, target);
      });
  mol1RingBondCountVect.setMCSBondBitsPass2();
  mol2RingBondCountVect.setMCSBondBitsPass2();
  bool mol1Honored = mol1RingBondCountVect.isRingFusionHonored();
  bool mol2Honored = mol2RingBondCountVect.isRingFusionHonored();
  if (p.BondCompareParameters.MatchFusedRingsStrict) {
    res = mol1Honored && mol2Honored;
  } else {
    res = mol1Honored || mol2Honored;
  }
  return res;
}

bool FinalMatchCheckFunction(const std::uint32_t c1[], const std::uint32_t c2[],
                             const ROMol& mol1, const FMCS::Graph& query,
                             const ROMol& mol2, const FMCS::Graph& target,
                             const MCSParameters* p) {
  PRECONDITION(p, "p must not be NULL");
  if ((p->BondCompareParameters.MatchFusedRings ||
       p->BondCompareParameters.MatchFusedRingsStrict) &&
      !ringFusionCheck(c1, c2, mol1, query, mol2, target, *p)) {
    return false;
  }
  if (p->AtomCompareParameters.MatchChiralTag &&
      !FinalChiralityCheckFunction(c1, c2, mol1, query, mol2, target, p)) {
    return false;
  }
  const auto ip = dynamic_cast<const detail::MCSParametersInternal*>(p);
  if (ip && ip->UserFinalMatchChecker) {
    return ip->UserFinalMatchChecker(c1, c2, mol1, query, mol2, target, p);
  }
  return true;
}

bool FinalChiralityCheckFunction(const std::uint32_t c1[],
                                 const std::uint32_t c2[], const ROMol& mol1,
                                 const FMCS::Graph& query, const ROMol& mol2,
                                 const FMCS::Graph& target,
                                 const MCSParameters* /*unused*/) {
  const unsigned int qna = boost::num_vertices(query);  // getNumAtoms()
  // check chiral atoms only:
  for (unsigned int i = 0; i < qna; ++i) {
    const auto a1 = mol1.getAtomWithIdx(query[c1[i]]);
    const auto ac1 = a1->getChiralTag();

    const auto a2 = mol2.getAtomWithIdx(target[c2[i]]);
    const auto ac2 = a2->getChiralTag();

    ///*------------------ OLD Code :
    // ???: non chiral query atoms ARE ALLOWED TO MATCH to Chiral target atoms
    // (see test for issue 481)
    if (a1->getDegree() <
            3 ||  // #688: doesn't deal with "explicit" Hs properly
        !(ac1 == Atom::CHI_TETRAHEDRAL_CW ||
          ac1 == Atom::CHI_TETRAHEDRAL_CCW)) {
      continue;  // skip non chiral center QUERY atoms
    }
    if (!(ac2 == Atom::CHI_TETRAHEDRAL_CW ||
          ac2 == Atom::CHI_TETRAHEDRAL_CCW)) {
      return false;
    }
    //--------------------
    /* More accurate check:

            if( !(ac1 == Atom::CHI_TETRAHEDRAL_CW || ac1 ==
       Atom::CHI_TETRAHEDRAL_CCW)
             && !(ac2 == Atom::CHI_TETRAHEDRAL_CW || ac2 ==
       Atom::CHI_TETRAHEDRAL_CCW))
                continue; // skip check if both atoms are non chiral center

            if(!(   (ac1 == Atom::CHI_TETRAHEDRAL_CW || ac1 ==
       Atom::CHI_TETRAHEDRAL_CCW)
                 && (ac2 == Atom::CHI_TETRAHEDRAL_CW || ac2 ==
       Atom::CHI_TETRAHEDRAL_CCW)))//ac2 != ac1)
                 return false; // both atoms must be chiral or not without a
       query priority
    */
    const unsigned int a1Degree =
        boost::out_degree(c1[i], query);  // a1.getDegree();
    // number of all connected atoms in a seed
    if (a1Degree > a2->getDegree()) {  // #688 was != . // FIX issue 631
      // printf("atoms Degree (%u, %u) %u [%u], %u\n", query[c1[i]],
      // target[c2[i]], a1Degree, a1.getDegree(), a2.getDegree());
      if (1 == a1Degree && a1->getDegree() == a2->getDegree()) {
        continue;  // continue to grow the seed
      } else {
        return false;
      }
    }

    INT_LIST qOrder;
    for (unsigned int j = 0; j < qna && qOrder.size() != a1Degree; ++j) {
      const auto qB = mol1.getBondBetweenAtoms(query[c1[i]], query[c1[j]]);
      if (qB) {
        qOrder.push_back(qB->getIdx());
      }
    }

    // #688
    INT_LIST qmoOrder;
    {
      for (const auto& nbri :
           boost::make_iterator_range(mol1.getAtomBonds(a1))) {
        int dbidx = mol1[nbri]->getIdx();
        if (std::find(qOrder.begin(), qOrder.end(), dbidx) != qOrder.end()) {
          qmoOrder.push_back(dbidx);
        }
        //            else
        //                qmoOrder.push_back(-1);
      }
    }
    int qPermCount =  // was: a1.getPerturbationOrder(qOrder);
        static_cast<int>(countSwapsToInterconvert(qmoOrder, qOrder));

    INT_LIST mOrder;
    for (unsigned int j = 0; j < qna && mOrder.size() != a2->getDegree(); ++j) {
      const auto mB = mol2.getBondBetweenAtoms(target[c2[i]], target[c2[j]]);
      if (mB) {
        mOrder.push_back(mB->getIdx());
      }
    }

    // #688
    while (mOrder.size() < a2->getDegree()) {
      mOrder.push_back(-1);
    }
    INT_LIST moOrder;
    for (const auto& nbri : boost::make_iterator_range(mol2.getAtomBonds(a2))) {
      int dbidx = mol2[nbri]->getIdx();
      if (std::find(mOrder.begin(), mOrder.end(), dbidx) != mOrder.end()) {
        moOrder.push_back(dbidx);
      } else {
        moOrder.push_back(-1);
      }
    }

    int mPermCount =  // was: a2.getPerturbationOrder(mOrder);
        static_cast<int>(countSwapsToInterconvert(moOrder, mOrder));
    //----

    if ((qPermCount % 2 == mPermCount % 2 &&
         a1->getChiralTag() != a2->getChiralTag()) ||
        (qPermCount % 2 != mPermCount % 2 &&
         a1->getChiralTag() == a2->getChiralTag())) {
      return false;
    }
  }

  // check double bonds ONLY (why ???)
  std::map<unsigned int, unsigned int> qMap;
  for (unsigned int j = 0; j < qna; ++j) {
    qMap[query[c1[j]]] = j;
  }
  for (const auto& bondIdx : boost::make_iterator_range(boost::edges(query))) {
    const auto qBnd = mol1.getBondWithIdx(query[bondIdx]);
    if (qBnd->getBondType() != Bond::DOUBLE ||
        qBnd->getStereo() <= Bond::STEREOANY) {
      continue;
    }
    // don't think this can actually happen, but check to be sure:
    if (qBnd->getStereoAtoms().size() != 2) {  // MUST check it in the seed, not
                                               // in full query molecule, but
                                               // never happens !!!
      continue;
    }

    const auto mBnd =
        mol2.getBondBetweenAtoms(target[c2[qMap[qBnd->getBeginAtomIdx()]]],
                                 target[c2[qMap[qBnd->getEndAtomIdx()]]]);
    CHECK_INVARIANT(mBnd, "Matching bond not found");
    if (mBnd->getBondType() != Bond::DOUBLE ||
        mBnd->getStereo() <= Bond::STEREOANY) {
      continue;
    }
    // don't think this can actually happen, but check to be sure:
    if (mBnd->getStereoAtoms().size() != 2) {
      continue;
    }

    unsigned int end1Matches = 0;
    unsigned int end2Matches = 0;
    if (target[c2[qMap[qBnd->getBeginAtomIdx()]]] ==
        rdcast<unsigned int>(mBnd->getBeginAtomIdx())) {
      // query Begin == mol Begin
      if (target[c2[qMap[qBnd->getStereoAtoms()[0]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[0])) {
        end1Matches = 1;
      }
      if (target[c2[qMap[qBnd->getStereoAtoms()[1]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[1])) {
        end2Matches = 1;
      }
    } else {
      // query End == mol Begin
      if (target[c2[qMap[qBnd->getStereoAtoms()[0]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[1])) {
        end1Matches = 1;
      }
      if (target[c2[qMap[qBnd->getStereoAtoms()[1]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[0])) {
        end2Matches = 1;
      }
    }
    // std::cerr<<"  bnd: "<<qBnd->getIdx()<<":"<<qBnd->getStereo()<<" -
    // "<<mBnd->getIdx()<<":"<<mBnd->getStereo()<<"  --  "<<end1Matches<<"
    // "<<end2Matches<<std::endl;
    if (mBnd->getStereo() == qBnd->getStereo() &&
        (end1Matches + end2Matches) == 1) {
      return false;
    }
    if (mBnd->getStereo() != qBnd->getStereo() &&
        (end1Matches + end2Matches) != 1) {
      return false;
    }
  }
  return true;
}

bool FinalChiralityCheckFunction_1(const short unsigned int c1[],
                                   const short unsigned int c2[],
                                   const ROMol& mol1, const FMCS::Graph& query,
                                   const ROMol& mol2, const FMCS::Graph& target,
                                   const MCSParameters*) {
  const unsigned int qna = boost::num_vertices(query);  // getNumAtoms()
  // check chiral atoms:
  for (unsigned int i = 0; i < qna; ++i) {
    const auto a1 = mol1.getAtomWithIdx(query[c1[i]]);
    const auto ac1 = a1->getChiralTag();
    if (!(ac1 == Atom::CHI_TETRAHEDRAL_CW ||
          ac1 == Atom::CHI_TETRAHEDRAL_CCW)) {
      continue;  // skip non chiral center query atoms
    }
    const auto a2 = mol2.getAtomWithIdx(target[c2[i]]);
    const auto ac2 = a2->getChiralTag();
    if (!(ac2 == Atom::CHI_TETRAHEDRAL_CW ||
          ac2 == Atom::CHI_TETRAHEDRAL_CCW)) {
      continue;  // skip non chiral center TARGET atoms even if query atom is
    }
    // chiral
    ////                return false;
    // both atoms are chiral:
    const unsigned int a1Degree =
        boost::out_degree(c1[i], query);  // a1.getDegree();
    if (a1Degree != a2->getDegree()) {  // number of all connected atoms in seed
      return false;                     // ???
    }
    INT_LIST qOrder;
    for (unsigned int j = 0; j < qna && qOrder.size() != a1Degree; ++j) {
      const auto qB = mol1.getBondBetweenAtoms(query[c1[i]], query[c1[j]]);
      if (qB) {
        qOrder.push_back(qB->getIdx());
      }
    }

    int qPermCount = a1->getPerturbationOrder(qOrder);
    INT_LIST mOrder;
    for (unsigned int j = 0; j < qna && mOrder.size() != a2->getDegree(); ++j) {
      const auto mB = mol2.getBondBetweenAtoms(target[c2[i]], target[c2[j]]);
      if (mB) {
        mOrder.push_back(mB->getIdx());
      }
    }
    int mPermCount = a2->getPerturbationOrder(mOrder);

    if ((qPermCount % 2 == mPermCount % 2 &&
         a1->getChiralTag() != a2->getChiralTag()) ||
        (qPermCount % 2 != mPermCount % 2 &&
         a1->getChiralTag() == a2->getChiralTag())) {
      return false;
    }
  }

  // check double bonds ONLY (why ???)
  std::map<unsigned int, unsigned int> qMap;
  for (unsigned int j = 0; j < qna; ++j) {
    qMap[query[c1[j]]] = j;
  }
  for (const auto& bondIdx : boost::make_iterator_range(boost::edges(query))) {
    const auto qBnd = mol1.getBondWithIdx(query[bondIdx]);
    if (qBnd->getBondType() != Bond::DOUBLE ||
        qBnd->getStereo() <= Bond::STEREOANY) {
      continue;
    }
    // don't think this can actually happen, but check to be sure:
    if (qBnd->getStereoAtoms().size() != 2) {  // MUST check it in the seed, not
                                               // in full query molecule, but
                                               // never happens !!!
      continue;
    }

    const Bond* mBnd =
        mol2.getBondBetweenAtoms(target[c2[qMap[qBnd->getBeginAtomIdx()]]],
                                 target[c2[qMap[qBnd->getEndAtomIdx()]]]);
    CHECK_INVARIANT(mBnd, "Matching bond not found");
    if (mBnd->getBondType() != Bond::DOUBLE ||
        mBnd->getStereo() <= Bond::STEREOANY) {
      continue;
    }
    // don't think this can actually happen, but check to be sure:
    if (mBnd->getStereoAtoms().size() != 2) {
      continue;
    }

    unsigned int end1Matches = 0;
    unsigned int end2Matches = 0;
    if (target[c2[qMap[qBnd->getBeginAtomIdx()]]] == mBnd->getBeginAtomIdx()) {
      // query Begin == mol Begin
      if (target[c2[qMap[qBnd->getStereoAtoms()[0]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[0])) {
        end1Matches = 1;
      }
      if (target[c2[qMap[qBnd->getStereoAtoms()[1]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[1])) {
        end2Matches = 1;
      }
    } else {
      // query End == mol Begin
      if (target[c2[qMap[qBnd->getStereoAtoms()[0]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[1])) {
        end1Matches = 1;
      }
      if (target[c2[qMap[qBnd->getStereoAtoms()[1]]]] ==
          rdcast<unsigned int>(mBnd->getStereoAtoms()[0])) {
        end2Matches = 1;
      }
    }
    // std::cerr<<"  bnd: "<<qBnd->getIdx()<<":"<<qBnd->getStereo()<<" -
    // "<<mBnd->getIdx()<<":"<<mBnd->getStereo()<<"  --  "<<end1Matches<<"
    // "<<end2Matches<<std::endl;
    if (mBnd->getStereo() == qBnd->getStereo() &&
        (end1Matches + end2Matches) == 1) {
      return false;
    }
    if (mBnd->getStereo() != qBnd->getStereo() &&
        (end1Matches + end2Matches) != 1) {
      return false;
    }
  }
  return true;
}

namespace detail {
MCSParametersInternal::MCSParametersInternal(const MCSParameters* params)
    : MCSParameters(params) {
  UserFinalMatchChecker = FinalMatchChecker;
  FinalMatchChecker = FinalMatchCheckFunction;
}
}  // end namespace detail

}  // namespace RDKit
