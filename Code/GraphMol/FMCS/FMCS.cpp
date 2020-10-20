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

void MCSParameters::setMCSAtomTyperFromConstChar(const char *atomComp) {
  PRECONDITION(atomComp, "atomComp must not be NULL");
  static const std::map<const char *, AtomComparator> atomCompStringToEnum = {
    { "Any", AtomCompareAny },
    { "Elements", AtomCompareElements },
    { "Isotopes", AtomCompareIsotopes },
    { "AnyHeavy", AtomCompareAnyHeavyAtom }
  };
  try {
    setMCSAtomTyperFromEnum(atomCompStringToEnum.at(atomComp));
  }
  catch(const std::out_of_range&) {
    // we accept "def" as a no-op
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

void MCSParameters::setMCSBondTyperFromConstChar(const char *bondComp) {
  PRECONDITION(bondComp, "bondComp must not be NULL");
  static const std::map<const char *, BondComparator> bondCompStringToEnum = {
    { "Any", BondCompareAny },
    { "Order", BondCompareOrder },
    { "OrderExact", BondCompareOrderExact }
  };
  try {
    setMCSBondTyperFromEnum(bondCompStringToEnum.at(bondComp));
  }
  catch(const std::out_of_range&) {
    // we accept "def" as a no-op
  }
}

void parseMCSParametersJSON(const char* json, MCSParameters* params) {
  if (params && json && 0 != strlen(json)) {
    std::istringstream ss;
    ss.str(json);
    boost::property_tree::ptree pt;
    boost::property_tree::read_json(ss, pt);

    RDKit::MCSParameters& p = *params;
    p.MaximizeBonds = pt.get<bool>("MaximizeBonds", p.MaximizeBonds);
    p.Threshold = pt.get<double>("Threshold", p.Threshold);
    p.Timeout = pt.get<unsigned>("Timeout", p.Timeout);
    p.AtomCompareParameters.MatchValences =
        pt.get<bool>("MatchValences", p.AtomCompareParameters.MatchValences);
    p.AtomCompareParameters.MatchChiralTag =
        pt.get<bool>("MatchChiralTag", p.AtomCompareParameters.MatchChiralTag);
    p.AtomCompareParameters.MatchFormalCharge = pt.get<bool>(
        "MatchFormalCharge", p.AtomCompareParameters.MatchFormalCharge);
    p.AtomCompareParameters.RingMatchesRingOnly = pt.get<bool>(
        "RingMatchesRingOnly", p.AtomCompareParameters.RingMatchesRingOnly);
    p.BondCompareParameters.RingMatchesRingOnly = pt.get<bool>(
        "RingMatchesRingOnly", p.BondCompareParameters.RingMatchesRingOnly);
    p.BondCompareParameters.CompleteRingsOnly = pt.get<bool>(
        "CompleteRingsOnly", p.BondCompareParameters.CompleteRingsOnly);
    p.BondCompareParameters.MatchFusedRings = pt.get<bool>(
        "MatchFusedRings", p.BondCompareParameters.MatchFusedRings);
    p.BondCompareParameters.MatchFusedRingsStrict = pt.get<bool>(
        "MatchFusedRingsStrict", p.BondCompareParameters.MatchFusedRingsStrict);
    p.BondCompareParameters.MatchStereo =
        pt.get<bool>("MatchStereo", p.BondCompareParameters.MatchStereo);

    p.setMCSAtomTyperFromConstChar(pt.get<std::string>("AtomCompare", "def").c_str());
    p.setMCSBondTyperFromConstChar(pt.get<std::string>("BondCompare", "def").c_str());

    p.InitialSeed = pt.get<std::string>("InitialSeed", "");
  }
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
                  double threshold, unsigned timeout, bool verbose,
                  bool matchValences, bool ringMatchesRingOnly,
                  bool completeRingsOnly, bool matchChiralTag,
                  AtomComparator atomComp, BondComparator bondComp) {
  return findMCS(mols, maximizeBonds, threshold, timeout, verbose,
                 matchValences, ringMatchesRingOnly, completeRingsOnly,
                 matchChiralTag, atomComp, bondComp, IgnoreRingFusion);
}

MCSResult findMCS(const std::vector<ROMOL_SPTR>& mols, bool maximizeBonds,
                  double threshold, unsigned timeout, bool verbose,
                  bool matchValences, bool ringMatchesRingOnly,
                  bool completeRingsOnly, bool matchChiralTag,
                  AtomComparator atomComp, BondComparator bondComp,
                  RingComparator ringComp) {
  auto* ps = new MCSParameters();
  ps->MaximizeBonds = maximizeBonds;
  ps->Threshold = threshold;
  ps->Timeout = timeout;
  ps->Verbose = verbose;
  ps->setMCSAtomTyperFromEnum(atomComp);
  ps->AtomCompareParameters.MatchValences = matchValences;
  ps->AtomCompareParameters.MatchChiralTag = matchChiralTag;
  ps->AtomCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  ps->setMCSBondTyperFromEnum(bondComp);
  ps->BondCompareParameters.RingMatchesRingOnly = ringMatchesRingOnly;
  ps->BondCompareParameters.CompleteRingsOnly = completeRingsOnly;
  ps->BondCompareParameters.MatchFusedRings = (ringComp != IgnoreRingFusion);
  ps->BondCompareParameters.MatchFusedRingsStrict =
      (ringComp == StrictRingFusion);
  MCSResult res = findMCS(mols, ps);
  delete ps;
  return res;
}

bool MCSProgressCallbackTimeout(const MCSProgressData& stat,
                                const MCSParameters& params, void* userData) {
  RDUNUSED_PARAM(stat);
  PRECONDITION(userData, "userData must not be NULL");
  auto* t0 = (unsigned long long*)userData;
  unsigned long long t = nanoClock();
  return t - *t0 <= params.Timeout * 1000000ULL;
}

// PREDEFINED FUNCTORS:

//=== ATOM COMPARE ========================================================
bool checkAtomRingMatch(const MCSAtomCompareParameters& p, const ROMol& mol1,
                        unsigned int atom1, const ROMol& mol2,
                        unsigned int atom2) {
  if (p.RingMatchesRingOnly) {
    bool atom1inRing = queryIsAtomInRing(mol1.getAtomWithIdx(atom1));
    bool atom2inRing = queryIsAtomInRing(mol2.getAtomWithIdx(atom2));
    return atom1inRing == atom2inRing;
  } else {
    return true;
  }
}

bool checkAtomCharge(const MCSAtomCompareParameters& p,
                     const ROMol& mol1, unsigned int atom1,
                     const ROMol& mol2, unsigned int atom2) {
  RDUNUSED_PARAM(p);
  const Atom& a1 = *mol1.getAtomWithIdx(atom1);
  const Atom& a2 = *mol2.getAtomWithIdx(atom2);
  return a1.getFormalCharge() == a2.getFormalCharge();
}

bool checkAtomChirality(const MCSAtomCompareParameters& p,
                        const ROMol& mol1, unsigned int atom1,
                        const ROMol& mol2, unsigned int atom2) {
  RDUNUSED_PARAM(p);
  const Atom& a1 = *mol1.getAtomWithIdx(atom1);
  const Atom& a2 = *mol2.getAtomWithIdx(atom2);
  Atom::ChiralType ac1 = a1.getChiralTag();
  Atom::ChiralType ac2 = a2.getChiralTag();
  if (ac1 == Atom::CHI_TETRAHEDRAL_CW || ac1 == Atom::CHI_TETRAHEDRAL_CCW) {
    return (ac2 == Atom::CHI_TETRAHEDRAL_CW ||
            ac2 == Atom::CHI_TETRAHEDRAL_CCW);
  }
  return true;
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
  if (p.RingMatchesRingOnly) {
    return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
  }

  return true;
}

bool MCSAtomCompareElements(const MCSAtomCompareParameters& p,
                            const ROMol& mol1, unsigned int atom1,
                            const ROMol& mol2, unsigned int atom2, void*) {
  const Atom& a1 = *mol1.getAtomWithIdx(atom1);
  const Atom& a2 = *mol2.getAtomWithIdx(atom2);
  if (a1.getAtomicNum() != a2.getAtomicNum()) {
    return false;
  }
  if (p.MatchValences && a1.getTotalValence() != a2.getTotalValence()) {
    return false;
  }
  if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.RingMatchesRingOnly) {
    return checkAtomRingMatch(p, mol1, atom1, mol2, atom2);
  }
  return true;
}

bool MCSAtomCompareIsotopes(const MCSAtomCompareParameters& p,
                            const ROMol& mol1, unsigned int atom1,
                            const ROMol& mol2, unsigned int atom2, void* ud) {
  RDUNUSED_PARAM(ud);
  // ignore everything except isotope information:
  // if( ! MCSAtomCompareElements (p, mol1, atom1, mol2, atom2, ud))
  //    return false;
  const Atom& a1 = *mol1.getAtomWithIdx(atom1);
  const Atom& a2 = *mol2.getAtomWithIdx(atom2);
  if (a1.getIsotope() != a2.getIsotope()) {
    return false;
  }
  if (p.MatchChiralTag && !checkAtomChirality(p, mol1, atom1, mol2, atom2)) {
    return false;
  }
  if (p.MatchFormalCharge && !checkAtomCharge(p, mol1, atom1, mol2, atom2)) {
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
  const Atom& a1 = *mol1.getAtomWithIdx(atom1);
  const Atom& a2 = *mol2.getAtomWithIdx(atom2);
  // Any atom, including H, matches another atom of the same type,  according to
  // the other flags
  if (a1.getAtomicNum() == a2.getAtomicNum() ||
      (a1.getAtomicNum() > 1 && a2.getAtomicNum() > 1)) {
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
    for (size_t i = 0; i <= Bond::ZERO;
         i++) {  // fill cells of the same and unspecified type
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
  inline bool isEqual(unsigned i, unsigned j) const {
    return MatchMatrix[i][j];
  }
};

bool checkBondStereo(const MCSBondCompareParameters& p,
                     const ROMol& mol1, unsigned int bond1,
                     const ROMol& mol2, unsigned int bond2) {
  RDUNUSED_PARAM(p);
  const Bond* b1 = mol1.getBondWithIdx(bond1);
  const Bond* b2 = mol2.getBondWithIdx(bond2);
  Bond::BondStereo bs1 = b1->getStereo();
  Bond::BondStereo bs2 = b2->getStereo();
  if (b1->getBondType() == Bond::DOUBLE && b2->getBondType() == Bond::DOUBLE) {
    if (bs1 > Bond::STEREOANY && !(bs2 > Bond::STEREOANY)) {
      return false;
    }
  }
  return true;
}

bool checkBondRingMatch(const MCSBondCompareParameters&, const ROMol&,
                               unsigned int bond1, const ROMol& mol2,
                               unsigned int bond2, void* v_ringMatchMatrixSet) {
  if (!v_ringMatchMatrixSet) {
    throw "v_ringMatchMatrixSet is NULL";  // never
  }
  auto* ringMatchMatrixSet =
      static_cast<FMCS::RingMatchTableSet*>(v_ringMatchMatrixSet);

  const std::vector<size_t>& ringsIdx1 =
      ringMatchMatrixSet->getQueryBondRings(bond1);  // indices of rings
  const std::vector<size_t>& ringsIdx2 =
      ringMatchMatrixSet->getTargetBondRings(&mol2, bond2);  // indices of rings
  bool bond1inRing = !ringsIdx1.empty();
  bool bond2inRing = !ringsIdx2.empty();

  // bond are both either in a ring or not
  return (bond1inRing == bond2inRing);
}

bool MCSBondCompareAny(const MCSBondCompareParameters& p, const ROMol& mol1,
                       unsigned int bond1, const ROMol& mol2,
                       unsigned int bond2, void* ud) {
  if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
    return false;
  }
  if (p.RingMatchesRingOnly) {
    return checkBondRingMatch(p, mol1, bond1, mol2, bond2, ud);
  }
  return true;
}

bool MCSBondCompareOrder(const MCSBondCompareParameters& p, const ROMol& mol1,
                         unsigned int bond1, const ROMol& mol2,
                         unsigned int bond2, void* ud) {
  static const BondMatchOrderMatrix match(true);  // ignore Aromatization
  const Bond* b1 = mol1.getBondWithIdx(bond1);
  const Bond* b2 = mol2.getBondWithIdx(bond2);
  Bond::BondType t1 = b1->getBondType();
  Bond::BondType t2 = b2->getBondType();
  if (match.isEqual(t1, t2)) {
    if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
      return false;
    }
    if (p.RingMatchesRingOnly) {
      return checkBondRingMatch(p, mol1, bond1, mol2, bond2, ud);
    }
    return true;
  }
  return false;
}

bool MCSBondCompareOrderExact(const MCSBondCompareParameters& p,
                              const ROMol& mol1, unsigned int bond1,
                              const ROMol& mol2, unsigned int bond2, void* ud) {
  static const BondMatchOrderMatrix match(false);  // AROMATIC != SINGLE
  const Bond* b1 = mol1.getBondWithIdx(bond1);
  const Bond* b2 = mol2.getBondWithIdx(bond2);
  Bond::BondType t1 = b1->getBondType();
  Bond::BondType t2 = b2->getBondType();
  if (match.isEqual(t1, t2)) {
    if (p.MatchStereo && !checkBondStereo(p, mol1, bond1, mol2, bond2)) {
      return false;
    }
    if (p.RingMatchesRingOnly) {
      return checkBondRingMatch(p, mol1, bond1, mol2, bond2, ud);
    }
    return true;
  }
  return false;
}

//=== RING COMPARE ========================================================
inline bool ringFusionCheck(const std::uint32_t c1[],
                            const std::uint32_t c2[], const ROMol& mol1,
                            const FMCS::Graph& query, const ROMol& mol2,
                            const FMCS::Graph& target,
                            const MCSParameters* p) {
  PRECONDITION(p, "p must not be NULL");
  const RingInfo* ri2 = mol2.getRingInfo();
  const VECT_INT_VECT& br2 = ri2->bondRings();
  std::vector<size_t> nonFusedBonds(br2.size(), 0);
  std::vector<size_t> numMcsBondRings(br2.size(), 0);
  RDKit::FMCS::Graph::BOND_ITER_PAIR bpIter = boost::edges(query);
  size_t i = 0;
  // numMcsBondRings stores the number of bonds which
  // are part of the MCS for each ring
  for (auto it = bpIter.first; it != bpIter.second; ++it) {
    const Bond* b =
        mol2.getBondBetweenAtoms(target[c2[boost::source(*it, query)]],
                                 target[c2[boost::target(*it, query)]]);
    if (!b) {
      continue;
    }
    unsigned int bi = b->getIdx();
    if (!ri2->numBondRings(bi)) {
      continue;
    }
    for (i = 0; i < br2.size(); ++i) {
      if (std::find(br2[i].begin(), br2[i].end(), bi) != br2[i].end()) {
        ++numMcsBondRings[i];
      }
    }
  }
  // nonFusedBonds stores the number of non-fused bonds
  // (i.e., which belong to a single ring) for each ring
  for (i = 0; i < br2.size(); ++i) {
    for (auto bi : br2[i]) {
      if (ri2->numBondRings(bi) == 1) {
        ++nonFusedBonds[i];
      }
    }
  }
  /* if a ring has at least one bond which is part of the MCS,
     we need to check how many bonds are actually part of the MCS:
     if they are greater than the number of fused bonds but lower
     than the number of non-fused bonds, then the ring is not
     complete and we should return true straightaway
     We need to check that the number of bonds in a ring
     is greater than the number of fused bonds because that
     guarantees that this ring is actually involved in the
     MCS and not just adjacent to another ring which is
     indeed part of the MCS. */
  bool missingFusedBond = false;
  for (i = 0; i < br2.size(); ++i) {
    if (numMcsBondRings[i] &&
        numMcsBondRings[i] > (br2[i].size() - nonFusedBonds[i]) &&
        numMcsBondRings[i] < nonFusedBonds[i]) {
      if (p->BondCompareParameters.CompleteRingsOnly) {
        return true;
      } else {
        continue;
      }
    }
    if (numMcsBondRings[i] == nonFusedBonds[i] &&
        numMcsBondRings[i] < br2[i].size()) {
      missingFusedBond = true;
    }
  }
  /* If we found that the MCS is missing a fused bond, we may need to
     check against the smaller molecule as well.
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
     where there is indeed a missing fused bond. */
  bool missingFusedBond2 = false;
  if (missingFusedBond ^ p->BondCompareParameters.MatchFusedRingsStrict) {
    // if we are in permissive mode we allow one of the molecules to miss
    // fused bonds, but not both. In strict mode we allow neither.
    if (!p->BondCompareParameters.MatchFusedRingsStrict) {
      missingFusedBond = false;
    }
    const RingInfo* ri1 = mol1.getRingInfo();
    const VECT_INT_VECT& br1 = ri1->bondRings();
    std::set<unsigned int> mcsRingBondIdxSet;
    // put all MCS bond indices which belong to one or more rings in a set
    for (auto it = bpIter.first; it != bpIter.second; ++it) {
      const Bond* b =
          mol1.getBondBetweenAtoms(query[c1[boost::source(*it, query)]],
                                   query[c1[boost::target(*it, query)]]);
      if (b && ri1->numBondRings(b->getIdx())) {
        mcsRingBondIdxSet.insert(b->getIdx());
      }
    }
    for (i = 0; i < br1.size(); ++i) {
      size_t mcsBonds = 0;
      size_t mcsFusedBonds = 0;
      size_t nonFusedBonds = 0;
      // for each ring, check how many bonds belong to MCS (mcsBonds),
      // how many bonds in the MCS belong to multiple rings (mcsFusedBonds).
      // and how many bonds in each ring belong to a single ring (nonFusedBonds)
      for (auto bi : br1[i]) {
        bool isFused = (ri1->numBondRings(bi) > 1);
        if (!isFused) {
          ++nonFusedBonds;
        }
        if (std::find(mcsRingBondIdxSet.begin(), mcsRingBondIdxSet.end(), bi) !=
            mcsRingBondIdxSet.end()) {
          ++mcsBonds;
          if (isFused) {
            ++mcsFusedBonds;
          }
        }
      }
      if (!p->BondCompareParameters.CompleteRingsOnly &&
          mcsBonds > mcsFusedBonds &&
          mcsBonds - mcsFusedBonds < nonFusedBonds) {
        continue;
      }
      // if the ring is part of the MCS, and we found more bonds than the fused
      // ones (i.e., the ring is not simply adjacent to an MCS ring) but less
      // than the bonds which are part of this ring, then some fused bond is
      // missing.
      if (mcsBonds && mcsBonds > mcsFusedBonds && mcsBonds < br1[i].size()) {
        missingFusedBond2 = true;
        break;
      }
    }
  }
  // if both missingFusedBond and missingFusedBond2 are false, return true
  return (!missingFusedBond && !missingFusedBond2);
}

bool FinalMatchCheckFunction(const std::uint32_t c1[], const std::uint32_t c2[],
                             const ROMol& mol1, const FMCS::Graph& query,
                             const ROMol& mol2, const FMCS::Graph& target,
                             const MCSParameters* p) {
  PRECONDITION(p, "p must not be NULL");
  if ((p->BondCompareParameters.MatchFusedRings ||
       p->BondCompareParameters.MatchFusedRingsStrict) &&
      !ringFusionCheck(c1, c2, mol1, query, mol2, target, p)) {
    return false;
  }
  if (p->AtomCompareParameters.MatchChiralTag &&
      !FinalChiralityCheckFunction(c1, c2, mol1, query, mol2, target, p)) {
    return false;
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
    const Atom& a1 = *mol1.getAtomWithIdx(query[c1[i]]);
    Atom::ChiralType ac1 = a1.getChiralTag();

    const Atom& a2 = *mol2.getAtomWithIdx(target[c2[i]]);
    Atom::ChiralType ac2 = a2.getChiralTag();

    ///*------------------ OLD Code :
    // ???: non chiral query atoms ARE ALLOWED TO MATCH to Chiral target atoms
    // (see test for issue 481)
    if (a1.getDegree() < 3 ||  //#688: doesn't deal with "explicit" Hs properly
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
    const unsigned a1Degree =
        boost::out_degree(c1[i], query);  // a1.getDegree();
    // number of all connected atoms in a seed
    if (a1Degree > a2.getDegree()) {  //#688 was != . // FIX issue 631
      // printf("atoms Degree (%u, %u) %u [%u], %u\n", query[c1[i]],
      // target[c2[i]], a1Degree, a1.getDegree(), a2.getDegree());
      if (1 == a1Degree && a1.getDegree() == a2.getDegree()) {
        continue;  // continue to grow the seed
      } else {
        return false;
      }
    }

    INT_LIST qOrder;
    for (unsigned int j = 0; j < qna && qOrder.size() != a1Degree; ++j) {
      const Bond* qB = mol1.getBondBetweenAtoms(query[c1[i]], query[c1[j]]);
      if (qB) {
        qOrder.push_back(qB->getIdx());
      }
    }

    //#688
    INT_LIST qmoOrder;
    {
      ROMol::OEDGE_ITER dbeg, dend;
      boost::tie(dbeg, dend) = mol1.getAtomBonds(&a1);
      for (; dbeg != dend; dbeg++) {
        int dbidx = (*dbeg)->getIdx();
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
    for (unsigned int j = 0; j < qna && mOrder.size() != a2.getDegree(); ++j) {
      const Bond* mB = mol2.getBondBetweenAtoms(target[c2[i]], target[c2[j]]);
      if (mB) {
        mOrder.push_back(mB->getIdx());
      }
    }

    //#688
    while (mOrder.size() < a2.getDegree()) {
      mOrder.push_back(-1);
    }
    INT_LIST moOrder;
    ROMol::OEDGE_ITER dbeg, dend;
    boost::tie(dbeg, dend) = mol2.getAtomBonds(&a2);
    for (; dbeg != dend; dbeg++) {
      int dbidx = (*dbeg)->getIdx();
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
         a1.getChiralTag() != a2.getChiralTag()) ||
        (qPermCount % 2 != mPermCount % 2 &&
         a1.getChiralTag() == a2.getChiralTag())) {
      return false;
    }
  }

  // check double bonds ONLY (why ???)
  const unsigned int qnb = boost::num_edges(query);
  std::map<unsigned int, unsigned int> qMap;
  for (unsigned int j = 0; j < qna; ++j) {
    qMap[query[c1[j]]] = j;
  }
  RDKit::FMCS::Graph::BOND_ITER_PAIR bpIter = boost::edges(query);
  RDKit::FMCS::Graph::EDGE_ITER bIter = bpIter.first;
  for (unsigned int i = 0; i < qnb; i++, ++bIter) {
    const Bond* qBnd = mol1.getBondWithIdx(query[*bIter]);
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

bool FinalChiralityCheckFunction_1(const short unsigned c1[],
                                   const short unsigned c2[], const ROMol& mol1,
                                   const FMCS::Graph& query, const ROMol& mol2,
                                   const FMCS::Graph& target,
                                   const MCSParameters* p) {
  RDUNUSED_PARAM(p);
  const unsigned int qna = boost::num_vertices(query);  // getNumAtoms()
  // check chiral atoms:
  for (unsigned int i = 0; i < qna; ++i) {
    const Atom& a1 = *mol1.getAtomWithIdx(query[c1[i]]);
    Atom::ChiralType ac1 = a1.getChiralTag();
    if (!(ac1 == Atom::CHI_TETRAHEDRAL_CW ||
          ac1 == Atom::CHI_TETRAHEDRAL_CCW)) {
      continue;  // skip non chiral center query atoms
    }
    const Atom& a2 = *mol2.getAtomWithIdx(target[c2[i]]);
    Atom::ChiralType ac2 = a2.getChiralTag();
    if (!(ac2 == Atom::CHI_TETRAHEDRAL_CW ||
          ac2 == Atom::CHI_TETRAHEDRAL_CCW)) {
      continue;  // skip non chiral center TARGET atoms even if query atom is
    }
                 // chiral
                 ////                return false;
    // both atoms are chiral:
    const unsigned a1Degree =
        boost::out_degree(c1[i], query);  // a1.getDegree();
    if (a1Degree != a2.getDegree()) {  // number of all connected atoms in seed
      return false;                  // ???
    }
    INT_LIST qOrder;
    for (unsigned int j = 0; j < qna && qOrder.size() != a1Degree; ++j) {
      const Bond* qB = mol1.getBondBetweenAtoms(query[c1[i]], query[c1[j]]);
      if (qB) {
        qOrder.push_back(qB->getIdx());
      }
    }

    int qPermCount = a1.getPerturbationOrder(qOrder);
    INT_LIST mOrder;
    for (unsigned int j = 0; j < qna && mOrder.size() != a2.getDegree(); ++j) {
      const Bond* mB = mol2.getBondBetweenAtoms(target[c2[i]], target[c2[j]]);
      if (mB) {
        mOrder.push_back(mB->getIdx());
      }
    }
    int mPermCount = a2.getPerturbationOrder(mOrder);

    if ((qPermCount % 2 == mPermCount % 2 &&
         a1.getChiralTag() != a2.getChiralTag()) ||
        (qPermCount % 2 != mPermCount % 2 &&
         a1.getChiralTag() == a2.getChiralTag())) {
      return false;
    }
  }

  // check double bonds ONLY (why ???)
  const unsigned int qnb = boost::num_edges(query);
  std::map<unsigned int, unsigned int> qMap;
  for (unsigned int j = 0; j < qna; ++j) {
    qMap[query[c1[j]]] = j;
  }
  RDKit::FMCS::Graph::BOND_ITER_PAIR bpIter = boost::edges(query);
  RDKit::FMCS::Graph::EDGE_ITER bIter = bpIter.first;
  for (unsigned int i = 0; i < qnb; i++, ++bIter) {
    const Bond* qBnd = mol1.getBondWithIdx(query[*bIter]);
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

}  // namespace RDKit
