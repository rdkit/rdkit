//
//  Copyright (C) 2001-2023 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "AlignMolecules.h"
#include <Geometry/Transform3D.h>
#include <Numerics/Vector.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/QueryOps.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <Numerics/Alignment/AlignPoints.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <RDGeneral/RDThreads.h>

namespace RDKit {
namespace MolAlign {

namespace details {
void symmetrizeTerminalAtoms(RWMol &mol) {
  // clang-format off
  static const std::string qsmarts =
      "[{atomPattern};$([{atomPattern}]-[*]=[{atomPattern}]),$([{atomPattern}]=[*]-[{atomPattern}])]~[*]";
  static std::map<std::string, std::string> replacements = {
      {"{atomPattern}", "O,N;D1"}};
  // clang-format on
  static SmartsParserParams ps;
  ps.replacements = &replacements;
  static const std::unique_ptr<RWMol> qry{SmartsToMol(qsmarts, ps)};
  CHECK_INVARIANT(qry, "bad query pattern");

  auto matches = SubstructMatch(mol, *qry);
  if (matches.empty()) {
    return;
  }

  QueryBond qb;
  qb.setQuery(makeSingleOrDoubleBondQuery());
  for (const auto &match : matches) {
    mol.getAtomWithIdx(match[0].second)->setFormalCharge(0);
    auto obond = mol.getBondBetweenAtoms(match[0].second, match[1].second);
    CHECK_INVARIANT(obond, "could not find expected bond");
    mol.replaceBond(obond->getIdx(), &qb);
  }
}
}  // namespace details
namespace {
double alignConfsOnAtomMap(const Conformer &prbCnf, const Conformer &refCnf,
                           const MatchVectType &atomMap,
                           RDGeom::Transform3D &trans,
                           const RDNumeric::DoubleVector *weights, bool reflect,
                           unsigned int maxIterations) {
  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  for (const auto &mi : atomMap) {
    prbPoints.push_back(&prbCnf.getAtomPos(mi.first));
    refPoints.push_back(&refCnf.getAtomPos(mi.second));
  }
  double ssr = RDNumeric::Alignments::AlignPoints(
      refPoints, prbPoints, trans, weights, reflect, maxIterations);
  return ssr / static_cast<double>(prbPoints.size());
}

void getAllMatchesPrbRef(const ROMol &prbMol, const ROMol &refMol,
                         std::vector<MatchVectType> &matches, int maxMatches,
                         bool symmetrizeConjugatedTerminalGroups) {
  bool uniquify = false;
  bool recursionPossible = true;
  bool useChirality = false;
  bool useQueryQueryMatches = false;

  std::unique_ptr<RWMol> prbMolSymm;
  if (symmetrizeConjugatedTerminalGroups) {
    prbMolSymm.reset(new RWMol(prbMol));
    details::symmetrizeTerminalAtoms(*prbMolSymm);
  }
  const auto &prbMolForMatch = prbMolSymm ? *prbMolSymm : prbMol;
  SubstructMatch(refMol, prbMolForMatch, matches, uniquify, recursionPossible,
                 useChirality, useQueryQueryMatches, maxMatches);

  if (matches.empty()) {
    throw MolAlignException(
        "No sub-structure match found between the reference and probe mol");
  }
  if (matches.size() > 1e6) {
    std::string name;
    prbMol.getPropIfPresent(common_properties::_Name, name);
    std::cerr << "Warning in " << __FUNCTION__ << ": " << matches.size()
              << " matches detected for molecule " << name << ", this may "
              << "lead to a performance slowdown.\n";
  }
}

double calcMSDInternal(const Conformer &prbCnf, const Conformer &refCnf,
                       const MatchVectType &atomMap,
                       const RDNumeric::DoubleVector *weights) {
  unsigned int npt = atomMap.size();
  std::unique_ptr<RDNumeric::DoubleVector> unitWeights;
  if (!weights) {
    unitWeights.reset(new RDNumeric::DoubleVector(npt, 1.0));
    weights = unitWeights.get();
  } else {
    PRECONDITION(npt == weights->size(), "Mismatch in number of weights");
  }
  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  for (const auto &mi : atomMap) {
    prbPoints.push_back(&prbCnf.getAtomPos(mi.first));
    refPoints.push_back(&refCnf.getAtomPos(mi.second));
  }
  double ssr = 0.;
  const RDGeom::Point3D *rpt;
  const RDGeom::Point3D *ppt;
  for (unsigned int i = 0; i < npt; ++i) {
    rpt = refPoints[i];
    ppt = prbPoints[i];
    ssr += (*weights)[i] * (*ppt - *rpt).lengthSq();
  }
  return ssr / static_cast<double>(npt);
}

double getBestRMSInternal(const ROMol &prbMol, const ROMol &refMol, int prbCid,
                          int refCid, const std::vector<MatchVectType> &matches,
                          RDGeom::Transform3D *trans, MatchVectType *bestMatch,
                          const RDNumeric::DoubleVector *weights, bool reflect,
                          unsigned int maxIters, unsigned int numThreads) {
  PRECONDITION(!matches.empty(), "matches must not be empty");
#ifndef RDK_BUILD_THREADSAFE_SSS
  numThreads = 1;
#endif
  double msdBest = std::numeric_limits<double>::max();
  const Conformer &prbCnf = prbMol.getConformer(prbCid);
  const Conformer &refCnf = refMol.getConformer(refCid);
  const MatchVectType *bestMatchPtr = &matches[0];

  if (numThreads == 1) {
    for (const auto &matche : matches) {
      RDGeom::Transform3D tmpTrans;
      double msd = trans ? alignConfsOnAtomMap(prbCnf, refCnf, matche, tmpTrans,
                                               weights, reflect, maxIters)
                         : calcMSDInternal(prbCnf, refCnf, matche, weights);
      if (msd < msdBest) {
        msdBest = msd;
        bestMatchPtr = &matche;
        if (trans) {
          trans->assign(tmpTrans);
        }
      }
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    std::vector<std::thread> tg;
    std::vector<
        std::vector<std::tuple<double, unsigned int, RDGeom::Transform3D>>>
        rmsds(numThreads);
    for (auto ti = 0u; ti < numThreads; ++ti) {
      auto func = [&](unsigned int tidx) {
        for (auto midx = tidx; midx < matches.size(); midx += numThreads) {
          auto matche = matches[midx];
          RDGeom::Transform3D tmpTrans;
          auto msd = trans
                         ? alignConfsOnAtomMap(prbCnf, refCnf, matche, tmpTrans,
                                               weights, reflect, maxIters)
                         : calcMSDInternal(prbCnf, refCnf, matche, weights);
          rmsds[tidx].emplace_back(msd, midx, tmpTrans);
        }
      };
      tg.emplace_back(std::thread(func, ti));
    }
    for (auto &thread : tg) {
      if (thread.joinable()) {
        thread.join();
      }
    }
    for (const auto &rv : rmsds) {
      for (const auto &res : rv) {
        const auto &[msd, midx, tf] = res;
        if (msd < msdBest) {
          msdBest = msd;
          bestMatchPtr = &matches[midx];
          if (trans) {
            trans->assign(tf);
          }
        }
      }
    }
  }
#endif
  if (bestMatch) {
    *bestMatch = *bestMatchPtr;
  }
  return sqrt(msdBest);
}
}  // namespace

double getAlignmentTransform(const ROMol &prbMol, const ROMol &refMol,
                             RDGeom::Transform3D &trans, int prbCid, int refCid,
                             const MatchVectType *atomMap,
                             const RDNumeric::DoubleVector *weights,
                             bool reflect, unsigned int maxIterations) {
  const Conformer &prbCnf = prbMol.getConformer(prbCid);
  const Conformer &refCnf = refMol.getConformer(refCid);
  MatchVectType match;
  if (!atomMap) {
    // we have to figure out the mapping between the two molecule
    const bool recursionPossible = true;
    const bool useChirality = false;
    const bool useQueryQueryMatches = true;
    if (SubstructMatch(refMol, prbMol, match, recursionPossible, useChirality,
                       useQueryQueryMatches)) {
      atomMap = &match;
    } else {
      throw MolAlignException(
          "No sub-structure match found between the probe and query mol");
    }
  }
  double msd = alignConfsOnAtomMap(prbCnf, refCnf, *atomMap, trans, weights,
                                   reflect, maxIterations);
  return sqrt(msd);
}

double alignMol(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
                const MatchVectType *atomMap,
                const RDNumeric::DoubleVector *weights, bool reflect,
                unsigned int maxIterations) {
  RDGeom::Transform3D trans;
  double res = getAlignmentTransform(prbMol, refMol, trans, prbCid, refCid,
                                     atomMap, weights, reflect, maxIterations);
  // now transform the relevant conformation on prbMol
  Conformer &conf = prbMol.getConformer(prbCid);
  MolTransforms::transformConformer(conf, trans);
  return res;
}

double getBestAlignmentTransform(
    const ROMol &prbMol, const ROMol &refMol, RDGeom::Transform3D &bestTrans,
    MatchVectType &bestMatch, int prbCid, int refCid,
    const std::vector<MatchVectType> &map, int maxMatches,
    bool symmetrizeConjugatedTerminalGroups,
    const RDNumeric::DoubleVector *weights, bool reflect, unsigned int maxIters,
    int numThreads) {
  std::vector<MatchVectType> allMatches;
  if (map.empty()) {
    getAllMatchesPrbRef(prbMol, refMol, allMatches, maxMatches,
                        symmetrizeConjugatedTerminalGroups);
  }
  const auto &matches = map.empty() ? allMatches : map;
  auto bestRMS = getBestRMSInternal(prbMol, refMol, prbCid, refCid, matches,
                                    &bestTrans, &bestMatch, weights, reflect,
                                    maxIters, getNumThreadsToUse(numThreads));
  return bestRMS;
}

double getBestRMS(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
                  const std::vector<MatchVectType> &map, int maxMatches,
                  bool symmetrizeConjugatedTerminalGroups,
                  const RDNumeric::DoubleVector *weights, int numThreads) {
  std::vector<MatchVectType> allMatches;
  if (map.empty()) {
    getAllMatchesPrbRef(prbMol, refMol, allMatches, maxMatches,
                        symmetrizeConjugatedTerminalGroups);
  }
  const auto &matches = map.empty() ? allMatches : map;
  RDGeom::Transform3D trans;
  bool reflect = false;
  unsigned int maxIters = 50;
  auto bestRMS = getBestRMSInternal(prbMol, refMol, prbCid, refCid, matches,
                                    &trans, nullptr, weights, reflect, maxIters,
                                    getNumThreadsToUse(numThreads));

  // Perform a final alignment to the best alignment...
  MolTransforms::transformConformer(prbMol.getConformer(prbCid), trans);
  return bestRMS;
}

std::vector<double> getAllConformerBestRMS(
    const ROMol &mol, int numThreads, const std::vector<MatchVectType> &map,
    int maxMatches, bool symmetrizeConjugatedTerminalGroups,
    const RDNumeric::DoubleVector *weights) {
  numThreads = getNumThreadsToUse(numThreads);
  std::vector<MatchVectType> allMatches;
  if (map.empty()) {
    getAllMatchesPrbRef(mol, mol, allMatches, maxMatches,
                        symmetrizeConjugatedTerminalGroups);
  }
  const auto &matches = map.empty() ? allMatches : map;
  std::vector<double> res;
  RDGeom::Transform3D trans;
  bool reflect = false;
  unsigned int maxIters = 50;
  std::vector<int> cids;
  for (auto cit = mol.beginConformers(); cit != mol.endConformers(); ++cit) {
    cids.push_back((*cit)->getId());
  }
  if (numThreads == 1) {
    for (auto ci = 0u; ci < mol.getNumConformers(); ++ci) {
      for (auto cj = 0u; cj < ci; ++cj) {
        res.push_back(getBestRMSInternal(mol, mol, cids[ci], cids[cj], matches,
                                         &trans, nullptr, weights, reflect,
                                         maxIters, 1));
      }
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    std::vector<std::pair<unsigned int, unsigned int>> pairs;
    for (auto ci = 0u; ci < mol.getNumConformers(); ++ci) {
      for (auto cj = 0u; cj < ci; ++cj) {
        pairs.emplace_back(cids[ci], cids[cj]);
      }
    }
    std::vector<std::vector<std::pair<unsigned int, double>>> rmsds(numThreads);
    auto func = [&](unsigned int tidx) {
      RDGeom::Transform3D trans;
      bool reflect = false;
      unsigned int maxIters = 50;
      for (auto i = tidx; i < pairs.size(); i += numThreads) {
        auto rms = getBestRMSInternal(mol, mol, pairs[i].first, pairs[i].second,
                                      matches, &trans, nullptr, weights,
                                      reflect, maxIters, 1);
        rmsds[tidx].emplace_back(i, rms);
      }
    };
    std::vector<std::thread> tg;
    for (auto ti = 0; ti < numThreads; ++ti) {
      tg.emplace_back(std::thread(func, ti));
    }
    for (auto &thread : tg) {
      if (thread.joinable()) {
        thread.join();
      }
    }
    res.resize(pairs.size());
    for (const auto &tres : rmsds) {
      for (const auto &v : tres) {
        res[v.first] = v.second;
      }
    }
  }
#endif
  return res;
}

double CalcRMS(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
               const std::vector<MatchVectType> &map, int maxMatches,
               bool symmetrizeConjugatedTerminalGroups,
               const RDNumeric::DoubleVector *weights) {
  std::vector<MatchVectType> allMatches;
  if (map.empty()) {
    getAllMatchesPrbRef(prbMol, refMol, allMatches, maxMatches,
                        symmetrizeConjugatedTerminalGroups);
  }
  const auto &matches = map.empty() ? allMatches : map;
  bool reflect = false;
  unsigned int maxIters = 50;
  unsigned int numThreads = 1;
  return getBestRMSInternal(prbMol, refMol, prbCid, refCid, matches, nullptr,
                            nullptr, weights, reflect, maxIters, numThreads);
}

double CalcRMS(ROMol &prbMol, const ROMol &refMol, int prbCid, int refCid,
               const std::vector<MatchVectType> &map, int maxMatches,
               const RDNumeric::DoubleVector *weights) {
  return CalcRMS(prbMol, refMol, prbCid, refCid, map, maxMatches, false,
                 weights);
}

void _fillAtomPositions(RDGeom::Point3DConstPtrVect &pts, const Conformer &conf,
                        const std::vector<unsigned int> *atomIds = nullptr) {
  unsigned int na = conf.getNumAtoms();
  pts.clear();
  if (atomIds == nullptr) {
    unsigned int ai;
    pts.reserve(na);
    for (ai = 0; ai < na; ++ai) {
      pts.push_back(&conf.getAtomPos(ai));
    }
  } else {
    pts.reserve(atomIds->size());
    std::vector<unsigned int>::const_iterator cai;
    for (cai = atomIds->begin(); cai != atomIds->end(); cai++) {
      pts.push_back(&conf.getAtomPos(*cai));
    }
  }
}

void alignMolConformers(ROMol &mol, const std::vector<unsigned int> *atomIds,
                        const std::vector<unsigned int> *confIds,
                        const RDNumeric::DoubleVector *weights, bool reflect,
                        unsigned int maxIters, std::vector<double> *RMSlist) {
  if (mol.getNumConformers() == 0) {
    // nothing to be done ;
    return;
  }

  RDGeom::Point3DConstPtrVect refPoints, prbPoints;
  int cid = -1;
  if ((confIds != nullptr) && (confIds->size() > 0)) {
    cid = confIds->front();
  }
  const Conformer &refCnf = mol.getConformer(cid);
  _fillAtomPositions(refPoints, refCnf, atomIds);

  // now loop throught the remaininf conformations and transform them
  RDGeom::Transform3D trans;
  double ssd;
  if (confIds == nullptr) {
    unsigned int i = 0;
    ROMol::ConformerIterator cnfi;
    // Conformer *conf;
    for (cnfi = mol.beginConformers(); cnfi != mol.endConformers(); cnfi++) {
      // conf = (*cnfi);
      i += 1;
      if (i == 1) {
        continue;
      }
      _fillAtomPositions(prbPoints, *(*cnfi), atomIds);
      ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
                                               weights, reflect, maxIters);
      if (RMSlist) {
        ssd /= (prbPoints.size());
        RMSlist->push_back(sqrt(ssd));
      }
      MolTransforms::transformConformer(*(*cnfi), trans);
    }
  } else {
    std::vector<unsigned int>::const_iterator cai;
    unsigned int i = 0;
    for (cai = confIds->begin(); cai != confIds->end(); cai++) {
      i += 1;
      if (i == 1) {
        continue;
      }
      Conformer &conf = mol.getConformer(*cai);
      _fillAtomPositions(prbPoints, conf, atomIds);
      ssd = RDNumeric::Alignments::AlignPoints(refPoints, prbPoints, trans,
                                               weights, reflect, maxIters);
      if (RMSlist) {
        ssd /= (prbPoints.size());
        RMSlist->push_back(sqrt(ssd));
      }
      MolTransforms::transformConformer(conf, trans);
    }
  }
}
}  // namespace MolAlign
}  // namespace RDKit
