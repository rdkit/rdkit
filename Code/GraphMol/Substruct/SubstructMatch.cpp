//
//  Copyright (C) 2001-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/utils.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDThreads.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Resonance.h>
#include <GraphMol/MolBundle.h>
#include <GraphMol/Chirality.h>

#include "SubstructMatch.h"
#include "SubstructUtils.h"
#include <GraphMol/GenericGroups/GenericGroups.h>
#include <boost/smart_ptr.hpp>
#include <map>

#if BOOST_VERSION == 106400
#include <boost/serialization/array_wrapper.hpp>
#endif

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <mutex>
#include <thread>
#include <future>
#endif

#include "vf2.hpp"

using boost::make_iterator_range;

namespace RDKit {
namespace detail {

namespace {
bool hasChiralLabel(const Atom *at) {
  PRECONDITION(at, "bad atom");
  return at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
         at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW;
}

bool enhancedStereoIsOK(
    const ROMol &mol, const ROMol &query,
    std::unordered_map<unsigned int, unsigned int> &q_to_mol,
    const std::unordered_map<unsigned int, StereoGroup const *>
        &molStereoGroups,
    const std::unordered_map<unsigned int, bool> &matches) {
  std::unordered_map<unsigned int, StereoGroup const *> molAtomsToQueryGroups;

  // If the query has stereo groups:
  // * OR only matches AND or OR (not absolute)
  // * AND only matches OR
  for (auto &&sg : query.getStereoGroups()) {
    if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
      continue;
    }
    // StereoGroup const* matched_mol_group = nullptr;
    const bool is_and = sg.getGroupType() == StereoGroupType::STEREO_AND;
    for (auto &&a : sg.getAtoms()) {
      auto mol_group = molStereoGroups.find(q_to_mol[a->getIdx()]);
      if (mol_group == molStereoGroups.end()) {
        // group matching absolute. not ok.
        return false;
      } else if (is_and && mol_group->second->getGroupType() !=
                               StereoGroupType::STEREO_AND) {
        // AND matching OR. not ok.
        return false;
      }

      molAtomsToQueryGroups[q_to_mol[a->getIdx()]] = &sg;
    }
  }

  // If the mol has stereo groups:
  // * All atoms must either be the same or opposite, you can't mix
  // * Only one stereogroup must cover all matched atoms in the mol stereo group
  for (auto &&sg : mol.getStereoGroups()) {
    if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
      continue;
    }
    bool doesMatch;
    bool seen = false;
    StereoGroup const *QGroup = nullptr;

    for (auto &&a : sg.getAtoms()) {
      auto thisDoesMatch = matches.find(a->getIdx());
      if (thisDoesMatch == matches.end()) {
        // not matched
        continue;
      }

      auto pos = molAtomsToQueryGroups.find(a->getIdx());
      auto thisQGroup =
          pos == molAtomsToQueryGroups.end() ? nullptr : pos->second;
      if (!seen) {
        doesMatch = thisDoesMatch->second;
        QGroup = thisQGroup;
        seen = true;
      } else if (doesMatch != thisDoesMatch->second) {
        // diastereomer. not ok.
        return false;
      } else if (thisQGroup != QGroup) {
        // mix of groups in query. not ok.
        return false;
      }
    }
  }

  return true;
}

}  // namespace

typedef std::map<unsigned int, QueryAtom::QUERYATOM_QUERY *> SUBQUERY_MAP;

typedef struct {
  ResonanceMolSupplier &resMolSupplier;
  const ROMol &query;
  const SubstructMatchParameters &params;
} ResSubstructMatchHelperArgs_;

void MatchSubqueries(const ROMol &mol, QueryAtom::QUERYATOM_QUERY *q,
                     const SubstructMatchParameters &params,
                     SUBQUERY_MAP &subqueryMap,
                     std::vector<RecursiveStructureQuery *> &locked);

bool insertIfNeeded(std::set<MatchVectType> &matches, const MatchVectType &m) {
  bool shouldInsert = true;
  std::unordered_set<int> matchAsSet;
  std::transform(m.begin(), m.end(),
                 std::inserter(matchAsSet, matchAsSet.begin()),
                 [](const std::pair<int, int> &p) { return p.second; });
  for (auto it = matches.begin(); it != matches.end(); ++it) {
    std::unordered_set<int> existingMatchAsSet;
    std::transform(
        it->begin(), it->end(),
        std::inserter(existingMatchAsSet, existingMatchAsSet.begin()),
        [](const std::pair<int, int> &p) { return p.second; });
    if (matchAsSet == existingMatchAsSet) {
      if (m < *it) {
        matches.erase(it);
      } else {
        shouldInsert = false;
      }
      break;
    }
  }
  if (shouldInsert) {
    matches.insert(m);
  }
  return shouldInsert;
}

bool tryToInsert(std::set<MatchVectType> &matches, const MatchVectType &match,
                 const SubstructMatchParameters &params) {
  if (matches.size() == params.maxMatches) {
    return false;
  }
  if (!params.uniquify) {
    matches.insert(match);
  } else {
    insertIfNeeded(matches, match);
  }
  return true;
}

void ResSubstructMatchHelper_(const ResSubstructMatchHelperArgs_ &args,
                              std::set<MatchVectType> *matches, unsigned int bi,
                              unsigned int ei);

typedef std::list<
    std::pair<MolGraph::vertex_descriptor, MolGraph::vertex_descriptor>>
    ssPairType;

}  // namespace detail

MolMatchFinalCheckFunctor::MolMatchFinalCheckFunctor(
    const ROMol &query, const ROMol &mol, const SubstructMatchParameters &ps)
    : d_query(query), d_mol(mol), d_params(ps) {
  if (d_params.useEnhancedStereo) {
    for (const auto &sg : d_mol.getStereoGroups()) {
      if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
        continue;
      }
      for (const auto a : sg.getAtoms()) {
        d_molStereoGroups[a->getIdx()] = &sg;
      }
    }
  }
}

bool MolMatchFinalCheckFunctor::operator()(const std::uint32_t q_c[],
                                           const std::uint32_t m_c[]) {
  if (d_params.extraFinalCheck || d_params.useGenericMatchers) {
    // EFF: we can no-doubt do better than this
    std::vector<unsigned int> aids(m_c, m_c + d_query.getNumAtoms());
    for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
      aids[i] = m_c[i];
    }
    if (d_params.useGenericMatchers &&
        !GenericGroups::genericAtomMatcher(d_mol, d_query, aids)) {
      return false;
    }
    if (d_params.extraFinalCheck && !d_params.extraFinalCheck(d_mol, aids)) {
      return false;
    }
  }

  HashedStorageType match;
  if (d_params.uniquify) {
    match.resize(d_mol.getNumAtoms());
#ifdef RDK_INTERNAL_BITSET_HAS_HASH
    match.reset();
#else
    std::fill(match.begin(), match.end(), 0);
#endif
    for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
      match[m_c[i]] = 1;
    }
    if (matchesSeen.find(match) != matchesSeen.end()) {
      return false;
    }
  }

  if (!d_params.useChirality) {
    if (d_params.uniquify) {
      matchesSeen.insert(match);
    }
    return true;
  }

  std::unordered_map<unsigned int, bool> matches;

  // check chiral atoms:
  for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
    const Atom *qAt = d_query.getAtomWithIdx(q_c[i]);

    // With less than 3 neighbors we can't establish CW/CCW parity,
    // so query will be a match if it has any kind of chirality.
    if (qAt->getDegree() < 3 || !detail::hasChiralLabel(qAt)) {
      continue;
    }
    const Atom *mAt = d_mol.getAtomWithIdx(m_c[i]);
    if (!detail::hasChiralLabel(mAt)) {
      return false;
    }
    if (qAt->getDegree() > mAt->getDegree()) {
      return false;
    }

    INT_LIST qOrder;
    INT_LIST mOrder;
    for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
      const Bond *qB = d_query.getBondBetweenAtoms(q_c[i], q_c[j]);
      const Bond *mB = d_mol.getBondBetweenAtoms(m_c[i], m_c[j]);
      if (qB && mB) {
        mOrder.push_back(mB->getIdx());
        qOrder.push_back(qB->getIdx());
        if (mOrder.size() == qAt->getDegree()) {
          break;
        }
      }
    }
    CHECK_INVARIANT(qOrder.size() == qAt->getDegree(), "missing matches");
    CHECK_INVARIANT(qOrder.size() == mOrder.size(), "bad matches");
    int qPermCount = qAt->getPerturbationOrder(qOrder);

    unsigned unmatchedNeighbors = mAt->getDegree() - mOrder.size();
    mOrder.insert(mOrder.end(), unmatchedNeighbors, -1);

    INT_LIST moOrder;
    for (const auto &bond : make_iterator_range(d_mol.getAtomBonds(mAt))) {
      int dbidx = d_mol[bond]->getIdx();
      if (std::find(mOrder.begin(), mOrder.end(), dbidx) != mOrder.end()) {
        moOrder.push_back(dbidx);
      } else {
        moOrder.push_back(-1);
      }
    }

    int mPermCount =
        static_cast<int>(countSwapsToInterconvert(moOrder, mOrder));

    const bool requireMatch = qPermCount % 2 == mPermCount % 2;
    const bool labelsMatch = qAt->getChiralTag() == mAt->getChiralTag();
    const bool matchOK = requireMatch == labelsMatch;

    // if this is not part of a stereogroup and doesn't match, return false
    auto msg = d_molStereoGroups.find(m_c[i]);
    if (msg == d_molStereoGroups.end()) {
      if (!matchOK) {
        return false;
      }
    } else {
      matches[m_c[i]] = matchOK;
    }
  }

  std::unordered_map<unsigned int, unsigned int> q_to_mol;
  for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
    q_to_mol[q_c[j]] = m_c[j];
  }

  if (d_params.useEnhancedStereo) {
    if (!detail::enhancedStereoIsOK(d_mol, d_query, q_to_mol, d_molStereoGroups,
                                    matches)) {
      return false;
    }
  }

  // now check double bonds
  for (const auto &qBnd : d_query.bonds()) {
    if (qBnd->getBondType() != Bond::DOUBLE ||
        qBnd->getStereo() <= Bond::STEREOANY) {
      continue;
    }

    // don't think this can actually happen, but check to be sure:
    if (qBnd->getStereoAtoms().size() != 2) {
      continue;
    }

    const Bond *mBnd = d_mol.getBondBetweenAtoms(
        q_to_mol[qBnd->getBeginAtomIdx()], q_to_mol[qBnd->getEndAtomIdx()]);
    CHECK_INVARIANT(mBnd, "Matching bond not found");
    if (mBnd->getBondType() != Bond::DOUBLE ||
        qBnd->getStereo() <= Bond::STEREOANY) {
      continue;
    }
    // don't think this can actually happen, but check to be sure:
    if (mBnd->getStereoAtoms().size() != 2) {
      continue;
    }

    unsigned int end1Matches = 0;
    unsigned int end2Matches = 0;
    if (q_to_mol[qBnd->getBeginAtomIdx()] == mBnd->getBeginAtomIdx()) {
      // query Begin == mol Begin
      if (q_to_mol[qBnd->getStereoAtoms()[0]] ==
          static_cast<unsigned>(mBnd->getStereoAtoms()[0])) {
        end1Matches = 1;
      }
      if (q_to_mol[qBnd->getStereoAtoms()[1]] ==
          static_cast<unsigned>(mBnd->getStereoAtoms()[1])) {
        end2Matches = 1;
      }
    } else {
      // query End == mol Begin
      if (q_to_mol[qBnd->getStereoAtoms()[0]] ==
          static_cast<unsigned>(mBnd->getStereoAtoms()[1])) {
        end1Matches = 1;
      }
      if (q_to_mol[qBnd->getStereoAtoms()[1]] ==
          static_cast<unsigned>(mBnd->getStereoAtoms()[0])) {
        end2Matches = 1;
      }
    }

    const unsigned totalMatches = end1Matches + end2Matches;
    const auto mStereo =
        Chirality::translateEZLabelToCisTrans(mBnd->getStereo());
    const auto qStereo =
        Chirality::translateEZLabelToCisTrans(qBnd->getStereo());

    if (mStereo == qStereo && totalMatches == 1) {
      return false;
    }
    if (mStereo != qStereo && totalMatches != 1) {
      return false;
    }
  }
  if (d_params.uniquify) {
    matchesSeen.insert(match);
  }
  return true;
}

namespace detail {

class AtomLabelFunctor {
 public:
  AtomLabelFunctor(const ROMol &query, const ROMol &mol,
                   const SubstructMatchParameters &ps)
      : d_query(query), d_mol(mol), d_params(ps){};
  bool operator()(unsigned int i, unsigned int j) const {
    bool res = false;
    if (d_params.useChirality) {
      const Atom *qAt = d_query.getAtomWithIdx(i);
      if (qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
          qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
        const Atom *mAt = d_mol.getAtomWithIdx(j);
        if (mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
            mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW) {
          return false;
        }
      }
    }
    res = atomCompat(d_query[i], d_mol[j], d_params);
    return res;
  }

 private:
  const ROMol &d_query;
  const ROMol &d_mol;
  const SubstructMatchParameters &d_params;
};
class BondLabelFunctor {
 public:
  BondLabelFunctor(const ROMol &query, const ROMol &mol,
                   const SubstructMatchParameters &ps)
      : d_query(query), d_mol(mol), d_params(ps){};
  bool operator()(MolGraph::edge_descriptor i,
                  MolGraph::edge_descriptor j) const {
    if (d_params.useChirality) {
      const Bond *qBnd = d_query[i];
      if (qBnd->getBondType() == Bond::DOUBLE &&
          qBnd->getStereo() > Bond::STEREOANY) {
        const Bond *mBnd = d_mol[j];
        if (mBnd->getBondType() == Bond::DOUBLE &&
            mBnd->getStereo() <= Bond::STEREOANY) {
          return false;
        }
      }
    }
    bool res = bondCompat(d_query[i], d_mol[j], d_params);
    return res;
  }

 private:
  const ROMol &d_query;
  const ROMol &d_mol;
  const SubstructMatchParameters &d_params;
};
void ResSubstructMatchHelper_(const ResSubstructMatchHelperArgs_ &args,
                              std::set<MatchVectType> *matches, unsigned int bi,
                              unsigned int ei) {
  for (unsigned int i = bi;
       (matches->size() < args.params.maxMatches) && (i < ei); ++i) {
    ROMol *mol = args.resMolSupplier[i];
    std::vector<MatchVectType> matchesTmp =
        SubstructMatch(*mol, args.query, args.params);
    for (const auto &match : matchesTmp) {
      if (!tryToInsert(*matches, match, args.params)) {
        break;
      }
    }
    delete mol;
  }
};

struct RecursiveLocker {
  std::vector<RecursiveStructureQuery *> locked;
  RecursiveLocker(const ROMol &query, const bool recursionPossible) {
    if (recursionPossible) {
      locked.reserve(query.getNumAtoms());
    }
  }

  ~RecursiveLocker() {
    for (auto v : locked) {
      v->clear();
#ifdef RDK_BUILD_THREADSAFE_SSS
      v->d_mutex.unlock();
#endif
    }
  }
};
}  // namespace detail

// ----------------------------------------------
//
// find all matches
std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const ROMol &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> matches;
  if (!mol.getNumAtoms() || !query.getNumAtoms()) {
    return matches;
  }

  detail::RecursiveLocker locker(query, params.recursionPossible);

  if (params.recursionPossible) {
    detail::SUBQUERY_MAP subqueryMap;
    ROMol::ConstAtomIterator atIt;
    for (const auto atom : query.atoms()) {
      if (atom->hasQuery()) {
        // std::cerr<<"recurse from atom "<<(*atIt)->getIdx()<<std::endl;
        detail::MatchSubqueries(mol, atom->getQuery(), params, subqueryMap,
                                locker.locked);
      }
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, params);
  detail::BondLabelFunctor bondLabeler(query, mol, params);
  MolMatchFinalCheckFunctor matchChecker(query, mol, params);

  std::list<detail::ssPairType> pms;
#if 0
  bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
				atomLabeler,bondLabeler,pms);
#else
  bool found =
      boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
                     bondLabeler, matchChecker, pms, params.maxMatches);
#endif
  if (found) {
    unsigned int nQueryAtoms = query.getNumAtoms();
    matches.reserve(pms.size());
    MatchVectType matchVect(nQueryAtoms);
    for (const auto &pairs : pms) {
      for (const auto &pair : pairs) {
        matchVect[pair.first] = pair;
      }
      matches.push_back(matchVect);
    }
  }
  return matches;
}

std::vector<MatchVectType> SubstructMatch(
    const MolBundle &bundle, const ROMol &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> res;
  for (unsigned int i = 0; i < bundle.size() && !res.size(); ++i) {
    res = SubstructMatch(*bundle[i], query, params);
  }
  return res;
}

std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const MolBundle &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> res;
  for (unsigned int i = 0; i < query.size() && !res.size(); ++i) {
    res = SubstructMatch(mol, *query[i], params);
  }
  return res;
}

std::vector<MatchVectType> SubstructMatch(
    const MolBundle &mol, const MolBundle &query,
    const SubstructMatchParameters &params) {
  std::vector<MatchVectType> res;
  for (unsigned int i = 0; i < mol.size() && !res.size(); ++i) {
    for (unsigned int j = 0; j < query.size() && !res.size(); ++j) {
      res = SubstructMatch(*mol[i], *query[j], params);
    }
  }
  return res;
}

// ----------------------------------------------
//
// find all matches in a ResonanceMolSupplier object
//
//
std::vector<MatchVectType> SubstructMatch(
    ResonanceMolSupplier &resMolSupplier, const ROMol &query,
    const SubstructMatchParameters &params) {
  std::set<MatchVectType> matches;
  detail::ResSubstructMatchHelperArgs_ args = {resMolSupplier, query, params};
  unsigned int nt =
      std::min(resMolSupplier.length(), getNumThreadsToUse(params.numThreads));
  if (nt == 1) {
    detail::ResSubstructMatchHelper_(args, &matches, 0,
                                     resMolSupplier.length());
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    std::vector<std::future<void>> tg;
    std::vector<std::set<MatchVectType> *> matchesThread(nt);
    unsigned int ei = 0;
    double dpt =
        static_cast<double>(resMolSupplier.length()) / static_cast<double>(nt);
    double dc = 0.0;
    for (unsigned int ti = 0; ti < nt; ++ti) {
      matchesThread[ti] = new std::set<MatchVectType>();
      unsigned int bi = ei;
      dc += dpt;
      ei = static_cast<unsigned int>(floor(dc));
      tg.emplace_back(std::async(std::launch::async,
                                 detail::ResSubstructMatchHelper_, args,
                                 matchesThread[ti], bi, ei));
    }
    for (auto &fut : tg) {
      fut.get();
    }

    for (unsigned int ti = 0; ti < nt; ++ti) {
      for (const auto &match : *matchesThread[ti]) {
        if (!detail::tryToInsert(matches, match, args.params)) {
          break;
        }
      }
      delete matchesThread[ti];
    }
  }
#endif
  return std::vector<MatchVectType>(matches.begin(), matches.end());
}

namespace detail {
unsigned int RecursiveMatcher(const ROMol &mol, const ROMol &query,
                              std::vector<int> &matches,
                              SUBQUERY_MAP &subqueryMap,
                              const SubstructMatchParameters &params,
                              std::vector<RecursiveStructureQuery *> &locked) {
  SubstructMatchParameters lparams = params;
  lparams.maxMatches = std::max(params.maxRecursiveMatches, params.maxMatches);
  lparams.uniquify = false;
  for (auto qAtom : query.atoms()) {
    if (qAtom->hasQuery()) {
      MatchSubqueries(mol, qAtom->getQuery(), lparams, subqueryMap, locked);
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, lparams);
  detail::BondLabelFunctor bondLabeler(query, mol, lparams);
  MolMatchFinalCheckFunctor matchChecker(query, mol, lparams);

  matches.clear();
  matches.resize(0);
  std::list<detail::ssPairType> pms;
#if 0
      bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
				    atomLabeler,bondLabeler,pms);
#else
  bool found = boost::vf2_all(query.getTopology(), mol.getTopology(),
                              atomLabeler, bondLabeler, matchChecker, pms);
#endif
  unsigned int res = 0;
  if (found) {
    matches.reserve(pms.size());
    for (const auto &pairs : pms) {
      if (!query.hasProp(common_properties::_queryRootAtom)) {
        matches.push_back(pairs.begin()->second);
      } else {
        int rootIdx;
        query.getProp(common_properties::_queryRootAtom, rootIdx);
        bool found = false;
        for (const auto &pairIter : pairs) {
          if (pairIter.first == static_cast<unsigned int>(rootIdx)) {
            matches.push_back(pairIter.second);
            found = true;
            break;
          }
        }
        if (!found) {
          BOOST_LOG(rdErrorLog)
              << "no match found for queryRootAtom" << std::endl;
        }
      }
      if (matches.size() == lparams.maxMatches) {
        break;
      }
    }
    res = matches.size();
  }
  // std::cout << " <<< RecursiveMatcher: " << int(query) << std::endl;
  return res;
}

void MatchSubqueries(const ROMol &mol, QueryAtom::QUERYATOM_QUERY *query,
                     const SubstructMatchParameters &params,
                     SUBQUERY_MAP &subqueryMap,
                     std::vector<RecursiveStructureQuery *> &locked) {
  PRECONDITION(query, "bad query");
  // std::cout << "*-*-* MS: " << query << std::endl;
  // std::cout << "\t\t" << typeid(*query).name() << std::endl;
  if (query->getDescription() == "RecursiveStructure") {
    auto *rsq = (RecursiveStructureQuery *)query;
#ifdef RDK_BUILD_THREADSAFE_SSS
    rsq->d_mutex.lock();
#endif
    locked.push_back(rsq);
    rsq->clear();
    bool matchDone = false;
    if (rsq->getSerialNumber() &&
        subqueryMap.find(rsq->getSerialNumber()) != subqueryMap.end()) {
      // we've matched an equivalent serial number before, just
      // copy in the matches:
      matchDone = true;
      auto orsq =
          (const RecursiveStructureQuery *)subqueryMap[rsq->getSerialNumber()];
      for (auto setIter = orsq->beginSet(); setIter != orsq->endSet();
           ++setIter) {
        rsq->insert(*setIter);
      }
      // std::cerr<<" copying results for query serial number:
      // "<<rsq->getSerialNumber()<<std::endl;
    }

    if (!matchDone) {
      ROMol const *queryMol = rsq->getQueryMol();
      // in case we are reusing this query, clear its contents now.
      if (queryMol) {
        std::vector<int> matchStarts;
        unsigned int res = RecursiveMatcher(mol, *queryMol, matchStarts,
                                            subqueryMap, params, locked);
        if (res) {
          for (int &matchStart : matchStarts) {
            rsq->insert(matchStart);
          }
        }
      }
      if (rsq->getSerialNumber()) {
        subqueryMap[rsq->getSerialNumber()] = query;
        // std::cerr << " storing results for query serial number: "
        //           << rsq->getSerialNumber() << " " << rsq->size() <<
        //           std::endl;
      }
    }
  } else {
    // std::cout << "\tmsq1: ";
  }

  // now recurse over our children (these things can be nested)
  for (auto childIt = query->beginChildren(); childIt != query->endChildren();
       ++childIt) {
    MatchSubqueries(mol, childIt->get(), params, subqueryMap, locked);
  }
  // std::cout << "<<- back " << (int)query << std::endl;
}

}  // end of namespace detail
}  // namespace RDKit
