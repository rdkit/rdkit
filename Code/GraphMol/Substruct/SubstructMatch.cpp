//
//  Copyright (C) 2001-2019 Greg Landrum and Rational Discovery LLC
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
#include "GraphMol/Chirality.h"

#include "SubstructMatch.h"
#include "SubstructUtils.h"
#include <boost/smart_ptr.hpp>
#include <map>

#if BOOST_VERSION == 106400
#include <boost/serialization/array_wrapper.hpp>
#endif

#ifdef RDK_THREADSAFE_SSS
#include <mutex>
#include <thread>
#include <future>
#endif

#include "ullmann.hpp"
#include "vf2.hpp"

namespace RDKit {
namespace detail {
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

bool matchCompare(const std::pair<int, int> &a, const std::pair<int, int> &b);
bool matchVectCompare(const MatchVectType &a, const MatchVectType &b);
bool isToBeAddedToVector(std::vector<MatchVectType> &matches,
                         const MatchVectType &m);
void mergeMatchVect(std::vector<MatchVectType> &matches,
                    const std::vector<MatchVectType> &matchesTmp,
                    const ResSubstructMatchHelperArgs_ &args);
void ResSubstructMatchHelper_(const ResSubstructMatchHelperArgs_ &args,
                              std::vector<MatchVectType> *matches,
                              unsigned int bi, unsigned int ei);

typedef std::list<
    std::pair<MolGraph::vertex_descriptor, MolGraph::vertex_descriptor>>
    ssPairType;

class MolMatchFinalCheckFunctor {
 public:
  MolMatchFinalCheckFunctor(const ROMol &query, const ROMol &mol,
                            const SubstructMatchParameters &ps)
      : d_query(query), d_mol(mol), d_params(ps){};
  bool operator()(const boost::detail::node_id c1[],
                  const boost::detail::node_id c2[]) const {
    // std::cerr << "  check! " << df_useChirality << std::endl;
    if (!d_params.useChirality) return true;
    // for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
    //   std::cerr << "    " << c1[i] << " " << c2[i] << std::endl;
    // }

    // check chiral atoms:
    for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
      const Atom *qAt = d_query.getAtomWithIdx(c1[i]);
      if (qAt->getDegree() <
              3 ||  // FIX: doesn't deal with "explicit" Hs properly
          (qAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
           qAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW))
        continue;
      const Atom *mAt = d_mol.getAtomWithIdx(c2[i]);
      if (mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
          mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW)
        return false;
      if (qAt->getDegree() > mAt->getDegree()) return false;
      INT_LIST qOrder;
      for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
        const Bond *qB = d_query.getBondBetweenAtoms(c1[i], c1[j]);
        if (qB) {
          qOrder.push_back(qB->getIdx());
          if (qOrder.size() == qAt->getDegree()) break;
        }
      }
      int qPermCount = qAt->getPerturbationOrder(qOrder);

      INT_LIST mOrder;
      for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
        const Bond *mB = d_mol.getBondBetweenAtoms(c2[i], c2[j]);
        if (mB) {
          mOrder.push_back(mB->getIdx());
          if (mOrder.size() == mAt->getDegree()) break;
        }
      }
      while (mOrder.size() < mAt->getDegree()) {
        mOrder.push_back(-1);
      }
      INT_LIST moOrder;
      ROMol::OEDGE_ITER dbeg, dend;
      boost::tie(dbeg, dend) = d_mol.getAtomBonds(mAt);
      while (dbeg != dend) {
        int dbidx = d_mol[*dbeg]->getIdx();
        if (std::find(mOrder.begin(), mOrder.end(), dbidx) != mOrder.end())
          moOrder.push_back(dbidx);
        else
          moOrder.push_back(-1);
        ++dbeg;
      }
      int mPermCount =
          static_cast<int>(countSwapsToInterconvert(moOrder, mOrder));
      // std::cerr << "qorder: ";
      // std::copy(qOrder.begin(), qOrder.end(),
      //           std::ostream_iterator<int>(std::cerr, ", "));
      // std::cerr << std::endl;
      // std::cerr << "moOrder: ";
      // std::copy(moOrder.begin(), moOrder.end(),
      //           std::ostream_iterator<int>(std::cerr, ", "));
      // std::cerr << std::endl;
      // std::cerr << "morder: ";
      // std::copy(mOrder.begin(), mOrder.end(),
      //           std::ostream_iterator<int>(std::cerr, ", "));
      // std::cerr << std::endl;
      // std::cerr << "qPerm: " << qPermCount << " mPerm: " << mPermCount
      //           << " qtag: " << qAt->getChiralTag()
      //           << " mtag: " << mAt->getChiralTag() << std::endl;
      if ((qPermCount % 2 == mPermCount % 2 &&
           qAt->getChiralTag() != mAt->getChiralTag()) ||
          (qPermCount % 2 != mPermCount % 2 &&
           qAt->getChiralTag() == mAt->getChiralTag()))
        return false;
    }

    // now check double bonds
    for (unsigned int i = 0; i < d_query.getNumBonds(); ++i) {
      const Bond *qBnd = d_query.getBondWithIdx(i);
      if (qBnd->getBondType() != Bond::DOUBLE ||
          qBnd->getStereo() <= Bond::STEREOANY)
        continue;

      // don't think this can actually happen, but check to be sure:
      if (qBnd->getStereoAtoms().size() != 2) continue;

      std::map<unsigned int, unsigned int> qMap;
      for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
        qMap[c1[j]] = j;
      }
      const Bond *mBnd = d_mol.getBondBetweenAtoms(
          c2[qMap[qBnd->getBeginAtomIdx()]], c2[qMap[qBnd->getEndAtomIdx()]]);
      CHECK_INVARIANT(mBnd, "Matching bond not found");
      if (mBnd->getBondType() != Bond::DOUBLE ||
          qBnd->getStereo() <= Bond::STEREOANY)
        continue;
      // don't think this can actually happen, but check to be sure:
      if (mBnd->getStereoAtoms().size() != 2) continue;

      unsigned int end1Matches = 0;
      unsigned int end2Matches = 0;
      if (c2[qMap[qBnd->getBeginAtomIdx()]] == mBnd->getBeginAtomIdx()) {
        // query Begin == mol Begin
        if (c2[qMap[qBnd->getStereoAtoms()[0]]] == mBnd->getStereoAtoms()[0])
          end1Matches = 1;
        if (c2[qMap[qBnd->getStereoAtoms()[1]]] == mBnd->getStereoAtoms()[1])
          end2Matches = 1;
      } else {
        // query End == mol Begin
        if (c2[qMap[qBnd->getStereoAtoms()[0]]] == mBnd->getStereoAtoms()[1])
          end1Matches = 1;
        if (c2[qMap[qBnd->getStereoAtoms()[1]]] == mBnd->getStereoAtoms()[0])
          end2Matches = 1;
      }

      const unsigned totalMatches = end1Matches + end2Matches;
      const auto mStereo =
          Chirality::translateEZLabelToCisTrans(mBnd->getStereo());
      const auto qStereo =
          Chirality::translateEZLabelToCisTrans(qBnd->getStereo());

      if (mStereo == qStereo && totalMatches == 1) return false;
      if (mStereo != qStereo && totalMatches != 1) return false;
    }

    return true;
  }

 private:
  const ROMol &d_query;
  const ROMol &d_mol;
  const SubstructMatchParameters &d_params;
};

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
            mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW)
          return false;
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
            mBnd->getStereo() <= Bond::STEREOANY)
          return false;
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
void mergeMatchVect(std::vector<MatchVectType> &matches,
                    const std::vector<MatchVectType> &matchesTmp,
                    const ResSubstructMatchHelperArgs_ &args) {
  for (auto it = matchesTmp.begin();
       (matches.size() < args.params.maxMatches) && (it != matchesTmp.end());
       ++it) {
    if ((std::find(matches.begin(), matches.end(), *it) == matches.end()) &&
        (!args.params.uniquify || isToBeAddedToVector(matches, *it)))
      matches.push_back(*it);
  }
};
void ResSubstructMatchHelper_(const ResSubstructMatchHelperArgs_ &args,
                              std::vector<MatchVectType> *matches,
                              unsigned int bi, unsigned int ei) {
  for (unsigned int i = bi;
       (matches->size() < args.params.maxMatches) && (i < ei); ++i) {
    ROMol *mol = args.resMolSupplier[i];
    std::vector<MatchVectType> matchesTmp =
        SubstructMatch(*mol, args.query, args.params);
    mergeMatchVect(*matches, matchesTmp, args);
    delete mol;
  }
};
}  // namespace detail

// ----------------------------------------------
//
// find all matches
std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const ROMol &query,
    const SubstructMatchParameters &params) {
  std::vector<RecursiveStructureQuery *> locked;
#ifdef RDK_THREADSAFE_SSS
  locked.reserve(query.getNumAtoms());
#endif
  if (params.recursionPossible) {
    detail::SUBQUERY_MAP subqueryMap;
    ROMol::ConstAtomIterator atIt;
    for (atIt = query.beginAtoms(); atIt != query.endAtoms(); atIt++) {
      if ((*atIt)->getQuery()) {
        // std::cerr<<"recurse from atom "<<(*atIt)->getIdx()<<std::endl;
        detail::MatchSubqueries(mol, (*atIt)->getQuery(), params, subqueryMap,
                                locked);
      }
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, params);
  detail::BondLabelFunctor bondLabeler(query, mol, params);
  detail::MolMatchFinalCheckFunctor matchChecker(query, mol, params);

  std::list<detail::ssPairType> pms;
#if 0
    bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
                                  atomLabeler,bondLabeler,pms);
#else
  bool found =
      boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
                     bondLabeler, matchChecker, pms, params.maxMatches);
#endif
  std::vector<MatchVectType> matches;
  if (found) {
    unsigned int nQueryAtoms = query.getNumAtoms();
    matches.reserve(pms.size());
    for (std::list<detail::ssPairType>::const_iterator iter1 = pms.begin();
         iter1 != pms.end(); ++iter1) {
      MatchVectType matchVect;
      matchVect.resize(nQueryAtoms);
      for (const auto &iter2 : *iter1) {
        matchVect[iter2.first] = std::pair<int, int>(iter2.first, iter2.second);
      }
      matches.push_back(matchVect);
    }
    if (params.uniquify) {
      removeDuplicates(matches, mol.getNumAtoms());
    }
  }

#ifdef RDK_THREADSAFE_SSS
  if (params.recursionPossible) {
    BOOST_FOREACH (RecursiveStructureQuery *v, locked)
      v->d_mutex.unlock();
  }
#endif
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
  std::vector<MatchVectType> matches;
  detail::ResSubstructMatchHelperArgs_ args = {resMolSupplier, query, params};
  unsigned int nt =
      std::min(resMolSupplier.length(), getNumThreadsToUse(params.numThreads));
  if (nt == 1)
    detail::ResSubstructMatchHelper_(args, &matches, 0,
                                     resMolSupplier.length());
#ifdef RDK_THREADSAFE_SSS
  else {
    std::vector<std::future<void>> tg;
    std::vector<std::vector<MatchVectType> *> matchesThread(nt);
    unsigned int ei = 0;
    double dpt =
        static_cast<double>(resMolSupplier.length()) / static_cast<double>(nt);
    double dc = 0.0;
    for (unsigned int ti = 0; ti < nt; ++ti) {
      matchesThread[ti] = new std::vector<MatchVectType>();
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

    unsigned int matchSize = 0;
    for (unsigned int ti = 0; ti < nt; ++ti)
      matchSize += matchesThread[ti]->size();
    matches.reserve(matchSize);
    for (unsigned int ti = 0; ti < nt; ++ti) {
      mergeMatchVect(matches, *(matchesThread[ti]), args);
      delete matchesThread[ti];
    }
  }
#endif
  std::sort(matches.begin(), matches.end(), detail::matchVectCompare);
  return matches;
}

namespace detail {
unsigned int RecursiveMatcher(const ROMol &mol, const ROMol &query,
                              std::vector<int> &matches,
                              SUBQUERY_MAP &subqueryMap,
                              const SubstructMatchParameters &params,
                              std::vector<RecursiveStructureQuery *> &locked) {
  ROMol::ConstAtomIterator atIt;
  for (atIt = query.beginAtoms(); atIt != query.endAtoms(); atIt++) {
    if ((*atIt)->getQuery()) {
      MatchSubqueries(mol, (*atIt)->getQuery(), params, subqueryMap, locked);
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, params);
  detail::BondLabelFunctor bondLabeler(query, mol, params);
  detail::MolMatchFinalCheckFunctor matchChecker(query, mol, params);

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
    for (std::list<detail::ssPairType>::const_iterator iter1 = pms.begin();
         iter1 != pms.end(); ++iter1) {
      if (!query.hasProp(common_properties::_queryRootAtom)) {
        matches.push_back(iter1->begin()->second);
      } else {
        int rootIdx;
        query.getProp(common_properties::_queryRootAtom, rootIdx);
        bool found = false;
        for (const auto &pairIter : *iter1) {
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
  // std::cout << "*-*-* MS: " << (int)query << std::endl;
  // std::cout << "\t\t" << typeid(*query).name() << std::endl;
  if (query->getDescription() == "RecursiveStructure") {
    RecursiveStructureQuery *rsq = (RecursiveStructureQuery *)query;
#ifdef RDK_THREADSAFE_SSS
    rsq->d_mutex.lock();
    locked.push_back(rsq);
#endif
    rsq->clear();
    bool matchDone = false;
    if (rsq->getSerialNumber() &&
        subqueryMap.find(rsq->getSerialNumber()) != subqueryMap.end()) {
      // we've matched an equivalent serial number before, just
      // copy in the matches:
      matchDone = true;
      const RecursiveStructureQuery *orsq =
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
        // std::cerr<<" storing results for query serial number:
        // "<<rsq->getSerialNumber()<<std::endl;
      }
    }
  } else {
    // std::cout << "\tmsq1: ";
  }

  // now recurse over our children (these things can be nested)
  Queries::Query<int, Atom const *, true>::CHILD_VECT_CI childIt;
  // std::cout << query << " " << query->endChildren()-query->beginChildren() <<
  // std::endl;
  for (childIt = query->beginChildren(); childIt != query->endChildren();
       childIt++) {
    MatchSubqueries(mol, childIt->get(), params, subqueryMap, locked);
  }
  // std::cout << "<<- back " << (int)query << std::endl;
}

bool matchCompare(const std::pair<int, int> &a, const std::pair<int, int> &b) {
  return (a.second < b.second);
}

bool matchVectCompare(const MatchVectType &a, const MatchVectType &b) {
  for (unsigned int i = 0; i < std::min(a.size(), b.size()); ++i) {
    if (a[i].second != b[i].second) return (a[i].second < b[i].second);
  }
  return (a.size() < b.size());
}

bool isToBeAddedToVector(std::vector<MatchVectType> &matches,
                         const MatchVectType &m) {
  bool isToBeAdded = true;
  MatchVectType mCopy = m;
  std::sort(mCopy.begin(), mCopy.end(), matchCompare);
  for (auto it = matches.end(); isToBeAdded && it != matches.begin();) {
    --it;
    isToBeAdded = (it->size() != mCopy.size());
    if (!isToBeAdded) {
      MatchVectType matchCopy = *it;
      std::sort(matchCopy.begin(), matchCopy.end(), matchCompare);
      for (unsigned int i = 0; !isToBeAdded && (i < matchCopy.size()); ++i)
        isToBeAdded = (mCopy[i].second != matchCopy[i].second);
      if (!isToBeAdded) {
        for (unsigned int i = 0; !isToBeAdded && (i < m.size()); ++i)
          isToBeAdded = (m[i].second < (*it)[i].second);
        if (isToBeAdded) {
          matches.erase(it);
          break;
        }
      }
    }
  }
  return isToBeAdded;
}
}  // end of namespace detail
}  // namespace RDKit
