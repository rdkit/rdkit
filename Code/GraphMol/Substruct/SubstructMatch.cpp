//
//  Copyright (C) 2001-2015 Greg Landrum and Rational Discovery LLC
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
#include "SubstructMatch.h"
#include "SubstructUtils.h"
#include <boost/smart_ptr.hpp>
#include <map>

#ifdef RDK_THREADSAFE_SSS
#include <boost/thread/mutex.hpp>
#endif

#include "ullmann.hpp"
#include "vf2.hpp"

namespace RDKit {
namespace detail {
typedef std::map<unsigned int, QueryAtom::QUERYATOM_QUERY *> SUBQUERY_MAP;

typedef struct {
  ResonanceMolSupplier &resMolSupplier;
  const ROMol &query;
  bool uniquify;
  bool recursionPossible;
  bool useChirality;
  bool useQueryQueryMatches;
  unsigned int maxMatches;
} ResSubstructMatchHelperArgs_;

void MatchSubqueries(const ROMol &mol, QueryAtom::QUERYATOM_QUERY *q,
                     bool useChirality, SUBQUERY_MAP &subqueryMap,
                     bool useQueryQueryMatches,
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

typedef std::list<std::pair<MolGraph::vertex_descriptor,
                            MolGraph::vertex_descriptor> > ssPairType;

class MolMatchFinalCheckFunctor {
 public:
  MolMatchFinalCheckFunctor(const ROMol &query, const ROMol &mol,
                            bool useChirality)
      : d_query(query), d_mol(mol), df_useChirality(useChirality){};
  bool operator()(const boost::detail::node_id c1[],
                  const boost::detail::node_id c2[]) const {
    // std::cerr << "  check! " << df_useChirality << std::endl;
    if (!df_useChirality) return true;
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
          (qBnd->getStereo() != Bond::STEREOZ &&
           qBnd->getStereo() != Bond::STEREOE))
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
          (mBnd->getStereo() != Bond::STEREOZ &&
           mBnd->getStereo() != Bond::STEREOE))
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
      // std::cerr<<"  bnd: "<<qBnd->getIdx()<<":"<<qBnd->getStereo()<<" -
      // "<<mBnd->getIdx()<<":"<<mBnd->getStereo()<<"  --  "<<end1Matches<<"
      // "<<end2Matches<<std::endl;
      if (mBnd->getStereo() == qBnd->getStereo() &&
          (end1Matches + end2Matches) == 1)
        return false;
      if (mBnd->getStereo() != qBnd->getStereo() &&
          (end1Matches + end2Matches) != 1)
        return false;
    }

    return true;
  }

 private:
  const ROMol &d_query;
  const ROMol &d_mol;
  bool df_useChirality;
};

class AtomLabelFunctor {
 public:
  AtomLabelFunctor(const ROMol &query, const ROMol &mol, bool useChirality,
                   bool useQueryQueryMatches)
      : d_query(query),
        d_mol(mol),
        df_useChirality(useChirality),
        df_useQueryQueryMatches(useQueryQueryMatches){};
  bool operator()(unsigned int i, unsigned int j) const {
    bool res = false;
    if (df_useChirality) {
      const Atom *qAt = d_query.getAtomWithIdx(i);
      if (qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
          qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
        const Atom *mAt = d_mol.getAtomWithIdx(j);
        if (mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
            mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW)
          return false;
      }
    }
    res = atomCompat(d_query[i], d_mol[j], df_useQueryQueryMatches);
    return res;
  }

 private:
  const ROMol &d_query;
  const ROMol &d_mol;
  bool df_useChirality;
  bool df_useQueryQueryMatches;
};
class BondLabelFunctor {
 public:
  BondLabelFunctor(const ROMol &query, const ROMol &mol, bool useChirality,
                   bool useQueryQueryMatches)
      : d_query(query), d_mol(mol), df_useChirality(useChirality) {
    RDUNUSED_PARAM(useQueryQueryMatches);
  };
  bool operator()(MolGraph::edge_descriptor i,
                  MolGraph::edge_descriptor j) const {
    if (df_useChirality) {
      const BOND_SPTR qBnd = d_query[i];
      if (qBnd->getBondType() == Bond::DOUBLE &&
          (qBnd->getStereo() == Bond::STEREOZ ||
           qBnd->getStereo() == Bond::STEREOE)) {
        const BOND_SPTR mBnd = d_mol[j];
        if (mBnd->getBondType() == Bond::DOUBLE &&
            !(mBnd->getStereo() == Bond::STEREOZ ||
              mBnd->getStereo() == Bond::STEREOE))
          return false;
      }
    }
    bool res = bondCompat(d_query[i], d_mol[j]);
    return res;
  }

 private:
  const ROMol &d_query;
  const ROMol &d_mol;
  bool df_useChirality;
  // bool df_useQueryQueryMatches;
};
void mergeMatchVect(std::vector<MatchVectType> &matches,
                    const std::vector<MatchVectType> &matchesTmp,
                    const ResSubstructMatchHelperArgs_ &args) {
  for (std::vector<MatchVectType>::const_iterator it = matchesTmp.begin();
       (matches.size() < args.maxMatches) && (it != matchesTmp.end()); ++it) {
    if ((std::find(matches.begin(), matches.end(), *it) == matches.end()) &&
        (!args.uniquify || isToBeAddedToVector(matches, *it)))
      matches.push_back(*it);
  }
};
void ResSubstructMatchHelper_(const ResSubstructMatchHelperArgs_ &args,
                              std::vector<MatchVectType> *matches,
                              unsigned int bi, unsigned int ei) {
  for (unsigned int i = bi; (matches->size() < args.maxMatches) && (i < ei);
       ++i) {
    ROMol *mol = args.resMolSupplier[i];
    std::vector<MatchVectType> matchesTmp;
    SubstructMatch(*mol, args.query, matchesTmp, args.uniquify,
                   args.recursionPossible, args.useChirality,
                   args.useQueryQueryMatches, args.maxMatches);
    mergeMatchVect(*matches, matchesTmp, args);
    delete mol;
  }
};
}

// ----------------------------------------------
//
// find one match
//
bool SubstructMatch(const ROMol &mol, const ROMol &query,
                    MatchVectType &matchVect, bool recursionPossible,
                    bool useChirality, bool useQueryQueryMatches) {
  // std::cerr<<"begin match"<<std::endl;
  std::vector<RecursiveStructureQuery *> locked;
#ifdef RDK_THREADSAFE_SSS
  locked.reserve(query.getNumAtoms());
#endif
  if (recursionPossible) {
    ROMol::ConstAtomIterator atIt;
    detail::SUBQUERY_MAP subqueryMap;
    for (atIt = query.beginAtoms(); atIt != query.endAtoms(); atIt++) {
      if ((*atIt)->getQuery()) {
        detail::MatchSubqueries(mol, (*atIt)->getQuery(), useChirality,
                                subqueryMap, useQueryQueryMatches, locked);
      }
    }
  }
  // std::cerr<<"main matching"<<std::endl;

  matchVect.clear();
  matchVect.resize(0);

  detail::MolMatchFinalCheckFunctor matchChecker(query, mol, useChirality);
  detail::AtomLabelFunctor atomLabeler(query, mol, useChirality,
                                       useQueryQueryMatches);
  detail::BondLabelFunctor bondLabeler(query, mol, useChirality,
                                       useQueryQueryMatches);

  detail::ssPairType match;
#if 0
    bool res=boost::ullmann(query.getTopology(),mol.getTopology(),
                            atomLabeler,bondLabeler,match);
#else
  bool res = boost::vf2(query.getTopology(), mol.getTopology(), atomLabeler,
                        bondLabeler, matchChecker, match);
#endif
  if (res) {
    matchVect.resize(query.getNumAtoms());
    for (detail::ssPairType::const_iterator iter = match.begin();
         iter != match.end(); ++iter) {
      matchVect[iter->first] = std::pair<int, int>(iter->first, iter->second);
    }
  }

#ifdef RDK_THREADSAFE_SSS
  if (recursionPossible) {
    BOOST_FOREACH (RecursiveStructureQuery *v, locked)
      v->d_mutex.unlock();
  }
#endif
  return res;
}

// ----------------------------------------------
//
// find one match in ResonanceMolSupplier object
//
bool SubstructMatch(ResonanceMolSupplier &resMolSupplier, const ROMol &query,
                    MatchVectType &matchVect, bool recursionPossible,
                    bool useChirality, bool useQueryQueryMatches) {
  bool match = false;
  for (unsigned int i = 0; !match && (i < resMolSupplier.length()); ++i) {
    ROMol *mol = resMolSupplier[i];
    match = SubstructMatch(*mol, query, matchVect, recursionPossible,
                           useChirality, useQueryQueryMatches);
    delete mol;
  }
  return match;
}

// ----------------------------------------------
//
// find all matches
//
//  NOTE: this blows out the contents of matches
//
unsigned int SubstructMatch(const ROMol &mol, const ROMol &query,
                            std::vector<MatchVectType> &matches, bool uniquify,
                            bool recursionPossible, bool useChirality,
                            bool useQueryQueryMatches,
                            unsigned int maxMatches) {
  std::vector<RecursiveStructureQuery *> locked;
#ifdef RDK_THREADSAFE_SSS
  locked.reserve(query.getNumAtoms());
#endif
  if (recursionPossible) {
    detail::SUBQUERY_MAP subqueryMap;
    ROMol::ConstAtomIterator atIt;
    for (atIt = query.beginAtoms(); atIt != query.endAtoms(); atIt++) {
      if ((*atIt)->getQuery()) {
        // std::cerr<<"recurse from atom "<<(*atIt)->getIdx()<<std::endl;
        detail::MatchSubqueries(mol, (*atIt)->getQuery(), useChirality,
                                subqueryMap, useQueryQueryMatches, locked);
      }
    }
  }

  matches.clear();
  matches.resize(0);

  detail::AtomLabelFunctor atomLabeler(query, mol, useChirality,
                                       useQueryQueryMatches);
  detail::BondLabelFunctor bondLabeler(query, mol, useChirality,
                                       useQueryQueryMatches);
  detail::MolMatchFinalCheckFunctor matchChecker(query, mol, useChirality);

  std::list<detail::ssPairType> pms;
#if 0
    bool found=boost::ullmann_all(query.getTopology(),mol.getTopology(),
                                  atomLabeler,bondLabeler,pms);
#else
  bool found =
      boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
                     bondLabeler, matchChecker, pms, maxMatches);
#endif
  unsigned int res = 0;
  if (found) {
    unsigned int nQueryAtoms = query.getNumAtoms();
    matches.reserve(pms.size());
    for (std::list<detail::ssPairType>::const_iterator iter1 = pms.begin();
         iter1 != pms.end(); ++iter1) {
      MatchVectType matchVect;
      matchVect.resize(nQueryAtoms);
      for (detail::ssPairType::const_iterator iter2 = iter1->begin();
           iter2 != iter1->end(); ++iter2) {
        matchVect[iter2->first] =
            std::pair<int, int>(iter2->first, iter2->second);
      }
      matches.push_back(matchVect);
    }
    if (uniquify) {
      removeDuplicates(matches, mol.getNumAtoms());
    }
    res = matches.size();
  }

#ifdef RDK_THREADSAFE_SSS
  if (recursionPossible) {
    BOOST_FOREACH (RecursiveStructureQuery *v, locked)
      v->d_mutex.unlock();
  }
#endif
  return res;
}

// ----------------------------------------------
//
// find all matches in a ResonanceMolSupplier object
//
//  NOTE: this blows out the contents of matches
//
unsigned int SubstructMatch(ResonanceMolSupplier &resMolSupplier,
                            const ROMol &query,
                            std::vector<MatchVectType> &matches, bool uniquify,
                            bool recursionPossible, bool useChirality,
                            bool useQueryQueryMatches, unsigned int maxMatches,
                            int numThreads) {
  matches.clear();
  detail::ResSubstructMatchHelperArgs_ args = {
      resMolSupplier, query, uniquify, recursionPossible, useChirality,
      useQueryQueryMatches, maxMatches};
  unsigned int nt =
      std::min(resMolSupplier.length(), getNumThreadsToUse(numThreads));
  if (nt == 1)
    detail::ResSubstructMatchHelper_(args, &matches, 0,
                                     resMolSupplier.length());
#ifdef RDK_THREADSAFE_SSS
  else {
    boost::thread_group tg;
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
      tg.add_thread(new boost::thread(detail::ResSubstructMatchHelper_, args,
                                      matchesThread[ti], bi, ei));
    }
    tg.join_all();
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
  return matches.size();
}

namespace detail {
unsigned int RecursiveMatcher(const ROMol &mol, const ROMol &query,
                              std::vector<int> &matches, bool useChirality,
                              SUBQUERY_MAP &subqueryMap,
                              bool useQueryQueryMatches,
                              std::vector<RecursiveStructureQuery *> &locked) {
  ROMol::ConstAtomIterator atIt;
  for (atIt = query.beginAtoms(); atIt != query.endAtoms(); atIt++) {
    if ((*atIt)->getQuery()) {
      MatchSubqueries(mol, (*atIt)->getQuery(), useChirality, subqueryMap,
                      useQueryQueryMatches, locked);
    }
  }

  detail::AtomLabelFunctor atomLabeler(query, mol, useChirality,
                                       useQueryQueryMatches);
  detail::BondLabelFunctor bondLabeler(query, mol, useChirality,
                                       useQueryQueryMatches);
  detail::MolMatchFinalCheckFunctor matchChecker(query, mol, useChirality);

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
        for (detail::ssPairType::const_iterator pairIter = iter1->begin();
             pairIter != iter1->end(); ++pairIter) {
          if (pairIter->first == static_cast<unsigned int>(rootIdx)) {
            matches.push_back(pairIter->second);
            found = true;
            break;
          }
        }
        if (!found) {
          BOOST_LOG(rdErrorLog) << "no match found for queryRootAtom"
                                << std::endl;
        }
      }
    }
    res = matches.size();
  }
  // std::cout << " <<< RecursiveMatcher: " << int(query) << std::endl;
  return res;
}

void MatchSubqueries(const ROMol &mol, QueryAtom::QUERYATOM_QUERY *query,
                     bool useChirality, SUBQUERY_MAP &subqueryMap,
                     bool useQueryQueryMatches,
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
      for (RecursiveStructureQuery::CONTAINER_TYPE::const_iterator setIter =
               orsq->beginSet();
           setIter != orsq->endSet(); ++setIter) {
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
        unsigned int res =
            RecursiveMatcher(mol, *queryMol, matchStarts, useChirality,
                             subqueryMap, useQueryQueryMatches, locked);
        if (res) {
          for (std::vector<int>::iterator i = matchStarts.begin();
               i != matchStarts.end(); i++) {
            rsq->insert(*i);
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
    MatchSubqueries(mol, childIt->get(), useChirality, subqueryMap,
                    useQueryQueryMatches, locked);
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
  for (std::vector<MatchVectType>::iterator it = matches.end();
       isToBeAdded && it != matches.begin();) {
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
}
