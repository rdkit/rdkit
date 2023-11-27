//  Copyright (c) 2017-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
#include "SubstructLibrary.h"
#include <RDGeneral/RDThreads.h>
#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/GeneralizedSubstruct/XQMol.h>
#include <boost/dynamic_bitset.hpp>

namespace RDKit {

using namespace GeneralizedSubstruct;

bool SubstructLibraryCanSerialize() {
#ifdef RDK_USE_BOOST_SERIALIZATION
  return true;
#else
  return false;
#endif
}

struct Bits {
  const ExplicitBitVect *queryBits;
  const FPHolderBase *fps;
  SubstructMatchParameters params;

  Bits(const FPHolderBase *fingerprints, const ROMol &m,
       const SubstructMatchParameters &ssparams)
      : fps(fingerprints), params(ssparams) {
    if (fps) {
      queryBits = fps->makeFingerprint(m);
    } else {
      queryBits = nullptr;
    }
  }

  Bits(const FPHolderBase *fingerprints, const TautomerQuery &m,
       const SubstructMatchParameters &ssparams)
      : fps(nullptr), params(ssparams) {
    if (fingerprints) {
      const TautomerPatternHolder *tp =
          dynamic_cast<const TautomerPatternHolder *>(fingerprints);
      if (!tp) {
        BOOST_LOG(rdWarningLog) << "Pattern fingerprints for tautomersearch "
                                   "aren't tautomer fingerprints, ignoring..."
                                << std::endl;
        queryBits = nullptr;
        fps = nullptr;
      } else {
        fps = fingerprints;
        queryBits = m.patternFingerprintTemplate(tp->getNumBits());
      }
    } else {
      queryBits = nullptr;
    }
  }

  // FIX complete this
  Bits(const FPHolderBase *fingerprints, const ExtendedQueryMol &xqm,
       const SubstructMatchParameters &ssparams)
      : fps(fingerprints), params(ssparams) {
    if (fps) {
      const auto *tph = dynamic_cast<const TautomerPatternHolder *>(fps);
      const auto *ph = dynamic_cast<const PatternHolder *>(fps);
      if (std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm.xqmol)) {
        queryBits = fps->makeFingerprint(
            *std::get<ExtendedQueryMol::RWMol_T>(xqm.xqmol));
      } else if (std::holds_alternative<ExtendedQueryMol::MolBundle_T>(
                     xqm.xqmol)) {
        auto &bndl = std::get<ExtendedQueryMol::MolBundle_T>(xqm.xqmol);
        auto tqb = new ExplicitBitVect(ph->getNumBits());
        queryBits = tqb;
        for (auto mol : bndl->getMols()) {
          auto tfp = fps->makeFingerprint(*mol);
          *tqb &= *tfp;
          delete tfp;
        }
      } else if (std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(
                     xqm.xqmol)) {
        auto &tq = std::get<ExtendedQueryMol::TautomerQuery_T>(xqm.xqmol);
        if (!tph) {
          BOOST_LOG(rdWarningLog) << "Pattern fingerprints for tautomersearch "
                                     "aren't tautomer fingerprints, ignoring..."
                                  << std::endl;
          queryBits = nullptr;
          fps = nullptr;
        } else {
          queryBits = tq->patternFingerprintTemplate(tph->getNumBits());
        }
      } else if (std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(
                     xqm.xqmol)) {
        if (!tph) {
          BOOST_LOG(rdWarningLog) << "Pattern fingerprints for tautomersearch "
                                     "aren't tautomer fingerprints, ignoring..."
                                  << std::endl;
          queryBits = nullptr;
          fps = nullptr;
        } else {
          auto &bndl = std::get<ExtendedQueryMol::TautomerBundle_T>(xqm.xqmol);
          auto tqb = new ExplicitBitVect(ph->getNumBits());
          queryBits = tqb;
          for (auto &tq : *bndl) {
            auto tfp = tq->patternFingerprintTemplate(tph->getNumBits());
            *tqb &= *tfp;
            delete tfp;
          }
        }
      } else {
        queryBits = nullptr;
      }
    } else {
      queryBits = nullptr;
    }
  }

  bool check(unsigned int idx) const {
    if (fps) {
      return fps->passesFilter(idx, *queryBits);
    }
    return true;
  }
};

unsigned int SubstructLibrary::addMol(const ROMol &m) {
  unsigned int idx = mols->addMol(m);
  if (fps) {
    unsigned int fpidx = fps->addMol(m);
    CHECK_INVARIANT(idx == fpidx,
                    "#mols different than #fingerprints in SubstructLibrary");
  }
  if (keyholder.get() != nullptr) {
    unsigned int keyidx = keyholder->addMol(m);
    CHECK_INVARIANT(idx == keyidx,
                    "#mols different than #keys in SubstructLibrary");
  }

  return idx;
}

namespace {
//  Return true if the pattern contains a ring query
bool query_needs_rings(const ROMol &in_query) {
  for (const auto &atom : in_query.atoms()) {
    if (atom->hasQuery()) {
      if (describeQuery(atom).find("Ring") != std::string::npos) {
        return true;
      } else if (atom->getQuery()->getDescription() == "RecursiveStructure") {
        auto *rsq = (RecursiveStructureQuery *)atom->getQuery();
        if (query_needs_rings(*rsq->getQueryMol())) {
          return true;
        }
      }
    }
  }
  for (const auto &bond : in_query.bonds()) {
    if (bond->hasQuery()) {
      if (describeQuery(bond).find("Ring") != std::string::npos) {
        return true;
      }
    }
  }
  return false;
}

bool query_needs_rings(const TautomerQuery &in_query) {
  return query_needs_rings(in_query.getTemplateMolecule());
}

bool query_needs_rings(const ExtendedQueryMol &xqm) {
  if (std::holds_alternative<ExtendedQueryMol::RWMol_T>(xqm.xqmol)) {
    return query_needs_rings(*std::get<ExtendedQueryMol::RWMol_T>(xqm.xqmol));
  } else if (std::holds_alternative<ExtendedQueryMol::TautomerQuery_T>(
                 xqm.xqmol)) {
    return query_needs_rings(
        std::get<ExtendedQueryMol::TautomerQuery_T>(xqm.xqmol)
            ->getTemplateMolecule());
  } else if (std::holds_alternative<ExtendedQueryMol::MolBundle_T>(xqm.xqmol)) {
    for (const auto &mol :
         std::get<ExtendedQueryMol::MolBundle_T>(xqm.xqmol)->getMols()) {
      if (query_needs_rings(*mol)) {
        return true;
      }
    }
    return false;
  } else if (std::holds_alternative<ExtendedQueryMol::TautomerBundle_T>(
                 xqm.xqmol)) {
    for (const auto &tq :
         *std::get<ExtendedQueryMol::TautomerBundle_T>(xqm.xqmol)) {
      if (query_needs_rings(tq->getTemplateMolecule())) {
        return true;
      }
    }
    return false;
  }
  return true;  // if we somehow get here, we better assume that rings are
                // necessary
}

template <class Query>
void SubSearcher(const Query &in_query, const Bits &bits,
                 const MolHolderBase &mols, unsigned int start,
                 unsigned int &end, unsigned int numThreads,
                 const bool needs_rings, int &counter, const int maxResults,
                 boost::dynamic_bitset<> &found,
                 const std::vector<unsigned int> &searchOrder,
                 std::vector<unsigned int> *idxs) {
  PRECONDITION(searchOrder.empty() || searchOrder.size() >= end,
               "bad searchOrder data");
  // we copy the query so that we don't end up with lock contention for
  // recursive matchers when using multiple threads
  Query query(in_query);
  for (unsigned int idx = start; idx < end; idx += numThreads) {
    unsigned int sidx = idx;
    if (!searchOrder.empty()) {
      sidx = searchOrder[idx];
    }
    if (!bits.check(sidx) || found[sidx]) {
      continue;
    }
    // need shared_ptr as it (may) control the lifespan of the
    //  returned molecule!
    const boost::shared_ptr<ROMol> &m = mols.getMol(sidx);
    ROMol *mol = m.get();
    if (!mol) {
      continue;
    }
    if (needs_rings &&
        (!mol->getRingInfo() || !mol->getRingInfo()->isSymmSssr())) {
      MolOps::symmetrizeSSSR(*mol);
    }

    if (!SubstructMatch(*mol, query, bits.params).empty()) {
      ++counter;
      found.set(sidx);
      if (idxs) {
        idxs->push_back(sidx);
        if (maxResults > 0 && counter == maxResults) {
          // if we reached maxResults, record the last idx we processed and
          // bail out
          end = idx;
          break;
        }
      }
    }
  }
}

template <class Query>
int internalGetMatches(const Query &query, MolHolderBase &mols,
                       const FPHolderBase *fps, unsigned int startIdx,
                       unsigned int endIdx,
                       const SubstructMatchParameters &params, int numThreads,
                       int maxResults, boost::dynamic_bitset<> &found,
                       const std::vector<unsigned int> &searchOrder,
                       std::vector<unsigned int> *idxs) {
  PRECONDITION(startIdx < mols.size(), "startIdx out of bounds");
  PRECONDITION(searchOrder.empty() || startIdx < searchOrder.size(),
               "startIdx out of bounds");
  PRECONDITION(endIdx > startIdx, "endIdx > startIdx");

  // do not do any work if no results were requested
  if (maxResults == 0) {
    return 0;
  }

  endIdx = std::min(mols.size(), endIdx);
  if (!searchOrder.empty()) {
    endIdx = std::min(static_cast<unsigned int>(searchOrder.size()), endIdx);
  }

  numThreads = static_cast<int>(getNumThreadsToUse(numThreads));
  numThreads = std::min(numThreads, static_cast<int>(endIdx));

  bool needs_rings = query_needs_rings(query);
  Bits bits(fps, query, params);
  int counter = 0;

#ifdef RDK_BUILD_THREADSAFE_SSS
  if (numThreads > 1) {
    std::vector<int> counterVect(numThreads, 0);
    int maxResultsPerThread = maxResults;
    if (maxResults > 0) {
      maxResultsPerThread /= numThreads;
    }
    std::vector<int> maxResultsVect(numThreads, maxResultsPerThread);
    std::vector<unsigned int> endIdxVect(numThreads, endIdx);
    if (maxResults > 0) {
      int excess = maxResults % numThreads;
      for (int i = 0; i < excess; ++i) {
        ++maxResultsVect[i];
      }
    }
    std::vector<std::future<void>> thread_group;
    std::vector<std::vector<unsigned int>> internal_results;
    std::vector<boost::dynamic_bitset<>> internal_found(numThreads);
    if (idxs) {
      internal_results.resize(numThreads);
    }
    int thread_group_idx;
    for (thread_group_idx = 0; thread_group_idx < numThreads;
         ++thread_group_idx) {
      internal_found[thread_group_idx] = found;
      // need to use boost::ref otherwise things are passed by value
      thread_group.emplace_back(std::async(
          std::launch::async, SubSearcher<Query>, std::ref(query), bits,
          std::ref(mols), startIdx + thread_group_idx,
          std::ref(endIdxVect[thread_group_idx]), numThreads, needs_rings,
          std::ref(counterVect[thread_group_idx]),
          maxResultsVect[thread_group_idx],
          std::ref(internal_found[thread_group_idx]), std::ref(searchOrder),
          idxs ? &internal_results[thread_group_idx] : nullptr));
    }
    unsigned int maxEndIdx;
    if (maxResults > 0) {
      // If we are running with maxResults in a multi-threaded settings,
      // some threads may have screened more molecules than others.
      // If maxResults was close to the theoretical maximum, some threads
      // might have even run out of molecules to screen without reaching
      // maxResults so we need to make sure that all threads have screened as
      // many molecules as the most productive thread if we want
      // multi-threaded runs to yield the same results independently from the
      // number of threads.
      thread_group_idx = 0;
      for (auto &fut : thread_group) {
        fut.get();
        counter += counterVect[thread_group_idx++];
      }
      thread_group.clear();
      // Find out out the max number of molecules that was screened by the
      // most productive thread and do the same in all other threads, unless
      // the max number of molecules was reached
      maxEndIdx = *std::max_element(endIdxVect.begin(), endIdxVect.end());
      for (thread_group_idx = 0; thread_group_idx < numThreads;
           ++thread_group_idx) {
        if (endIdxVect[thread_group_idx] >= maxEndIdx) {
          continue;
        }
        internal_found[thread_group_idx] = found;
        // need to use boost::ref otherwise things are passed by value
        thread_group.emplace_back(std::async(
            std::launch::async, SubSearcher<Query>, std::ref(query), bits,
            std::ref(mols), endIdxVect[thread_group_idx] + numThreads,
            std::ref(maxEndIdx), numThreads, needs_rings,
            std::ref(counterVect[thread_group_idx]), -1,
            std::ref(internal_found[thread_group_idx]), std::ref(searchOrder),
            &internal_results[thread_group_idx]));
      }
    }
    for (auto &fut : thread_group) {
      fut.get();
    }
    for (thread_group_idx = 0; thread_group_idx < numThreads;
         ++thread_group_idx) {
      if (idxs) {
        idxs->insert(idxs->end(), internal_results[thread_group_idx].begin(),
                     internal_results[thread_group_idx].end());
      }
      found |= internal_found[thread_group_idx];
      // If there was no maxResults, we still need to count, otherwise
      // this has already been done previously
      if (maxResults < 0) {
        counter += counterVect[thread_group_idx];
      }
    }
  } else {
    // if this is running single-threaded, no need to suffer the overhead of
    // std::async
    SubSearcher(query, bits, mols, startIdx, endIdx, 1, needs_rings, counter,
                maxResults, found, searchOrder, idxs);
  }
  if (idxs && numThreads > 1) {
    // the sort is necessary to ensure consistency across runs with different
    // numbers of threads
    if (!searchOrder.empty()) {
      std::transform(idxs->begin(), idxs->end(), idxs->begin(),
                     [searchOrder](unsigned int v) -> unsigned int {
                       return std::find(searchOrder.begin(), searchOrder.end(),
                                        v) -
                              searchOrder.begin();
                     });
    }
    std::sort(idxs->begin(), idxs->end());
    // we may have actually accumulated more results than maxResults due
    // to the top up above, so trim the results down if that's the case
    if (maxResults > 0 &&
        idxs->size() > static_cast<unsigned int>(maxResults)) {
      idxs->resize(maxResults);
    }
    if (!searchOrder.empty()) {
      std::transform(idxs->begin(), idxs->end(), idxs->begin(),
                     [searchOrder](unsigned int v) -> unsigned int {
                       return searchOrder[v];
                     });
    }
  }
#else
  SubSearcher(query, bits, mols, startIdx, endIdx, 1, needs_rings, counter,
              maxResults, found, searchOrder, idxs);
#endif

  delete bits.queryBits;

  return counter;
}

int molbundleGetMatches(const MolBundle &query, MolHolderBase &mols,
                        const FPHolderBase *fps, unsigned int startIdx,
                        unsigned int endIdx,
                        const SubstructMatchParameters &params, int numThreads,
                        int maxResults,
                        const std::vector<unsigned int> &searchOrder,
                        std::vector<unsigned int> *idxs) {
  int res = 0;
  boost::dynamic_bitset<> found(mols.size());
  for (const auto &qmol : query.getMols()) {
    maxResults -= res;
    res += internalGetMatches(*qmol, mols, fps, startIdx, endIdx, params,
                              numThreads, maxResults, found, searchOrder, idxs);
  }
  return res;
}

}  // namespace

std::vector<unsigned int> SubstructLibrary::getMatches(
    const ROMol &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads,
    int maxResults) const {
  std::vector<unsigned int> idxs;
  boost::dynamic_bitset<> found(mols->size());
  internalGetMatches(query, *mols, fps, startIdx, endIdx, params, numThreads,
                     maxResults, found, searchOrder, &idxs);
  return idxs;
}

std::vector<unsigned int> SubstructLibrary::getMatches(
    const TautomerQuery &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads,
    int maxResults) const {
  std::vector<unsigned int> idxs;
  boost::dynamic_bitset<> found(mols->size());
  internalGetMatches(query, *mols, fps, startIdx, endIdx, params, numThreads,
                     maxResults, found, searchOrder, &idxs);
  return idxs;
}

std::vector<unsigned int> SubstructLibrary::getMatches(
    const MolBundle &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads,
    int maxResults) const {
  std::vector<unsigned int> idxs;
  molbundleGetMatches(query, *mols, fps, startIdx, endIdx, params, numThreads,
                      maxResults, searchOrder, &idxs);
  return idxs;
}

std::vector<unsigned int> SubstructLibrary::getMatches(
    const ExtendedQueryMol &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads,
    int maxResults) const {
  std::vector<unsigned int> idxs;
  boost::dynamic_bitset<> found(mols->size());
  internalGetMatches(query, *mols, fps, startIdx, endIdx, params, numThreads,
                     maxResults, found, searchOrder, &idxs);
  return idxs;
}

unsigned int SubstructLibrary::countMatches(
    const ROMol &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads) const {
  boost::dynamic_bitset<> found(mols->size());
  return internalGetMatches(query, *mols, fps, startIdx, endIdx, params,
                            numThreads, -1, found, searchOrder, nullptr);
}

unsigned int SubstructLibrary::countMatches(
    const TautomerQuery &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads) const {
  boost::dynamic_bitset<> found(mols->size());
  return internalGetMatches(query, *mols, fps, startIdx, endIdx, params,
                            numThreads, -1, found, searchOrder, nullptr);
}
unsigned int SubstructLibrary::countMatches(
    const MolBundle &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads) const {
  return molbundleGetMatches(query, *mols, fps, startIdx, endIdx, params,
                             numThreads, -1, searchOrder, nullptr);
}

unsigned int SubstructLibrary::countMatches(
    const ExtendedQueryMol &query, unsigned int startIdx, unsigned int endIdx,
    const SubstructMatchParameters &params, int numThreads) const {
  boost::dynamic_bitset<> found(mols->size());
  return internalGetMatches(query, *mols, fps, startIdx, endIdx, params,
                            numThreads, -1, found, searchOrder, nullptr);
}

bool SubstructLibrary::hasMatch(const ROMol &query, unsigned int startIdx,
                                unsigned int endIdx,
                                const SubstructMatchParameters &params,
                                int numThreads) const {
  const int maxResults = 1;
  return getMatches(query, startIdx, endIdx, params, numThreads, maxResults)
             .size() > 0;
}

bool SubstructLibrary::hasMatch(const TautomerQuery &query,
                                unsigned int startIdx, unsigned int endIdx,
                                const SubstructMatchParameters &params,
                                int numThreads) const {
  const int maxResults = 1;
  return getMatches(query, startIdx, endIdx, params, numThreads, maxResults)
             .size() > 0;
}
bool SubstructLibrary::hasMatch(const MolBundle &query, unsigned int startIdx,
                                unsigned int endIdx,
                                const SubstructMatchParameters &params,
                                int numThreads) const {
  const int maxResults = 1;
  return getMatches(query, startIdx, endIdx, params, numThreads, maxResults)
             .size() > 0;
}
bool SubstructLibrary::hasMatch(const ExtendedQueryMol &query,
                                unsigned int startIdx, unsigned int endIdx,
                                const SubstructMatchParameters &params,
                                int numThreads) const {
  const int maxResults = 1;
  return getMatches(query, startIdx, endIdx, params, numThreads, maxResults)
             .size() > 0;
}

void SubstructLibrary::toStream(std::ostream &ss) const {
#ifndef RDK_USE_BOOST_SERIALIZATION
  RDUNUSED_PARAM(ss);
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_oarchive ar(ss);
  ar << *this;
#endif
}

std::string SubstructLibrary::Serialize() const {
  std::stringstream ss;
  toStream(ss);
  return ss.str();
}

void SubstructLibrary::initFromStream(std::istream &ss) {
#ifndef RDK_USE_BOOST_SERIALIZATION
  RDUNUSED_PARAM(ss);
  PRECONDITION(0, "Boost SERIALIZATION is not enabled")
#else
  boost::archive::text_iarchive ar(ss);
  ar >> *this;
#endif
}

void SubstructLibrary::initFromString(const std::string &text) {
  std::stringstream ss(text);
  initFromStream(ss);
}

}  // namespace RDKit
