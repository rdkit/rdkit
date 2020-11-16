//  Copyright (c) 2017-2019, Novartis Institutes for BioMedical Research Inc.
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
#ifdef RDK_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

#include <GraphMol/Substruct/SubstructMatch.h>


namespace RDKit {

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
  bool recursionPossible;
  bool useChirality;
  bool useQueryQueryMatches;

  Bits(const FPHolderBase *fps, const ROMol &m, bool recursionPossible,
       bool useChirality, bool useQueryQueryMatches)
      : fps(fps),
        recursionPossible(recursionPossible),
        useChirality(useChirality),
        useQueryQueryMatches(useQueryQueryMatches) {
    if (fps) {
      queryBits = fps->makeFingerprint(m);
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
  unsigned int size = mols->addMol(m);
  if (fps) {
    unsigned int fpsize = fps->addMol(m);
    CHECK_INVARIANT(size == fpsize,
                    "#mols different than #fingerprints in SubstructLibrary");
  }
  return size;
}

namespace {
//  Return true if the pattern contains a ring query
bool query_needs_rings(const ROMol &in_query) {
  for (auto &atom: in_query.atoms()) {
    if(atom->hasQuery()) {
      if (describeQuery(atom).find("Ring") != std::string::npos) {
        return true;
      }
    }
  }
  for (auto &bond: in_query.bonds()) {
    if(bond->hasQuery()) {
      if (describeQuery(bond).find("Ring") != std::string::npos) {
        return true;
      }
    }
  }
  return false;
}

void SubSearcher(const ROMol &in_query, const Bits &bits,
                 const MolHolderBase &mols, unsigned int start,
                 unsigned int &end, unsigned int numThreads,
                 const bool needs_rings, int &counter, const int maxResults,
                 std::vector<unsigned int> *idxs) {
  ROMol query(in_query);
  MatchVectType matchVect;
  for (unsigned int idx = start; idx < end; idx += numThreads) {
    if (!bits.check(idx)) {
      continue;
    }
    // need shared_ptr as it (may) control the lifespan of the
    //  returned molecule!
    const boost::shared_ptr<ROMol> &m = mols.getMol(idx);
    ROMol *mol = m.get();
    if (needs_rings && (!mol->getRingInfo() || !mol->getRingInfo()->isInitialized())) {
      MolOps::symmetrizeSSSR(*mol);
    }
    
    if (SubstructMatch(*mol, query, matchVect, bits.recursionPossible,
                       bits.useChirality, bits.useQueryQueryMatches)) {
      ++counter;
      if (idxs) {
        idxs->push_back(idx);
        if (maxResults > 0 && counter == maxResults) {
          // if we reached maxResults, record the last idx we processed and bail
          // out
          end = idx;
          break;
        }
      }
    }
  }
}

int internalGetMatches(
    const ROMol &query, MolHolderBase &mols, const FPHolderBase *fps,
    unsigned int startIdx, unsigned int endIdx, bool recursionPossible,
    bool useChirality, bool useQueryQueryMatches,
    int numThreads = -1, int maxResults = 1000,
    std::vector<unsigned int> *idxs = nullptr) {
  PRECONDITION(startIdx < mols.size(), "startIdx out of bounds");
  PRECONDITION(endIdx > startIdx, "endIdx > startIdx");

  // do not do any work if no results were requested
  if (maxResults == 0) {
    return 0;
  }

  endIdx = std::min(mols.size(), endIdx);

  numThreads = static_cast<int>(getNumThreadsToUse(numThreads));
  numThreads = std::min(numThreads, static_cast<int>(endIdx));

  bool needs_rings = query_needs_rings(query);
  Bits bits(fps, query, recursionPossible, useChirality, useQueryQueryMatches);
  int counter = 0;

#ifdef RDK_THREADSAFE_SSS
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
    if (idxs) {
      internal_results.resize(numThreads);
    }
    int thread_group_idx;
    for (thread_group_idx = 0; thread_group_idx < numThreads;
         ++thread_group_idx) {
      // need to use boost::ref otherwise things are passed by value
      thread_group.emplace_back(
          std::async(std::launch::async, SubSearcher, std::ref(query), bits,
                     std::ref(mols), startIdx + thread_group_idx,
                     std::ref(endIdxVect[thread_group_idx]), numThreads,
                     needs_rings, std::ref(counterVect[thread_group_idx]),
                     maxResultsVect[thread_group_idx],
                     idxs ? &internal_results[thread_group_idx] : nullptr));
    }
    if (maxResults > 0) {
      // If we are running with maxResults in a multi-threaded settings,
      // some threads may have screened more molecules than others.
      // If maxResults was close to the theoretical maximum, some threads
      // might have even run out of molecules to screen without reaching
      // maxResults so we need to make sure that all threads have screened as
      // many molecules as the most productive thread if we want multi-threaded
      // runs to yield the same results independently from the number of
      // threads.
      thread_group_idx = 0;
      for (auto &fut : thread_group) {
        fut.get();
        counter += counterVect[thread_group_idx++];
      }
      thread_group.clear();
      // Find out out the max number of molecules that was screened by the most
      // productive thread and do the same in all other threads, unless the
      // max number of molecules was reached
      unsigned int maxEndIdx =
          *std::max_element(endIdxVect.begin(), endIdxVect.end());
      for (thread_group_idx = 0; thread_group_idx < numThreads;
           ++thread_group_idx) {
        if (endIdxVect[thread_group_idx] >= maxEndIdx) {
          continue;
        }
        // need to use boost::ref otherwise things are passed by value
        thread_group.emplace_back(std::async(
            std::launch::async, SubSearcher, std::ref(query), bits,
            std::ref(mols), endIdxVect[thread_group_idx] + numThreads,
            std::ref(maxEndIdx), numThreads, needs_rings,
            std::ref(counterVect[thread_group_idx]), -1,
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
                maxResults, idxs);
  }
  if (idxs) {
    // the sort is necessary to ensure consistency across runs with different
    // numbers of threads
    std::sort(idxs->begin(), idxs->end());
    // we may have actually accumulated more results than maxResults due
    // to the top up above, so trim the results down if that's the case
    if (maxResults > 0 &&
        idxs->size() > static_cast<unsigned int>(maxResults)) {
      idxs->resize(maxResults);
    }
  }
#else
  SubSearcher(query, bits, mols, startIdx, endIdx, 1, needs_rings, counter, maxResults, idxs);
#endif

  delete bits.queryBits;

  return counter;
}

}

std::vector<unsigned int> SubstructLibrary::getMatches(
    const ROMol &query, bool recursionPossible, bool useChirality,
    bool useQueryQueryMatches, int numThreads, int maxResults) {
  return getMatches(query, 0, mols->size(), recursionPossible, useChirality,
                    useQueryQueryMatches, numThreads, maxResults);
}

std::vector<unsigned int> SubstructLibrary::getMatches(
    const ROMol &query, unsigned int startIdx, unsigned int endIdx,
    bool recursionPossible, bool useChirality, bool useQueryQueryMatches,
    int numThreads, int maxResults) {
  std::vector<unsigned int> idxs;
  internalGetMatches(query, *mols, fps, startIdx, endIdx,
                     recursionPossible, useChirality,
                     useQueryQueryMatches, numThreads, maxResults, &idxs);
  return idxs;
}

unsigned int SubstructLibrary::countMatches(const ROMol &query,
                                            bool recursionPossible,
                                            bool useChirality,
                                            bool useQueryQueryMatches,
                                            int numThreads) {
  return countMatches(query, 0, mols->size(), recursionPossible, useChirality,
                      useQueryQueryMatches, numThreads);
}

unsigned int SubstructLibrary::countMatches(
    const ROMol &query, unsigned int startIdx, unsigned int endIdx,
    bool recursionPossible, bool useChirality, bool useQueryQueryMatches,
    int numThreads) {
  return internalGetMatches(query, *mols, fps, startIdx, endIdx,
                            recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads, -1);
}

bool SubstructLibrary::hasMatch(const ROMol &query, bool recursionPossible,
                                bool useChirality, bool useQueryQueryMatches,
                                int numThreads) {
  const int maxResults = 1;
  return getMatches(query, recursionPossible, useChirality,
                    useQueryQueryMatches, numThreads, maxResults)
             .size() > 0;
}

bool SubstructLibrary::hasMatch(const ROMol &query, unsigned int startIdx,
                                unsigned int endIdx, bool recursionPossible,
                                bool useChirality, bool useQueryQueryMatches,
                                int numThreads) {
  const int maxResults = 1;
  return getMatches(query, startIdx, endIdx, recursionPossible, useChirality,
                    useQueryQueryMatches, numThreads, maxResults)
             .size() > 0;
}

void SubstructLibrary::toStream(std::ostream &ss) const {
#ifndef RDK_USE_BOOST_SERIALIZATION
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

}
