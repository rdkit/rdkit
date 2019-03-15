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

#include <atomic>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit {

bool SubstructLibraryCanSerialize() {
#ifdef RDK_USE_BOOST_SERIALIZATION
  return true;
#else
  return false;
#endif
}

Bits::Bits(const FPHolderBase *fps, const ROMol &m, bool recursionPossible,
       bool useChirality, bool useQueryQueryMatches)
      : fps(fps),
        recursionPossible(recursionPossible),
        useChirality(useChirality),
        useQueryQueryMatches(useQueryQueryMatches) {
    if (fps) {
      queryBits = fps->makeFingerprint(m);
    } else
      queryBits = nullptr;
  }

bool Bits::check(unsigned int idx) const {
    if (fps) {
      return fps->passesFilter(idx, *queryBits);
    }
    return true;
  }

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

// end is exclusive here
void SubSearcher(const ROMol &in_query, const Bits &bits,
                 const MolHolderBase &mols, std::vector<unsigned int> &idxs,
                 unsigned int start, unsigned int end, unsigned int numThreads,
                 std::atomic<int> &counter, const int maxResults) {
  ROMol query(in_query);
  MatchVectType matchVect;
  for (unsigned int idx = start;
       idx < end && (maxResults == -1 || counter < maxResults);
       idx += numThreads) {
    if (!bits.check(idx)) continue;
    // need shared_ptr as it (may) control the lifespan of the
    //  returned molecule!
    const boost::shared_ptr<ROMol> &m = mols.getMol(idx);
    const ROMol *mol = m.get();
    if (SubstructMatch(*mol, query, matchVect, bits.recursionPossible,
                       bits.useChirality, bits.useQueryQueryMatches)) {
      // this is squishy when updating the counter.  While incrementing is
      // atomic
      // several substructure runs can update the counter beyond the maxResults
      //  This okay: if we get one or two extra, we can fix it on the way out
      if (maxResults != -1 && counter >= maxResults) break;
      idxs.push_back(idx);
      if (maxResults != -1) counter++;
    }
  }
}

// end is inclusive here
void SubSearchMatchCounter(const ROMol &in_query, const Bits &bits,
                           const MolHolderBase &mols, unsigned int start,
                           unsigned int end, int numThreads,
                           std::atomic<int> &counter) {
  ROMol query(in_query);
  MatchVectType matchVect;
  for (unsigned int idx = start; idx < end; idx += numThreads) {
    if (!bits.check(idx)) continue;
    // need shared_ptr as it (may) controls the lifespan of the
    //  returned molecule!
    const boost::shared_ptr<ROMol> &m = mols.getMol(idx);
    const ROMol *mol = m.get();
    if (SubstructMatch(*mol, query, matchVect, bits.recursionPossible,
                       bits.useChirality, bits.useQueryQueryMatches)) {
      counter++;
    }
  }
}

std::vector<unsigned int> internalGetMatches(
    const ROMol &query, MolHolderBase &mols, const FPHolderBase *fps,
    unsigned int startIdx, unsigned int endIdx, bool recursionPossible,
    bool useChirality, bool useQueryQueryMatches, int numThreads = -1,
    int maxResults = 1000) {
  PRECONDITION(startIdx < mols.size(), "startIdx out of bounds");
  PRECONDITION(endIdx > startIdx, "endIdx > startIdx");
  if (numThreads == -1)
    numThreads = (int)getNumThreadsToUse(numThreads);
  else
    numThreads = std::min(numThreads, (int)getNumThreadsToUse(numThreads));

  endIdx = std::min(mols.size(), endIdx);
  if (endIdx < static_cast<unsigned int>(numThreads)) numThreads = endIdx;

  std::vector<std::future<void>> thread_group;
  std::atomic<int> counter(0);
  std::vector<std::vector<unsigned int>> internal_results(numThreads);

  Bits bits(fps, query, recursionPossible, useChirality, useQueryQueryMatches);

  for (int thread_group_idx = 0; thread_group_idx < numThreads;
       ++thread_group_idx) {
    // need to use boost::ref otherwise things are passed by value
    thread_group.emplace_back(
        std::async(std::launch::async, SubSearcher, std::ref(query), bits,
                   std::ref(mols), std::ref(internal_results[thread_group_idx]),
                   startIdx + thread_group_idx, endIdx, numThreads,
                   std::ref(counter), maxResults));
  }
  for (auto &fut : thread_group) {
    fut.get();
  }
  delete bits.queryBits;

  std::vector<unsigned int> results;
  for (int thread_group_idx = 0; thread_group_idx < numThreads;
       ++thread_group_idx) {
    results.insert(results.end(), internal_results[thread_group_idx].begin(),
                   internal_results[thread_group_idx].end());
  }

  // this is so we don't really have to do locking on the atomic counter...
  if (maxResults != -1 && rdcast<int>(results.size()) > maxResults)
    results.resize(maxResults);

  return results;
}

int internalMatchCounter(const ROMol &query, MolHolderBase &mols,
                         const FPHolderBase *fps, unsigned int startIdx,
                         unsigned int endIdx, bool recursionPossible,
                         bool useChirality, bool useQueryQueryMatches,
                         int numThreads = -1) {
  PRECONDITION(startIdx < mols.size(), "startIdx out of bounds");
  PRECONDITION(endIdx > startIdx, "endIdx > startIdx");

  endIdx = std::min(mols.size(), endIdx);

  if (numThreads == -1)
    numThreads = (int)getNumThreadsToUse(numThreads);
  else
    numThreads = std::min(numThreads, (int)getNumThreadsToUse(numThreads));

  if (endIdx < static_cast<unsigned int>(numThreads)) numThreads = endIdx;

  std::vector<std::future<void>> thread_group;
  std::atomic<int> counter(0);

  Bits bits(fps, query, recursionPossible, useChirality, useQueryQueryMatches);
  for (int thread_group_idx = 0; thread_group_idx < numThreads;
       ++thread_group_idx) {
    // need to use boost::ref otherwise things are passed by value
    thread_group.emplace_back(
        std::async(std::launch::async, SubSearchMatchCounter, std::ref(query),
                   bits, std::ref(mols), startIdx + thread_group_idx, endIdx,
                   numThreads, std::ref(counter)));
  }
  for (auto &thread : thread_group) {
    thread.get();
  }
  delete bits.queryBits;
  return (int)counter;
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
  return internalGetMatches(query, *mols, fps, startIdx, endIdx,
                            recursionPossible, useChirality,
                            useQueryQueryMatches, numThreads, maxResults);
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
  return internalMatchCounter(query, *mols, fps, startIdx, endIdx,
                              recursionPossible, useChirality,
                              useQueryQueryMatches, numThreads);
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
