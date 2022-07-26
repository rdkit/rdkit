//
//  Copyright (C) 2020 Brian P. Kelley
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SubstructLibrary.h"
#include "PatternFactory.h"
#include <RDGeneral/RDThreads.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#include <future>
#endif

namespace RDKit {
namespace {
void fillPatterns(const SubstructLibrary &slib, const FPHolderBase &fph,
                  std::vector<ExplicitBitVect *> &fps, unsigned int start,
                  unsigned int end, unsigned int numThreads) {
  for (unsigned int idx = start; idx < end; idx += numThreads) {
    auto mol = slib.getMol(idx);
    if (mol.get()) {
      fps[idx] = fph.makeFingerprint(*mol.get());
    } else {
      // Make an empty FP
      fps[idx] = fph.makeFingerprint(ROMol());
    }
  }
}
}  // namespace

namespace {
void internalAddPatterns(SubstructLibrary &sslib, int numThreads,
                         boost::shared_ptr<FPHolderBase> *patterns) {
  PRECONDITION(sslib.getFpHolder().get() == nullptr,
               "Substruct library already has fingerprints");
  numThreads = (int)getNumThreadsToUse(numThreads);

  boost::shared_ptr<FPHolderBase> ptr;
  if (patterns) {
    ptr = *patterns;
  } else {
    ptr = boost::shared_ptr<FPHolderBase>(new PatternHolder);
  }

  std::vector<ExplicitBitVect *> &fps = ptr->getFingerprints();
  unsigned int startIdx = 0;
  unsigned int endIdx = sslib.getMolecules().size();
  fps.resize(endIdx);

#ifdef RDK_BUILD_THREADSAFE_SSS
  std::vector<std::future<void>> thread_group;
  for (int thread_group_idx = 0; thread_group_idx < numThreads;
       ++thread_group_idx) {
    // need to use std::ref otherwise things are passed by value
    thread_group.emplace_back(std::async(
        std::launch::async, fillPatterns, std::ref(sslib), std::ref(*ptr),
        std::ref(fps), startIdx + thread_group_idx, endIdx, numThreads));
  }
  for (auto &fut : thread_group) {
    fut.get();
  }
#else
  fillPatterns(sslib, *ptr, fps, 0, sslib.size(), 1);
#endif
  if (ptr->size() != sslib.size()) {
    throw ValueErrorException(
        "Number of fingerprints generated not equal to current number of "
        "molecules");
  }

  sslib.getFpHolder() = ptr;
  sslib.resetHolders();
}
}  // namespace

void addPatterns(SubstructLibrary &sslib, int numThreads) {
  internalAddPatterns(sslib, numThreads, nullptr);
}

void addPatterns(SubstructLibrary &sslib,
                 boost::shared_ptr<FPHolderBase> patterns, int numThreads) {
  PRECONDITION(patterns.get() != nullptr, "Null PatternHolder");
  internalAddPatterns(sslib, numThreads, &patterns);
}

}  // namespace RDKit
