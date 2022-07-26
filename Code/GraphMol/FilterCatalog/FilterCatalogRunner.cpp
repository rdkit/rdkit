//  Copyright (c) 2019 Brian P Kelley
//  All rights reserved.
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "FilterCatalog.h"
#include "Filters.h"
#include "FilterMatchers.h"
#include <GraphMol/SmilesParse/SmilesParse.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <RDGeneral/RDThreads.h>
#include <thread>
#include <future>
#endif

namespace RDKit {
namespace {
boost::shared_ptr<FilterCatalogEntry> &makeBadSmilesEntry() {
  static boost::shared_ptr<FilterCatalogEntry> bad_smiles(
      new FilterCatalogEntry("no valid RDKit molecule",
                             boost::shared_ptr<FilterMatcherBase>()));
  return bad_smiles;
}
void CatalogSearcher(
    const FilterCatalog &fc, const std::vector<std::string> &smiles,
    std::vector<std::vector<FilterCatalog::CONST_SENTRY>> &results, int start,
    int numThreads) {
  for (unsigned int idx = start; idx < smiles.size(); idx += numThreads) {
    std::unique_ptr<ROMol> mol(SmilesToMol(smiles[idx]));
    if (mol.get()) {
      results[idx] = fc.getMatches(*mol);
    } else {
      results[idx].push_back(makeBadSmilesEntry());
    }
  }
}
}  // namespace

std::vector<std::vector<boost::shared_ptr<const FilterCatalogEntry>>>
RunFilterCatalog(const FilterCatalog &fc,
                 const std::vector<std::string> &smiles, int numThreads) {
  // preallocate results so the threads don't move the vector around in memory
  //  There is one result per input smiles
  std::vector<std::vector<FilterCatalog::CONST_SENTRY>> results(smiles.size());

#ifdef RDK_BUILD_THREADSAFE_SSS
  std::vector<std::future<void>> thread_group;
  numThreads = (int)getNumThreadsToUse(numThreads);
  for (int thread_group_idx = 0; thread_group_idx < numThreads;
       ++thread_group_idx) {
    // need to use std::ref otherwise things are passed by value
    thread_group.emplace_back(std::async(
        std::launch::async, CatalogSearcher, std::ref(fc), std::ref(smiles),
        std::ref(results), thread_group_idx, numThreads));
  }
  for (auto &fut : thread_group) {
    fut.get();
  }

#else
  int start = 0;
  numThreads = 1;
  CatalogSearcher(fc, smiles, results, start, numThreads);
#endif
  return results;
}

}  // namespace RDKit
