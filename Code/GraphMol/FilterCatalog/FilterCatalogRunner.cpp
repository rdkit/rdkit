//  Copyright (c) 2019 Brian P Kelley
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

#include "FilterCatalog.h"
#include "Filters.h"
#include "FilterMatchers.h"
#include <GraphMol/SmilesParse/SmilesParse.h>

#ifdef RDK_THREADSAFE_SSS
#include <RDGeneral/RDThreads.h>
#include <thread>
#include <future>
#endif

namespace RDKit {
namespace {
void CatalogSearcher(const FilterCatalog &fc,
		     const std::vector<std::string> &smiles,
		     std::vector<std::vector<FilterCatalog::CONST_SENTRY>> &results,
		     int start,
		     int numThreads) {
  for(unsigned int idx = start;
      idx < smiles.size();
      idx += numThreads) {
    std::unique_ptr<ROMol> mol(SmilesToMol(smiles[idx]));
    if(mol.get()) {
      results[idx] = fc.getMatches(*mol);
    }
  } 
}
}
  
std::vector<std::vector<FilterCatalog::CONST_SENTRY>> RunFilterCatalog(
              const FilterCatalog &fc,
	      const std::vector<std::string> &smiles,
	      int numThreads) {
  // preallocate results, hopefully this won't grow in the thread...
  std::vector<std::vector<FilterCatalog::CONST_SENTRY>> results(smiles.size());

#ifdef RDK_THREADSAFE_SSS  
  if (numThreads == -1)
    numThreads = (int)getNumThreadsToUse(numThreads);
  else
    numThreads = std::min(numThreads, (int)getNumThreadsToUse(numThreads));

  std::vector<std::future<void>> thread_group;  
  for (int thread_group_idx = 0; thread_group_idx < numThreads;
       ++thread_group_idx) {
    // need to use boost::ref otherwise things are passed by value
    thread_group.emplace_back(
		   std::async(std::launch::async, CatalogSearcher,
			      std::ref(fc),
			      std::ref(smiles),
			      std::ref(results),
			      thread_group_idx,
			      numThreads));
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
  
}
