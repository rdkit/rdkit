//  Copyright (c) 2019-2022 Brian P Kelley
//  All rights reserved.
//
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "ThreadedFilterCatalogRunner.h"

namespace RDKit {

namespace {
struct DoNothingThreadInitializer {
};
}
  
std::vector<std::vector<boost::shared_ptr<const FilterCatalogEntry>>>
RunFilterCatalog(const FilterCatalog &fc,
			 const std::vector<std::string> &smiles, int numThreads) {
  return ThreadedRunFilterCatalog<DoNothingThreadInitializer>(fc, smiles, numThreads);
}
  
}  // namespace RDKit
