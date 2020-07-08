//
//  Copyright (C) 2017 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_MATCH_DATA
#define RGROUP_MATCH_DATA
#include "RGroupData.h"

namespace RDKit {
typedef boost::shared_ptr<RGroupData> RData;
typedef std::map<int, RData> R_DECOMP;

//! RGroupMatch is the decomposition for a single molecule
struct RGroupMatch {
  size_t   core_idx; // index of the matching core
  R_DECOMP rgroups;  // rlabel->RGroupData mapping

  RGroupMatch(size_t core_index, R_DECOMP input_rgroups)
      : core_idx(core_index), rgroups(std::move(input_rgroups)) {}
};  

}
#endif
