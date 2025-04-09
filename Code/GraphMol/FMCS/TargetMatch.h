//
//  Copyright (C) 2014 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "FMCS.h"
#include "SubstructMatchCustom.h"
#include "MatchTable.h"

namespace RDKit {
namespace FMCS {
class Seed;
struct Target;

struct TargetMatch {
  bool Empty{true};
  size_t MatchedAtomSize{0};
  size_t MatchedBondSize{0};
  std::vector<unsigned int> TargetAtomIdx;
  std::vector<unsigned int> TargetBondIdx;
  boost::dynamic_bitset<> VisitedTargetBonds;
  boost::dynamic_bitset<> VisitedTargetAtoms;  // for checking rings
 public:
  TargetMatch() {}
  TargetMatch(const TargetMatch &src) { *this = src; }
  TargetMatch &operator=(const TargetMatch &src) {
    Empty = src.Empty;
    if (!Empty) {
      MatchedAtomSize = src.MatchedAtomSize;
      MatchedBondSize = src.MatchedBondSize;
      TargetAtomIdx = src.TargetAtomIdx;
      TargetBondIdx = src.TargetBondIdx;
      VisitedTargetBonds = src.VisitedTargetBonds;
      VisitedTargetAtoms = src.VisitedTargetAtoms;
    }
    return *this;
  }
  bool empty() const { return Empty; }
  void clear() {
    Empty = true;

    TargetAtomIdx.clear();
    TargetBondIdx.clear();
    VisitedTargetBonds.clear();
    VisitedTargetAtoms.clear();
  }
  void init(const Seed &seed, const match_V_t &match, const ROMol &query,
            const Target &target);
};
}  // namespace FMCS
}  // namespace RDKit
