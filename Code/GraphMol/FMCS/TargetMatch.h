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
#include "MatchTable.h"

namespace RDKit {
namespace FMCS {
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
            const Target &target) {
    TargetAtomIdx.clear();
    TargetAtomIdx.resize(query.getNumAtoms(), NotSet);
    TargetBondIdx.clear();
    TargetBondIdx.resize(query.getNumBonds(), NotSet);
    VisitedTargetBonds.resize(target.Molecule->getNumBonds());
    VisitedTargetAtoms.resize(target.Molecule->getNumAtoms());
    VisitedTargetBonds.reset();
    VisitedTargetAtoms.reset();

    MatchedAtomSize = match.size();
    for (const auto &m : match) {
      TargetAtomIdx[seed.MoleculeFragment.Atoms.at(m.first)->getIdx()] =
          m.second;
      VisitedTargetAtoms.set(m.second);
    }

    MatchedBondSize = 0;
    for (const auto bond : seed.MoleculeFragment.Bonds) {
      unsigned int i = bond->getBeginAtomIdx();
      unsigned int j = bond->getEndAtomIdx();
      unsigned int ti = TargetAtomIdx.at(i);
      unsigned int tj = TargetAtomIdx.at(j);
      const auto tb = target.Molecule->getBondBetweenAtoms(ti, tj);
      if (tb) {
        ++MatchedBondSize;
        TargetBondIdx[bond->getIdx()] = tb->getIdx();  // add
        VisitedTargetBonds.set(tb->getIdx());
      }
    }
    Empty = false;
  }
};
}  // namespace FMCS
}  // namespace RDKit
