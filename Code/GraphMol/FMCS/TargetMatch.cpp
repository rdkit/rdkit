//
//  Copyright (C) 2024 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include "TargetMatch.h"
#include "Target.h"
#include "Seed.h"
namespace RDKit {
namespace FMCS {

void TargetMatch::init(const Seed &seed, const match_V_t &match,
                       const ROMol &query, const Target &target) {
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
    TargetAtomIdx[seed.MoleculeFragment.Atoms.at(m.first)->getIdx()] = m.second;
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

}  // namespace FMCS
}  // namespace RDKit
