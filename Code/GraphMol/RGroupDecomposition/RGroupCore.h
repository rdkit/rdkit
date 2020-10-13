//
//  Copyright (C) 2020 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RGROUP_CORE
#define RGROUP_CORE

#include "../RDKitBase.h"
#include "RGroupUtils.h"

namespace RDKit {

//! RCore is the core common to a series of molecules
struct RCore {
  boost::shared_ptr<RWMol> core;
  boost::shared_ptr<RWMol> labelledCore;
  std::set<int> core_atoms_with_user_labels;
  RCore(){};
  RCore(const RWMol &c, bool onlyMatchAtRGroups = false) : core(new RWMol(c)) {
    if (onlyMatchAtRGroups) {
      findIndicesWithRLabel();
    }
  }
  void findIndicesWithRLabel() {
    // First find all the core atoms that have user
    //  label and but their indices into core_atoms_with_user_labels
    for (const auto atom : core->atoms()) {
      if (atom->hasProp(RLABEL)) {
        core_atoms_with_user_labels.insert(atom->getIdx());
      }
    }
  }
};

}  // namespace RDKit
#endif
