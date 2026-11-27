//
//  Copyright (C) 2026 ETH Zurich
//  Created by: Katharina Buchthal
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RDKIT_GRAPHMOL_JOINTRING_H
#define RDKIT_GRAPHMOL_JOINTRING_H
#include <RDGeneral/export.h>
#include <vector>

namespace RDKit {

class ROMol;

// a joint ring: indices of corresponding atoms
struct RDKIT_GRAPHMOL_EXPORT JointRingConstraint {
  std::vector<unsigned int> joint;
  std::vector<unsigned int> excluded;
};

// joint ring ids an atom is associated with
struct RDKIT_GRAPHMOL_EXPORT JointRingIds {
  std::vector<unsigned int> joint;
  std::vector<unsigned int> excluded;
};

RDKIT_GRAPHMOL_EXPORT std::vector<JointRingConstraint> collectJRConstraints(
    const ROMol &mol);
}  // namespace RDKit
#endif