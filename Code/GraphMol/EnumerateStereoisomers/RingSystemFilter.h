//
//  Copyright (C) 2004-2026 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <set>
#include <utility>

namespace RDKit {
class ROMol;
}

void getRingPatternsParityRelations(
    const RDKit::ROMol &mol,
    std::set<std::pair<unsigned int, unsigned int>> &same,
    std::set<std::pair<unsigned int, unsigned int>> &opposite);
