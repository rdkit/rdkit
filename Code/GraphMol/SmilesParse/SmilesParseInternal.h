//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include "SmilesParse.h"

#include <GraphMol/RDMol.h>

#include <vector>

namespace SmilesParseInternal {
bool parseSmiles(const char* text, RDKit::RDMol& mol,
                 RDKit::SmilesParseTemp& temp);
}
