//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include "StructChecker.h"

namespace RDKit {
namespace StructureCheck {
bool StripSmallFragments(RWMol &mol, bool verbose=false);
void AddMWMF(
    RWMol &mol,
    bool pre);  // set mol formula & mass properties "MW_PRE" or "MW_POST"
}
}
