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

    unsigned getAtomicNumber(const std::string symbol);
    bool AtomSymbolMatch(const std::string symbol, const std::string pattern);

    bool TransformAugmentedAtoms(RWMol &mol, const std::vector<std::pair<AugmentedAtom, AugmentedAtom> > &aapair);
    bool CheckAtoms(const ROMol &mol, const std::vector<AugmentedAtom> &good_atoms);
 }
}

