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
#include "Utilites.h"

namespace RDKit {
 namespace StructureCheck {

     RDKit::Bond::BondType convertBondType(AABondType bt);
     AABondType convertBondType(RDKit::Bond::BondType rdbt);

    unsigned getAtomicNumber(const std::string symbol);
    bool AtomSymbolMatch(const std::string symbol, const std::string pattern);
    bool LigandMatches(const Atom &a, const Bond &b, const Ligand &l, bool use_charge = false);
    bool RecMatch(const ROMol &mol, std::vector<unsigned> &match, unsigned int level,
                    const AugmentedAtom &aa, const std::vector<Neighbourhood> &nbp);
    bool AAMatch(const ROMol &mol, unsigned i,
                    std::vector<unsigned> &match,
                    const AugmentedAtom &aa,
                    const std::vector<unsigned> &atom_ring_status,
                    const std::vector<Neighbourhood> &nbp);

    bool TransformAugmentedAtoms(RWMol &mol, const std::vector<std::pair<AugmentedAtom, AugmentedAtom> > &aapair);
    bool CheckAtoms(const ROMol &mol, const std::vector<AugmentedAtom> &good_atoms);
 }
}

