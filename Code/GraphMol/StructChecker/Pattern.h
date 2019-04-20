//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include "StructChecker.h"
#include "Utilites.h"

namespace RDKit {
namespace StructureCheck {

RDKIT_STRUCTCHECKER_EXPORT RDKit::Bond::BondType convertBondType(AABondType bt);
RDKIT_STRUCTCHECKER_EXPORT AABondType
convertBondType(RDKit::Bond::BondType rdbt);

RDKIT_STRUCTCHECKER_EXPORT unsigned getAtomicNumber(const std::string symbol);
RDKIT_STRUCTCHECKER_EXPORT bool AtomSymbolMatch(const std::string symbol,
                                                const std::string pattern);
RDKIT_STRUCTCHECKER_EXPORT bool LigandMatches(const Atom &a, const Bond &b,
                                              const Ligand &l,
                                              bool use_charge = false);
RDKIT_STRUCTCHECKER_EXPORT bool isBondTypeMatch(const RDKit::Bond &b,
                                                AABondType lbt);
RDKIT_STRUCTCHECKER_EXPORT bool RecMatch(const ROMol &mol, unsigned atomIdx,
                                         const AugmentedAtom &aa,
                                         const std::vector<Neighbourhood> &nbp,
                                         bool verbose);
RDKIT_STRUCTCHECKER_EXPORT bool AAMatch(
    const ROMol &mol, unsigned i, const AugmentedAtom &aa,
    const std::vector<unsigned> &atom_ring_status,
    const std::vector<Neighbourhood> &nbp, bool verbose);

RDKIT_STRUCTCHECKER_EXPORT bool TransformAugmentedAtoms(
    RWMol &mol,
    const std::vector<std::pair<AugmentedAtom, AugmentedAtom>> &aapair,
    bool verbose);
RDKIT_STRUCTCHECKER_EXPORT bool CheckAtoms(
    const ROMol &mol, const std::vector<AugmentedAtom> &good_atoms,
    bool verbose);
}  // namespace StructureCheck
}  // namespace RDKit
