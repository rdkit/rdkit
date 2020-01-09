//
//  Copyright (C) 2001-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_SMILESPARSEOPS_H
#define _RD_SMILESPARSEOPS_H
#include <GraphMol/Bond.h>

namespace RDKit {
class RWMol;
class Atom;
}  // namespace RDKit
namespace SmilesParseOps {
void CheckRingClosureBranchStatus(RDKit::Atom *atom, RDKit::RWMol *mp);
void ReportParseError(const char *message, bool throwIt = true);
void CleanupAfterParseError(RDKit::RWMol *mol);
// This uses SMARTS semantics: unspecified bonds are treated as
// aromatic or single.
void AddFragToMol(RDKit::RWMol *mol, RDKit::RWMol *frag,
                  RDKit::Bond::BondType bondOrder = RDKit::Bond::UNSPECIFIED,
                  RDKit::Bond::BondDir bondDir = RDKit::Bond::NONE);
RDKit::Bond::BondType GetUnspecifiedBondType(const RDKit::RWMol *mol,
                                             const RDKit::Atom *atom1,
                                             const RDKit::Atom *atom2);
void CloseMolRings(RDKit::RWMol *mol, bool toleratePartials);
void SetUnspecifiedBondTypes(RDKit::RWMol *mol);
void AdjustAtomChiralityFlags(RDKit::RWMol *mol);
void CleanupAfterParsing(RDKit::RWMol *mol);
void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
                       std::string::const_iterator &pos);
//! removes formal charge, isotope, etc. Primarily useful for QueryAtoms
void ClearAtomChemicalProps(RDKit::Atom *atom);
};  // namespace SmilesParseOps

#endif
