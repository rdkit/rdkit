//
//  Copyright (C) 2001-2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <string_view>

#include <RDGeneral/export.h>
#ifndef RD_SMILESPARSEOPS_H
#define RD_SMILESPARSEOPS_H
#include <GraphMol/Bond.h>

namespace RDKit {
class RWMol;
class Atom;
class QueryBond;
}  // namespace RDKit
namespace SmilesParseOps {
RDKIT_SMILESPARSE_EXPORT void CheckRingClosureBranchStatus(RDKit::Atom *atom,
                                                           RDKit::RWMol *mp);
RDKIT_SMILESPARSE_EXPORT void ReportParseError(const char *message,
                                               bool throwIt = true);
RDKIT_SMILESPARSE_EXPORT void CleanupAfterParseError(RDKit::RWMol *mol);
// This uses SMARTS semantics: unspecified bonds are treated as
// aromatic or single.
RDKIT_SMILESPARSE_EXPORT void AddFragToMol(
    RDKit::RWMol *mol, RDKit::RWMol *frag,
    RDKit::Bond::BondType bondOrder = RDKit::Bond::UNSPECIFIED,
    RDKit::Bond::BondDir bondDir = RDKit::Bond::NONE);
RDKIT_SMILESPARSE_EXPORT RDKit::Bond::BondType GetUnspecifiedBondType(
    const RDKit::RWMol *mol, const RDKit::Atom *atom1,
    const RDKit::Atom *atom2);
RDKIT_SMILESPARSE_EXPORT void CheckChiralitySpecifications(RDKit::RWMol *mol,
                                                           bool strict);
RDKIT_SMILESPARSE_EXPORT void CloseMolRings(RDKit::RWMol *mol,
                                            bool toleratePartials);
RDKIT_SMILESPARSE_EXPORT void SetUnspecifiedBondTypes(RDKit::RWMol *mol);
RDKIT_SMILESPARSE_EXPORT void AdjustAtomChiralityFlags(RDKit::RWMol *mol);
RDKIT_SMILESPARSE_EXPORT void CleanupAfterParsing(RDKit::RWMol *mol);
RDKIT_SMILESPARSE_EXPORT void parseCXExtensions(
    RDKit::RWMol &mol, const std::string &extText,
    std::string::const_iterator &pos, unsigned int startAtomIdx = 0,
    unsigned int startBondIdx = 0);
inline void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
                              unsigned int startAtomIdx,
                              unsigned int startBondIdx) {
  auto iter = extText.begin();
  parseCXExtensions(mol, extText, iter, startAtomIdx, startBondIdx);
};
//! removes formal charge, isotope, etc. Primarily useful for QueryAtoms
RDKIT_SMILESPARSE_EXPORT void ClearAtomChemicalProps(RDKit::Atom *atom);

//! returns whether or not the combination of tag and permutation provided are
//! legal
RDKIT_SMILESPARSE_EXPORT bool checkChiralPermutation(int chiralTag,
                                                     int permutation);

//! this is a bit of a hack to try and get nicer "SMILES" from
//! a SMARTS molecule
RDKIT_SMILESPARSE_EXPORT RDKit::QueryBond *getUnspecifiedQueryBond(
    const RDKit::Atom *a1, const RDKit::Atom *a2);

namespace detail {
constexpr auto _needsDetectBondStereo = "_needsDetectBondStereo";
constexpr auto _needsDetectAtomStereo = "_needsDetectAtomStereo";

void printSyntaxErrorMessage(std::string_view input,
                             std::string_view err_message,
                             unsigned int bad_token_position);
}  // namespace detail
}  // namespace SmilesParseOps

#endif
