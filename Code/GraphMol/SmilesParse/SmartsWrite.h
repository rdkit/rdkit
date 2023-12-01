//
//  Copyright (C) 2004-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_SMARTSWRITE_H_012020
#define RD_SMARTSWRITE_H_012020

#include <string>
#include <vector>

namespace RDKit {
class Atom;
class Bond;
namespace SmartsWrite {
//! returns the SMARTS for an Atom
RDKIT_SMILESPARSE_EXPORT std::string GetAtomSmarts(const Atom *qatom);
//! returns the SMARTS for a Bond
RDKIT_SMILESPARSE_EXPORT std::string GetBondSmarts(const Bond *qbond,
                                                   int atomToLeftIdx = -1);
}  // namespace SmartsWrite

class ROMol;
//! returns the SMARTS for a molecule
RDKIT_SMILESPARSE_EXPORT std::string MolToSmarts(const ROMol &mol,
                                                 bool doIsomericSmarts = true,
                                                 int rootedAtAtom = -1);

RDKIT_SMILESPARSE_EXPORT std::string MolFragmentToSmarts(
    const ROMol &mol, const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr, bool doIsomericSmarts = true);

//! returns the CXSMARTS for a molecule
RDKIT_SMILESPARSE_EXPORT std::string MolToCXSmarts(
    const ROMol &mol, bool doIsomericSmarts = true);

RDKIT_SMILESPARSE_EXPORT std::string MolFragmentToCXSmarts(
    const ROMol &mol, const std::vector<int> &atomsToUse,
    const std::vector<int> *bondsToUse = nullptr, bool doIsomericSmarts = true);
};  // namespace RDKit

#endif
