//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_SMARTSWRITE_H
#define _RD_SMARTSWRITE_H

#include <string>

namespace RDKit {
class QueryAtom;
class QueryBond;
namespace SmartsWrite {
//! returns the SMARTS for a QueryAtom
RDKIT_SMILESPARSE_EXPORT std::string GetAtomSmarts(const QueryAtom *qatom);
//! returns the SMARTS for a QueryBond
RDKIT_SMILESPARSE_EXPORT std::string GetBondSmarts(const QueryBond *qbond,
                                                   int atomToLeftIdx = -1);
}  // namespace SmartsWrite

class ROMol;
//! returns the SMARTS for a molecule
RDKIT_SMILESPARSE_EXPORT std::string MolToSmarts(ROMol &mol,
                                                 bool doIsomericSmarts = true);
};  // namespace RDKit

#endif
