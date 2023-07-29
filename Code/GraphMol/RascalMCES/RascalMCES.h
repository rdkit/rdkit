//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDKIT_RASCAL_MCES_H
#define RDKIT_RASCAL_MCES_H

#include <vector>

#include <GraphMol/RascalMCES/RascalOptions.h>
#include <GraphMol/RascalMCES/RascalResult.h>
namespace RDKit {
class ROMol;

namespace RascalMCES {

// Find one or more MCESs between the two molecules.  The MCES is the
// Maximum Common Edge Substructure, and is the largest set of bonds
// common to the 2 molecules.
/*!
 *
 * @param mol1 : first molecule
 * @param mol2 : second molecule for MCES determination.
 * @param opts : (optional) set of options controlling the MCES determination
 * @return : vector of RascalResult objects.
 */
RDKIT_RASCALMCES_EXPORT std::vector<RascalResult> rascalMces(
    const ROMol &mol1, const ROMol &mol2,
    const RascalOptions &opts = RascalOptions());

}  // namespace RascalMCES
}  // namespace RDKit
#endif  // RDKIT_RASCAL_MCES_H
