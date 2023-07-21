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
namespace RDKit {
class ROMol;

namespace RascalMCES {

class RascalResult;

// Find one or more MCESs between the two molecules.
RDKIT_RASCAL_MCES_EXPORT std::vector<RascalResult> rascalMces(
    const ROMol &mol1, const ROMol &mol2, RascalOptions opts = RascalOptions());

}  // namespace RascalMCES
}  // namespace RDKit
#endif  // RDKIT_RASCAL_MCES_H
