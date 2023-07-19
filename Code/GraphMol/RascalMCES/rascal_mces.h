//
// Copyright (C) David Cosgrove 2023
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_RASCAL_MCES_H
#define RDKIT_RASCAL_MCES_H

#include <vector>

namespace RDKit {
class ROMol;

namespace RascalMCES {

class RascalOptions;
class RascalResult;

// Find one or more MCESs between the two molecules.
std::vector<RascalResult> rascal_mces(const ROMol &mol1, const ROMol &mol2,
                                      RascalOptions opts);

}  // namespace RascalMCES
}  // namespace RDKit
#endif  // RDKIT_RASCAL_MCES_H
