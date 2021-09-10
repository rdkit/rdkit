//
//  Copyright (C) 2018-2021 Susan H. Leung and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef RD_ACIDBASE_CATALOG_UTILS_H
#define RD_ACIDBASE_CATALOG_UTILS_H

#include <GraphMol/RDKitBase.h>
#include "AcidBaseCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class AcidBaseCatalogParams;

RDKIT_MOLSTANDARDIZE_EXPORT std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>>
readPairs(std::string fileName);
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>>
readPairs(std::istream &inStream, int nToRead = -1);
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>>
readPairs(
    const std::vector<std::tuple<std::string, std::string, std::string>> &data);

}  // namespace MolStandardize
}  // namespace RDKit

#endif
