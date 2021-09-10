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
#ifndef RD_TRANSFORM_CATALOG_UTILS_H
#define RD_TRANSFORM_CATALOG_UTILS_H

#include <GraphMol/RDKitBase.h>
#include "TransformCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class TransformCatalogParams;

RDKIT_MOLSTANDARDIZE_EXPORT std::vector<std::shared_ptr<ChemicalReaction>>
readTransformations(std::string fileName);
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<std::shared_ptr<ChemicalReaction>>
readTransformations(std::istream &inStream, int nToRead = -1);
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<std::shared_ptr<ChemicalReaction>>
readTransformations(
    const std::vector<std::pair<std::string, std::string>> &data);

}  // namespace MolStandardize
}  // namespace RDKit

#endif
