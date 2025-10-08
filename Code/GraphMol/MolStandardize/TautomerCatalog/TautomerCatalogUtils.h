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
#ifndef RD_TAUTOMER_CATALOG_UTILS_H
#define RD_TAUTOMER_CATALOG_UTILS_H

#include <GraphMol/RDKitBase.h>
#include "TautomerCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/Bond.h>
#include <utility>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class TautomerCatalogParams;

// typedef enum {
//	SINGLE,
//	DOUBLE,
//	TRIPLE,
//	AROMATIC,
//} BondType;

// typedef std::vector<ROMol*, std::string, std::string> tautomerTransform;

RDKIT_MOLSTANDARDIZE_EXPORT std::vector<Bond::BondType> stringToBondType(
    std::string bond_str);
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<int> stringToCharge(
    std::string charge_str);

RDKIT_MOLSTANDARDIZE_EXPORT std::vector<TautomerTransform> readTautomers(
    std::string fileName);
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<TautomerTransform> readTautomers(
    std::istream &inStream, int nToRead = -1);
RDKIT_MOLSTANDARDIZE_EXPORT std::vector<TautomerTransform> readTautomers(
    const std::vector<
        std::tuple<std::string, std::string, std::string, std::string>> &data);
}  // namespace MolStandardize
}  // namespace RDKit

#endif
