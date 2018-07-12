#ifndef __RD_ACIDBASE_CATALOG_UTILS_H__
#define __RD_ACIDBASE_CATALOG_UTILS_H__

#include <GraphMol/RDKitBase.h>
#include "AcidBaseCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class AcidBaseCatalogParams;

std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> readPairs(std::string fileName);
std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> readPairs(std::istream &inStream, int nToRead = -1);

} // namespace MolStandardize
} // namespace RDKit

#endif
