#ifndef __RD_TRANSFORM_CATALOG_UTILS_H__
#define __RD_TRANSFORM_CATALOG_UTILS_H__

#include <GraphMol/RDKitBase.h>
#include "TransformCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class TransformCatalogParams;

std::vector<std::shared_ptr<ChemicalReaction>> readTransformations(std::string fileName);
std::vector<std::shared_ptr<ChemicalReaction>> readTransformations(std::istream &inStream, int nToRead = -1);

} // namespace MolStandardize
} // namespace RDKit

#endif
