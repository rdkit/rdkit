#ifndef __RD_FRAGMENT_CATALOG_UTILS_H__
#define __RD_FRAGMENT_CATALOG_UTILS_H__

#include <GraphMol/RDKitBase.h>
#include "FragmentCatalogParams.h"
#include <GraphMol/Substruct/SubstructMatch.h>
#include <iostream>

namespace RDKit {
class ROMol;

namespace MolStandardize {
class FragmentCatalogParams;

std::vector<std::shared_ptr<ROMol>> readFuncGroups(std::string fileName);
std::vector<std::shared_ptr<ROMol>> readFuncGroups(std::istream &inStream, int nToRead = -1);
//ROMol* removeFrags(const ROMol &mol, const FragmentCatalogParams *params);

} // namespace MolStandardize
} // namespace RDKit

#endif
