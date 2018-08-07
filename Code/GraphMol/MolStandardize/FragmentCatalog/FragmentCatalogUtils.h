//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
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
std::vector<std::shared_ptr<ROMol>> readFuncGroups(std::istream &inStream,
                                                   int nToRead = -1);
// ROMol* removeFrags(const ROMol &mol, const FragmentCatalogParams *params);

}  // namespace MolStandardize
}  // namespace RDKit

#endif
