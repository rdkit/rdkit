//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#include <RDGeneral/export.h>
#ifndef RDKIT_HYPERSPACESUBSTRUCTURESEARCH_H
#define RDKIT_HYPERSPACESUBSTRUCTURESEARCH_H

#include <vector>

namespace RDKit {
class ROMol;
namespace HyperspaceSSSearch {

// Perform a substructure search with the given query molecule across
// the library defined in lib_name.
/*!
 *
 * @param query : query molecule
 * @param libName : name of library containing synthon-based molecule library
 * @param maxBondSplits : the maximum number of bonds to be used in each split.
 * If maxBondSplits is 1, all molecules are returned
 * where 1 bond is removed from the molecule resulting in 2 fragments (single
 * ring bonds will not be removed, as that will only give 1 fragment).  If
 * maxBondSplits is 2, all 1 and 2 bond splits are returned, and likewise for
 * maxBondSplits=3.  Values higher than that will give reset to 3 as we only
 * anticipate a maximum of 4 component libraries, and values lower than 1 will
 * be reset to 1.
 * @return : vector of unique pointers to ROMol - the hits
 */
std::vector<std::unique_ptr<ROMol>> SSSearch(const ROMol &query,
                                             unsigned int maxBondSplits,
                                             const std::string &libName);
}  // namespace HyperspaceSSSearch
}  // namespace RDKit

#endif  // RDKIT_HYPERSPACESUBSTRUCTURESEARCH_H
