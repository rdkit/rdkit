//
//  Copyright (C) 2003-2009 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#ifndef _RD_SUBGRAPHUTILS_H_
#define _RD_SUBGRAPHUTILS_H_

#include "Subgraphs.h"
#include <boost/tuple/tuple.hpp>
#include <cstdint>

namespace RDKit {
class ROMol;

namespace Subgraphs {
//! used to return path discriminators (three unsigned ints):
typedef boost::tuples::tuple<std::uint32_t, std::uint32_t, std::uint32_t>
    DiscrimTuple;

RDKIT_SUBGRAPHS_EXPORT DiscrimTuple calcPathDiscriminators(
    const ROMol &mol, const PATH_TYPE &path, bool useBO = true,
    std::vector<std::uint32_t> *extraInvars = nullptr);
RDKIT_SUBGRAPHS_EXPORT PATH_LIST uniquifyPaths(const ROMol &mol,
                                               const PATH_LIST &allPathsb,
                                               bool useBO = true);

// Return the list of bond that connect a list of atoms
// ASSUMPTION: the atoms specified in the list are connected
RDKIT_SUBGRAPHS_EXPORT PATH_TYPE bondListFromAtomList(const ROMol &mol,
                                                      const PATH_TYPE &atomIds);

// create a new molecule object from a part of molecule "mol". The part of
// of the molecule is specified as a list of bonds in "path".
// the optional argument "useQuery" will set all the bond and atoms in the
// the new molecule to "QueryAtoms" and "QueryBonds" instead of regular Atoms
// and Bonds
//  atomIdxMap provides a mapping between the atomsIds in mol to the atomIds in
// the newly created sub-molecule (the molecule that is returned)
RDKIT_SUBGRAPHS_EXPORT ROMol *pathToSubmol(const ROMol &mol,
                                           const PATH_TYPE &path, bool useQuery,
                                           std::map<int, int> &atomIdxMap);
RDKIT_SUBGRAPHS_EXPORT ROMol *pathToSubmol(const ROMol &mol,
                                           const PATH_TYPE &path,
                                           bool useQuery = false);
}  // end of namespace Subgraphs
}  // namespace RDKit

#endif
