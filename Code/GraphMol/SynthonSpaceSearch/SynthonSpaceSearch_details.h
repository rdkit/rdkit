//
// Copyright (C) David Cosgrove 2024.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDKIT_SYNTHONSPACESEARCHDETAILS_H
#define RDKIT_SYNTHONSPACESEARCHDETAILS_H

#include <vector>

#include <RDGeneral/export.h>

namespace RDKit {
class ROMol;
namespace SynthonSpaceSearch {

namespace details {
// Find all combinations of M things selected from N.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<unsigned int>>
combMFromN(unsigned int m, unsigned int n);
// Find all permutations of M things selected from N.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<unsigned int>>
permMFromN(unsigned int m, unsigned int n);
// Split the molecule into fragments.  maxBondSplits gives the maximum number
// of bonds to be used in each split.  There will a vector of vectors of
// molecules, 1 inner vector for each split i.e. maxBondSplits in total, the
// first with 1 split, the 2nd with 2 etc.  Each inner vector contains the
// fragments from a split molecule.  The maxBondSplits will be constrained to
// between 1 and 4 inclusive, so if it is supplied outside that range, it will
// be altered.  Also, you can't split a molecule on 3 bonds if it only contains
// 2.
RDKIT_SYNTHONSPACESEARCH_EXPORT std::vector<std::vector<std::unique_ptr<ROMol>>>
splitMolecule(const ROMol &query, unsigned int maxBondSplits);
// Counts the number of [1*], [2*]...[4*] in the string.
RDKIT_SYNTHONSPACESEARCH_EXPORT int countConnections(const std::string &smiles);
}  // namespace details
}  // namespace SynthonSpaceSearch
}  // namespace RDKit

#endif  // RDKIT_SYNTHONSPACESUBSTRUCTURESEARCH_H
