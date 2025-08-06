//
//  Copyright (C) 2025 David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

// Does a 2 molecule, single fragment MCSS using a vanilla Bron-Kerbosch
// algorithm with a constrained correspondence graph.  A pair of atoms
// in mol1 will only match a pair of atoms in mol2 if they are the same
// minimum number of bonds apart.

#ifndef TWOMOLMCSS_H
#define TWOMOLMCSS_H

#include <string>
#include <vector>

#include <RDGeneral/export.h>

namespace RDKit {
class ROMol;
// Make an approximate MCSS from the 2 molecules where 2 atoms in mol1
// match 2 atoms in mol2 only if they are the same topological distance
// apart.  This may be smaller than a full MCSS and will only be a
// single fragment.
RDKIT_FMCS_EXPORT void TwoMolMCSS(
    const ROMol &mol1, const ROMol &mol2,
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
        &maxCliques);

// Take a clique, where the first element in each pair refers to
// atoms in the given molecule and make a SMARTS string from it
// by deleting all atoms from mol that aren't in the clique.
RDKIT_FMCS_EXPORT std::string makeSMARTSFromMCSS(
    const ROMol &mol,
    const std::vector<std::pair<unsigned int, unsigned int>> &clique);
}  // namespace RDKit
#endif  // TWOMOLMCSS_H
