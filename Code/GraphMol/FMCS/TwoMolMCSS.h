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
#include <GraphMol/FMCS/MatchTable.h>

namespace RDKit {
class ROMol;

// Make an MCSS from the 2 molecules using Koch's modification of
// the Bron-Kerbosch algorithm so that it only returns a single
// fragment.  There may be more than 1 result.  If uniquify is true,
// symmetrically equivalent solutions involving the same atoms will
// be reduced to a single example.  Each element of a pair in the
// output vectors is the index of an atom in mol1 that matches
// an atom in mol2.  The MCSSs must be of at least size minMCCSSize.
RDKIT_FMCS_EXPORT void TwoMolMCSS(
    const ROMol &mol1, const ROMol &mol2, unsigned int minMCSSSize,
    const FMCS::MatchTable &atomMatchTable,
    const FMCS::MatchTable &bondMatchTable, bool uniquify,
    std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
        &maxCliques);

// Take a clique, where the first element in each pair refers to
// atoms in the given molecule and make a SMARTS string from it
// by deleting all atoms from mol that aren't in the clique.
RDKIT_FMCS_EXPORT std::string makeSMARTSFromMCSS(
    const ROMol &mol1,
    const std::vector<std::pair<unsigned int, unsigned int>> &clique);
}  // namespace RDKit
#endif  // TWOMOLMCSS_H
