//
//  Copyright (C) 2015 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include <vector>
#include <string>
#include <stdexcept>
#include "../RDKitBase.h"

namespace RDKit {

namespace MMPA {
//! fragments a Molecule for processing with the Matched Molecular Pairs
//!  MMPA algorithm (Hussain et al)
/*!
  \param mol           Molecule to fragment
  \param result        Vector of Core and Sidechain results from the various
 cuts
  \param maxCuts        Maximum number of times to cut the molecule to generate
                        fragments.  A max cut of 3 will fragment with 1,2 and 3
                        cuts.
  \param maxCutBonds  Set the bond limit for determining which molecules
                        to analyze.  If a molecule has more than
                        this number of cutabble bonds, ignore.

 \return true if the molecule was fragmented, false otherwise.
*/

RDKIT_MMPA_EXPORT bool fragmentMol(
    const ROMol& mol, std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>>& result,
    unsigned int maxCuts = 3, unsigned int maxCutBonds = 20,
    const std::string& pattern = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]");

//! fragments a Molecule for processing with the Matched Molecular Pairs
//!  MMPA algorithm (Hussain et al)
/*!
  \param mol           Molecule to fragment
  \param result        Vector of Core and Sidechain results from the various
 cuts
  \param minCuts        Minimum number of times to cut the molecule to generate
                        fragments.
  \param maxCuts        Maximum number of times to cut the molecule to generate
                        fragments.
  \param maxCutBonds  Set the bond limit for determining which molecules
                        to analyze.  If a molecule has more than
                        this number of cutabble bonds, ignore.

 \return true if the molecule was fragmented, false otherwise.
*/
RDKIT_MMPA_EXPORT bool fragmentMol(
    const ROMol& mol, std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>>& result,
    unsigned int minCuts, unsigned int maxCuts, unsigned int maxCutBonds,
    const std::string& pattern = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]");

//! fragments a Molecule for processing with the Matched Molecular Pairs
//!  MMPA algorithm (Hussain et al)
/*!
  \param mol           Molecule to fragment
  \param result        Vector of Core and Sidechain results from the various
                        cuts
  \param bondsToCut      Vector of bond indices to use as cut points
  \param minCuts        Minimum number of times to cut the molecule to generate
                        fragments.
  \param maxCuts        Maximum number of times to cut the molecule to generate
                        fragments.
 \return true if the molecule was fragmented, false otherwise.
*/
RDKIT_MMPA_EXPORT bool fragmentMol(
    const ROMol& mol, std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>>& result,
    const std::vector<unsigned int>& bondsToCut, unsigned int minCuts = 1,
    unsigned int maxCuts = 3);

}  // namespace MMPA
}  // namespace RDKit
