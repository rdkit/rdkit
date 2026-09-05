//
//
//  Copyright (C) 2020 Schrödinger, LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once

#include <RDGeneral/BoostStartInclude.h>
#include <boost/dynamic_bitset.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <RDGeneral/export.h>

namespace RDKit {

class ROMol;

namespace CIPLabeler_detail {
RDKIT_CIPLABELER_EXPORT bool decrementRemainingCallCountAndCheck();
}

namespace CIPLabeler {

/*
  Some very symmetrical mols can cause pseudo-infinite processing
  (e.g. dodecahedrane). To avoid this, the caller can limit the total number
  of recursive comparisons performed by assignCIPLabels. The limit applies
  across both the preliminary and full labeling passes. If the limit is
  exhausted, the following error is thrown.
*/

class RDKIT_CIPLABELER_EXPORT MaxIterationsExceeded
    : public std::runtime_error {
 public:
  explicit MaxIterationsExceeded()
      : std::runtime_error(
            "Max Iterations Exceeded in CIP label calculation"){};
};

/**
 * Calculate Stereochemical labels based on an accurate implementation
 * of the CIP rules.
 *
 * This is a C++ port of https://github.com/SiMolecule/centres, which was
 * originally written by John Mayfield in Java. The original algorithm was
 * described in:
 *
 * Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, A.,
 * Redkin, D. Algorithmic Analysis of Cahn--Ingold--Prelog Rules of
 * Stereochemistry: Proposals for Revised Rules and a Guide for Machine
 * Implementation. J. Chem. Inf. Model. 2018, 58, 1755-1765.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \note only atoms with chiral tags and double bonds with proper
 *          bond directions will be labelled.
 *   \note Labels will be stored under the common_properties::_CIPCode
 *          property of the relevant atoms/bonds.
 *   \param maxRecursiveIterations - maximum total number of recursive
 *          comparisons across all labeling passes; zero means no
 *          caller-specified limit.
 */
RDKIT_CIPLABELER_EXPORT void assignCIPLabels(
    ROMol &mol, unsigned int maxRecursiveIterations = 0);

/**
 * Overload that allows selecting which atoms and/or bonds will be labeled.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \param atoms - bitset with the atom indexes to be labeled.
 *      A nonempty bitset must have mol.getNumAtoms() bits. Unselected atom
 *      configurations may still be used as auxiliary stereochemical
 *      dependencies, but their primary labels are not written.
 *
 *   \param bonds - bitset with the bond indexes to be labeled.
 *      A nonempty bitset must have mol.getNumBonds() bits. Unselected bond
 *      configurations may still be used as auxiliary stereochemical
 *      dependencies, but their primary labels are not written.
 *
 *   \param maxRecursiveIterations - maximum total number of recursive
 *      comparisons across the preliminary and full labeling passes. A value
 *      of zero means no caller-specified limit. A value of 1,250,000 takes
 *      about 1 second. Most structures require fewer than 10,000 comparisons.
 *      A peptide with MW~3000 took about 100 comparisons, and a 20,000 MW
 *      protein took about 600 comparisons.
 *
 *   \note Partial labeling does not set common_properties::_CIPComputed,
 *      which denotes a complete molecule-wide calculation.
 *
 */
RDKIT_CIPLABELER_EXPORT void assignCIPLabels(
    ROMol &mol, const boost::dynamic_bitset<> &atoms,
    const boost::dynamic_bitset<> &bonds,
    unsigned int maxRecursiveIterations = 0);

}  // namespace CIPLabeler
}  // namespace RDKit
