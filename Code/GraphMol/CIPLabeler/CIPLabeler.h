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
#include <RDGeneral/export.h>

namespace RDKit {

class ROMol;

namespace CIPLabeler {

/**
 * Calculate Stereochemical labels based on an accurate implementation
 * of the CIP rules, described in:
 *
 *  Hanson, R. M., Musacchio, S., Mayfield, J. W., Vainio, M. J., Yerin, A.,
 *  Redkin, D. Algorithmic Analysis of Cahn−Ingold−Prelog Rules of
 *  Stereochemistry: Proposals for Revised Rules and a Guide for Machine
 *  Implementation. J. Chem. Inf. Model. 2018, 58, 1755-1765.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \note only atoms with chiral tags and double bonds with proper
 *          bond directions will be labelled.
 *   \note Labels will be stored under the common_properties::_CIPCode
 *          property of the relevant atoms/bonds.
 */
RDKIT_CIPLABELER_EXPORT void assignCIPLabels(ROMol &mol);
}
}