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
 * of the CIP rules.
 *
 *   \param mol - the molecule to be labelled.
 *
 *   \note only atoms with chiral tags and double bonds with proper
 *         bond directions will be labelled.
 */
RDKIT_CIPLABELER_EXPORT void assignCIPLabels(ROMol &mol);
}
}