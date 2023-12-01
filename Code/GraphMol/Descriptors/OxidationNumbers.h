//
//  Copyright (c) 2023, David Cosgrove, CozChemIx Limited
//  All rights reserved.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Calculate the oxidation numbers (states) of the atoms in a molecule.
// Based on the code at
// https://github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/src/linchemin/cheminfo/functions.py#L544
// and therefore also subject to the MIT licence as detailed at
// https://github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/LICENSE

#include <RDGeneral/export.h>

#ifndef RD_OXIDATION_NUMBERS_MAR2023
#define RD_OXIDATION_NUMBERS_MAR2023

namespace RDKit {
class Atom;
class ROMol;
namespace Descriptors {

/*!
 * Calculates the oxidation numbers (states) of the atoms in a molecule
 * and stores them in the property _OxidationNumber on the atoms.  Uses Pauling
 * electronegativies.
 * This is experimental code, still under development.
 *
 * @param mol the molecule of interest
 */
RDKIT_DESCRIPTORS_EXPORT void calcOxidationNumbers(const ROMol &mol);

}  // end of namespace Descriptors
}  // end of namespace RDKit
#endif
