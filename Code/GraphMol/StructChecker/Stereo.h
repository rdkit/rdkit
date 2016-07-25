//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#pragma once
#include "StructChecker.h"

namespace RDKit {
 namespace StructureCheck {

// ??? stereo parity
static const int ODD      = 1;
static const int EVEN     = 2;
static const int UNMARKED = 3;

// return codes for DubiousStereochemistry()
static const int EITHER_BOND_FOUND              = 1;
static const int STEREO_BOND_AT_NON_STEREO_ATOM = 2;
static const int ZEROED_Z_COORDINATES           = 4;
static const int CONVERTED_TO_2D                = 8;
    /* DubiousStereochemistry:
    * Checks if there is some ill-defined stereochemistry in the
    * molecule *mp. The function returns a bit set integer which defines
    * the problems encountered.
    */
int DubiousStereochemistry(RWMol &mol);

     /* FixDubious3DMolecule:
     * Checks if the structure has 3D coordinates and/or flat sp3-carbons with stereo-bonds and
     * converts the designation to 2D, clearing away any Z-component of the coordinates.
     * Real 3D structures without stereo designations go through untouched.
     */
int FixDubious3DMolecule  (RWMol &mol);

    // Removes ill-defined stereodescriptors.
void RemoveDubiousStereochemistry(RWMol& mol);

    /*
    * Checks if all potential stereocenters are either completely undefined
    * or attributed with hashes and wedges according to MDL rules.
    */
bool CheckStereo(const ROMol& mol);

    /*
    * Checks if any two atoms in *mp come closer than 10% of the
    * average bond length or if an atom is too close the line
    * between two bonded atoms.
    */
bool AtomClash(RWMol &mol, double clash_limit);

 }
}

