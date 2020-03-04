//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/export.h>
#pragma once
#include "StructChecker.h"
#include "Utilites.h"

namespace RDKit {
namespace StructureCheck {

// ??? stereo parity
static const int ODD = 1;
static const int EVEN = 2;
static const int UNMARKED = 3;

static const int ALLENE_PARITY = -2;
static const int ILLEGAL_REPRESENTATION = -1;
static const int UNDEFINED_PARITY = 0;
static const int ODD_PARITY = 1;
static const int EVEN_PARITY = 2;
static inline int INVERT_PARITY(int p) { return ((p) == 0 ? (0) : (3 - (p))); }

// bond color:
static const int CIS = 1;
static const int TRANS = 2;

// return codes for DubiousStereochemistry()
static const int EITHER_BOND_FOUND = 1;
static const int STEREO_BOND_AT_NON_STEREO_ATOM = 2;
static const int ZEROED_Z_COORDINATES = 4;
static const int CONVERTED_TO_2D = 8;
/* DubiousStereochemistry:
 * Checks if there is some ill-defined stereochemistry in the
 * molecule *mp. The function returns a bit set integer which defines
 * the problems encountered.
 */
RDKIT_STRUCTCHECKER_EXPORT int DubiousStereochemistry(RWMol &mol);

/* FixDubious3DMolecule:
 * Checks if the structure has 3D coordinates and/or flat sp3-carbons with
 * stereo-bonds and
 * converts the designation to 2D, clearing away any Z-component of the
 * coordinates.
 * Real 3D structures without stereo designations go through untouched.
 */
RDKIT_STRUCTCHECKER_EXPORT int FixDubious3DMolecule(RWMol &mol);

// Removes ill-defined stereodescriptors.
RDKIT_STRUCTCHECKER_EXPORT void RemoveDubiousStereochemistry(RWMol &mol);

/*
 * Checks if all potential stereocenters are either completely undefined
 * or attributed with hashes and wedges according to MDL rules.
 */
RDKIT_STRUCTCHECKER_EXPORT bool CheckStereo(const ROMol &mol);

/*
 * Checks if any two atoms in *mp come closer than 10% of the
 * average bond length or if an atom is too close the line
 * between two bonded atoms.
 */
RDKIT_STRUCTCHECKER_EXPORT bool AtomClash(RWMol &mol, double clash_limit);

/*
 * Computes the stereo parity of atom number iatom in *mp relative
 * to its numbering. The immediate neighbours are defined by *nbp
 * to speed up processing.
 */
RDKIT_STRUCTCHECKER_EXPORT int AtomParity(const ROMol &mol, unsigned iatom,
                                          const Neighbourhood &nbp);

/*
 * Sets the color field of the defined double bonds in *mp to CIS,
 * TRANS, or NONE depending on the ligands with the lowest numbering[].
 * It returns the number of defined double bonds found.
 */
RDKIT_STRUCTCHECKER_EXPORT int CisTransPerception(
    const ROMol &mol, const std::vector<RDGeom::Point3D> &points,
    const std::vector<unsigned> &numbering, std::vector<unsigned> &bondColor);
}  // namespace StructureCheck
}  // namespace RDKit
