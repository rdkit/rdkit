//
//  Copyright (c) 2010, Novartis Institutes for BioMedical Research Inc.
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

/************************************************************************/
/*                                                                      */
/*  File:     stereo.h                                                  */
/*                                                                      */
/*  Purpose:  Declares the funtions providing stereochemistry           */
/*            checking.                                                 */
/*                                                                      */
/*  Revisions:        1.0.0   11-Sep-92       Creation of module        */
/*                                                                      */
/************************************************************************/

#include "ssmatch.h"

#define CIS    1
#define TRANS  2

#define ALLENE_PARITY          (-2)
#define ILLEGAL_REPRESENTATION (-1)
#define UNDEFINED_PARITY       0
#define ODD_PARITY             1
#define EVEN_PARITY            2
#define INVERT_PARITY(p) ((p)==0 ? (0) : (3-(p)))

/* return codes for DubiousStereochemistry() */
#define EITHER_BOND_FOUND              1
#define STEREO_BOND_AT_NON_STEREO_ATOM 2
#define ZEROED_Z_COORDINATES           4
#define CONVERTED_TO_2D                8

struct npoint_t
   {
      double x, y, z;
      int number;
   };

extern
double Volume(struct npoint_t tetra[4]);
/*
 * Computes the signed volume of the tetrahedron defined by
 * the four points in tetra[].
 */

extern
int FixDubious3DMolecule(struct reaccs_molecule_t *mp);
/*
 * Checks if the structure has 3D coordinates and/or flat sp3-carbons with stereo-bonds and
 * converts the designation to 2D, clearing away any Z-component of the coordinates.
 * Real 3D structures without stereo designations go through untouched.
 */

extern
int CheckStereo(struct reaccs_molecule_t *mp);
/*
 * Checks if all potential stereocenters are either completely undefined
 * or attributed with hashes and wedges according to MDL rules.
 */

extern
int DubiousStereochemistry(struct reaccs_molecule_t *mp);
/*
 * Checks if there is some ill-defined stereochemistry in the
 * molecule *mp. Currently the function looks for EITHER bonds only.
 */

extern
void RemoveDubiousStereochemistry(struct reaccs_molecule_t *mp);
/*
 * Removes ill-defined stereodescriptors. Currently, only EITHER
 * bonds are removed.
 */

extern
void SetCollisionLimit(int percent);
/*
 * Sets the limit for atom/atom and atom/bond collision in percent
 * of the average bond length.
 */

extern
int AtomClash(struct reaccs_molecule_t *mp);
/*
 * Checks if any two atoms in *mp come closer than 10% of the
 * average bond length or if an atom is too close the line
 * between two bonded atoms.
 */

extern
int AtomParity(struct reaccs_molecule_t *mp,
               int iatom,
           neighbourhood_t *nbp);
/*
 * Computes the stereo parity of atom number iatom in *mp relative
 * to its original numbering. The immediate neighbours are defined by *nbp
 * to speed up processing.
 */

extern
int CisTransPerception(struct reaccs_molecule_t *mp,
                      int                       numbering[]);
/*
 * Sets the color field of the defined double bonds in *mp to CIS,
 * TRANS, or NONE depending on numbering[]. It returns the number
 * of defined double bonds found.
 */

extern
int NoParityDefined(struct reaccs_molecule_t *mp);
/*
 * Returns TRUE if the molecule *mp has no defined stereocenter.
 */

int AllCentersRefined(struct reaccs_molecule_t *mp,
               int                       numbering[]);
/*
 * Checks if all defined stereocenters in *mp are unambiguously
 * numbered by numbering, i.e. all four ligands get different
 * numbers. It returns TRUE if this is the case and FALSE otherwise.
 */

int IsStereoMatch(struct reaccs_molecule_t *mp, neighbourhood_t *nbp,
                  struct reaccs_molecule_t *qp, neighbourhood_t *qnbp,
                  ssmatch_t *matches, int invert_stereo);
/**
 * Checks if the match matches into mp obeys the center stereochemistry of qp.
 * Uses the mirror image of *qp if invert_stereo == TRUE.
 */
