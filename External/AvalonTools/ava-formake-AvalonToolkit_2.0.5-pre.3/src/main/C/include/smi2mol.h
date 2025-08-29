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
/*    File:           smi2mol.h                                         */
/*                                                                      */
/*    Purpose:        This file declares a function to convert a        */
/*                    SMILES string into the corresponding MDL          */
/*                    molecule data structure.                          */
/*                                                                      */
/*    History:        11-Jul-1994     Start of development.             */
/*                                                                      */
/************************************************************************/

#define DO_LAYOUT      	  	0x0001
#define SMARTS_PERCEPTION 	0x0002
#define EXPECT_SMARTS     	0x0002  // just to get a better name
#define ISOMERIC_SMILES   	0x0004
#define DROP_TRIVIAL_HYDROGENS  0x0008
#define CANONICAL_ORDER         0x0010
#define CLASSIC_AROMATICITY     0x0020
#define AROMATIC_SMILES         0x0040
#define TRUE_DB_STEREO          0x0080
#define TO_MONOMER              0x1000
#define DY_AROMATICITY          0x2000
#define USE_Z                   0x4000
#define PERCEIVE_SHORTCUTS      0x8000

extern
struct reaccs_molecule_t *SMIToMOL(const char *smiles, int flags);
/*
 * Converts the molecule described by smiles[] to the corresponding
 * MDL data structure and returns a pointer to it.
 *
 * flags is used to modify the behaviour of the conversion. Currently,
 * it is used to request a layout calculation.
 */

extern
char *MOLToSMI(struct reaccs_molecule_t *mp, int flags);
/*
 * Returns a pointer to the SMILES string representing the connectivity
 * of the molecule *mp. The memory must be deallocated by the caller.
 *
 * flags can be ISOMERIC_SMILES which retains isomeric SMILES information
 * and/or SMARTS_PERCEPTION which performs a SMARTS_PERCEPTION before
 * *mp is translated into SMILES.
 */

char *MOLToSMIExt(struct reaccs_molecule_t *mp, int iso,
                  int numbering[],
                  char **coordpp);
/*
 * Returns a pointer to the SMILES string representing the connectivity
 * of the molecule *mp. The memory must be deallocated by the caller.
 *
 * The stereochemistry of the structure is included in the SMILES if
 * iso was set to TRUE.
 *
 * The numbering[] array is used to break ties when generating the out SMILES.
 * It may be NULL which means an arbitrary order. Index origin of numbering[]
 * is 0.
 *
 * The coordinates of the atoms in the molecule are written as comma
 * separated values to a string and a pointer to this string is assigned
 * to (*coordpp).
 */

struct reaccs_reaction_t *SMIToRXN(char *smiles);
/*
 * Converts a reactions smiles to a REACCS reaction data structure.
 */
