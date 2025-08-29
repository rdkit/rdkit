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
/*    File:           shortcut.h                                        */
/*                                                                      */
/*    Purpose:        Defines the data structures and constants for     */
/*                    shortcut handling.                                */
/*                                                                      */
//    History:        10-May-2012     First version                     */
//                                                                      */
/************************************************************************/

#include "pattern.h"
#include "reaccs.h"
#include "ssmatch.h"

#define TYPE_MASK               0x00FF
#define AMINO_ACIDS             0x0001
#define PROTECTING_GROUP        0x0002

#define LEVEL_MASK              0xFF00
#define STANDARD_SHORTCUT       0x0100
#define EXTENDED_SHORTCUT       0x0200
#define NON_STANDARD_SHORTCUT   0x0400
#define N_METHYL_SHORTCUT       0x0800
#define CATCH_ALL               0x8000

struct shortcut_input_t
{
    char *smiles;              // R-group coded fragment SMILES. Gets converted into MOL query with attachment allowed only at R-groups
    char *single_letter_code;  // shortest (single letter) code
    char *shortcut;            // longer code
    int shortcut_class;        // flags indicating what kind of shortcut this is.
    char *r1class;             // Augmented Atom class the R1 connection atom must be in.
    char *r2class;             // Augmented Atom class the R2 connection atom must be in.
    char *r3class;             // Augmented Atom class the R3 connection atom must be in.
};

struct shortcut_t
{
    struct reaccs_molecule_t *query;
    int nattachments;         // number of defined attachments
    int attachment_idx1;      // index of first attachment atom
    int attachment_idx2;      // index of second attachment atom or -1 if none
    int attachment_idx3;      // index of third attachment atom or -1 if none
    char *short_code;
    char *long_code;
    int shortcut_class;       // flags indicating what kind of shortcut this is.
    int rank;
    unsigned long *aa_mask;   // required bit mask of augmented atom classes for each shortcut atom
    struct shortcut_t *next;
};

struct aa_class_t
{
    char name[20];                // Name of aa_class. Note that there can be multiple augmented atom patterns with the same name and bit_pattern.
    unsigned long bit_mask;       // Bit pattern for quickly checking preperceived atoms.
    augmented_atom_t aa_pattern;  // augmented atom to match
    int no_ring_atom;             // indicates if this atom must not be a ring atom
    struct reaccs_molecule_t *smiles_query;     // If not NULL, this is an additional query that atom number 1 must match
    struct aa_class_t *next;      // pointer to next augmented atom in list
};

struct shortcut_match_t
{
    ssmatch_t *match;
    int match_index1;      // index of atom to which the first attachment atom matches
    int match_index2;      // index of atom to which the second attachment atom matches or -1 if no second attachment
    int match_index3;      // index of atom to which the third attachment atom matches or -1 if no third attachment
    int color;
    int thread_color;      // color of the left-most node in the thread
    int left_color;        // color of atom connected to match_index1 if any
    int right_color;       // color of atom connected to match_index2 if any
    struct shortcut_t *shortcut;
    int mirror_match;      // TRUE if only the mirror image of shortcut actually matched
};

extern void InitShortcuts();
/*
 * Perceive standard shortcut definitions and initialize derived data structures.
 */

extern struct shortcut_t *MakeShortcut(char * smiles, char *shortCode, char *logCode,
                                       int shortcut_class,
                                       char *r1patterns, char *r2patterns, char *r3patterns,
                                       struct shortcut_t *old_list,
                                       struct aa_class_t *aa_classes,
                                       int rank);
/*
 * Prepends a newly created shorcut descriptor to old_list.
 */

struct reaccs_molecule_t *ApplyShortcuts(struct reaccs_molecule_t *mp, int flags);
/*
 * Apply the currently registered shortcuts to mp and return the updated molecules.
 */

