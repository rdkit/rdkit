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

/**********************************************************************/
/*                                                                    */
/*    File:           canonizer.h                                     */
/*                                                                    */
/*    Purpose:        This file declares the functions and constants  */
/*                    to generate canonical SMILES of various types   */
/*                    from an MDL molecule data structure.            */
/*                                                                    */
/*    History:        07-Feb-2002     Start of development.           */
/*                                                                    */
/**********************************************************************/

#define DB_STEREO       0x01    /* use double-bond stereochemistry        */
#define CENTER_STEREO   0x02    /* use tetrahedral center stereochemistry */
#define SCRAMBLE        0x04    /* scramble CT before processing          */
#define COMPONENT_SET   0x08    /* list every component only once         */
#define MAIN_PART       0x10    /* only keep non-trivial pieces           */
#define FORCE_OCTET     0x20    /* Make it comply with Octet rule         */
#define REMOVE_ISOTOPE  0x40    /* Remove isotopic labelling              */
#define REMOVE_MAPPING  0x80    /* clear away labelling by mapping field  */

/*
 * Return a pointer to a newly allocated character string that contains
 * the canonicalized form of insmiles[].
 *
 * It will return NULL if anything goes wrong.
 *
 * flags is used to modify processing as to which stereo features shall
 * be used and which components to retain.
 */
extern
char *CanSmiles(char *insmiles, int flags);

/**
 * This function removes the SMILES of any remaining standard counter ions
 * or solvent molecules from smiles[] and returnes the modified string.
 * The table of 'salt_array' is provides as a parameter and assumed to
 * be already canonicalized.
 * Note: The function modifies smiles in situ.
 */
extern
char *RemoveStandardComponents(char *smiles, char **salt_array);

#define MAXLINE 20000

/**
 * Needed to reproducibly scramble for testing.
 */
void SetSeed(int seed_value);
