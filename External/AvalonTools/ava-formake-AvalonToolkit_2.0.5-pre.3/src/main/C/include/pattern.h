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
/*    File:           pattern.h                                         */
/*                                                                      */
/*    Purpose:        This file declares the functions to translate     */
/*                    one augmented atom pattern to another. This       */
/*                    service is mainly used for structure              */
/*                    standardization e.g. by struchk.exe.              */
/*                                                                      */
/************************************************************************/

#ifndef _PATTERN_H_
#define _PATTERN_H_      1

#include "utilities.h"

typedef struct LIGAND_T
   {
      char *atom_symbol;
      short charge;
      short radical;
      short bond_type;
      short s_count;    /* substitution count, 0 means "don't care" */
   } ligand_t;

typedef struct AUGMENTED_ATOM_T
   {
      short n_ligands;
      short charge;
      short radical;
      char *atom_symbol;
      char *short_name;
      short topography; /* 1 means ring, 2 means chain, 0 (default) means "don't care" */
      ligand_t ligands[MAXNEIGHBOURS];
   } augmented_atom_t;

typedef augmented_atom_t aa_pair[2];

extern
int AAMatch(struct reaccs_molecule_t *mp,
            unsigned int i,
            unsigned int *match,
            augmented_atom_t *aap,
            int atom_ring_status[],
            neighbourhood_t *nbp); 

extern int StringToAugmentedAtom(char *string, augmented_atom_t *aap);
/*
 * Reads string[] and constructs the corresponding augmented atom *aap.
 * Returns TRUE if scan was successful and FALSE otherwise.
 * The syntax of a naugmented atom string is as follows:
 *
 *   <AAString> ::= <Atom Symbol> [<Charge>]
 *               {'(' <Bond Symbol> <Atom Symbol> [<Charge>] ')'}.
 * Bond symbols and charge descriptors are defined in the symbol tables
 * 'bond_to_string' and 'charge_to_string'.
 */

extern char *AAToString(char buffer[],
                        struct reaccs_molecule_t *mp,
                        int i,
                        neighbourhood_t *np);
/*
 * Prints a string representation of the neighbourhood of atom
 * i in molecule *mp to the string buffer.
 * It uses the neighbourhood information stored in *np.
 */

extern int AAPrint(FILE *fp, struct reaccs_molecule_t *mp,
                   int i, neighbourhood_t *np);
/*
 * Prints a string representation of the neighbourhood of atom
 * i in molecule *mp to the file *fp.
 * It uses the neighbourhood information stored in *np.
 */

extern char aa_check_version[];        /* is filled by ReadAugmentedAtoms */

extern augmented_atom_t *ReadAugmentedAtoms(FILE *fp, int *natomp);
/*
 * Reads file fp and constructs an array of augmented atom
 * descriptions.
 * The function expects the number of augmented atom description
 * on the first line and the strings corresponding to the data structures
 * on the *natomp following lines.
 * Only the portion of the strings between the two '"' characters is used,
 * other characters are regarded as comments.
 */

extern char aa_trans_version[];        /* is filled by ReadAAPairs */

extern aa_pair *ReadAAPairs(FILE *fp, int *ntransp);
/*
 * Reads file fp and constructs an array of pairs of augmented atom
 * descriptions. The first one being the search pattern and the second
 * one the target pattern.
 * The function expects the number of augmented atom description pairs
 * on the first line and the strings corresponding to the data structures
 * on the *ntransp following lines.
 * Only the portion of the strings between the two '"' characters is used,
 * other characters are regarded as comments.
 */

extern
struct reaccs_molecule_t * AlternativeAATransformation(
                            struct reaccs_molecule_t *mp,
                            int                       next_atom,
                            aa_pair                  *tfm_table,
                            int                       ntfm,
                            struct reaccs_molecule_t *old_list,
                            int                       use_old);
/*
 * Searches the atoms im *mp starting from next_atom for an atom for which there
 * is at least one transformation in tfm_table[0..ntfm]. It then generates
 * all molecules with the different transformations for that atom and
 * recursively calls itself for all the resulting molecules with next_atom
 * set to the atom following the currently transformed one. It returns the
 * list of all transformed molecules prepended to old_list. If no
 * transformation is possible *mp follwed by old_list is returned.
 * nbp[] is the neighbouring information for *mp used to speed-up access.
 * If "use_old" is set, then, in addition to the transformed molecule, the
 * original one is also kept.
 */

extern int ConvertMoleculeToString(struct reaccs_molecule_t *mp,
                                   char                     *string,
                                   int                       nbuf);
/*
 * Computes a string representation of molecule *mp. The result is
 * put into string[0..nbuf-1]. The function returns TRUE if it
 * successfully translated the structure and FALSE otherwise.
 */

extern
int TransformAugmentedAtoms(struct reaccs_molecule_t *mp,
                       augmented_atom_t table[][2],
                            int nentries);
/*
 * Scans the molecule *mp for occurrences of the augmented atoms in
 * table[0..nentries-1][0] and replaces the matches with the corresponding
 * atom in table[0..nentries-1][1].
 */

#endif
