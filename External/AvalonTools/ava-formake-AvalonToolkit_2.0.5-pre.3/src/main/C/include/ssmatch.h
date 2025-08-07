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
/*      File:           ssmatch.h                                       */
/*                                                                      */
/*      purpose:        implements data structures for the substructure */
/*                      matching of reaccs molecules                    */
/*                                                                      */
//      history:        creation:       23-Feb-1990                     */
//                                                                      */
/************************************************************************/
 
#ifndef SSMATCH_H
#define SSMATCH_H

#define NOT_ASSIGNED            (-1)
#define SINGLE_MATCH            TRUE
#define MULTIPLE_MATCHES        FALSE

// bond type flags used in aromaticity matching
#define BTF_SINGLE	0x01
#define BTF_DOUBLE	0x02
#define BTF_TRIPLE	0x04
#define BTF_AROMATIC	0x08
#define BTF_DYAROMATIC	0x10

typedef struct SSMATCH_T ssmatch_t;

struct SSMATCH_T
   {
      int match_atoms[MAXATOMS];
      int n_match;
      ssmatch_t *next;
   };

extern ssmatch_t *NewSSMatch(void);
/*
 * Creates an empty new substructure mapping structure. It either allocates
 * it or takes it from the free list.
 */

extern void FreeSSMatch(ssmatch_t *matchp);
/*
 * Places the argument into the list of free substructure match records.
 */

extern void FreeSSMatchHeap(void);
/*
 * Frees all the ssmatch_t structures in the list free_ss_matches.
 */
 

extern 
ssmatch_t *SSMatch(struct reaccs_molecule_t *mp,
                   struct reaccs_molecule_t *ssp,
                   int                       single_match);
/*
 * Returns the list of matchings of the substructure *ssp in structure *mp.
 * 'single_match' is TRUE if only the existance of a substructure is to be
 * checked.
 *
 * Actually a wrapper aroung SSMatchExt for compatibility.
 */
 
extern 
ssmatch_t *SSMatchExt(struct reaccs_molecule_t *mp,
                      unsigned long aa_mask[], // perceived set if augnmented atom classes for shortcut matching or NULL
                      struct reaccs_molecule_t *ssp,
                      unsigned long saa_mask[], // set of rquired augmented atom classes for the SS atom or NULL
                      int                       single_match);
/*
 * Returns the list of matchings of the substructure *ssp in structure *mp.
 * 'single_match' is TRUE if only the existance of a substructure is to be
 * checked.
 */

extern
ssmatch_t *FGMatch(struct reaccs_molecule_t *mp,
                   struct reaccs_molecule_t *ssp,
                   int                       nlimit);
/*
 * Returns the list of matchings of the substructure *ssp in structure *mp.
 * 'nlimit' atoms form the match, the remaining atoms define the allowed
 * neighbourhood of the matching atoms (optional neighbours).
 */

extern
ssmatch_t *RSSMatch(struct reaccs_molecule_t *mp,struct reaccs_molecule_t *ssp);
/*
 * Returns the list of matchings of the substructure *ssp in structure *mp.
 * Atoms with a zero mapping field must not have any connections other than
 * than to also mapped atoms.
 */

extern char *fp_prefixes[];

#define USE_RING_PATTERN      0x000001
#define USE_RING_PATH         0x000002
#define USE_ATOM_SYMBOL_PATH  0x000004
#define USE_ATOM_CLASS_PATH   0x000008
#define USE_ATOM_COUNT        0x000010
#define USE_AUGMENTED_ATOM    0x000020
#define USE_HCOUNT_PATH       0x000040
#define USE_HCOUNT_CLASS_PATH 0x000080
#define USE_HCOUNT_PAIR       0x000100
#define USE_BOND_PATH         0x000200
#define USE_AUGMENTED_BOND    0x000400
#define USE_RING_SIZE_COUNTS  0x000800
#define USE_DEGREE_PATH       0x001000
#define USE_CLASS_SPIDERS     0x002000
#define USE_FEATURE_PAIRS     0x004000
#define USE_ALL_FEATURES      0x007FFF
#define USE_FG_FEATURES       (USE_ATOM_SYMBOL_PATH|\
                               USE_HCOUNT_PATH|\
                               USE_HCOUNT_PAIR|\
                               USE_AUGMENTED_BOND|\
                               USE_AUGMENTED_ATOM)

#define USE_SCAFFOLD_IDS      0x100000
#define USE_SCAFFOLD_COLORS   0x200000
#define USE_SCAFFOLD_LINKS    0x400000
#define USE_SHORTCUT_LABELS   0x800000

#define USE_NON_SSS_BITS      0xF00000

#define USE_DY_AROMATICITY    0x0001
#define ACCUMULATE_BITS       0x0002

/*
 * Sets the bits in fingerprint[0..nbytes-1] that correspond to the paths
 * with up to maxatoms nodes. as_query must be set to TRUE if the fingerprint
 * is to be computed for a query.
 */
extern
int SetFingerprintBits(struct reaccs_molecule_t *mp,
                       char *fingerprint, int nbytes,
                       int which_bits,
                       int as_query,
		               int fpflags);

/*
 * Counts the bit multiplicity in fp_counts[0..ncounts-1] that correspond to the paths
 * with up to some defined size. which_bits is a set of flags that defines
 * which algorithms are triggered. as_query must be set to TRUE if the
 * fingerprint is to be computed for a query, which is mostly for hydrogen
 * counts.
 *
 * exclude_atom is the number of the atom that must not be touched by each
 * fragment. exclude_atom == 0 means don't exclude any atoms.
 */
extern
int SetFingerprintCountsWithFocus(struct reaccs_molecule_t *mp,
                                int *fp_counts, int ncounts,
                                int which_bits,
                                int as_query,
                                int fpflags,
                                int exclude_atom);

/**
 * Removes entries from matches when a previous match hits the same
 * set of atoms. This is a short-cut for symmetry equivalent matches,
 * which is not exact!
 */
ssmatch_t *PruneDuplicateAtomSets(ssmatch_t *matches);

/*
 * Counts the number of substructure matches listed in smp.
 */
int CountMatches(ssmatch_t *smp);

/*
 * Prints a substructure match onto file *fp.
 */
void PrintMatch(FILE *fp, ssmatch_t *match);

/*
 * Sets the bits in fingerprint[0..nbytes-1] that correspond to the paths
 * with up to some defined size. which_bits is a set of flags that defines
 * which algorithms are triggered. as_query must be set to TRUE if the
 * fingerprint is to be computed for a query, which is mostly for hydrogen
 * counts.
 *
 * Collects only those bits that require focus_atom in the molecule.
 */
int AccumulateFingerprintBitsWithFocus(struct reaccs_molecule_t *mp,
                                       char *fingerprint, int nbytes,
                                       int which_bits,
                                       int as_query,
		                               int fpflags,
                                       int focus_atom);
/*
 * Counts the number of bits set in fingerprint[0..nbytes-1].
 */
int CountBits(char *fingerprint, int nbytes);

#endif
