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
/*    File:           shortcut.c                                        */
/*                                                                      */
/*    Purpose:        Implements the functions to represent shortcut    */
/*                    rules for abbreviated structure display.          */
/*                                                                      */
//    History:        10-May-2012     First version                     */
//                                                                      */
/************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "shortcut.h"
#include "smi2mol.h"
#include "utilities.h"
#include "perceive.h"
#include "ssmatch.h"
#include "stereo.h"

struct aa_class_string_t
{
    char *class_string;         // identifier of this atom environment
    char *aadef;                // AugmentedAtom pattern for first test
    int no_ring_atom;           // TRUE if only chain atom will match
    char *smiles_pattern;       // SMILES used to define query pattern(at it's atom 1) more precisely,
                                // e.g. to force S-S bridges
};

static struct aa_class_string_t aa_class_strings[] =
{
    {"URETHANE",       "C(=O)(-O)(-N)",   FALSE, NULL},
    {"AMINO_N",        "N(-C)",           FALSE, NULL},
    {"AMINO_N",        "N(-C)(-C)",       FALSE, NULL},
    {"AMINO_N_ME",     "N(-C)(-C)(-C)",   FALSE, NULL},
    {"AMINO_PRO",      "N(-C)(-C)(-C)",   FALSE, NULL},
    {"AMINO_PRO",      "N(-C)(-C)",       FALSE, NULL},
    {"AMINO_C",        "C(=O)(-C)(-N,O)", FALSE, NULL},
    {"SS_BRIDGE",      "C(-S)(-C)",       FALSE, "CSSC"},
    {"S_ETHER",        "C(-S)(-C)",       FALSE, "CSC"},
    {"S_ETHER",        "C(-C)(-S)(-C)",   FALSE, "C(C)SC"},
    {"AA_SUBST",       "C(-C)(-C)",       TRUE,  NULL},
};

#define CAPPING_RANK 10000

static struct shortcut_input_t shortcut_input[] =
    {
        // regular aminoacids (ordered by start rang of cyclic peptide perception)
        {"C(CN[*+])([*+2])=O",	                             "G",	"Gly",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@@H]1(N(CCC1)[*+])C([*+2])=O",	             "P",	"Pro",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_PRO","AMINO_C",	""},
        {"[C@H]([C@@H](C)O)(C([*+2])=O)N[*+]",	             "T",	"Thr",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"C(CCNC(N)=N)[C@@H](C([*+2])=O)N[*+]",	             "R",	"Arg",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H](CC(N)=O)(C([*+2])=O)N[*+]",	             "N",	"Asn",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H](CC(O)=O)(C([*+2])=O)N[*+]",	             "D",	"Asp",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CS)N[*+]",                      "C",	"Cys",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(C[SeH])N[*+]",                  "U",	"Sec",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H](CCC(O)=O)(C([*+2])=O)N[*+]",	             "E",	"Glu",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H](CCC(N)=O)(C([*+2])=O)N[*+]",	             "Q",	"Gln",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"C1(N=CNC=1)C[C@@H](C([*+2])=O)N[*+]",	             "H",	"His",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"C1(NC=NC=1)C[C@@H](C([*+2])=O)N[*+]",	             "H",	"His",	AMINO_ACIDS|STANDARD_SHORTCUT,  // tautomer
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H]([C@H](CC)C)(C([*+2])=O)N[*+]",	             "I",	"Ile",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H](CCCCN)(C([*+2])=O)N[*+]",	             "K",	"Lys",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"C1C=N[C@H]([C@@H]1C)C(=O)NCCCCC(C(=O)[*+2])N[*+]", "O",	"Pyl",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H](CC(C)C)(C([*+2])=O)N[*+]",	             "L",	"Leu",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"C1(C=CC=CC=1)C[C@@H](C([*+2])=O)N[*+]",	     "F",	"Phe",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CCSC)N[*+]",	             "M",	"Met",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CO)N[*+]",	                     "S",	"Ser",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CC1C=CC(=CC=1)O)N[*+]",	     "Y",	"Tyr",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"C12=C(NC=C1C[C@@H](C([*+2])=O)N[*+])C=CC=C2",	     "W",	"Trp",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(C)N[*+]",	                     "A",	"Ala",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        {"[C@H](C(C)C)(C([*+2])=O)N[*+]",	             "V",	"Val",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},

        // regular N-methylated aminoacids
        {"C(CN(C)[*+])([*+2])=O",	                     "MeG",	"MeGly",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H]([C@@H](C)O)(C([*+2])=O)N(C)[*+]",	     "MeT",	"MeThr",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"C(CCNC(N)=N)[C@@H](C([*+2])=O)N(C)[*+]",	     "MeR",	"MeArg",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H](CC(N)=O)(C([*+2])=O)N(C)[*+]",	             "MeN",	"MeAsn",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H](CC(O)=O)(C([*+2])=O)N(C)[*+]",	             "MeD",	"MeAsp",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CS)N(C)[*+]",                   "MeC",	"MeCys",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(C[SeH])N(C)[*+]",               "MeU",	"MeSec",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H](CCC(O)=O)(C([*+2])=O)N(C)[*+]",	             "MeE",	"MeGlu",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H](CCC(N)=O)(C([*+2])=O)N(C)[*+]",	             "MeQ",	"MeGln",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"C1(N=CNC=1)C[C@@H](C([*+2])=O)N(C)[*+]",	     "MeH",	"MeHis",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"C1(NC=NC=1)C[C@@H](C([*+2])=O)N(C)[*+]",	     "MeH",	"MeHis",	AMINO_ACIDS|N_METHYL_SHORTCUT,  // tautomer
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H]([C@H](CC)C)(C([*+2])=O)N(C)[*+]",	     "MeI",	"MeIle",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H](CCCCN)(C([*+2])=O)N(C)[*+]",	             "MeK",	"MeLys",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H](CC(C)C)(C([*+2])=O)N(C)[*+]",	             "MeL",	"MeLeu",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"C1(C=CC=CC=1)C[C@@H](C([*+2])=O)N(C)[*+]",	     "MeF",	"MePhe",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CCSC)N(C)[*+]",	             "MeM",	"MeMet",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CO)N(C)[*+]",	             "MeS",	"MeSer",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(CC1C=CC(=CC=1)O)N(C)[*+]",	     "MeY",	"MeTyr",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"C12=C(NC=C1C[C@@H](C([*+2])=O)N(C)[*+])C=CC=C2",   "MeW",	"MeTrp",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@@H](C([*+2])=O)(C)N(C)[*+]",                    "MeA",	"MeAla",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},
        {"[C@H](C(C)C)(C([*+2])=O)N(C)[*+]",	             "MeV",	"MeVal",	AMINO_ACIDS|N_METHYL_SHORTCUT,
         "AMINO_N_ME",	"AMINO_C",	""},

        // http://www.americanpeptide.com/commerce/support/technical/references/aaa.jsp
        // http://www.americanpeptide.com/corp/aac_101708_p1.pdf
        // http://www.stn-international.com/uploads/tx_ptgsarelatedfiles/protseq.pdf
        // Pyroglutamic acid
        {"[C@@H]1(N(C(=O)CC1)[*+])C([*+2])=O",	             "Pyr",	"Pyr",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_PRO","AMINO_C",	""},
         // Homocysteine
        {"[C@@H](C([*+2])=O)(CCS)N[*+]",                     "Hcy",	"Hcy",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
         // Ornithine
        {"[C@@H](C([*+2])=O)(CCCN)N[*+]",                    "Orn",	"Orn",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // Norvaline, note the non-standard use of 'nor' prefix
        {"[C@H](CCC)(C([*+2])=O)N[*+]",	                     "Nva",	"Nva",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // Norleucine, note the non-standard use of 'nor' prefix
        {"[C@H](CCCC)(C([*+2])=O)N[*+]",	             "Nle",	"Nle",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // Nortyrosine, note the standard use of 'nor' for 'remove one CH2'
        {"[C@H](c1ccc(O)cc1)(C([*+2])=O)N[*+]",	             "Nty",	"Nty",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // Homoserine
        {"[C@@H](C([*+2])=O)(CCO)N[*+]",	             "Hse",	"Hse",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},

        // Extended list of shortcuts
        // http://www.thelabrat.com/protocols/aminoacidtable.shtml
        // 2-Aminobutyric acid:
        {"[C@@H](C([*+2])=O)(CC)N[*+]",	                     "Abu",	"Abu",	AMINO_ACIDS|EXTENDED_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},

        // fragments with three attachments
        // cystein with three attachment points forming cystine bridges
        {"[C@@H](C([*+2])=O)(C[*+3])N[*+]",	             "C",	"Cys",	AMINO_ACIDS|STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	"SS_BRIDGE"},
// TODO: the following two rules may be enabled if layout of DURAMYCIN gets fixed!
        // {"[C@@H](C([*+2])=O)(C(C)[*+3])N[*+]",	             "Abu",	"Abu",	AMINO_ACIDS|STANDARD_SHORTCUT,
        // "AMINO_N",	"AMINO_C",	"S_ETHER"},
        // {"[C@@H](C([*+2])=O)(C[*+3])N[*+]",	             "A",	"Ala",	AMINO_ACIDS|STANDARD_SHORTCUT,
        //   "AMINO_N",	"AMINO_C",	"S_ETHER"},
        // Lysin with 3 attachment points
        // {"[C@H](CCCCN[*+3])(C([*+2])=O)N[*+]",	             "K",	"Lys",	AMINO_ACIDS|STANDARD_SHORTCUT,
        //  "AMINO_N",	"AMINO_C",	"AMINO_N"},

        // Non-standard shortcuts
        // Dehydroalanine
        {"C(C([*+2])=O)(=C)N[*+]",	                     "Dha",	"Dha",	AMINO_ACIDS|NON_STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // cyclohexylalanine
        {"C1(CCCCC1)C[C@@H](C([*+2])=O)N[*+]",	             "Cha",	"Cha",	AMINO_ACIDS|NON_STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // biphenylalanine
        {"[C@@H](C([*+2])=O)(Cc1cc(c2ccccc2)ccc1)N[*+]",     "BPhe",    "BPhe",	AMINO_ACIDS|NON_STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // alpha-naphtylalanine
        {"[C@@H](C([*+2])=O)(Cc1cccc2ccccc12)N[*+]",	     "aNal",    "aNal",	AMINO_ACIDS|NON_STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},
        // beta-naphtylalanine
        {"[C@@H](C([*+2])=O)(Cc1cc(cccc2)c2cc1)N[*+]",	     "bNal",    "bNal",	AMINO_ACIDS|NON_STANDARD_SHORTCUT,
         "AMINO_N",	"AMINO_C",	""},

        // catch-all aminoacid shortcut
        {"[C@@H](C([*+2])=O)(C[*+3])N[*+]",	             "AA",	"AA",	AMINO_ACIDS|CATCH_ALL,
         "AMINO_N",	"AMINO_C",	"AA_SUBST"},

        // catch-all methylated aminoacid shortcut
        {"[C@@H](C([*+2])=O)(C[*+3])N(C)[*+]",	             "MeAA",	"MeAA",	AMINO_ACIDS|CATCH_ALL,
         "AMINO_N_ME",	"AMINO_C",	"AA_SUBST"},

        // protecting groups
        {"c1ccccc1COC(=O)[*+]",	                             "Z",	"Cbz",	PROTECTING_GROUP|STANDARD_SHORTCUT,
         "URETHANE",	"",	        ""},
        {"CC(=O)[*+]",	                                     "Ac",	"Ac",	PROTECTING_GROUP|STANDARD_SHORTCUT,
         "AMINO_C",	"",	        ""},
    };

static void StripRAtoms(struct reaccs_molecule_t *mp, char *code)
/*
 * Removes the R-atoms from the query molecule mp and color the attached atoms with their mass if any.
 * This trick is used to remember the attachment position.
 *
 * The method also set query_H_count and sub_desc to represent the fact that only R-atoms indicate possible
 * substitutions.
 */
{
    struct reaccs_atom_t *ap, *ap1, *ap2;
    struct reaccs_bond_t *bp;
    int i, j, next_map, tmp_color;
    int rmap[10];
    neighbourhood_t *nbp;
    int rgnum;
    struct prop_line_t *pl_pruned, *pl_tmp;
    int *H_count;

    MakeHydrogensImplicit(mp);
    /* Set up hydrogen count fields in structure for matching */
    H_count = TypeAlloc(mp->n_atoms+1, int);
    ComputeImplicitH(mp, H_count);

    nbp = TypeAlloc(mp->n_atoms, neighbourhood_t);
    SetupNeighbourhood(mp, nbp, mp->n_atoms);

    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    {
        ap1 = mp->atom_array+(bp->atoms[0]-1);
        if (0 == strcmp(ap1->atom_symbol, "R#"))    // make RGP labels look like mass labels
        {
           if (!GetNumProperty(mp->prop_lines, "M  RGP", 1+(ap1-mp->atom_array), &rgnum)) rgnum = 1;
           ap1->mass_difference = rgnum;
           strcpy(ap1->atom_symbol, "R");
        }
        if (0 == strcmp(ap1->atom_symbol, "R"))
        {
            ap2 = mp->atom_array+(bp->atoms[1]-1);
            ap2->sub_desc = nbp[bp->atoms[1]-1].n_ligands;
            if (ap1->charge > 0)
                ap2->color |= (1<<(ap1->charge-1));
            else if (ap1->mass_difference > 0)
                ap2->color |= (1<<(ap1->mass_difference-1));
            else
                ap2->color |= 1;
// fprintf(stderr, "marking atom %d with color %d/%d\n", (ap2-mp->atom_array)+1, ap2->color, ap1->mass_difference);
           ap2->mapping = ap1->mass_difference;
        }
        ap1 = mp->atom_array+(bp->atoms[1]-1);
        if (0 == strcmp(ap1->atom_symbol, "R#"))    // make RGP labels look like mass labels
        {
           if (!GetNumProperty(mp->prop_lines, "M  RGP", 1+(ap1-mp->atom_array), &rgnum)) rgnum = 1;
           ap1->mass_difference = rgnum;
           strcpy(ap1->atom_symbol, "R");
        }
        if (0 == strcmp(ap1->atom_symbol, "R"))
        {
            ap2 = mp->atom_array+(bp->atoms[0]-1);
            ap2->sub_desc = nbp[bp->atoms[0]-1].n_ligands;
            if (ap1->charge > 0)
                ap2->color |= (1<<(ap1->charge-1));
            else if (ap1->mass_difference > 0)
                ap2->color |= (1<<(ap1->mass_difference-1));
            else
                ap2->color |= 1;
// fprintf(stderr, "marking atom %d with color %d/%d\n", (ap2-mp->atom_array)+1, ap2->color, ap1->mass_difference);
            ap2->mapping = ap1->mass_difference;
        }
    }
    // Remove RGP properties if any
    pl_pruned = (struct prop_line_t *)NULL;
    while (mp->prop_lines)
    {
        if (0 == strncmp(mp->prop_lines->text, "M  RGP", 6)) // found an R-group line => remove it
        {
            pl_tmp = mp->prop_lines->next;
            MyFree((char *)mp->prop_lines);
            mp->prop_lines = pl_tmp;
        }
        else
        {
            pl_tmp = mp->prop_lines;
            mp->prop_lines = mp->prop_lines->next;
            pl_tmp->next = pl_pruned;
            pl_pruned = pl_tmp;
        }
    }
    mp->prop_lines = pl_pruned;
    // define query_H_count and sub_desc fields
    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    {
       if (ap->color != NONE)   // atom with attachment => fix hydrogen count only
       {
          ap->sub_desc = NONE;
          ap->query_H_count = ZERO_COUNT+H_count[i+1];
       }
       else                     // other atom => fix substitution count
       {
          ap->sub_desc = SUB_AS_IS;
          ap->query_H_count = NONE;
       }
    }
    mp->n_props = CountPropLines(mp->prop_lines);
    MyFree((char *)H_count);
    /* Now, we remove the R-atoms from the connection table */
    while (1)
    {
       for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
          if (0 == strcmp(ap->atom_symbol, "R")) break;

       /* Return if there is no R atom */
       if (i == mp->n_atoms)
           break;
       else
           RemoveAtomFromMolecule(mp, i+1);
    }

// PrintREACCSMolecule(stderr,mp,code);

    MyFree((char *)nbp);
    return;
}

static char *ltrim(char *s) 
{ 
    while(isspace(*s)) s++; 
    return s; 
} 
 
static char *rtrim(char *s) 
{ 
    char* back = s + strlen(s); 
    while(isspace(*--back)); 
    *(back+1) = '\0'; 
    return s; 
} 
 
static char *trim(char *s) 
{ 
    return rtrim(ltrim(s));  
} 

struct shortcut_t *MakeShortcut(char * smiles, char *short_code, char *long_code,
                                int shortcut_class,
                                char *r1class, char *r2class, char *r3class,
                                struct shortcut_t *old_list,
                                struct aa_class_t *aa_classes,
                                int input_rank)
{
    struct shortcut_t *result;
    struct reaccs_molecule_t *mp;
    neighbourhood_t *nbp;
    char buffer[100], *tokp;
    struct aa_class_t *aacp;
    unsigned long bit_mask;
    int i;

    mp = SMIToMOL(smiles, DO_LAYOUT);
    if (IsNULL(mp)) return old_list;    // molecule did not parse
    StripRAtoms(mp, long_code);
    result = TypeAlloc(1, struct shortcut_t);
    result->next = old_list;
    result->query = mp;
    result->aa_mask = TypeAlloc(mp->n_atoms, unsigned long);
    result->short_code = short_code;
    result->long_code = long_code;
    result->shortcut_class = shortcut_class;
    result->nattachments = 0;
    result->attachment_idx1 = -1;
    result->attachment_idx2 = -1;
    result->attachment_idx3 = -1;
    result->rank = input_rank;
    nbp = TypeAlloc(mp->n_atoms, neighbourhood_t);
    SetupNeighbourhood(mp, nbp, mp->n_atoms);
    for (i=0; i<mp->n_atoms; i++)
        if (3 <= nbp[i].n_ligands  &&  nbp[i].n_ligands <= 4)
        {
            mp->atom_array[i].stereo_parity = AtomParity(mp, i+1, nbp+i);
// fprintf(stderr, "stereo_parity of atom %d with %d ligands in shortcut '%s' is %d\n", i+1, nbp[i].n_ligands, result->long_code, mp->atom_array[i].stereo_parity);
        }
    SetRingSizeFlags(mp, 14, nbp);
    MyFree((char *)nbp);
    if (!IsNULL(r1class)  &&  strlen(r1class) > 0)
    {
        strncpy(buffer, r1class, 99); buffer[99] = '\0';
        for (aacp=aa_classes; !IsNULL(aacp); aacp=aacp->next)
            if (0 == strcmp(aacp->name, r1class))
                bit_mask = aacp->bit_mask;
        for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].mapping == 1)
            {
                result->nattachments = 1;
                result->attachment_idx1 = i;
                result->aa_mask[i] |= bit_mask;
                mp->atom_array[i].mapping = NONE;
                break;
            }
    }
    if (!IsNULL(r2class)  &&  strlen(r2class) > 0)
    {
        strncpy(buffer, r2class, 99); buffer[99] = '\0';
        for (aacp=aa_classes; !IsNULL(aacp); aacp=aacp->next)
            if (0 == strcmp(aacp->name, r2class))
                bit_mask = aacp->bit_mask;
        for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].mapping == 2)
            {
                result->nattachments = 2;
                result->attachment_idx2 = i;
                result->aa_mask[i] |= bit_mask;
                mp->atom_array[i].mapping = NONE;
                break;
            }
    }
    if (!IsNULL(r3class)  &&  strlen(r3class) > 0)
    {
        strncpy(buffer, r3class, 99); buffer[99] = '\0';
        for (aacp=aa_classes; !IsNULL(aacp); aacp=aacp->next)
            if (0 == strcmp(aacp->name, r3class))
                bit_mask = aacp->bit_mask;
        for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].mapping == 3)
            {
                result->nattachments = 3;
                result->attachment_idx3 = i;
                result->aa_mask[i] |= bit_mask;
                mp->atom_array[i].mapping = NONE;
                break;
            }
    }
    // make sure that capping groups are not used to start shortcut sequences
    if (result->attachment_idx2 < 0) result->rank += CAPPING_RANK;
    return result;
}

// list of augmented atom classes to be used to constrain attachment atoms of shortcuts
static struct aa_class_t *aa_classes = (struct aa_class_t *)NULL;
static struct shortcut_t *shortcuts = (struct shortcut_t *)NULL;


// static void ErrorLog(char *message)
// {
//    FILE *fp;
// 
//    fp = fopen("error.log", "a+");
//    fprintf(fp, "shortcut.c: '%s'\n", message);
//    fclose(fp);
// }

void InitShortcuts()
/*
 * Perceive standard shortcut definitions and initialize derived data structures.
 */
{
    int i;
    struct aa_class_t *aacp;
    char *class_string, *last_class, *aadef;
    unsigned long bit_mask;
    struct shortcut_t *sc1, *sc2;
    struct reaccs_molecule_t *qp;

    if (!IsNULL(shortcuts)) return;  // already intialized
// ErrorLog("Initializing shortcut table");
// fprintf(stderr, "Initializing shortcut table\n");

    bit_mask = 0;
    last_class = "";
    // first load the augmented atom classes
    for (i=0; i<sizeof(aa_class_strings)/sizeof(struct aa_class_string_t); i++)
    {
        aacp = TypeAlloc(1, struct aa_class_t);
        aacp->next = aa_classes; aa_classes = aacp;
        class_string = aa_class_strings[i].class_string;
        if (bit_mask == 0)
            bit_mask = 1;
        else if (0 != strcmp(class_string, last_class)) // new class label found => step up mask
        {
            bit_mask <<= 1;
        }
        last_class = class_string;
        aadef        = aa_class_strings[i].aadef;
        strcpy(aacp->name, class_string);
        aacp->bit_mask = bit_mask;
        StringToAugmentedAtom(aadef, &aacp->aa_pattern);
        aacp->no_ring_atom = aa_class_strings[i].no_ring_atom;
        qp = (struct reaccs_molecule_t *)NULL;
        if (!IsNULL(aa_class_strings[i].smiles_pattern))
        {
           qp = SMIToMOL(aa_class_strings[i].smiles_pattern, DO_LAYOUT);
// PrintREACCSMolecule(stderr, qp, aa_class_strings[i].smiles_pattern);
        }
        aacp->smiles_query = qp;
// fprintf(stderr, "class_string(%d) = '%s', pattern = '%s', bit_mask = %lX\n",  i, class_string, aadef, bit_mask);
    }

    for (i=0; i<sizeof(shortcut_input)/sizeof(struct shortcut_input_t); i++)
    {
        // fprintf(stderr, "Perceiving '%s'\n", shortcut_input[i]);
        shortcuts = MakeShortcut(shortcut_input[i].smiles,
                                 shortcut_input[i].single_letter_code,
                                 shortcut_input[i].shortcut,
                                 shortcut_input[i].shortcut_class,
                                 shortcut_input[i].r1class,
                                 shortcut_input[i].r2class,
                                 shortcut_input[i].r3class,
                                 shortcuts, aa_classes, i+1);
// if (i == 21) PrintREACCSMolecule(stderr,shortcuts->query,shortcuts->long_code);
        // if (i == 0) break;
    }
    // invert list of shortcuts

    sc1 = shortcuts;
    shortcuts = (struct shortcut_t *)NULL;
    while (!IsNULL(sc1))
    {
        sc2 = sc1->next;
        sc1->next =shortcuts;
        shortcuts = sc1;
        sc1 = sc2;
    }
}

#define SEED_ATOM_MASK (1<<31)

void SetShortcutMasks(unsigned long *aa_mask,
                      struct reaccs_molecule_t *mp,
                      neighbourhood_t *nbp,
                      struct aa_class_t *aa_classes)
/*
 * Perceive the augmented atom classes to be used in matching.
 */
{
    int i;
    unsigned long old_mask;
    unsigned int match[MAXNEIGHBOURS+1];
    ssmatch_t *matches, *h;
    unsigned long *qry_mask;

    // make plenty of room for query FIRST_ATOM classes
    qry_mask = TypeAlloc(mp->n_atoms, unsigned long);
    for (i=0; i<mp->n_atoms; i++) qry_mask[i] = 0;
    qry_mask[0] = SEED_ATOM_MASK;
    while (aa_classes)
    {
        for (i=0; i<mp->n_atoms; i++)
            if (AAMatch(mp, i, match, &aa_classes->aa_pattern, (int *)NULL, nbp)  &&
                ((!aa_classes->no_ring_atom)  ||  (0 == mp->atom_array[i].rsize_flags)))
            {
                if (!IsNULL(aa_classes->smiles_query)  &&               // there is a more detailed query
                    mp->n_atoms >= aa_classes->smiles_query->n_atoms)   // this query is not too big
                {
                   old_mask = aa_mask[i];
                   aa_mask[i] = SEED_ATOM_MASK; // mark current atom as the one to match the first query atom
                   matches = SSMatchExt(mp,                       aa_mask, 
                                        aa_classes->smiles_query, qry_mask,
                                        TRUE);
                   aa_mask[i] = old_mask;       // clear the temporary mark
                   if (CountMatches(matches) <= 0) continue;    // don't mark this class
// fprintf(stderr, "atom %d matches AA-class %s\n", i+1, aa_classes->name);
                    // free storage allocated by match
                    while (!IsNULL(matches))
                    {
                       h = matches->next;
                       matches->next = (ssmatch_t *)NULL;
                       FreeSSMatch(matches);
                       matches = h;
                    }
                }
                aa_mask[i] |= aa_classes->bit_mask;
            }
        aa_classes = aa_classes->next;
    }
    MyFree((char *)qry_mask);
}

static int **ComputeUncoloredDistance(struct reaccs_molecule_t *mp)
   /*
    * Computes the graph distances of colored atoms along bonds with at least one uncolored atom.
    *
    * The caller needs to deallocate the 2D array returned from this call.
    */
{
   int **result;
   int i, j, k, tmp, ai0, ai1;

   result = TypeAlloc(mp->n_atoms, int *);
   for (i=0; i<mp->n_atoms; i++) result[i] = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++)
      for (j=0; j<mp->n_atoms; j++)
         if (i==j) result[i][j] = 0;
         else      result[i][j] = -1;   // INFINITY
   // initialize with connecting bonds
   for (i=0; i<mp->n_bonds; i++)
   {
      ai0 = mp->bond_array[i].atoms[0]-1;
      ai1 = mp->bond_array[i].atoms[1]-1;
      if (mp->atom_array[ai0].color != NONE  &&  mp->atom_array[ai1].color != NONE) continue;
      result[ai0][ai1] = 1; result[ai1][ai0] = 1;
   }
   for (k=0; k<mp->n_atoms; k++)
      for (i=0; i<mp->n_atoms; i++)
         for (j=0; j<mp->n_atoms; j++)
         {
            if (result[i][k] < 0  ||  result[k][j] < 0) continue;       // infinity => no improvement
            tmp = result[i][k]+result[k][j];
            if (result[i][j] < 0) result[i][j] = tmp;
            else result[i][j] = result[i][j]<tmp ? result[i][j] : tmp;
         }
#if 0
   for (i=0; i<mp->n_atoms; i++)
   {
      if (mp->atom_array[i].color == NONE) continue;
      for (j=i+1; j<mp->n_atoms; j++)
      {
         if (mp->atom_array[j].color == NONE) continue;
         if (result[i][j] >= 0  ||  result[j][i] >= 0)
            fprintf(stderr, "result[%d][%d] = %d, result[%d][%d] = %d\n", mp->atom_array[i].color, mp->atom_array[j].color, result[i][j], mp->atom_array[j].color, mp->atom_array[i].color, result[j][i]);
      }
   }
#endif

   return (result);
}

struct reaccs_molecule_t *ApplyShortcuts(struct reaccs_molecule_t *mp, int flags)
/*
 * Apply the currently registered shortcuts to mp and return the updated molecules.
 */
{
    ssmatch_t *matches, *h;
    int i, ii, imatch, j, thread_color, nmatches, nreal_atoms;
    int ai1, ai2, closest_start;
    int is_regular_match, is_mirror_match;
    int *good_atoms, *good_bonds;
    struct shortcut_t *sc_list;
    int *H_count;
    struct reaccs_bond_t *bp;
    neighbourhood_t *nbp, *qnbp;
    unsigned long *aa_mask;
    int next_match, next_new_match, best_start;
    struct shortcut_match_t match_array[MAXATOMS];
    int match_used[MAXATOMS];
    struct shortcut_match_t new_match_array[MAXATOMS];
    struct shortcut_match_t matcht;
    struct reaccs_atom_t *ap, atom;

    int **distance;

    char buffer[200];

    MakeHydrogensImplicit(mp);
    /* Set up hydrogen count fields in structure for matching */
    H_count = TypeAlloc(mp->n_atoms+1, int);
    ComputeImplicitH(mp, H_count);
    /* Add the explicit hydrogens to the implicit counts */
    ResetColors(mp);
    for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
    {
        if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
            H_count[bp->atoms[1]]++;
        else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
            H_count[bp->atoms[0]]++;
    }
    /* set the 'query_H_count' field to the correct value */
    for (j=0; j<mp->n_atoms; j++)
        if (H_count[j+1] >= 0) 
        {
            mp->atom_array[j].query_H_count = ZERO_COUNT+H_count[j+1];
        }
    MyFree((char *)H_count);
    nbp = TypeAlloc(mp->n_atoms, neighbourhood_t);
    SetupNeighbourhood(mp, nbp, mp->n_atoms);
    // Color the atoms by the first complete shortcut match
    next_match = 0;
    SetRingSizeFlags(mp, 14, nbp);
    // perceive shortcut related augmented atoms
    aa_mask = TypeAlloc(mp->n_atoms, unsigned long);
    SetShortcutMasks(aa_mask, mp, nbp, aa_classes);
    for (sc_list = shortcuts; !IsNULL(sc_list); sc_list = sc_list->next)
    {
        // skip if not of requested type
        if (0 == (flags & sc_list->shortcut_class & TYPE_MASK)) continue;
        // skip if not of requested level
        if (0 == (flags & sc_list->shortcut_class & LEVEL_MASK)) continue;
        matches = SSMatchExt(mp,              aa_mask, 
                             sc_list->query,  sc_list->aa_mask, 
                             FALSE);
        qnbp = TypeAlloc(sc_list->query->n_atoms, neighbourhood_t);
        SetupNeighbourhood(sc_list->query, qnbp, sc_list->query->n_atoms);
        nmatches = CountMatches(matches);
        if (nmatches != 0)
        {
// sprintf(buffer, "%d matches for '%s'", nmatches, sc_list->long_code); ErrorLog(buffer);
// fprintf(stderr, "%d matches for '%s'\n", nmatches, sc_list->long_code);
// PrintREACCSMolecule(stderr,sc_list->query,sc_list->long_code);
        }
        while (!IsNULL(matches))
        {
            h = matches->next;
            matches->next = (ssmatch_t *)NULL;
// if (0 == strcmp(sc_list->long_code, "Lys"))
// {
// PrintREACCSMolecule(stderr,mp,"peptide");
// PrintREACCSMolecule(stderr,sc_list->query,sc_list->long_code);
// for (i=0; i<matches->n_match; i++) fprintf(stderr, "%d=>%d, ", i+1, matches->match_atoms[i]+1);
// fprintf(stderr, "\n");
// }
            for (i=0; i<matches->n_match; i++)
                if (mp->atom_array[matches->match_atoms[i]].color != NONE) break;
            if (i == matches->n_match)  // completely independant new match found => color it in
            {
                is_regular_match = IsStereoMatch(mp, nbp, sc_list->query, qnbp, matches, FALSE);
                if (!is_regular_match)
                    is_mirror_match = IsStereoMatch(mp, nbp, sc_list->query, qnbp, matches, TRUE);
                else
                    is_mirror_match = FALSE;
// fprintf(stderr, "match candidate of shortcut %s, is_regular=%d, is_mirror=%d\n", sc_list->long_code, is_regular_match, is_mirror_match);
                if (!is_regular_match  &&  !is_mirror_match)
                {
                    FreeSSMatch(matches);
                }
                else
                {
                    match_array[next_match].match = matches;    // needs to be deallocated after uses
                    match_array[next_match].mirror_match = !is_regular_match && is_mirror_match;
                    match_array[next_match].shortcut = sc_list; // points to global data structure so not to be deallocated
                    match_array[next_match].left_color = match_array[next_match].right_color = match_array[next_match].color = 0;
                    if (sc_list->attachment_idx1 >= 0)
                        match_array[next_match].match_index1 = matches->match_atoms[sc_list->attachment_idx1];
                    else
                        match_array[next_match].match_index1 = -1;
                    if (sc_list->attachment_idx2 >= 0)
                        match_array[next_match].match_index2 = matches->match_atoms[sc_list->attachment_idx2];
                    else
                        match_array[next_match].match_index2 = -1;
                    if (sc_list->attachment_idx3 >= 0)
                        match_array[next_match].match_index3 = matches->match_atoms[sc_list->attachment_idx3];
                    else
                        match_array[next_match].match_index3 = -1;
                    match_array[next_match].color = 1+next_match;
                    for (i=0; i<matches->n_match; i++)
                        mp->atom_array[matches->match_atoms[i]].color = 1+next_match;
                    next_match++;
                }
            }
            else
            {
                FreeSSMatch(matches);
            }
            matches = h;
        }
        MyFree((char *)qnbp);
    }
    // clear query_H_count used for matching
    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
        ap->query_H_count = NONE;
    // collect left_color and right_color of shortcut fragments
    for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
    {
        if (mp->atom_array[bp->atoms[0]-1].color == mp->atom_array[bp->atoms[1]-1].color)   // bond does not cross fragments
            continue;
        for (i=0; i<next_match; i++)
            if (mp->atom_array[bp->atoms[0]-1].color == match_array[i].color)   // bond starts at this fragment
            {
                if (bp->atoms[0]-1 == match_array[i].match_index1  &&  match_array[i].match_index2 >= 0)
                    match_array[i].left_color = mp->atom_array[bp->atoms[1]-1].color;
                if (bp->atoms[0]-1 == match_array[i].match_index2)
                    match_array[i].right_color = mp->atom_array[bp->atoms[1]-1].color;
            }
            else if (mp->atom_array[bp->atoms[1]-1].color == match_array[i].color)   // bond ends at this fragment
            {
                if (bp->atoms[1]-1 == match_array[i].match_index1  &&  match_array[i].match_index2 >= 0)
                    match_array[i].left_color = mp->atom_array[bp->atoms[0]-1].color;
                if (bp->atoms[1]-1 == match_array[i].match_index2)
                    match_array[i].right_color = mp->atom_array[bp->atoms[0]-1].color;
            }
    }
    // remove orphan matches, i.e. matches that don't have any neighbouring non-capping shortcut matches
    for (i=0; i<next_match; i++)
    {
        if (match_array[i].left_color > 0  &&  match_array[i].right_color > 0) continue;    // central shortcut of sequence of at least 3
        // right shortcut of sequence of at least three
        if (match_array[i].left_color > 0  &&  match_array[i].right_color == 0  &&  match_array[match_array[i].left_color-1].left_color > 0) continue;
        // right shortcut of sequence of at least three
        if (match_array[i].right_color > 0  &&  match_array[i].left_color == 0  &&  match_array[match_array[i].right_color-1].right_color > 0) continue;
        // for now, keep copping shortcuts
        if (match_array[i].match_index2 == 0) continue;
        // clear away match color
        for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
        {
            if (ap->color == match_array[i].color) ap->color = NONE;
        }
        match_array[i].color = NONE;
    }
    // clear away colors of capping neighbours
    for (i=0; i<next_match; i++)
    {
        if (match_array[i].left_color > 0  &&  match_array[match_array[i].left_color-1].match_index2 < 0)
            match_array[i].left_color = 0;
        if (match_array[i].right_color > 0  &&  match_array[match_array[i].right_color-1].match_index2 < 0)
            match_array[i].right_color = 0;
    }
    for (i=0; i<next_match; i++)
    {
// sprintf(buffer, "%d <= %d(%s) => %d", match_array[i].left_color, match_array[i].color, match_array[i].shortcut->long_code, match_array[i].right_color);
// ErrorLog(buffer);
// fprintf(stderr, "%d <= %d(%s) => %d\n", match_array[i].left_color, match_array[i].color, match_array[i].shortcut->long_code, match_array[i].right_color);
    }
    // reorder matches to be from left to right
    for (i=0; i<next_match; i++)
    {
        match_used[i] = FALSE;
    }
    next_new_match = 0;
    for (;;)
    {
        // search new start of sequence
        best_start = -1;
        for (i=0; i<next_match; i++)
        {
            if (match_used[i]) continue;
            if (match_array[i].left_color == 0)     // leftmost shortcut is best
            {
                best_start = i;
// sprintf(buffer, "1: best_start = %d (%s)", best_start, match_array[best_start].shortcut->long_code); ErrorLog(buffer);
// fprintf(stderr, "1: best_start = %d (%s)\n", best_start, match_array[best_start].shortcut->long_code);
                break;
            }
            if (best_start < 0  ||  match_array[i].shortcut->rank < match_array[best_start].shortcut->rank)
            {
                best_start = i;
// sprintf(buffer, "2: best_start = %d (%s)", best_start, match_array[best_start].shortcut->long_code); ErrorLog(buffer);
// fprintf(stderr, "2: best_start = %d (%s)\n", best_start, match_array[best_start].shortcut->long_code);
            }
        }
        if (best_start == -1) break;    // no new sequence start found
        // trace the sequence and copy to new_match_array
        new_match_array[next_new_match++] = match_array[best_start];
        match_used[best_start] = TRUE;
        for (;;)
        {
            if (match_array[best_start].right_color == 0) break;    // end of sequence
            for (i=0; i<next_match; i++)
            {
                if (match_used[i]) continue;
                if (match_array[i].color == match_array[best_start].right_color) break;
            }
            if (i == next_match) break;
            best_start = i;
            new_match_array[next_new_match++] = match_array[best_start];
            match_used[best_start] = TRUE;
// fprintf(stderr, "next_shortcut = %d (%s)\n", best_start, match_array[best_start].shortcut->long_code);
        }
    }
    nreal_atoms = mp->n_atoms;
    // generate shortcut atoms
    for (i=0; i<next_match; i++)
    {
        if (new_match_array[i].color == NONE) continue;
// fprintf(stderr, "NEW_ATOMS: %d <= %d(%s) => %d\n",
// new_match_array[i].left_color, new_match_array[i].color, new_match_array[i].shortcut->long_code, new_match_array[i].right_color);
        // if (new_match_array[i].color == NONE  &&  new_match_array[i].match_index2 >= 0) continue;
        mp->atom_array = TypeRealloc(mp->atom_array, mp->n_atoms, mp->n_atoms+1, struct reaccs_atom_t);
        ap = mp->atom_array+(mp->n_atoms++);
        strcpy(ap->atom_symbol, "R");
        if (new_match_array[i].mirror_match)
        {
            strcpy(ap->atext, "D");
            strcat(ap->atext, new_match_array[i].shortcut->long_code);
        }
        else
            strcpy(ap->atext, new_match_array[i].shortcut->long_code);
        ap->color = new_match_array[i].color;
        h = new_match_array[i].match;
        ap->x = 0.0; ap->y = 0.0; ap->z = 0.0;
        for (j=0; j<h->n_match; j++)
        {
            ap->x += mp->atom_array[h->match_atoms[j]].x;
            ap->y += mp->atom_array[h->match_atoms[j]].y;
        }
        ap->x /= h->n_match; ap->y /= h->n_match;
        if (new_match_array[i].left_color   == NONE) thread_color = new_match_array[i].color;
        new_match_array[i].thread_color = thread_color;
    }

    // strip away covered atoms
    good_atoms = TypeAlloc(mp->n_atoms+1, int);
    for (i=1; i<=mp->n_atoms; i++) good_atoms[i] = TRUE;
    good_bonds = TypeAlloc(mp->n_bonds, int);
    for (i=0; i<mp->n_bonds; i++) good_bonds[i] = TRUE;
    // reassign connecting bonds
    for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
    {
        if (mp->atom_array[bp->atoms[0]-1].color == mp->atom_array[bp->atoms[1]-1].color)   // bond does not cross fragments
        {
            if (mp->atom_array[bp->atoms[0]-1].color != NONE) good_bonds[j] = FALSE;
            continue;
        }
        if (mp->atom_array[bp->atoms[1]-1].color != NONE)
        {
            for (i=nreal_atoms; i<mp->n_atoms; i++)
                if (mp->atom_array[bp->atoms[1]-1].color == mp->atom_array[i].color)
                {
                    bp->atoms[1] = i+1;
                    break;
                }
        }
        if (mp->atom_array[bp->atoms[0]-1].color != NONE)
        {
            for (i=nreal_atoms; i<mp->n_atoms; i++)
                if (mp->atom_array[bp->atoms[0]-1].color == mp->atom_array[i].color)
                {
                    bp->atoms[0] = i+1;
                    break;
                }
        }
    }
    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    {
        if (i < nreal_atoms  &&  ap->color != NONE) good_atoms[i+1] = FALSE;
ap->mapping = 0*ap->color;
    }
    StripMolecule(mp, good_atoms, good_bonds);

// PrintREACCSMolecule(stderr,mp,"After Stripping");

    // compute atom distance matrix
    distance = ComputeUncoloredDistance(mp);
    for (i=0; i<next_match; i++)
    {
       if (new_match_array[i].color == NONE) continue;
       if (new_match_array[i].right_color == NONE)  // end of thread => find closest other thread start
       {
          closest_start = -1;
          // search for atom to connect
          for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
             if (ap->color == new_match_array[i].color) break;
          if (j == mp->n_atoms) continue;        // not found
          ai1 = j+1;
          imatch = -1;
          // search for closest thread start in different thread
          for (ii=0; ii<next_match; ii++)
          {
             if (new_match_array[ii].thread_color == new_match_array[i].thread_color) // => ignore same thread
                continue;
             if (new_match_array[ii].color == NONE) continue;
             if (new_match_array[ii].left_color == NONE)  // start of thread => check if closest start
             {
                // search for candidate atom
                for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
                   if (mp->atom_array[j].color == new_match_array[ii].color) break;
                if (j == mp->n_atoms) continue;        // not found
                ai2 = j+1;
                if (distance[ai1-1][ai2-1] == -1) continue;     // not connected
                if (closest_start > 0  &&  distance[ai1-1][ai2-1] > distance[ai1-1][closest_start-1]) continue;        // poorer distance
                closest_start = ai2;    // better (or first) distance => save it
                imatch = ii;
// fprintf(stderr, "possibly connecting shortcut '%s' %d[%d] with shortcut '%s' %d[%d] by distance %d\n",
// new_match_array[i].shortcut->long_code,  mp->atom_array[ai1-1].color,           ai1,
// new_match_array[imatch].shortcut->long_code, mp->atom_array[closest_start-1].color, closest_start,
// distance[ai1-1][closest_start-1]);
             }
          }
          if (closest_start > 0)
          {
// fprintf(stderr, "connecting shortcut '%s' %d[%d] with shortcut '%s' %d[%d] by distance %d\n",
// new_match_array[i].shortcut->long_code,  mp->atom_array[ai1-1].color,           ai1,
// new_match_array[imatch].shortcut->long_code, mp->atom_array[closest_start-1].color, closest_start,
// distance[ai1-1][closest_start-1]);
             // saturating chain reference colors
             new_match_array[i].right_color = new_match_array[imatch].color;
             new_match_array[imatch].left_color = new_match_array[i].color;
             // adding NONE bond
             mp->bond_array = TypeRealloc(mp->bond_array, mp->n_bonds, mp->n_bonds+1, struct reaccs_bond_t);
             bp = mp->bond_array+(mp->n_bonds++);
             bp->bond_type = NONE;
             bp->atoms[0] = ai1; bp->atoms[1] = closest_start;
          }
       }
    }
    // for non-cicular sequences, find leftmost shortcut and swap it to the front
    for (i=0; i<next_match; i++)
    {
       if (new_match_array[i].color == NONE) continue;
       if (new_match_array[i].left_color == NONE)  // end of thread => find closest other thread start
          break;
    }
    if (i < next_match) // found a leftmost color
    {
       imatch = i;
       for (ii=0; ii<mp->n_atoms; ii++)
          if (mp->atom_array[ii].atext[0] != '\0'  &&  0 == strcmp("R", mp->atom_array[ii].atom_symbol))
             break;
       for (i=ii; i<mp->n_atoms; i++)
          if (mp->atom_array[i].atext[0] != '\0'  &&  0 == strcmp("R", mp->atom_array[i].atom_symbol)  &&  mp->atom_array[i].color == new_match_array[imatch].color)
             break;
       if (ii < mp->n_atoms  &&  i < mp->n_atoms)
       {
          atom = mp->atom_array[ii];
          mp->atom_array[ii] = mp->atom_array[i];
          mp->atom_array[i] = atom;
          for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
          {
             if (bp->atoms[0] == i+1)       bp->atoms[0] = ii+1;
             else if (bp->atoms[0] == ii+1) bp->atoms[0] = i+1;
             if (bp->atoms[1] == i+1)       bp->atoms[1] = ii+1;
             else if (bp->atoms[1] == ii+1) bp->atoms[1] = i+1;
          }
       }
    }
    // free distance matrix
    for (i=0; i<mp->n_atoms; i++) MyFree((char *)distance[i]);
    MyFree((char *)distance);

    // Free (return to SS heap) the utilized precise matches
    for (i=0; i<next_match; i++)
    {
       if (!IsNULL(match_array[i].match))
       {
          FreeSSMatch(match_array[i].match);
       }
    }
// PrintREACCSMolecule(stderr,mp,"COLORED");
    FreeSSMatchHeap();
    MyFree((char *)aa_mask);
    MyFree((char *)good_bonds);
    MyFree((char *)good_atoms);
    MyFree((char *)nbp);
    return mp;
}
