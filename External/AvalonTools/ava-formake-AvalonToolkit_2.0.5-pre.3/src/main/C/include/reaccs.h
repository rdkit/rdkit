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

/*************************************************************************/
/*                                                                       */
/*    File:           reaccs.h                                           */
/*                                                                       */
/*    Purpose:        This file contains the declarations of the         */
/*                    types used to handle REACCS data structures.       */
/*                    The 'reaccs' flavor is the format that was cabable */
/*                    of representing most features w/o resorting to     */
/*                    version 2 property lists, so these data structure  */
/*                    are used through-out the toolkit even for          */
/*                    non-reactions. Naming of fields is mainly historic.*/
/*                                                                       */
/*************************************************************************/

#ifndef _REACCS_H_
#define _REACCS_H_        1

#include "mdl.h"
#include "symbol_lists.h"

// number of atoms and bonds limited to be able to use version 2 MOL file format
#ifdef _Windows
#ifdef __WIN32__
#define MAXATOMS       999
#define MAXBONDS       999
#else
#define MAXATOMS       256
#define MAXBONDS       256
#endif
#else
#define MAXATOMS       999
#define MAXBONDS       999
#endif

// struct that represents a MOL file atom
struct reaccs_atom_t
   {
      float x, y, z;          /* atom coordinates in angstrom  */
      char atom_symbol[4];    /* chemical symbol of atom       */
      int mass_difference;    /* difference to isotope of mass */
                              /* closest to natural mixture    */
      /* int xcharge_radical;      obsolete, replaced by the following fields */
                              /* 0=none, 1=+3, 2=+2, 3=+1, 4=radical, */
                              /*         5=-1, 6=-2, 7=-3             */
      int charge;
      int radical;
      int stereo_parity;      /* 0=none, 1=odd, 2=even, 3=unmarked    */
                              /* Remark:                              */
                              /* this appears to be the SEMA-parity   */
                              /* this value is impossible to compute  */
                              /* w/o performing the SEMA canonization.*/
                              /* We use this field as a scratch field */
                              /* for numbered parity.                 */
      int query_H_count;
      int query_stereo_box;
      int dummy1;             /* valence                      */
      int dummy2;             /* H0 designator (CPSS) */
      int dummy3;             /* reaction component type (CPSS) */
      int dummy4;             /* reaction component number (CPSS) */
      int mapping;
      int second_stereo_parity; /* inversion (1) / retention (2) flag */
      int dummy5;               /* exact change flag */
      /* fortran input format: '(F10.4,F10.4,F10.4,X,A3,I2,I3,I3,I3,I3,4I3,I3,I3,I3)' */
      int sub_desc;           /* substitution count descriptor */
                              /* 0 == don't care, >0 == subst count */
                              /* -1 == no subst, -2 = as is */
      float value;            /* field used to store values during    */
                              /* processing, no chemical meaning */
      unsigned int rsize_flags; /* used by PerceiveRingSizes to flag the size */
                                /* of the rings this bond participates in */
      int color;              /* dummy field used for graph coloring  */
                              /* algorithms. No chemical meaning!     */
      char atext[MDL_MAXLINE+1];
                             /* atom text if any. Should only be used if symbol == "R" */
   };

                 /* MDL's codes for substitution counts */
#define SUB_ZERO      -1        /* no substitution */
#define SUB_AS_IS     -2        /* usually indicated by s* in molecule editors */
#define SUB_ONE        1        /* this and the subsequent values name the substitution count */
#define SUB_MORE       6        /* 6 and more substitutions are treated the same */

// struct that represents a MOL file bond
struct reaccs_bond_t
   {
      int atoms[2];             /* bond is from atom[0] to atom[1]      */
      int bond_type;            /* 1=single, 2=double, 3=triple;
                                   queries: 4=aromatic, 5=S/D, 6=S/A, 7=D/A */
      int stereo_symbol;        /* 0=none, 1=up, 6=down, 4=either,
                                   3=cis/trans either */
      int dummy;                /* used to store ring class 3,4 or 5,6
                                   or >6 ring */
      unsigned int rsize_flags; /* used by PerceiveRingSizes to flag the size */
                                /* of the rings this bond participates in */
      int topography;           /* 0=normal, 1=ring, 2=chain            */
      int reaction_mark;        /* 1=center, 0=unmarked, -1=not center  */
                                /* '(I3,I3,I3,I3,I3,I3,I3)'             */
      float value;              /* field used to store values during    */
                                /* processing                           */
      int bond_type_flags;	/* mainly used in MDL-style aromaticity */
                                /* handling for substructure testing    */
      int color;                /* dummy field used for graph coloring  */
                                /* algorithms. No chemical meaning!     */
   };

// reaction mark constants
#define  NOT_CENTER      (-1)
#define  CENTER            1
#define  UNCHANGED         2
#define  MAKE_BREAK        4
#define  UNCERTAIN         6
#define  CHANGE            8

// size of molname
#define MAXNAME         80
// size of name of writing program
#define NPROGNAME        8

struct stext_line_t
   {
      struct stext_line_t *next;  /* pointer to next line in list */
      float x, y;                 /* coordinates */
      char text[MDL_MAXLINE+1];   /* actual text                  */
   };

struct prop_line_t                /* property line                */
   {
      struct prop_line_t *next;   /* pointer to next line in list */
      char text[MDL_MAXLINE+1];   /* actual text                  */
   };

struct reaccs_molecule_t
   {
      char name[MAXNAME+1];     /* simple name of the molecule          */
                                /* '(A80)'                              */
      char user_initials[2+1];  /* first and last user initial          */
      char program_name[NPROGNAME+1];   /* name of creating program     */
      char date[6+1];           /* date of creation "MMDDYY"            */
      char time[4+1];           /* time of creation "hhmm"              */
      char dimensionality[2+1]; /* "2D" or "3D"                         */
      int  scale1;              /* first scaling factor (unknown use)   */
      float scale2;             /* second scaling factor (unknown use)  */
      float energy;             /* energy if "3D"-structure             */
      long registry_number;     /* internal registry number of molecule */
                                /* '(A2,A8,A6,A4,A2,I2,F10.5,F12.5,I6)' */
      char comment[MDL_MAXLINE+1];  /* one line of comment              */
                                /* '(A80)'                              */
      unsigned int n_atoms;     /* number of atoms                      */
      unsigned int n_bonds;     /* number of bonds                      */
      int n_atom_lists;         /* number of atom type lists for queries*/
      int dummy1;               /* not used; was reserved for nfragments*/
      int chiral_flag;          /* 1 if compound is chiral, 0 if not    */
      struct stext_line_t
         *stext_lines;                /* stext lines                          */
      int n_props;              /* number of additional properties lines*/
      struct prop_line_t
         *prop_lines;               /* uninterpreted property lines         */
      char version[8];
                                /* '(I3,I3,I3,I3,I3,I3,12X,I3,A6)'      */
      struct reaccs_atom_t
         *atom_array;           /* pointer to array of atom descriptors */
      struct reaccs_bond_t
         *bond_array;           /* pointer to array of bond descriptors */
      struct symbol_list_t
         *symbol_lists;         /* lists of allowed atom type strings */
      struct reaccs_molecule_t
         *next;                 /* pointer to next molecule in a list   */
      int color;                /* used in coloring algorithms testing  */
                                /* for connectedness or as a general    */
                                /* attribute of the molecule.           */
   };

// defines for the infamous 'chiral flag'
#define CHIRAL     1
#define ACHIRAL    0

struct reaccs_reaction_t
   {
      char name[MAXNAME+1];     /* simple name of the reaction          */
                                /* '(A80)'                              */
      char user_initials[4+1];  /* first, second, third, and last initial */
      char program_name[NPROGNAME+1];   /* name of creating program     */
      char date[6+1];           /* date of creation "MMDDYY"            */
      char time[4+1];           /* time of creation "hhmm"              */
      long registry_number;     /* internal registry number of reaction */
                                /* '(A4,2X,A8,2X,A6,A4,2X,I6)' */
      char comment[MDL_MAXLINE+1];  /* one line of comment              */
                                /* '(A80)'                              */
      unsigned int n_reactants; /* number of reactants                  */
      unsigned int n_products;  /* number of products                   */
                                /* '(I3,I3)'                            */
      struct reaccs_molecule_t
         *reactants;            /* list of reactant structures          */
      struct reaccs_molecule_t
         *products;             /* list of product structures             */

      int molecularity;         /* Stores the combined molecularity of  */
                                /* the strategic bonds in the reaction  */
      struct reaccs_reaction_t
         *next;                     /* pointer to next reaction in a list   */
   };

#define INTRAMOLECULAR     1
#define INTERMOLECULAR     2
#define REAGENT_ADDITION   4

#define       MAXDATA         255

// used to represent the data line section of an SDF file
struct data_line_t
   {
      struct data_line_t *next;
      char data[MAXDATA+1];
   };

#endif
