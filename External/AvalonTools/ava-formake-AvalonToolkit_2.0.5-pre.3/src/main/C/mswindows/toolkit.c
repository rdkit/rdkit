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
/************************************************************************
 *
 *   File:         toolkit.c
 *
 *   Purpose:      Dynamic link library to manipulate chemical structure
 *                 through a Windows VB application.
 *                 Note: The file depict.c provides imaging tools
 *                 for structures.
 *
 *************************************************************************/

#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

/* This set of macros is necessary to make a source file that exports */
/* the DLL entry points in uppercase on both VC++ and Borland C++     */
#if defined(__BORLANDC__)
#else
#define _export __declspec(dllexport) __stdcall
#endif

#define DEBUG 1
#define DEBUG_HEAP FALSE

static int checkHeap(int where)
{
   HGLOBAL heap;
   int result = 0;
   if (!DEBUG_HEAP) return (result);

   heap = GetProcessHeap();
   result = HeapValidate(heap,0,NULL);
   if (!result)
      fprintf(stderr, "HeapValidate(%d, %d, %d) returns %d\n",
         heap, 0, NULL, HeapValidate(heap,0,NULL));
   return (result);
}

static void WMFLog(char *message, int natoms, int nbonds, int length)
{
   FILE *fp;

   fp = fopen("wmf.log", "a+");
   fprintf(fp, "%15s %3d %3d %5d\n", message, natoms, nbonds, length);
   fclose(fp);
}

#include <math.h>
#include <string.h>

#include "local.h"
#include "forio.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "utilities.h"
#include "perceive.h"
#include "shortcut.h"
#include "layout.h"
#include "smi2mol.h"
#include "canonizer.h"
#include "ssmatch.h"

#include "mymacros.h"

#include "toolkit.h"


int _export FAR PASCAL CanonicalizeSmiles(LPSTR smiles,
                                          LPSTR cansmi,
                                          int nsmi,
                                          unsigned int level)
/*
 * Computes the canonical form of smiles and places it in cansmi[0..nsmi-1].
 * It returns the length of the canonical smiles if it worked, the negative
 * length if the size of the buffer was too small, and 0 in case of error.
 *
 * grouping can be NONE, COMPONENT_SET, MAIN_PARTS, or NON_STEREO_MAIN
 * to indicate what the level of precision for identity is.
 */
{
   char *tmpsmiles;
   int len;

   ShowMessageS("smiles = '%s'", "CanonicalizeSmiles", smiles);
   if (level == NON_STEREO_MAIN_PARTS)
      tmpsmiles = CanSmiles(smiles, COMPONENT_SET|MAIN_PART);
   else if (level == MAIN_PARTS)
      tmpsmiles = CanSmiles(smiles,
                            COMPONENT_SET|DB_STEREO|CENTER_STEREO|MAIN_PART);
   else if (level == COMPONENTS)
      tmpsmiles = CanSmiles(smiles, COMPONENT_SET|DB_STEREO|CENTER_STEREO);
   else
      tmpsmiles = CanSmiles(smiles, DB_STEREO|CENTER_STEREO);
   if (tmpsmiles == NULL) return (0);
   len = strlen(tmpsmiles);
   if (nsmi < len+1)
   {
      MyFree((char *)tmpsmiles);
      return (-len);
   }
   strcpy(cansmi, tmpsmiles);
   MyFree((char *)tmpsmiles);
   return (len);
}

int _export FAR PASCAL TransformSmiles(LPSTR insmiles,
                                       LPSTR outsmiles,
                                       int nsmi,
                                       unsigned int flags)
/*
 * Transforms the structure denoted by insmiles into an outsmiles applying the methods indicated by flags.
 *
 * The flag APPLY_SHORTCUTS abbreviates the common shortcuts and returns an extended SMILES format with atom texts
 * enclised in {} pairs.
 */
{
   char *tmpsmiles;
   int len;
   struct reaccs_molecule_t *mp;
   char *coordp;

   ShowMessageS("smiles = '%s'", "TransformSmiles", insmiles);
   mp = SMIToMOL(insmiles, DROP_TRIVIAL_HYDROGENS|DO_LAYOUT);
   if (flags & APPLY_SHORTCUTS)
   {
       InitShortcuts();
       mp = ApplyShortcuts(mp, AMINO_ACIDS|STANDARD_SHORTCUT|EXTENDED_SHORTCUT|NON_STANDARD_SHORTCUT);
   }
   tmpsmiles = MOLToSMIExt(mp, ISOMERIC_SMILES, (int *)NULL, &coordp);
   if (coordp)
   {
      MyFree((char *)coordp);
      coordp = NULL;
   }
   if (mp != NULL)
   {
       FreeMolecule(mp);
   }
   if (tmpsmiles == NULL) return (0);
   len = strlen(tmpsmiles);
   if (nsmi < len+1)
   {
      MyFree((char *)tmpsmiles);
      return (-len);
   }
   strcpy(outsmiles, tmpsmiles);
   MyFree((char *)tmpsmiles);
   return (len);
}

int _export FAR PASCAL HexFingerprintFromSmiles(LPSTR smiles,
                                                LPSTR hexfp,
                                                int nchar)
/*
 * Computes a hexadecimal representation of the Avalon fingerprint of smiles.
 * The function assumes that at least nchar+1 bytes are allocated.
 */
{
   struct reaccs_molecule_t *mp;
   char *fp;
   int i, nibble0, nibble1;
   mp = SMIToMOL(smiles, DROP_TRIVIAL_HYDROGENS);
   fp = TypeAlloc(nchar/2, char);


   if (mp != NULL)
   {
       SetFingerprintBits(mp, fp, nchar/2, USE_ALL_FEATURES, 0, 0);
       SetFingerprintBits(mp, fp, nchar/2, USE_ALL_FEATURES, 0, ACCUMULATE_BITS|USE_DY_AROMATICITY);
       FreeMolecule(mp);
   }
   for (i=0; i<nchar/2; i++)
   {
      nibble0 = 0x0F & (int)fp[i];
      if (nibble0 < 10) hexfp[2*i] = '0'+nibble0;
      else              hexfp[2*i] = 'A'+nibble0-10;

      nibble1 = 0x0F & ((0xF0 & (int)fp[i]) >> 4);
      if (nibble1 < 10) hexfp[2*i+1] = '0'+nibble1;
      else              hexfp[2*i+1] = 'A'+nibble1-10;
   }
   MyFree(fp);
}


int _export FAR PASCAL HexFingerprintFromMolString(LPSTR mol_string,
                                                   LPSTR hexfp,
                                                   int nchar,
                                                   int asQuery)
/*
 * Computes a hexadecimal representation of the Avalon fingerprint of mol_string.
 * as_query indicates if the FP should interpret the CT as a query.
 * The function assumes that at least nchar+1 bytes are allocated.
 */
{
   struct reaccs_molecule_t *mp;
   char *fp;
   int i, nibble0, nibble1;
   Fortran_FILE *ffp;
   int has_R = FALSE;
   neighbourhood_t *nbp;

   if (strchr(mol_string,'\n') != NULL)
   {    // MOL file string
      ffp = FortranStringOpen((char *)mol_string);
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
      {
         fprintf(stderr, "failed to read molecule\n");
         FreeMolecule(mp);
         mp = (struct reaccs_molecule_t *)NULL;
      }
      FortranClose(ffp);
   }
   else
   {   // single line ==> SMILES format
       mp = SMIToMOL(mol_string, DROP_TRIVIAL_HYDROGENS);
   }

   // quick check to avoid somewhat expensive processing
   if (mp != NULL)
       for (i=0; i<mp->n_atoms; i++)
          if (0 == strcmp(mp->atom_array[i].atom_symbol, "R")) has_R = TRUE;
   if (has_R)
   {
          MakeHydrogensImplicit(mp);
          nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
          SetupNeighbourhood(mp,nbp,mp->n_atoms);
          PerceiveMarkush(mp, nbp);
          SetRingSizeFlags(mp, 14, nbp);
          MyFree((char *)nbp);
   }

   fp = TypeAlloc(nchar/2, char);
   if (mp != NULL)
   {
       SetFingerprintBits(mp, fp, nchar/2, USE_ALL_FEATURES, asQuery, 0);
       if (asQuery == 0)
           SetFingerprintBits(mp, fp, nchar/2, USE_ALL_FEATURES, asQuery, ACCUMULATE_BITS|USE_DY_AROMATICITY);
       FreeMolecule(mp);
   }
   for (i=0; i<nchar/2; i++)
   {
      nibble0 = 0x0F & (int)fp[i];
      if (nibble0 < 10) hexfp[2*i] = '0'+nibble0;
      else              hexfp[2*i] = 'A'+nibble0-10;

      nibble1 = 0x0F & ((0xF0 & (int)fp[i]) >> 4);
      if (nibble1 < 10) hexfp[2*i+1] = '0'+nibble1;
      else              hexfp[2*i+1] = 'A'+nibble1-10;
   }
   // hexfp[nchar] = '\0';
   MyFree(fp);
}

int _export FAR PASCAL HeavyAtomCount(LPSTR smiles)
/*
 * Computes the number of heavy atoms in the molecule represented by smiles.
 */
{
   struct reaccs_molecule_t *mp;
   int result = 0;
   int i;
   mp = SMIToMOL(smiles, DROP_TRIVIAL_HYDROGENS);
   if (mp != NULL)
   {
       for (i=0; i<mp->n_atoms; i++)
           if (0 != strcmp(mp->atom_array[i].atom_symbol, "H")  &&
               0 != strcmp(mp->atom_array[i].atom_symbol, "D")  &&
               0 != strcmp(mp->atom_array[i].atom_symbol, "T"))
               result++;
       FreeMolecule(mp);
   }
   return result;
}

#define FPSIZE 128

int _export FAR PASCAL SmilesSimilarity(LPSTR smiles1,
                                        LPSTR smiles2)
/*
 * Computes the fingerprint-based similarity of the structures represented by
 * the two SMILES smiles1 and smiles2 and returns an integer percentage.
 */
{
   struct reaccs_molecule_t *mp;
   char *fp1, *fp2;
   int i, j, nand, nor;
   int and_byte, or_byte;
   mp = SMIToMOL(smiles1, DROP_TRIVIAL_HYDROGENS);
   fp1 = TypeAlloc(FPSIZE, char);
   if (mp != NULL)
   {
       SetFingerprintBits(mp, fp1, FPSIZE, USE_ALL_FEATURES, 0, USE_DY_AROMATICITY);
       FreeMolecule(mp);
   }
   mp = SMIToMOL(smiles2, DROP_TRIVIAL_HYDROGENS);
   fp2 = TypeAlloc(FPSIZE, char);
   if (mp != NULL)
   {
       SetFingerprintBits(mp, fp2, FPSIZE, USE_ALL_FEATURES, 0, USE_DY_AROMATICITY);
       FreeMolecule(mp);
   }
   nand = nor = 0;
   for (i=0; i<FPSIZE; i++)
   {
       and_byte = ((0xFF&fp1[i])&(0xFF&fp2[i]));
       or_byte  = ((0xFF&fp1[i])|(0xFF&fp2[i]));
       for (j=0; j<8; j++)
       {
           int mask = (1<<j);
           if (mask&and_byte) nand++;
           if (mask&or_byte)  nor++;
       }
   }
   MyFree(fp1);
   MyFree(fp2);

   if (nor > 0) return (int)(100.0*nand/(double)nor);
   else         return 0;
}

FILE *log_file=NULL;

int _export FAR PASCAL SmilesMatchesQueryCT(LPSTR smiles,
                                            LPSTR mol_string)
/*
 * Returns TRUE if the input smiles matches the query mol_string and false
 * otherwise.
 */
{
   struct reaccs_molecule_t *mp;
   struct reaccs_molecule_t *qp;
   int result = FALSE;
   Fortran_FILE *ffp;
   ssmatch_t * matches;
   struct reaccs_bond_t *bp;
   unsigned i;
   int *H_count;
   neighbourhood_t *nbp;

checkHeap(1002);

   mp = SMIToMOL(smiles, DROP_TRIVIAL_HYDROGENS);
   /* Illegal molecules don't match at all */
   if (mp == (struct reaccs_molecule_t *)NULL)
      return (result);

   // Read query
   ffp = FortranStringOpen((char *)mol_string);
   qp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,qp,""))
   {
      fprintf(stderr, "failed to read template\n");
      FreeMolecule(qp);
      qp = (struct reaccs_molecule_t *)NULL;
   }
   FortranClose(ffp);

checkHeap(1012);

   if (qp == (struct reaccs_molecule_t *)NULL)  /* illegal query => no match */
   {
      FreeMolecule(mp);
      return (result);
   }
   nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   SetRingSizeFlags(mp, 14, nbp);
   MyFree((char *)nbp);

   /* Compute query_H_count to match the required explicit hydrogens */
   MakeHydrogensImplicit(qp);
   nbp     = TypeAlloc(qp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(qp,nbp,qp->n_atoms);
   PerceiveMarkush(qp, nbp);
   SetRingSizeFlags(qp, 14, nbp);
   MyFree((char *)nbp);
   /* PerceiveRingSizes(qp); */

   /* Set up hydrogen count fields in structure for matching */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);

   /* Add the explicit hydrogens to the implicit counts */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
         H_count[bp->atoms[1]]++;
      else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
         H_count[bp->atoms[0]]++;
   }
   /* set the 'query_H_count' field to the correct value */
   for (i=0; i<mp->n_atoms; i++)
      if (H_count[i+1] >= 0) 
         mp->atom_array[i].query_H_count = ZERO_COUNT+H_count[i+1];


checkHeap(1022);

   matches = SSMatch(mp, qp, TRUE);
   if (matches != (ssmatch_t *)NULL)
   {
      result = TRUE;
      FreeSSMatch(matches);
   }

   FreeMolecule(mp); FreeMolecule(qp);

checkHeap(1032);

   MyFree((char *)H_count);
   return (result);
}
