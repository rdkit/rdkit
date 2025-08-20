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
/*   File:           depictutil.c                                       */
/*                                                                      */
/*   Purpose:        Utility functions used by depict.c. This set of    */
/*                   functions is designed to be non-Windows specific.  */
/*                   This module does not itself export any functions.  */
/*                                                                      */
/************************************************************************/

#include "depictutil.h"

#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <string.h>

#include "local.h"
#include "forio.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "utilities.h"
#include "perceive.h"
#include "layout.h"
#include "smi2mol.h"

#include "mymacros.h"

void SMILESStringToMOLFile(char* smiles, char* fname)
/*
 * Converts the SMILES string *smiles to a MOL-file named *fname.
 *
 * The file is cleared in case of error.
 */
{
   FILE *fp;
   struct reaccs_molecule_t *mp;
   int i;

   mp = SMIToMOL(smiles, DO_LAYOUT);
   /* The following is a patch to get the correct sizes into ISIS */
   fp = fopen(fname,"w");
   if (mp)
   {
      for (i=0; i<mp->n_atoms; i++)
      {
         mp->atom_array[i].x *= 0.5;
         mp->atom_array[i].y *= 0.5;
      }

      PrintREACCSMolecule(fp,mp,"");
      FreeMolecule(mp);
   }
   fclose(fp);
   return;
}

int CTStringToSmiles(char* ct_string, char* smiles, int nbuf)
/*
 * Write the isomeric SMILES corresponding to the MOL file in ct_string into
 * the buffer smiles[0..nbuf-1].
 *
 * The function returns the size of the SMILES written, 0 on error, and
 * the negative required size of the buffer if not enough space was provided.
 *
 * Note: The buffer must also provide space for the terminal '\0'.
 */
{
   struct reaccs_molecule_t *mp;
   char *tmp_smiles;
   Fortran_FILE *fp;
   int result;

   if (!ct_string  ||  !smiles) return (FALSE);

   result = FALSE;

   fp = FortranStringOpen(ct_string);
   if (fp)
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL == ReadREACCSMolecule(fp,mp,"") && mp->n_atoms > 0)
      {
	 tmp_smiles = MOLToSMI(mp, ISOMERIC_SMILES);
	 if (nbuf >= strlen(tmp_smiles)+1)
	 {
	    strcpy(smiles, tmp_smiles);
	    result = strlen(tmp_smiles);
	 }
	 else
	 {
	    result = -strlen(tmp_smiles)-1;
	 }
         MyFree(tmp_smiles);
      }
      else
      {
         /* NOP */
      }
   }
   else
   {
      /* NOP */
   }

   FreeMolecule(mp);
   FortranClose(fp);
   return (result);
}

void SmilesToMWMF(char* smiles, double *mwp, char* mfbuffer, int bufsize)
/*
 * Computes the molecular weight of the molecule defined
 * by smiles and puts it into *mwp. The buffer mfbuffer is filled with
 * a '\0' delimited string of the molecular formula. The caller is
 * responsible for providing enough space for this string.
 */
{
   struct reaccs_molecule_t *mp;

   (*mwp) = 0;	mfbuffer[0]= '\0';
   if (!smiles  ||  !mfbuffer) return;

   mp = SMIToMOL(smiles, 1);	/* Don't do a layout */

   if (mp)
   {
      (*mwp) = MolecularWeight(mp);
      MolecularFormula(mp, mfbuffer, bufsize);
      FreeMolecule(mp);
      return;
   }
   else
      return;
}

char* MOLFileToSMILESString(int* sizep, char* fname)
/*
 * Converts the MOL file stored in *fname to the corresponding
 * SMILES string and returns a pointer to it. The data has been
 * allocated for the caller which needs to deallocate it after copying
 * the SMILES.
 *
 * In case of error, the function returns a NULL pointer.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;

   char *smiles = NULL;

   ffp = FortranOpen(fname,"r");  /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
   {
      if (mp) FreeMolecule(mp);
      if (ffp) FortranClose(ffp);
      (*sizep) = -1;
      return (NULL);
   }
   FortranClose(ffp);

   smiles = MOLToSMI(mp, ISOMERIC_SMILES);          /* Convert MOL to SMILES */
   FreeMolecule(mp);
   if (smiles == NULL)
   {
      smiles = TypeAlloc(1,char);
      smiles[0] = '\0';
   }
   (*sizep) = strlen(smiles);
   return (smiles);
}

#define UNIQUESMILES	1
#define ADDCOORDINATES	2

char* MOLFileToSMILESExt(int* sizep, char* fname, int flags)
/*
 * Converts the MOL file stored in *fname to the corresponding
 * SMILES string and returns a pointer to it. The storage is allocated
 * for the caller who shall deallocate it after use/copy.
 *
 * The parameter flags defines the kind processing required. It
 * can request that stereochemistry be ignored (UNIQUESMILES) and
 * that a comma separated list of coordinates for each atom in the
 * SMILES be appended.
 *
 * In case of error, the function returns a NULL pointer.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;

   char *smiles    = NULL;
   char *smibuffer = NULL;

   char *coords = NULL;

   ffp = FortranOpen(fname,"r");  /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
   {
      if (mp) FreeMolecule(mp);
      if (ffp) FortranClose(ffp);
      (*sizep) = -1;
      return (NULL);
   }
   FortranClose(ffp);

   if (flags & UNIQUESMILES)
      smiles = MOLToSMIExt(mp, 0, (int *)NULL, &coords); /* ignore stereo */
   else
      smiles = MOLToSMIExt(mp, ISOMERIC_SMILES, (int *)NULL, &coords);
   FreeMolecule(mp);

   /* allocate persitent memory */
   if (coords && smiles)
   {
      smibuffer = TypeAlloc(strlen(smiles)+1 + strlen(coords)+1, char);

      strcpy(smibuffer, smiles);            /* copy SMILES */
      (*sizep) = strlen(smiles);

      if (flags & ADDCOORDINATES)
      {
	 strcat(smibuffer, " ");
	 strcat(smibuffer, coords);            /* append coordinates */
	 (*sizep) += 1+strlen(coords);
      }

      MyFree((char *) smiles);               /* free temporay storage */
      MyFree((char *) coords);               /* free temporay storage */
   }
   else
   {
      smibuffer = NULL;
      (*sizep) = -1;
   }

   return (smibuffer);
}

char* MOLFileToSMARTSString(int* sizep, char* fname)
/*
 * Converts the MOL file query stored in *fname to the corresponding
 * SMARTS string and returns a pointer to it. The storage is allocated
 * for the caller who shall deallocate it after use/copy.
 *
 * In case of error, the function returns a NULL pointer.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;
   char *smiles = NULL;

   ffp = FortranOpen(fname,"r");    /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,"")  ||
       mp->n_atoms == 0)
   {
      if (mp) FreeMolecule(mp);
      if (ffp) FortranClose(ffp);
      (*sizep) = -1;
      return (NULL);
   }
   FortranClose(ffp);

   /* Convert MOL to SMARTS */
   smiles = MOLToSMI(mp, ISOMERIC_SMILES | SMARTS_PERCEPTION);
   FreeMolecule(mp);
   if (smiles == NULL  ||  strlen(smiles) == 0)
   {
      (*sizep) = 0;
      return (NULL);
   }
   else
      (*sizep) = strlen(smiles);
   return (smiles);
}

int MOLFileToSMARTSBuffer(int* sizep, char* buffer, char* fname)
/*
 * Converts the MOL file query stored in *fname to the corresponding
 * SMARTS string and writes the result to buffer[]. Temporary storage is
 * removed. (*sizep) shall contain the buffer size on input and will
 * become the actual size of the copied SMARTS on output.
 *
 * In case of error, the function sets (*sizep) to -1 or -required_buffer_size
 * if the size of the buffer wasn't sufficient.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;
   char *smiles = NULL;


   ffp = FortranOpen(fname,"r");    /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,"")  ||  mp->n_atoms == 0)
   {
      if (mp)
      {
         FreeMolecule(mp);
         (*sizep) = -1;
      }
      else
         (*sizep) = -1;
      return (0);
   }
   FortranClose(ffp);

   /* Convert MOL to SMARTS */
   smiles = MOLToSMI(mp, ISOMERIC_SMILES | SMARTS_PERCEPTION);
   FreeMolecule(mp);
   if (smiles == NULL  ||  strlen(smiles) == 0)
   {
      (*sizep) = -1;
      return (0);
   }

   if ((*sizep) <= strlen(smiles))
   {
      (*sizep) = -strlen(smiles);
   }
   else
   {
      strcpy(buffer, smiles);            /* copy SMILES */
      (*sizep) = strlen(smiles);
   }

   MyFree((char *) smiles);               /* free temporay storage */

   return (1);
}
