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
/*                                                                   	*/
/*    File:           reaccsio.c                                        */
/*                                                                      */
/*    Purpose:        This file implements the routines to read RD-,    */
/*                    SD-, and RXN-files.                               */
//                                                                      */
//    History:        23-Dec-1992     Switched to ANSI prototypes.      */
//                                                                      */
//		              26-Jul-98	      Added definitions of routines     */
//				                      reading/writing molecules from	*/
//				                      and to strings (code by AG)       */
/*                                                                      */
/************************************************************************/

#ifdef _WIN32
#include <objbase.h>
#endif

#include "reaccsio.h"

#include <stdio.h>
#include <string.h>

#include "forio.h"
#include "local.h"
#include "reaccs.h"
#include "symbol_lists.h"
#include "utilities.h"

/**
 * Reads an ATOM description frm a version 3 MOL file.
 * mp is needed to save atom properties not stored with the atom.
 */
int ReadV30Atom(Fortran_FILE *fp,
                struct reaccs_atom_t *ap,
                struct reaccs_molecule_t *mp)
{
   int nitems;
   char buffer[MAX_BUFFER+1];
   int charge_radical;
   char prop[10][20];
   int seq;
   int idummy;
   float mass;
   int i;
   char symbol[100];
   struct symbol_list_t *list;

   if (fp->status != FORTRAN_NORMAL) return(fp->status);

   strncpy(buffer,fp->buffer,MAX_BUFFER);

   nitems = sscanf(buffer+strlen("M  V30"),
                   "%d %s %f %f %f %d %s %s %s %s %s %s %s %s %s %s",
                   &seq, symbol,
                   &ap->x,  &ap->y,  &ap->z,
                   &ap->mapping,
                   prop[0], prop[1], prop[2], prop[3], prop[4],
                   prop[5], prop[6], prop[7], prop[8], prop[9]);

   if (nitems < 6)
   {
      ShowMessageI("incorrect # (%d) of arguments on atom line",
                   "ReadV30Atom:",
                   nitems);
      ShowMessageS("buffer ==\n%s\n","ReadV30Atom",fp->buffer);
      GetBuffer(fp);
      return(FORTRAN_ERROR);
   }

   if (strlen(symbol) <= 3)
      strcpy(ap->atom_symbol, symbol);
   else     // need to parse symbol
   {
      mp->symbol_lists =
          ParseV30SymbolList(symbol, ap-mp->atom_array, mp, mp->symbol_lists);
      mp->n_atom_lists = 0;
      for (list=mp->symbol_lists; !IsNULL(list); list = list->next) mp->n_atom_lists++;
   }

   if (nitems >= 16)
      ShowMessageI("Warning %d properties for atom.\n",
                   "ReadV30Atom:", (nitems-6));

   /* Processing properties */
   for (i=0; i<nitems-6; i++)
   {
      if (1 == sscanf(prop[i], "CHG=%d", &ap->charge)) continue;
      if (1 == sscanf(prop[i], "RAD=%d", &ap->radical)) continue;
      if (1 == sscanf(prop[i], "HCOUNT=%d", &ap->query_H_count)) continue;
      if (1 == sscanf(prop[i], "CFG=%d", &ap->stereo_parity)) continue;
      if (1 == sscanf(prop[i], "VAL=%d", &ap->dummy1)) continue;
      if (1 == sscanf(prop[i], "SUBST=%d", &ap->sub_desc)) continue;
      if (1 == sscanf(prop[i], "MASS=%f", &mass))
      {
         ap->mass_difference = LookupMassDifference(mass, ap->atom_symbol);
         continue;
      }
      if (1 == sscanf(prop[i], "RBCNT=%d", &idummy)) continue;
      if (1 == sscanf(prop[i], "UNSAT=%d", &idummy)) continue;
      fprintf(stderr, "ignoring '%s' for atom %d\n", prop[i], seq);
   }
   GetBuffer(fp);
   return(FORTRAN_NORMAL);
}

/**
 * Read a V3.0 MOL file bond record.
 */
int ReadV30Bond(Fortran_FILE *fp, struct reaccs_bond_t *bp)
{
   int nitems;
   char buffer[MAX_BUFFER+1];
   char prop[10][20];
   int seq;
   int i;
   char symbol[10];

   if (fp->status != FORTRAN_NORMAL) return(fp->status);

   strncpy(buffer,fp->buffer,MAX_BUFFER);

   nitems = sscanf(buffer+strlen("M  V30"),
                   "%d %d %d %d %s %s %s %s %s %s %s %s %s %s",
                   &seq, &bp->bond_type,
                   &bp->atoms[0],  &bp->atoms[1],
                   prop[0], prop[1], prop[2], prop[3], prop[4],
                   prop[5], prop[6], prop[7], prop[8], prop[9]);

   if (nitems < 4)
   {
      ShowMessageI("incorrect # (%d) of arguments on bond line",
                   "ReadV30Atom:",
                   nitems);
      ShowMessageS("buffer ==\n%s\n","ReadV30Atom",fp->buffer);
      GetBuffer(fp);
      return(FORTRAN_ERROR);
   }

   if (nitems >= 14)
      ShowMessageI("Warning %d properties for bond.\n",
                   "ReadV30Bond:", (nitems-6));

   /* Processing properties */
   for (i=0; i<nitems-6; i++)
   {
      fprintf(stderr, "ignoring '%s' for bond %d\n", prop[i], seq);
   }
   GetBuffer(fp);
   return(FORTRAN_NORMAL);
}

int ReadREACCSAtom(Fortran_FILE *fp, struct reaccs_atom_t *ap)
{
   int nitems;
   char buffer[MAX_BUFFER+1];
   int charge_radical;

   if (fp->status != FORTRAN_NORMAL) return(fp->status);

   strncpy(buffer,fp->buffer,MAX_BUFFER);

   if (buffer[31] == ' ') buffer[32] = '?'; /* to cope with empty atom types */
   nitems = sscanf(buffer,
                  "%10f%10f%10f %s",
                  &ap->x,  &ap->y,  &ap->z,
                   ap->atom_symbol);
   BlankToZero(buffer+34);
   ap->mass_difference = 0; charge_radical = 0;
   ap->stereo_parity = NONE;
   ap->query_H_count = 0;
   ap->query_stereo_box = 0;
   ap->dummy1 = 0; ap->dummy2 = 0; ap->dummy3 = 0; ap->dummy4 = 0;
   ap->sub_desc = NONE;
   ap->mapping = NONE;
   ap->second_stereo_parity = NONE;
   ap->dummy5 = 0;
   ap->atext[0]='\0';
   sscanf(buffer+34,
      "%2d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d",
      &ap->mass_difference,
      &charge_radical, &ap->stereo_parity,
      &ap->query_H_count,  &ap->query_stereo_box,
      &ap->dummy1, &ap->dummy2, &ap->dummy3, &ap->dummy4,
      &ap->mapping,  &ap->second_stereo_parity,
      &ap->dummy5);

   if (nitems >= 4)
   {
      SplitChargeRadical(charge_radical, &ap->charge, &ap->radical);
      GetBuffer(fp);
      return(FORTRAN_NORMAL);
   }
   else
   {
      ShowMessageI("incorrect # (%d) of arguments on atom line",
                   "ReadREACCSAtom:",
                   nitems);
      ShowMessageS("buffer ==\n%s\n","ReadREACCSAtom",fp->buffer);
      GetBuffer(fp);
      return(FORTRAN_ERROR);
   }
}

/*
 * This section manages the flag that tells the printing functions to
 * strip any trailing zeros on MOL file output to save space.
 */
static int strip_zeroes = FALSE;

int SetStripZeroes(int do_strip)
{
   int result;
   result = strip_zeroes;
   strip_zeroes = do_strip;
   return (result);
}

void PrintMACCSAtom(FILE *fp, struct reaccs_atom_t *ap)
{
   fprintf(fp,"%10.4f%10.4f%10.4f",ap->x,ap->y,ap->z);
   fprintf(fp," %-3s%2d%3d",
              ap->atom_symbol, ap->mass_difference,
              (ANY_CHARGE == CombineChargeRadical(ap->charge,ap->radical)) ?
              0 : CombineChargeRadical(ap->charge,ap->radical));
   if (!strip_zeroes           ||
       ap->stereo_parity != 0  ||
       ap->query_H_count != 0)
      fprintf(fp,"%3d%3d\n", ap->stereo_parity, ap->query_H_count);
   else
      fprintf(fp,"\n");
}

// debugging strings
static char *color_strings[] =
{
   "NO_COLOR",
   "STRATEGIC_COLOR",
   "MODIFIED_COLOR",
   "ALPHA_COLOR",
   "ALLYL_COLOR",
   "CLOSURE_COLOR",
   "RING_COLOR",
   "SINGLE_COLOR",
   "MULTIPLE_COLOR",
   "PROTECTION_COLOR",
   "INTER_COLOR",
   "INTRA_COLOR",
   "REAGENT_COLOR",
   "ADD_REMOVE_COLOR",
   "DOUBLE_COLOR",
};

void PrintREACCSAtom(FILE *fp,
              struct reaccs_atom_t *ap)
{
   int i;
   fprintf(fp,"%10.4f%10.4f%10.4f",ap->x,ap->y,ap->z);
   fprintf(fp," %-3s%2d%3d",
              ap->atom_symbol, ap->mass_difference,
              (ANY_CHARGE == CombineChargeRadical(ap->charge,ap->radical)) ?
               0 : CombineChargeRadical(ap->charge,ap->radical));
   if (!strip_zeroes			||
       ap->stereo_parity != 0		||
       ap->query_H_count != 0		||
       ap->query_stereo_box != 0	||
       ap->dummy1 != 0			||
       ap->dummy2 != 0			||
       ap->dummy3 != 0			||
       ap->dummy4 != 0			||
       ap->dummy5 != 0			||
       ap->mapping != 0			||
       ap->second_stereo_parity != 0)
   {
      fprintf(fp,"%3d%3d%3d",
              ap->stereo_parity,
              ap->query_H_count,  ap->query_stereo_box);
      fprintf(fp,"%3d%3d%3d%3d",
              ap->dummy1, ap->dummy2, ap->dummy3, ap->dummy4);
      fprintf(fp,"%3d%3d%3d",
              ap->mapping, ap->second_stereo_parity, ap->dummy5);
      if (ap->color != NONE  && 0)
      {
          for (i=1; i<14; i++)
              if (ap->color & (1<<(i-1))) fprintf(fp," %s", color_strings[i]);
      }
      fprintf(fp,"\n");
   }
   else
      fprintf(fp,"\n");
}

#define MAX_BONDLINE_FIELDS 7
#define BONDLINE_FIELD_LEN 3

int ReadREACCSBond(Fortran_FILE *fp, struct reaccs_bond_t *bp)
{
   int nitems, i, j, k;
   int bond_line_len, n_chars, pos;
   int *ptrarray[MAX_BONDLINE_FIELDS];
   char c;
   char buffer[BONDLINE_FIELD_LEN+1];

   if (fp->status != FORTRAN_NORMAL) return(fp->status);

   bp->stereo_symbol = 0;
   bp->dummy = 0;
   bp->topography = 0;
   bp->reaction_mark = NONE;
   ptrarray[0] = &bp->atoms[0];
   ptrarray[1] = &bp->atoms[1];
   ptrarray[2] = &bp->bond_type;
   ptrarray[3] = &bp->stereo_symbol;
   ptrarray[4] = &bp->dummy;
   ptrarray[5] = &bp->topography;
   ptrarray[6] = &bp->reaction_mark;
   bond_line_len = strlen(fp->buffer);
   nitems = bond_line_len ? (bond_line_len - 1) / BONDLINE_FIELD_LEN + 1 : 0;
   if (nitems > MAX_BONDLINE_FIELDS)
      nitems = MAX_BONDLINE_FIELDS;
   for (i = 0; i < nitems; ++i)
   {
      pos = i * BONDLINE_FIELD_LEN;
      memset(buffer, 0, BONDLINE_FIELD_LEN + 1);
      n_chars = bond_line_len - pos;
      if (n_chars > BONDLINE_FIELD_LEN)
         n_chars = BONDLINE_FIELD_LEN;
      for (j = 0, k = 0; j < n_chars; ++j)
      {
         c = fp->buffer[pos + j];
         if (c != ' ')
            buffer[k++] = c;
      }
      sscanf(buffer, "%3d", ptrarray[i]);
   }
   if (nitems >= 3)
   {
      GetBuffer(fp);
      return(FORTRAN_NORMAL);
   }
   else
   {
      ShowMessageI("incorrect # (%d) of arguments on bond line",
                   "ReadREACCSBond",
                   nitems);
      ShowMessageS("buffer ==\n%s\n","ReadREACCSBond",buffer);
      return(FORTRAN_ERROR);
   }
}

void PrintMACCSBond(FILE *fp, struct reaccs_bond_t *bp)
{
   int stereo_symbol = bp->stereo_symbol;
   /* CIS_TRANS_SWAPPED is not leagal in a MOL file */
   if (stereo_symbol == CIS_TRANS_SWAPPED) stereo_symbol = CIS_TRANS_EITHER;

   fprintf(fp,"%3d%3d%3d%3d\n",
              bp->atoms[0],   bp->atoms[1],
              bp->bond_type,  stereo_symbol);
}

void PrintREACCSBond(FILE *fp,
                     struct reaccs_bond_t *bp)
{
   int stereo_symbol = bp->stereo_symbol;
   /* CIS_TRANS_SWAPPED is not leagal in a MOL file */
   if (stereo_symbol == CIS_TRANS_SWAPPED) stereo_symbol = CIS_TRANS_EITHER;

   fprintf(fp,"%3d%3d%3d%3d",
              bp->atoms[0],   bp->atoms[1],
              bp->bond_type,  stereo_symbol);
   if (!strip_zeroes		||
       bp->topography != 0	||
       bp->reaction_mark != 0)
      fprintf(fp,"%3d%3d%3d\n",
              bp->dummy,
              bp->topography, bp->reaction_mark);
   else
      fprintf(fp,"\n");
}

struct prop_line_t * ReadProperties(Fortran_FILE *fp,
                                   struct reaccs_molecule_t *mp,
                                   int nprops)
/*
 * Reads the nprops property lines off the file *fp. Interpretable lines
 * are stored in *mp while a pointer to the other lines is returned.
 */
{
   int nentries, n;
   int tmp_ats[8];
   int tmp_vals[8];
   struct prop_line_t *hp, *hhp, *result;
   int atom;
   float value;

   int old_nprops;

   old_nprops = nprops;
   result = (struct prop_line_t *)NULL;
   while (nprops-- > 0)
   {
      if (STRING_BEGINS(fp->buffer,"M  CHG"))    /* charge property */
      {
         sscanf(fp->buffer+strlen("M  CHG"), "%3d", &nentries);
         n = sscanf(fp->buffer+strlen("M  CHG")+3,
                    " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                    &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
                    &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
                    &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
                    &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
         if (n != nentries*2)
         {
            ShowMessageI("n = %d","ReadProperties", n);
            ShowMessageI("nentries = %d","ReadProperties", nentries);
            ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
         }
         for (n=0; n<nentries; n++)
            mp->atom_array[tmp_ats[n]-1].charge = tmp_vals[n];
      }
      else if (STRING_BEGINS(fp->buffer,"M  RAD")) /* radical property */
      {
         sscanf(fp->buffer+strlen("M  RAD"), "%3d", &nentries);
         n = sscanf(fp->buffer+strlen("M  RAD")+3,
                    " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                    &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
                    &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
                    &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
                    &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
         if (n != nentries*2)
         {
            ShowMessageI("n = %d","ReadProperties", n);
            ShowMessageI("nentries = %d","ReadProperties", nentries);
            ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
         }
         for (n=0; n<nentries; n++)
            mp->atom_array[tmp_ats[n]-1].radical = tmp_vals[n];
      }
      else if (STRING_BEGINS(fp->buffer,"M  SUB")) /* subst. degree property */
      {
         sscanf(fp->buffer+strlen("M  SUB"), "%3d", &nentries);
         n = sscanf(fp->buffer+strlen("M  SUB")+3,
                    " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                    &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
                    &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
                    &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
                    &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
         if (n != nentries*2)
         {
            ShowMessageI("n = %d","ReadProperties", n);
            ShowMessageI("nentries = %d","ReadProperties", nentries);
            ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
         }
         for (n=0; n<nentries; n++)
            mp->atom_array[tmp_ats[n]-1].sub_desc = tmp_vals[n];
      }
      else if (STRING_BEGINS(fp->buffer,"V  ")) /* value for PC files */
      {
         n = sscanf(fp->buffer+strlen("V  "), " %d %f", &atom,  &value);
         if (n == 2)
            mp->atom_array[atom-1].value = value;
         else
         {
            ShowMessage("value line ignored","ReadProperties");
            ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
         }
      }
      else if (STRING_BEGINS(fp->buffer,"M  END"))
      {
         if (nprops != 0)
         {
            if (old_nprops != 999)
               ShowMessage("M  END line was not last property line",
                           "ReadProperties");
            nprops = 0;
         }
         // GetBuffer(fp);
         break;
      }
      else if (STRING_BEGINS(fp->buffer,"M  RGP"))      /* R-Group Lines */
      {
         hp = TypeAlloc(1,struct prop_line_t);
         strncpy(hp->text, fp->buffer, MDL_MAXLINE);
         hp->text[MDL_MAXLINE] = '\0';
         hp->next = result; result = hp;
      }
      else if (STRING_BEGINS(fp->buffer,"G  "))       /* Group Lines */
      {
         hp = TypeAlloc(1,struct prop_line_t);
         strncpy(hp->text, fp->buffer, MDL_MAXLINE);
         hp->text[MDL_MAXLINE] = '\0';
         hp->next = result; result = hp;
         GetBuffer(fp); nprops--;
         hp = TypeAlloc(1,struct prop_line_t);
         strncpy(hp->text, fp->buffer, MDL_MAXLINE);
         hp->text[MDL_MAXLINE] = '\0';
         hp->next = result; result = hp;
      }
      else if (STRING_BEGINS(fp->buffer,"A  "))       /* atom text Lines */
      {
         n = sscanf(fp->buffer+strlen("A  "), " %d", &atom);
         if (n==1  &&  atom <= mp->n_atoms  && 0 == strcmp(mp->atom_array[atom-1].atom_symbol, "R"))
         {                                      // interpret atom text as atom label
             GetBuffer(fp); nprops--;
             strncpy(mp->atom_array[atom-1].atext, fp->buffer, MDL_MAXLINE);
             mp->atom_array[atom-1].atext[MDL_MAXLINE] = '\0';
         }
         else                                   // just store as standard property
         {
             hp = TypeAlloc(1,struct prop_line_t);
             strncpy(hp->text, fp->buffer, MDL_MAXLINE);
             hp->text[MDL_MAXLINE] = '\0';
             hp->next = result; result = hp;
             GetBuffer(fp); nprops--;
             hp = TypeAlloc(1,struct prop_line_t);
             strncpy(hp->text, fp->buffer, MDL_MAXLINE);
             hp->text[MDL_MAXLINE] = '\0';
             hp->next = result; result = hp;
          }
      }
      else if (STRING_BEGINS(fp->buffer,"M  ISO"))      /* isotopic labelling */
      {
         /* Just ignore these, because they are taken care of by      */
         /* the mass difference field. It would require a major       */
         /* change to modify the semantics now. Room for improvement! */
         sscanf(fp->buffer+strlen("M  ISO"), "%3d", &nentries);
         n = sscanf(fp->buffer+strlen("M  ISO")+3,
                    " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                    &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
                    &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
                    &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
                    &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
         if (n != nentries*2)
         {
            ShowMessageI("n = %d","ReadProperties", n);
            ShowMessageI("nentries = %d","ReadProperties", nentries);
            ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
         }
         /* Just handle mass difference for substitution point labels and some well-known other atoms, other labelings are captured as property lines */
         for (n=0; n<nentries; n++)
            if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "R"))
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n];
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "H")  &&    1 <= tmp_vals[n]  &&  tmp_vals[n] <=   3)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-1;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Li")  &&    6 <= tmp_vals[n]  &&  tmp_vals[n] <=   7)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-7;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "B")  &&   10 <= tmp_vals[n]  &&  tmp_vals[n] <=  11)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-11;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "C")  &&   12 <= tmp_vals[n]  &&  tmp_vals[n] <=  14)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-12;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "N")  &&   13 <= tmp_vals[n]  &&  tmp_vals[n] <=  15)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-14;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "O")  &&   16 <= tmp_vals[n]  &&  tmp_vals[n] <=  18)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-16;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "F")  &&   18 <= tmp_vals[n]  &&  tmp_vals[n] <=  19)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-19;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "P")  &&   31 <= tmp_vals[n]  &&  tmp_vals[n] <=  33)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-31;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "S")  &&   32 <= tmp_vals[n]  &&  tmp_vals[n] <=  36)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-32;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Cl")  &&   35 <= tmp_vals[n]  &&  tmp_vals[n] <=  37)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-35;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Co")  &&   56 <= tmp_vals[n]  &&  tmp_vals[n] <=  61)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-59;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol,"Br")  &&   79 <= tmp_vals[n]  &&  tmp_vals[n] <=  81)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-80;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "Y")  &&   88 <= tmp_vals[n]  &&  tmp_vals[n] <=  90)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-89;
            else if (0 == strcmp(mp->atom_array[tmp_ats[n]-1].atom_symbol, "I")  &&  123 <= tmp_vals[n]  &&  tmp_vals[n] <= 131)
               mp->atom_array[tmp_ats[n]-1].mass_difference = tmp_vals[n]-127;
            else { // make sure unknown isotopic labelings get represented as properties
                hp = TypeAlloc(1,struct prop_line_t);
                sprintf(hp->text, "M  ISO  1 %3d %3d", tmp_ats[n], tmp_vals[n]);
                hp->text[MDL_MAXLINE] = '\0';
                hp->next = result; result = hp;
            }
      }
      else if (STRING_BEGINS(fp->buffer,"M  ALS"))      /* atom type lists */
      {
         /* Just ignore these, because they are taken care of by      */
         /* the atom list lines read before.			      */
      }
      /* ignore short cut lines */
      else if (STRING_BEGINS(fp->buffer,"M  SLB") ||    // Sgroup Labels
               STRING_BEGINS(fp->buffer,"M  STY") ||    // Sgroup Type
               STRING_BEGINS(fp->buffer,"M  SAL") ||    // Sgroup Atom List
               STRING_BEGINS(fp->buffer,"M  SBL") ||    // Sgroup Bond List
               STRING_BEGINS(fp->buffer,"M  SDI") ||    // Sgroup Display Information, e.g. bracket positions
               STRING_BEGINS(fp->buffer,"M  SMT") ||    // Sgroup Subscript text
               STRING_BEGINS(fp->buffer,"M  SBV") ||    // Abbreviation SGroup bond and vector information
               STRING_BEGINS(fp->buffer,"M  SCL") ||    // Abbreviation Sgroup class
               STRING_BEGINS(fp->buffer,"M  SAP") ||    // Abbreviation Sgroup Attachment Point
               STRING_BEGINS(fp->buffer,"M  SDT") ||    // Data Sgroup Field Description
               STRING_BEGINS(fp->buffer,"M  SDD") ||    // Data Sgroup Display Information
               STRING_BEGINS(fp->buffer,"M  SED") ||    // Data Sgroup Data
               STRING_BEGINS(fp->buffer,"M  SDS"))      // Sgroup Expansion
      {
         /* Just ignore those line because they are display-only. */
      }
      else
      {
         ShowMessage("found unknown property","ReadProperties");
         ShowMessageS("buffer = '%s'","ReadProperties",fp->buffer);
         hp = TypeAlloc(1,struct prop_line_t);
         strncpy(hp->text, fp->buffer, MDL_MAXLINE);
         hp->text[MDL_MAXLINE] = '\0';
         hp->next = result; result = hp;
      }
      GetBuffer(fp);
   }

   /* restore order of list */
   hhp = result; result = (struct prop_line_t *)NULL;
   while (hhp)
   {
      hp = hhp; hhp = hp->next;
      hp->next = result; result = hp;
   }

   if (fp->buffer[0] != 'M'  &&  fp->status == FORTRAN_NORMAL)
   {
      ShowMessage("possible format error","ReadProperties");
      ShowMessageS("current line: '%s'","ReadProperties",fp->buffer);
   }
   else
       GetBuffer(fp);
   return (result);
}

void PrintPropLines(FILE *fp,
                    struct reaccs_molecule_t *mp,
                    struct prop_line_t *prop_lines)
/*
 * Prints the property lines pointed to by prop_lines to the file *fp.
 * Lines derived from the details of the molecule are also added.
 */
{
   int i;
   struct reaccs_atom_t *ap;
   int use_charge_radical;

   while (prop_lines)
   {
      fprintf(fp,"%s\n",prop_lines->text);
      prop_lines = prop_lines->next;
   }

   use_charge_radical = FALSE;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (ap->charge > 3  ||  ap->charge < -3)
         use_charge_radical = TRUE;
      if (ap->radical != NONE)
        use_charge_radical = TRUE;
   }

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (use_charge_radical && ap->charge != 0)
         fprintf(fp, "M  CHG  1 %3d %3d\n", i+1, ap->charge);
      if (ap->radical != NONE)
         fprintf(fp, "M  RAD  1 %3d %3d\n", i+1, ap->radical);
      if (ap->sub_desc != NONE)
         fprintf(fp, "M  SUB  1 %3d %3d\n", i+1, ap->sub_desc);
      if (ap->value != 0)
         fprintf(fp, "V  %3d %g\n", i+1, ap->value);
      if (ap->atext[0] != '\0'  &&  0 == strcmp(ap->atom_symbol, "R"))
      {
         fprintf(fp, "A  %3d\n", i+1);
         fprintf(fp, "%s\n", ap->atext);
      }
   }

   fprintf(fp,"M  END\n");
}

int ReadREACCSMolecule(Fortran_FILE *fp,
                       struct reaccs_molecule_t *mp,
                       const char *label)
/*
 * Reads the next molecule description in Fortran_FILE *fp starting
 * with a line that begins with label[]. The molecule is put into *mp
 * which must be preallocated by the caller.
 * The molecule must start at the current line if label == "".
 */
{
   int nitems;
   int n_stexts;
   struct stext_line_t *stp, *stph;
   int i;
   char buffer[MAX_BUFFER+1];
   int status;
   long regno=0;
   /* additional temps for V3000 format */
   int n_sgroups, n_3d;
   char regtext[20];

   if (IsNULL(mp))
   {
fprintf(stderr, "ReadREACCSMolecule: null molecule pointer\n");
       return(FORTRAN_ERROR);
   }

   // set some safe defaults
   mp->name[0] = '\0';
   mp->scale1 = 0; mp->scale2 = 0; mp->energy = 0; mp->registry_number = 0;
   mp->user_initials[0] = '\0';
   strcpy(mp->program_name,"DUMMY");
   mp->date[0]              = '\0';
   mp->time[0]                = '\0';
   mp->dimensionality[0]= '\0';
   mp->comment[0] = '\0';
   mp->n_atom_lists = 0;
   mp->dummy1       = 0;
   mp->chiral_flag  = 0;
   n_stexts         = 0;
   mp->stext_lines = (struct stext_line_t *)NULL;
   mp->n_props      = 0;
   mp->version[0]   = '\0';
   mp->atom_array = (struct reaccs_atom_t *)NULL;
   mp->bond_array = (struct reaccs_bond_t *)NULL;

   if (IsNULL(fp))
   {
fprintf(stderr, "ReadREACCSMolecule: null file pointer\n");
       return(FORTRAN_ERROR);
   }

   if (IsNULL(label))
   {
fprintf(stderr, "ReadREACCSMolecule: null label pointer\n");
       return(FORTRAN_ERROR);
   }

   if (fp->status != FORTRAN_NORMAL)
   {
// fprintf(stderr, "ReadREACCSMolecule: fp->status(1) = %d\n", fp->status);
       return(fp->status);
   }

   if (0 != strcmp(label,""))
   {
      if (!SearchString(fp,label,"$MFMT")) return (fp->status);
      if (STRING_BEGINS(fp->buffer,"$MFMT $MIREG"))
         sscanf(fp->buffer+strlen("$MFMT $MIREG"),"%ld",&regno);
      else
         regno = 0;
      GetBuffer(fp);                            /* skip header line */
      if (fp->status != FORTRAN_NORMAL)
       {
fprintf(stderr, "ReadREACCSMolecule: fp->status(2) = %d\n", fp->status);
          return(fp->status);
       }
   }

           /* patch for SD-files w/o molecule records */
   if (fp->buffer[0] == '>'  &&  strchr(fp->buffer,'<'))
   {
      mp->atom_array = TypeAlloc(0, struct reaccs_atom_t);
      mp->bond_array = TypeAlloc(0, struct reaccs_bond_t);
      return(fp->status);
   }

   strncpy(mp->name,fp->buffer,MAXNAME);
   mp->name[MAXNAME] = '\0';
   GetBuffer(fp);
   if (fp->status != FORTRAN_NORMAL) return(fp->status);

   strncpy(buffer,fp->buffer,MAX_BUFFER);
   BlankToZero(buffer+22);
   strcpy(mp->dimensionality,"2D");
   mp->user_initials[0] = '\0'; mp->program_name[0] = '\0';
   mp->date[0] = '\0'; mp->time[0] = '\0';
   mp->scale1 = 0; mp->scale2 = 0; mp->energy = 0; mp->registry_number = 0;

   nitems = sscanf(buffer,
                   "%2c%8c%6c%4c%2c%2d%10f%12f%6ld",
                   mp->user_initials,
                   mp->program_name,
                   mp->date, mp->time,
                   mp->dimensionality,
                  &mp->scale1, &mp->scale2,
                  &mp->energy,
                  &mp->registry_number);
   if (regno != 0) mp->registry_number = regno;

   mp->user_initials[2]  = '\0'; RemoveTrailingBlanks(mp->user_initials);
   mp->program_name[8]   = '\0'; RemoveTrailingBlanks(mp->program_name);
   mp->date[6]           = '\0'; RemoveTrailingBlanks(mp->date);
   mp->time[4]           = '\0'; RemoveTrailingBlanks(mp->time);
   mp->dimensionality[2] = '\0'; RemoveTrailingBlanks(mp->dimensionality);

   GetBuffer(fp);
   if (fp->status != FORTRAN_NORMAL) return(fp->status);

   strncpy(mp->comment,fp->buffer,MDL_MAXLINE);
   GetBuffer(fp);
   if (fp->status != FORTRAN_NORMAL) return(fp->status);
   RemoveTrailingBlanks(mp->comment);

   strncpy(buffer,fp->buffer,MAX_BUFFER);
   mp->version[0]   = '\0';
   if (strlen(buffer) > 34)
   {
      strncpy(mp->version, buffer+34, 6);
      mp->version[5] = '\0';
      // fix a common glitch
      if (NULL != strstr(mp->version, "V200")) strcpy(mp->version, "V2000");
   }

   if (0 == strcmp(mp->version, "V3000"))   // new MOL format
   {
       GetBuffer(fp);
       mp->n_atom_lists = 0;
       mp->dummy1       = 0;
       mp->chiral_flag  = 0;
       n_stexts         = 0;
       mp->n_props      = 0;
       mp->n_atoms      = 0;
       mp->n_atoms      = 0;
       while (fp->status == FORTRAN_NORMAL  &&
              0 != strcmp("M  END", fp->buffer))
       {
          if (0 == strncmp("M  V30 BEGIN CTAB", fp->buffer,
                                 strlen("M  V30 BEGIN CTAB")))
             GetBuffer(fp);
          else if (0 == strncmp("M  V30 COUNTS ", fp->buffer,
                                strlen("M  V30 COUNTS ")))
          {
             regtext[0] = '\0';
             nitems = sscanf(fp->buffer + strlen("M  V30 COUNTS "),
                             "%d %d %d %d %d %s",
                                        &mp->n_atoms, &mp->n_bonds,
                                        &n_sgroups,
                                        &n_3d,
                                        &mp->chiral_flag,
                                        &regtext);
             mp->atom_array = TypeAlloc(mp->n_atoms,struct reaccs_atom_t);
             mp->bond_array = TypeAlloc(mp->n_bonds,struct reaccs_bond_t);
             GetBuffer(fp);
          }
          else if (0 == strncmp("M  V30 BEGIN ATOM", fp->buffer,
                                strlen("M  V30 BEGIN ATOM")))
          {
             GetBuffer(fp);
             i=0;
             while (fp->status == FORTRAN_NORMAL  &&
                    0 != strcmp("M  V30 END ATOM", fp->buffer))
             {
                status = ReadV30Atom(fp,&mp->atom_array[i++], mp);
                if (status != FORTRAN_NORMAL) return(status);
             }
             if (fp->status == FORTRAN_NORMAL)
                 GetBuffer(fp); // skip 'M  V30 END ATOM'
          }
          else if (0 == strncmp("M  V30 BEGIN BOND", fp->buffer,
                                strlen("M  V30 BEGIN BOND")))
          {
             GetBuffer(fp);
             i=0;
             while (fp->status == FORTRAN_NORMAL  &&
                    0 != strcmp("M  V30 END BOND", fp->buffer))
             {
                status = ReadV30Bond(fp,&mp->bond_array[i++]);
                if (status != FORTRAN_NORMAL) return(status);
             }
             if (fp->status == FORTRAN_NORMAL)
                 GetBuffer(fp); // skip 'M  V30 END BOND'
          }
          else if (0 == strncmp("M  V30 BEGIN COLLECTION", fp->buffer,
                                strlen("M  V30 BEGIN COLLECTION")))
          {
             while (fp->status == FORTRAN_NORMAL  &&
                    0 != strcmp("M  V30 END COLLECTION", fp->buffer))
             {
                 GetBuffer(fp);
             }
             if (fp->status == FORTRAN_NORMAL)
                 GetBuffer(fp); // skip 'M  V30 END COLLECTION'
          }
          else if (0 == strncmp("M  V30 END CTAB", fp->buffer,
                                 strlen("M  V30 END CTAB")))
             GetBuffer(fp);
          else
          {
// fprintf(stderr, "skipping line:\t%s\n", fp->buffer);
             GetBuffer(fp);
          }
       }
       if (fp->status != FORTRAN_NORMAL) return(fp->status);
       // parsing was in V3000 but still using old representation
       strcpy(mp->version, "V2000");
   }
   else                                     // old MOL format
   {
      BlankToZero(buffer);
       mp->n_atom_lists = 0;
       mp->dummy1       = 0;
       mp->chiral_flag  = 0;
       n_stexts         = 0;
       mp->n_props      = 0;
       mp->n_atoms      = 0;
       mp->n_atoms      = 0;
       nitems = sscanf(buffer,"%3d%3d%3d%3d%3d%3d%*12c%3d",
                               &mp->n_atoms, &mp->n_bonds,
                               &mp->n_atom_lists,
                               &mp->dummy1,
                               &mp->chiral_flag,
                               &n_stexts,
                               &mp->n_props);
       if (mp->version[0] == '\0') mp->version[0] = ' ';

       if (nitems >= 2  &&
           mp->n_atoms >=0  &&  mp->n_atoms <= MAXATOMS &&
           mp->n_bonds >=0  &&  mp->n_bonds <= MAXBONDS)
       {
          GetBuffer(fp);
          if (fp->status != FORTRAN_NORMAL) return(fp->status);
       }
       else
       {
          ShowMessage("incorrect syntax of arguments on atom/bond number line\n", "ReadREACCSMolecule");
          fprintf(stderr,"buffer = '%s'\n", fp->buffer);
          return(FORTRAN_ERROR);
       }

       mp->atom_array = TypeAlloc(mp->n_atoms,struct reaccs_atom_t);
       for (i=0; i<mp->n_atoms; i++)
       {
          status = ReadREACCSAtom(fp,&mp->atom_array[i]);
          if (status != FORTRAN_NORMAL) return(status);
       }

       mp->bond_array = TypeAlloc(mp->n_bonds,struct reaccs_bond_t);
       for (i=0; i<mp->n_bonds; i++)
       {
          status = ReadREACCSBond(fp,&mp->bond_array[i]);
          if (status != FORTRAN_NORMAL) return(status);
       }

       mp->symbol_lists = ReadSymbolLists(fp,mp->n_atom_lists);

       mp->stext_lines = (struct stext_line_t *)NULL;
       for (i=0; i<n_stexts; i++)
       {
          stp = TypeAlloc(1,struct stext_line_t);
          sscanf(fp->buffer,"%f %f",&stp->x, &stp->y);
          GetBuffer(fp);
          strncpy(stp->text,fp->buffer,MDL_MAXLINE);
          GetBuffer(fp);
          stp->next = mp->stext_lines;
          mp->stext_lines = stp;
       }
       stp = mp->stext_lines;
       mp->stext_lines = (struct stext_line_t *)NULL;
       while (stp)
       {
          stph = stp->next;
          stp->next = mp->stext_lines;
          mp->stext_lines = stp;
          stp = stph;
       }

       mp->prop_lines = ReadProperties(fp,mp,mp->n_props);
       mp->n_props = CountPropLines(mp->prop_lines);
   }

   return(FORTRAN_NORMAL);
}

struct data_line_t *ReadMACCSDataLines(Fortran_FILE *fp)
/*
 * Reads data lines from FORTRAN file *fp, constructs a data list,
 * and returns this list.
 */
{
   struct data_line_t  head;
   struct data_line_t *tailp;
   struct data_line_t *hp;

   head.next = (struct data_line_t *)NULL;
   tailp = &head;
   while(fp->status == FORTRAN_NORMAL  &&
         !STRING_BEGINS(fp->buffer,"$$$$"))
   {
      hp = TypeAlloc(1,struct data_line_t);
      strncpy(hp->data,fp->buffer,MAXDATA);
      hp->next = (struct data_line_t *)NULL;
      tailp->next = hp;
      tailp = hp;
      GetBuffer(fp);
   }
   return (head.next);
}

struct data_line_t *ReadREACCSDataLines(Fortran_FILE *fp)
/*
 * Reads data lines from FORTRAN file *fp, constructs a data list,
 * and returns this list.
 */
{
   struct data_line_t  head;
   struct data_line_t *tailp;
   struct data_line_t *hp;

   head.next = (struct data_line_t *)NULL;
   tailp = &head;
   while(fp->status == FORTRAN_NORMAL  &&
     (STRING_BEGINS(fp->buffer,"$DTYPE") ||
   STRING_BEGINS(fp->buffer,"$DATUM")))
   {
      hp = TypeAlloc(1,struct data_line_t);
      strncpy(hp->data,fp->buffer,MAXDATA);
      hp->next = (struct data_line_t *)NULL;
      tailp->next = hp;
      tailp = hp;
      GetBuffer(fp);
   }
   return (head.next);
}

int NeededPropLines(struct reaccs_molecule_t *mp)
/*
 * Computes how many property lines are needed to represent:
 * the float values associated with the atoms,
 * the uncommon charges,
 * the radical descriptions, and
 * the general property lines.
 */
{
   int i, result;
   struct reaccs_atom_t *ap;
   int use_charge_radical;

   use_charge_radical = FALSE;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (ap->charge > 3  ||  ap->charge < -3)
         use_charge_radical = TRUE;
      if (ap->radical != NONE)
         use_charge_radical = TRUE;
   }

   for (i=result=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (use_charge_radical && ap->charge != 0) result++;
      if (ap->radical != NONE)                   result++;
      if (ap->sub_desc != NONE)                  result++;
      if (ap->value != 0)                        result++;
      // atom text properties
      if (ap->atom_symbol[0] == 'R'  &&  ap->atom_symbol[1] == '\0'  &&  ap->atext[0] != '\0') result+=2;
   }

   return (result+mp->n_props);
}

void PrintMACCSMolecule(FILE *fp,
                        struct reaccs_molecule_t *mp,
                        char *header)
{
   int i;
   struct stext_line_t *stp;

   if (header[0] != '\0') fprintf(fp,"%s\n",header);

   fprintf(fp,"%s\n",mp->name);

   fprintf(fp,"%-2s%-8s%-6s%-4s%-2s%2d%10.5f%12.5f%6ld\n",
              mp->user_initials,
              mp->program_name,
              mp->date, mp->time,
              mp->dimensionality,
              mp->scale1, mp->scale2,
              mp->energy,
              mp->registry_number);

   fprintf(fp,"%s\n",mp->comment);

   fprintf(fp,"%3d%3d%3d%3d%3d%3d            %3d %-6s\n",
           mp->n_atoms, mp->n_bonds,
           mp->n_atom_lists,
           mp->dummy1,
           mp->chiral_flag,
           CountSTextLines(mp->stext_lines),
           NeededPropLines(mp)+
           1,                        /* 'M  END' line */
           mp->version);

   for (i=0; i<mp->n_atoms; i++)
      PrintMACCSAtom(fp,&mp->atom_array[i]);

   for (i=0; i<mp->n_bonds; i++)
      PrintMACCSBond(fp,&mp->bond_array[i]);

   PrintSymbolLists(fp, mp->symbol_lists);

   for (stp=mp->stext_lines; stp; stp=stp->next)
   {
      fprintf(fp,"%10.4f%10.4f\n",stp->x,stp->y);
      fprintf(fp,"%s\n",stp->text);
   }

   PrintPropLines(fp, mp, mp->prop_lines);
}

void PrintREACCSMolecule(FILE *fp,
                         struct reaccs_molecule_t *mp,
                         const char *header)
{
   int i;
   struct stext_line_t *stp;

   if (header[0] != '\0') fprintf(fp,"%s\n",header);
   fprintf(fp,"%s\n",mp->name);

   fprintf(fp,"%-2s%-8s%-6s%-4s%-2s%2d%10.5f%12.5f%6ld\n",
       mp->user_initials,
           mp->program_name,
           mp->date, mp->time,
           mp->dimensionality,
           mp->scale1, mp->scale2,
           mp->energy,
           mp->registry_number);

   fprintf(fp,"%s\n",mp->comment);

   fprintf(fp,"%3d%3d%3d%3d%3d%3d            %3d %-6s\n",
           mp->n_atoms, mp->n_bonds,
           mp->n_atom_lists,
           mp->dummy1,
           mp->chiral_flag,
           CountSTextLines(mp->stext_lines),
           NeededPropLines(mp)+
           1,                        /* 'M  END' line */
           mp->version);

   for (i=0; i<mp->n_atoms; i++)
   {
      PrintREACCSAtom(fp,&mp->atom_array[i]);
   }

   for (i=0; i<mp->n_bonds; i++)
      PrintREACCSBond(fp,&mp->bond_array[i]);

   PrintSymbolLists(fp, mp->symbol_lists);

   for (stp=mp->stext_lines; stp; stp=stp->next)
   {
      fprintf(fp,"%10.4f%10.4f\n",stp->x,stp->y);
      fprintf(fp,"%s\n",stp->text);
   }

   PrintPropLines(fp, mp, mp->prop_lines);
}

struct reaccs_reaction_t *ReadREACCSReaction(Fortran_FILE *fp)
{
   int nitems;
   int i;
   long regno;
   struct reaccs_reaction_t *rp;
   struct reaccs_molecule_t *mp, *mph;

   regno = 0;

   while ((fp->status != FORTRAN_EOF) &&  /* search for header line */
      !STRING_BEGINS(fp->buffer,"$RXN"))
   {
      if (STRING_BEGINS(fp->buffer,"$RFMT $RIREG"))
         sscanf(fp->buffer+strlen("$RFMT $RIREG"),"%ld",&regno);
      GetBuffer(fp);
   }
   if (fp->status != FORTRAN_NORMAL)
      return((struct reaccs_reaction_t *)NULL);

   GetBuffer(fp);                         /* skip header line */
   if (fp->status != FORTRAN_NORMAL)
      return((struct reaccs_reaction_t *)NULL);

   rp = TypeAlloc(1,struct reaccs_reaction_t);
   rp->n_reactants = rp->n_products = 0;
   rp->reactants = rp->products = (struct reaccs_molecule_t *)NULL;
   rp->next = (struct reaccs_reaction_t *)NULL;

   strncpy(rp->name,fp->buffer,MAXNAME);
   rp->name[MAXNAME] = '\0';
   GetBuffer(fp);
   if (fp->status != FORTRAN_NORMAL)
   {
      FreeReaction(rp);
      return((struct reaccs_reaction_t *)NULL);
   }

   rp->user_initials[0] = '\0'; rp->program_name[0] = '\0'; rp->date[0] = '\0'; rp->time[0] = '\0';
   nitems = sscanf(fp->buffer,
                   "%4c%*2c%8c%*2c%6c%4c%*2c%6ld",
                   rp->user_initials,
                   rp->program_name,
                   rp->date, rp->time,
                  &rp->registry_number);

   /* There is a bug in the definitions of registry numbers. */
   /* Therefore we always take the one from the $RFMT line!! */
   rp->registry_number = regno;

   rp->user_initials[4] = '\0'; RemoveTrailingBlanks(rp->user_initials);
   rp->program_name[8]  = '\0'; RemoveTrailingBlanks(rp->program_name);
   rp->date[6]          = '\0'; RemoveTrailingBlanks(rp->date);
   rp->time[4]          = '\0'; RemoveTrailingBlanks(rp->time);

   if (nitems >= 2)
   {
      GetBuffer(fp);
      if (fp->status != FORTRAN_NORMAL)
      {
         FreeReaction(rp);
         return((struct reaccs_reaction_t *)NULL);
      }
   }
   else
   {
      ShowMessageI("incorrect # (%d) of arguments on reaction-id line\n",
                   "ReadREACCSReaction(1)",
                   nitems);
      ShowMessageS("buffer = '%s'","ReadREACCSReaction",fp->buffer);
      return((struct reaccs_reaction_t *)NULL);
   }

   strncpy(rp->comment,fp->buffer,MDL_MAXLINE); GetBuffer(fp);
   if (fp->status != FORTRAN_NORMAL)
   {
      FreeReaction(rp);
      return((struct reaccs_reaction_t *)NULL);
   }
   RemoveTrailingBlanks(rp->comment);

   nitems = sscanf(fp->buffer,
              "%3d%3d",
               &rp->n_reactants, &rp->n_products);

   if (nitems == 2)
   {
      GetBuffer(fp);
      if (fp->status != FORTRAN_NORMAL)
      {
         FreeReaction(rp);
         return((struct reaccs_reaction_t *)NULL);
      }
   }
   else
   {
      ShowMessageI("incorrect # (%d) of arguments on reactant/product line\n",
                   "ReadREACCSReaction(2)",
                   nitems);
      ShowMessageS("buffer = '%s'","ReadREACCSReaction",fp->buffer);
      if (fp->status != FORTRAN_NORMAL)
      {
         FreeReaction(rp);
         return((struct reaccs_reaction_t *)NULL);
      }
   }

   for (i=0; i<rp->n_reactants; i++)
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (IsNULL(mp))
      {
         ShowMessage("not enough storage!", "ReadREACCSReaction");
         exit(1);
      }
      mp->next = rp->reactants; rp->reactants = mp;
      if (ReadREACCSMolecule(fp,mp,"$MOL") != FORTRAN_NORMAL)
      {
         FreeReaction(rp);
         return((struct reaccs_reaction_t *)NULL);
      }
   }
   mp = (struct reaccs_molecule_t *)NULL; /* fix order of reactants */
   while (!IsNULL(rp->reactants))
   {
      mph = rp->reactants; rp->reactants = mph->next;
      mph->next = mp; mp = mph;
   }
   rp->reactants = mp;

   for (i=0; i<rp->n_products; i++)
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (IsNULL(mp))
      {
         ShowMessage("not enough storage!", "ReadREACCSReaction");
         exit(1);
      }
      mp->next = rp->products; rp->products = mp;
      if (ReadREACCSMolecule(fp,mp,"$MOL") != FORTRAN_NORMAL)
      {
         FreeReaction(rp);
         return((struct reaccs_reaction_t *)NULL);
      }
   }
   mp = (struct reaccs_molecule_t *)NULL; /* fix order of products */
   while (!IsNULL(rp->products))
   {
      mph = rp->products; rp->products = mph->next;
      mph->next = mp; mp = mph;
   }
   rp->products = mp;

   return(rp);
}

void PrintREACCSReaction(FILE *fp, struct reaccs_reaction_t *rp)
{
   struct reaccs_molecule_t *mp;

   fprintf(fp,"$RXN\n");
   fprintf(fp,"%s\n",rp->name);

   fprintf(fp,"%-4s  %-8s  %-6s%-4s  %6ld\n",
              rp->user_initials,
              rp->program_name,
              rp->date, rp->time,
              rp->registry_number);

   fprintf(fp,"%s\n",rp->comment);
   fprintf(fp,"%3d%3d\n",rp->n_reactants,rp->n_products);

   for (mp=rp->reactants; !IsNULL(mp); mp=mp->next)
      PrintREACCSMolecule(fp,mp,"$MOL");

   for (mp=rp->products; !IsNULL(mp); mp=mp->next)
      PrintREACCSMolecule(fp,mp,"$MOL");
}

struct reaccs_molecule_t *ReadReactionAgents(Fortran_FILE *fp)
/*
 * Reads the reagtion agents (solvents and catalysts) for a
 * reaction entry. It returns a list of those molecules.
 */
{
   struct reaccs_molecule_t *result, *mhp;

   result = (struct reaccs_molecule_t *)NULL;

   while (SearchString(fp,"$DATUM $MFMT","$RFMT"))
   {
      mhp = TypeAlloc(1,struct reaccs_molecule_t);
      ReadREACCSMolecule(fp,mhp,"$DATUM $MFMT");
      mhp->next = result; result = mhp;
   }
   return (result);
}

#define MAXLINE 300

static char * _ReadFile(FILE *fp)
/*
 * Read the file into a string
 * the returned string should be free()ed if not NULL
 */
{
   char * Buf;
   char Line[MAXLINE+1];
   int  Size;
   int  Len;

   Buf = (char *)malloc(1);
   Size = 1;
   if (Buf == NULL)
   {
      fprintf(stderr,"Error allocating memory (%d)\n", Size);
      abort();
   }

   *Buf = '\0';

   while (fgets(Line,MAXLINE,fp) != NULL)
   {
      Len = strlen(Line);
      if (!Len) break;
      Size += Len;
      Buf = (char *)realloc(Buf,Size);
      if (Buf == NULL)
      {
         fprintf(stderr,"Error allocating memory (%d)\n", Size);
         abort();
      }
      strcat(Buf,Line);
   }

   return Buf;
}



struct reaccs_molecule_t * MolStr2Mol(char * MolStr)
/*
 * Convert the MolFile in MolStr to a molecule struct
 * this should be deallocated with: FreeMolecule()
 * returns NULL on ERROR
 */
{
   Fortran_FILE *fp;
   struct reaccs_molecule_t *mp;

   fp = FortranStringOpen(MolStr);
   mp = TypeAlloc(1, struct reaccs_molecule_t);
   if (ReadREACCSMolecule(fp, mp, "") != FORTRAN_NORMAL)
   {
      FreeMolecule(mp);
      mp = NULL;
   }
   FortranClose(fp);

   return mp;
}

char * MolToMolStr(struct reaccs_molecule_t * mp)
/*
 * Convert a molecule struct to a char * MolFile
 * the returned string should be free()ed
 * In case of problems NULL will be returned and an ErrorMessage appended
 * to the message list.
 */
{
   FILE *fp;
   char * MolStr;
#ifndef _WIN32
   size_t bufsize;
#else
   wchar_t *tmpdir = NULL;
   int guid_string_buf_size = 0;
   size_t tmpfilename_buf_size = 0;
   wchar_t *tmpfilename_buf = NULL;
   int guid_string_len = 0;
   int swprintf_res;
   int remove_res;
   const int GUID_STRING_SIZE_INCR = 64;
   const int GUID_STRING_MAX_SIZE = GUID_STRING_SIZE_INCR * 10;
   const wchar_t *AVATMP_PREFIX = L"avatmp_";
   GUID guid;
   wchar_t *guid_string = NULL;
   HRESULT guid_res;
#endif

   if (IsNULL(mp)) return NULL;
#ifndef _WIN32
   bufsize = 5*80 + mp->n_atoms*80 + mp->n_bonds*80 + mp->n_props*80 + 1;
   fp = fmemopen(NULL, bufsize, "w+");
   /* File could not be created => log an error and return NULL */
   if (IsNULL(fp))
   {
      sprintf(msg_buffer, "Could not open memory mapped file of size %d", bufsize);
      AddMsgToList(msg_buffer);
      return NULL;
   }
#else
   tmpdir = _wgetenv(L"TMP");
   if (!tmpdir)
   {
      AddMsgToList("Could not retrieve TMP environment variable");
      return NULL;
   }
   guid_res = CoCreateGuid(&guid);
   if (guid_res != S_OK)
   {
      sprintf(msg_buffer, "Could not create tempfile GUID; error code: %d", guid_res);
      AddMsgToList(msg_buffer);
      return NULL;
   }
   while (guid_string_buf_size < GUID_STRING_MAX_SIZE
      && !(guid_string_len = StringFromGUID2(&guid, guid_string, guid_string_buf_size)))
   {
      guid_string_buf_size += GUID_STRING_SIZE_INCR;
      guid_string = (wchar_t *)realloc(guid_string, guid_string_buf_size * sizeof(wchar_t));
      if (!guid_string)
      {
         AddMsgToList("Could not allocate guid_string");
         return NULL;
      }
   }
   if (!guid_string_len)
   {
      AddMsgToList("Could not convert GUID to string");
      return NULL;
   }
   tmpfilename_buf_size = guid_string_len + wcslen(tmpdir) + wcslen(AVATMP_PREFIX) + 1;
   tmpfilename_buf = (wchar_t *)malloc((tmpfilename_buf_size + 2) * sizeof(wchar_t));
   if (!tmpfilename_buf)
   {
      free(guid_string);
      guid_string = NULL;
      AddMsgToList("Could not allocate tmpfilename_buf");
      return NULL;
   }
   swprintf_res = swprintf_s(tmpfilename_buf, tmpfilename_buf_size, L"%s\\%s%s", tmpdir, AVATMP_PREFIX, guid_string);
   free(guid_string);
   guid_string = NULL;
   if (swprintf_res == -1)
   {
      free(tmpfilename_buf);
      tmpfilename_buf = NULL;
      AddMsgToList("Could not create temporary filename");
      return NULL;
   }
   fp = _wfopen(tmpfilename_buf, L"w+");
   /* File could not be created => log an error and return NULL */
   if (IsNULL(fp))
   {
      free(tmpfilename_buf);
      tmpfilename_buf = NULL;
      AddMsgToList("Could not open temporary file");
      return NULL;
   }
#endif
   PrintREACCSMolecule(fp, mp,"");

   fputc('\0', fp);
   fflush(fp);
   rewind(fp);

   MolStr = _ReadFile(fp);
   fclose(fp);

#ifdef _WIN32
   remove_res = _wremove(tmpfilename_buf);
   free(tmpfilename_buf);
   tmpfilename_buf = NULL;
   if (remove_res == -1)
   {
      AddMsgToList("Could not delete temporary file");
      return NULL;
   }
#endif

   if (MolStr == NULL)
      AddMsgToList("PrintREACCSMolecule did return NULL");

   return MolStr;
}

char * MolToMolStrOld(struct reaccs_molecule_t * mp)
/*
 * Convert a molecule struct to a char * MolFile
 * the returned string should be free()ed
 * In case of problems NULL will be returned and an ErrorMessage appended
 * to the message list.
 */
{
   FILE *fp;
   const char *tempdir;
   char *tempfile;
   int idir;
   char * MolStr;

   tempfile = NULL;
   fp = tmpfile();
   /* File could not be created => log an error and return NULL */
   if (IsNULL(fp))
   {
      sprintf(msg_buffer,"Error opening tmpfile() for writing");
      AddMsgToList(msg_buffer);
      return NULL;
   }

   PrintREACCSMolecule(fp, mp,"");

   fputc('\0', fp);
   fflush(fp);
   rewind(fp);

   MolStr = _ReadFile(fp);
   fclose(fp);
   if (!IsNULL(tempfile))   // tmpfile() did work => remove the file after use
   {
      remove(tempfile);
      MyFree((char *)tempfile);
   }

   if (MolStr == NULL)
      AddMsgToList("PrintREACCSMolecule did return NULL");

   return MolStr;
}

struct symbol_list_t *ParseV30SymbolList(char *symbol,
                                         int iatom,
                                         struct reaccs_molecule_t *mp,
                                         struct symbol_list_t *old_list)
/*
 * Parses the contents of symbol into a new symbol_list_t structure and
 * prepends it to old_list. atom_symbol is set to 'L' if parsing was OK
 * and to 'Unk' otherwise.
 */
{
    struct symbol_list_t *list;
    char* bra_from;
    char* bra_to;
    int len;
    bra_from = strchr(symbol, '[');
    bra_to = strchr(symbol, ']');
    if (bra_from == NULL  ||  bra_to == NULL  ||  bra_to < bra_from)
    {
        list = old_list;
        strcpy(mp->atom_array[iatom].atom_symbol, "Unk");
        fprintf(stderr, "ParseV30SymbolList: Could not parse symbol '%s'\n",
                symbol);
    }
    else
    {
        len = (bra_to-bra_from)-1;
        list = TypeAlloc(1,struct symbol_list_t);
        list->next = old_list;
        list->atom = iatom+1;
        strncpy(list->string, bra_from+1, len);
        list->string[len] = '\0';
        if (0 == strncmp(symbol,"NOT",3)) list->logic = EXCLUSIVE;
        else                              list->logic = INCLUSIVE;
        strcpy(mp->atom_array[iatom].atom_symbol, "L");
    }
    return list;
}
