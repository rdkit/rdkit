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

/**************************************************************************/
/*                                                                        */
/*    File:           smi2mol.c                                           */
/*                                                                        */
/*    Purpose:        This file implements the functions needed to        */
/*                    convert a SMILES string into the corresponding      */
/*                    MDL molecule data structure. It uses layout.c       */
/*                    to compute coordinates and denormal.c to fix        */
/*                    bond normalization. It also implements conversion   */
/*                    from MOL format to SMILES and various other         */
/*                    related conversion, when called as an executable.   */
/*                                                                        */
/*    History:        11-Jul-1994     Start of development.               */
/*                    14-Feb-2020     Moved main method to smi2mol_main.c */
/*                                                                        */
/**************************************************************************/
#include "smi2mol.h"

#include <stdlib.h>

#include <string.h>
#include <ctype.h>

#include "denormal.h"
#include "graph.h"
#include "reaccs.h"
#include "rtutils.h"
#include "shortcut.h"
#include "stereo.h"
#include "symbol_lists.h"
#include "utilities.h"

#include "perceive.h"
#include "layout.h"

#define DY_ABOVE_FIRST        1
#define DY_BELOW_FIRST        2
#define DY_ABOVE_SECOND       4
#define DY_BELOW_SECOND       8

/*
 * This function changes a DY_ABOVE* color to the corresponding DY_BELOW*
 * color and vice-versa.
 * The coloring algorithm is used in bond stereo perception.
 */
int SwapFirstSecond(int color)
{
   int new_color=0;

   new_color = color;

// if (color != new_color) fprintf(stderr, "swapped color %d to %d\n", color, new_color);
   return (new_color);
}

#define DY_ABOVE_BOTH        (DY_ABOVE_FIRST|DY_ABOVE_SECOND)
#define DY_BELOW_BOTH        (DY_BELOW_FIRST|DY_BELOW_SECOND)

// #define MAXOPEN 999
#define MAXOPEN 99

#define VISITED  1

static char *StrAppend(char *result, char *cp, char *sp)
/*
 * Appends the string pointed to by sp and cp to result. The function
 * deallocates the space of the cp-string and reallocates result
 * to contain the new characters.
 * The string *sp is NOT deallocated.
 */
{
   int cplen, splen;

   // nothing to add
   if (!cp && !sp)
   {
      return (result);
   }

   if (cp) cplen = strlen(cp);
   else    cplen = 0;

   if (sp) splen = strlen(sp);
   else    splen = 0;

   if (result)   /* new fragment */
   {
      result = (char *) realloc((void *) result,
                                strlen(result)+  /* old string */
                                splen+           /* dot */
                                cplen+           /* new string */
                                1);              /* NUL */
      if (sp) strcat(result, sp);
      if (cp)
      {
         strcat(result, cp);
         MyFree(cp);
      }
   }
   else if (cp != NULL  &&  sp == NULL)
      result = cp;
   else
   {
      result = TypeAlloc(splen+cplen+1, char);
      if (sp) strcpy(result, sp);
      if (cp)
      {
         strcat(result, cp);
         MyFree(cp);
      }
   }

   return (result);
}

int PerceiveSMILESAtomStereo(struct reaccs_molecule_t *mp,
                             int                       icurr,
                             int                       lig_ind[4],
                             int                       nlig,
                             neighbourhood_t          *nbp,
                             int use_z)
/*
* Compute the stereo symbol to be placed at the atom with index icurr in
* *mp. The previous SMILES atom has index iprev and the neighbourhood
* of the current atom.
 *
 * lig_ind[0..nlig-1] contain the list of neighbour atoms as
 * cited in the emerging SMILES string.
 *
*/
{
   struct npoint_t tetra[4];
   struct reaccs_bond_t *bp;
   int i, j;
   double vol;

// fprintf(stderr,"%d: nlig = %d\n", icurr+1, nlig);
   /* non-stereogenic center? */
   if (nbp->n_ligands >  4  ||
       nbp->n_ligands <  3  ||
       nlig           != 4) return (0);

   /* no stereo bond attached? */
   for (i=0; i<nbp->n_ligands; i++)
   {
      bp = &mp->bond_array[nbp->bonds[i]];
      if (bp->atoms[0] == icurr+1  &&
          (bp->stereo_symbol == UP  ||
           bp->stereo_symbol == DOWN  ||
           (use_z  &&  mp->atom_array[bp->atoms[1]-1].z  != 0.0)))
      break;
   }
   if (i == nbp->n_ligands) return (0);

   for (i=0; i<nlig; i++)     /* set up tetrahedron */
      if (lig_ind[i] == icurr) /* implicit H */
      {
         tetra[i].x = mp->atom_array[icurr].x;
         tetra[i].y = mp->atom_array[icurr].y;
         tetra[i].z = 0;
// fprintf(stderr,"(%.1g, %.1g, %.1g) ", tetra[i].x, tetra[i].y, tetra[i].z);
      }
      else
      {
         tetra[i].x = mp->atom_array[lig_ind[i]].x;
         tetra[i].y = mp->atom_array[lig_ind[i]].y;
         for (j=0; j<nbp->n_ligands; j++)
            if (nbp->atoms[j] == lig_ind[i])
            {
               bp = &mp->bond_array[nbp->bonds[j]];
               if (bp->atoms[0] != icurr+1)
                  tetra[i].z = 0;
               else if (bp->stereo_symbol == UP)
                  tetra[i].z = 1;
               else if (bp->stereo_symbol == DOWN)
                  tetra[i].z = (-1);
               else if (use_z)
                  tetra[i].z = mp->atom_array[lig_ind[i]].z;
               else
                  tetra[i].z = 0;
// fprintf(stderr,"(%.1g, %.1g, %.1g) ", tetra[i].x, tetra[i].y, tetra[i].z);
               break;
            }
         if (j == nbp->n_ligands)
            ShowMessage("could not find all ligands from lig_ind[]",
                        "PerceiveSMILESAtomStereo");
      }

   vol = Volume(tetra);
// fprintf(stderr,"vol = %.2g\n", vol);
   if      (vol < -0.01) return (-1);
   else if (vol >  0.01) return (1);
   else                  return (0);
}

char *SmilesAtom(struct reaccs_molecule_t *mp,
                 int index,
                 int lig_ind[4], int nlig,
                 neighbourhood_t *nbp,
                 int H_count[],
                 int expect_smarts,
                 int aromatic_smiles,
                 int use_z)
/*
 * Returns a pointer to a string which describes the atom at index
 * in the molecule *mp in SMILES notation.
 * The space is allocated from the heap.
 *
 * lig_ind[0..nlig-1] contain the list of neighbour atoms as
 * cited in the emerging SMILES string.
 *
 * nbp[] describes the neighbourhood of the current atom.
 *
 * The function assumes that hydrogen counts have been perceived
 * and stored in H_count[1..mp->n_atoms] (note: index origin == 1!!!)
 */
{
   char *result, buffer[60], tmp[10];
   struct reaccs_atom_t *ap;
   struct symbol_list_t *slp;
   int i;
   int rgnum;
   int stereo;
   int start_of_symbol;
   char *cp;
   int aromatic_level;
   int exo_methyl;

   struct prop_line_t *pl;
   char *atom_text;
   int n;

   /* find relevant atom text, if any */
   atom_text = (char *)NULL;
   for (pl = mp->prop_lines; pl; pl=pl->next)
   {
      if (STRING_BEGINS(pl->text, "A  "))
      {
         n = 0;
         sscanf(pl->text+3, "%d", &n);
         if (n == index+1) atom_text = pl->next->text;
      }
   }

   ap = mp->atom_array+index;

   stereo = PerceiveSMILESAtomStereo(mp, index, lig_ind, nlig, nbp, use_z);

   if (AtomSymbolMatch(ap->atom_symbol, "A,B,C,N,O,P,S,F,Cl,Br,I,R")  &&
       ap->value == 0.0                                             &&
       ap->charge == NONE                                           &&
       ap->radical == NONE                                          &&
       ap->mass_difference == NONE                                  &&
       ap->query_H_count == NONE                                    &&
       (ap->mapping == NONE || 0 == strcmp(ap->atom_symbol,"R"))    &&
       (ap->sub_desc == NONE  ||  !expect_smarts)                   &&
       atom_text == (char *)NULL                                    &&
       !(AtomSymbolMatch(ap->atom_symbol, "F,Cl,Br,I") && nlig > 1) &&  // shortcut atom symbol not used for non-trivial halogenes
       stereo == 0)
   {                                                   /* Easy cases */
      if (ap->atom_symbol[0] == 'R'  || /* use R-atoms for '*'s in SMILES */
          ap->atom_symbol[0] == 'A')    /* convert A atoms to SMARTS symbol */
      {
         if (ap->atom_symbol[0] == 'R'  &&  ap->atom_symbol[1] == '\0'  &&  ap->atext[0] != '\0')
         {
             result = TypeAlloc(strlen(ap->atext)+3, char);
             strcpy(result, "{"); strcat(result, ap->atext); strcat(result, "}");                
         }
         else
         {
             result = TypeAlloc(strlen(ap->atom_symbol)+1, char);
             strcpy(result, "*");                
         }
      }
      else
      {
         result = TypeAlloc(strlen(ap->atom_symbol)+1, char);
         strcpy(result, ap->atom_symbol);
      }

      /* check for aromatic atom and patch result[0] accordingly */
      if (AtomSymbolMatch(ap->atom_symbol, "C,N,O,P,S,B"))
      {
         aromatic_level = 0;
         exo_methyl = FALSE;
         for (i=0; i<nbp->n_ligands; i++)
         {
            if (mp->bond_array[nbp->bonds[i]].bond_type == SINGLE)
               exo_methyl = TRUE;
            if (mp->bond_array[nbp->bonds[i]].bond_type == AROMATIC)
               aromatic_level |= 1;
            else if (mp->bond_array[nbp->bonds[i]].bond_type == SINGLE_DOUBLE)
               aromatic_level |= 2;
            else if (mp->bond_array[nbp->bonds[i]].bond_type == SINGLE_AROMATIC)
               aromatic_level |= 2;
            else if (mp->bond_array[nbp->bonds[i]].bond_type == DOUBLE_AROMATIC)
               aromatic_level |= 2;
            else if (mp->bond_array[nbp->bonds[i]].bond_type == ANY_BOND)
               aromatic_level |= 2;
            if (mp->bond_array[nbp->bonds[i]].topography == CHAIN)
               aromatic_level |= 4;
         }
         if (nbp->n_ligands != 1                     ||
             !AtomSymbolMatch(ap->atom_symbol, "C")  ||
             !expect_smarts)
            exo_methyl = FALSE;
// fprintf(stderr, "arom = %d, sma = %d\n", aromatic_level, expect_smarts);
         if (aromatic_level == 1)        /* aromatic bond was found */
         {
            if (expect_smarts  || aromatic_smiles)
               result[0] += 'c'-'C';
            else
               strcpy(result, "*");     /* aromatic when not expected */
         }
         else if (aromatic_level & 1  &&
                  (expect_smarts  ||  aromatic_smiles))
            result[0] += 'c'-'C';
         else if ((aromatic_level == 5  ||  (aromatic_level & 2))  &&
                  !(expect_smarts  ||  aromatic_smiles))
            strcpy(result, "*");        /* aromatic when not expected */
         else if (aromatic_level & 2  ||  exo_methyl)
         {                              /* could be aromatic or non-aromatic */
            MyFree(result);
            result = TypeAlloc(strlen("[c,C]")+1, char);
            sprintf(result, "[%c,%c]",
                    ap->atom_symbol[0], ap->atom_symbol[0]+('c'-'C'));
         }
      }
      return (result);
   }
   else if (atom_text)
   {
      if (expect_smarts)
      {
         result = TypeAlloc(strlen(atom_text)+3, char);
         if (atom_text[0] != '[') strcpy(result, "[");
         strcat(result, atom_text);
         if (atom_text[0] != '[') strcat(result, "]");
         return (result);
      }
      else
      {
         result = TypeAlloc(strlen("*")+1, char);
         if (0 == strcmp(atom_text, "C")) /* handle a bug in structure macros */
            strcpy(result, "C");
         else
            strcpy(result, "*");
         return (result);
      }
   }

   else if (ap->atom_symbol[0] == 'R'  &&
            ap->atom_symbol[1] == '#')              /* attachment points */
   {
      result = TypeAlloc(strlen("[*+10]")+1, char);
      if (!GetNumProperty(mp->prop_lines, "M  RGP", index+1, &rgnum))
         if (ap->mass_difference > 0)
            rgnum = ap->mass_difference;
         else
            rgnum = 0;
      if (rgnum == 1)
         strcpy(result, "[*+]");
      else
      {
        sprintf(result,"[*+%d]", rgnum);
      }
      return (result);
   }
   else if (ap->atom_symbol[0] == 'R'  &&
            '0' < ap->atom_symbol[1]   &&  ap->atom_symbol[1] <= '9') /* ChemDraw attachment points */
   {
      result = TypeAlloc(strlen("[*+10]")+1, char);
      rgnum = ap->atom_symbol[1]-'0';
      if (rgnum == 1)
         strcpy(result, "[*+]");
      else
      {
        sprintf(result,"[*+%d]", rgnum);
      }
      return (result);
   }
   else if (ap->atom_symbol[0] == 'R'  &&
            ap->atom_symbol[1] == '\0') /* attachment points (value notation) */
   {
      result = TypeAlloc(strlen("[*+10]")+1, char);
      rgnum = (int)ap->value;
      if (ap->mass_difference > 0  &&  rgnum == 0)
         rgnum = ap->mass_difference;
      if (ap->charge > 0  &&  rgnum == 0)
         rgnum = ap->charge;

      if (rgnum == 1)
         strcpy(result, "[*+]");
      else
      {
        sprintf(result,"[*+%d]", rgnum);
      }
      return (result);
   }
   else if (0 == strcmp(ap->atom_symbol, "L"))        /* atom type list */
   {
      if (expect_smarts)
      {
         buffer[0] = '\0';
         strcat(buffer, "[");
         /* now do real processing */
         for (slp = mp->symbol_lists; !IsNULL(slp); slp=slp->next)
         {
            if (slp->atom == index+1) break;
         }
         if (!IsNULL(slp))
         {
            start_of_symbol = TRUE;
            for (cp=slp->string; (*cp); )
            {
               if (*cp == ' '  ||  *cp == ',')
               {
                  start_of_symbol = TRUE;
                  cp++;
                  continue;
               }

               if (start_of_symbol  &&  cp != slp->string)
               {
                  if (slp->logic)
                     strcat(buffer, ",");
                  else
                     strcat(buffer, ";");
               }

               if (!slp->logic && start_of_symbol)
                  strcat(buffer, "!");

               if (start_of_symbol)
               {
                  start_of_symbol = FALSE;
                  if (cp[0] == 'C'  &&  !isalpha(cp[1]))
                  {
                     strcat(buffer, "#6"); /* possibly aromatic carbon */
                     cp++;
                     continue;
                  }
                  if (cp[0] == 'N'  &&  !isalpha(cp[1]))
                  {
                     strcat(buffer, "#7"); /* possibly aromatic nitrogen */
                     cp++;
                     continue;
                  }
                  if (cp[0] == 'O'  &&  !isalpha(cp[1]))
                  {
                     strcat(buffer, "#8"); /* possibly aromatic oxigen */
                     cp++;
                     continue;
                  }
                  if (cp[0] == 'S'  &&  !isalpha(cp[1]))
                  {
                     strcat(buffer, "#16"); /* possibly aromatic sulfur */
                     cp++;
                     continue;
                  }
               }

               buffer[strlen(buffer)+1] = '\0';
               buffer[strlen(buffer)]   = *cp;
               cp++;
            }
         }
         else
            strcat(buffer, "*");  /* ignore unexpected atom list symbol */
                                  /* This may happen when hydrogens are */
                                  /* folded in. Needs a fix in utilities.c */

         if (ap->sub_desc != NONE)
         {
            if (ap->sub_desc > 0)                /* required substitution */
            {
               sprintf(tmp,";D%d", ap->sub_desc);
               strcat(buffer, tmp);
            }
            else if (ap->sub_desc == (-1))        /* no substitution */
            {
               strcat(buffer, ";D0");
            }
            else if (ap->sub_desc == (-2))        /* substitution as is */
            {
               sprintf(tmp,";D%d", nbp->n_ligands);
               strcat(buffer, tmp);
            }
         }

         strcat(buffer, "]");
         result = TypeAlloc(strlen(buffer)+1, char);
         strcpy(result, buffer);
         return (result);
      }
      else
      {
         result = TypeAlloc(strlen("*")+1, char);
         strcpy(result, "*");
         return (result);
      }
   }
   else if (0 == strcmp(ap->atom_symbol, "Q")  ||
            0 == strcmp(ap->atom_symbol, "A"))
   {
      if (expect_smarts)
      {
         if (0 == strcmp(ap->atom_symbol, "Q"))
         {
            result = TypeAlloc(strlen("[!#6;!#1]")+1, char);
            strcpy(result, "[!#6;!#1]");
         }
         else
         {
            result = TypeAlloc(strlen("[!#1]")+1, char);
            strcpy(result, "[!#1]");
         }
         return (result);
      }
      else
      {
         result = TypeAlloc(strlen("*")+1, char);
         strcpy(result, "*");
         return (result);
      }
   }

   buffer[0] = '\0';                                 /* not so easy cases */
   strcpy(buffer, "[");
   if (ap->mass_difference != NONE)
   {
      for (i=0; ptable[i].symbol; i++)
         if (0 == strcmp(ptable[i].symbol, ap->atom_symbol)) break;
      if (ptable[i].symbol == NULL)  /* element not found */
         sprintf(tmp,"%d",ap->mass_difference);
      else
         sprintf(tmp,"%d",(int)(ptable[i].mass + 0.5) + ap->mass_difference);
      strcat(buffer, tmp);
   }
   else
   {
      if (strcmp(ap->atom_symbol, "D") == 0) /* deuterium is special */
      {
         sprintf(tmp,"2");
         strcat(buffer, tmp);
      }
      else if (strcmp(ap->atom_symbol, "T") == 0) /* tritium is special */
      {
         sprintf(tmp,"3");
         strcat(buffer, tmp);
      }
   }

   if (strcmp(ap->atom_symbol, "D") == 0)
      strcat(buffer, "H");
   else if (strcmp(ap->atom_symbol, "T") == 0)
      strcat(buffer, "H");
   else
   {
      if (0 == strcmp(ap->atom_symbol, "A"))        /* SMARTS 'any' */
         strcat(buffer, "*");
      else                                        /* all other atom types */
         strcat(buffer, ap->atom_symbol);
   }

   /* check for aromatic atom and patch buffer[0] accordingly */
   if (AtomSymbolMatch(ap->atom_symbol, "C,N,O,P,S,B"))
   {
      aromatic_level = 0;
      for (i=0; i<nbp->n_ligands; i++)
         if (mp->bond_array[nbp->bonds[i]].bond_type == AROMATIC)
            aromatic_level |= 1;
         else if (mp->bond_array[nbp->bonds[i]].bond_type == SINGLE_AROMATIC)
            aromatic_level |= 2;
         else if (mp->bond_array[nbp->bonds[i]].bond_type == DOUBLE_AROMATIC)
            aromatic_level |= 2;
         else if (mp->bond_array[nbp->bonds[i]].bond_type == ANY_BOND)
            aromatic_level |= 2;
      if (aromatic_level == 1)        /* aromatic bond was found */
         buffer[strlen(buffer)-1] += 'c'-'C';
      else if (aromatic_level & 2)        /* could be aromatic or non-aromatic */
      {
         sprintf(buffer+(strlen(buffer)-1), "%c,%c",
                 ap->atom_symbol[0], ap->atom_symbol[0]+('c'-'C'));
      }
   }

   if (stereo == 1)         strcat(buffer, "@@");
   else if (stereo == (-1)) strcat(buffer, "@");

   /* Don't use computed H counts for query atoms */
   if (H_count[index+1] != 0                            && 
       !expect_smarts                                   &&
       (ap->sub_desc == NONE  ||  !expect_smarts) &&
       ap->query_H_count == 0)
   {
      if (H_count[index+1] == 1)
         strcat(buffer,"H");
      else
      {
         strcat(buffer,"H0");
         buffer[strlen(buffer)-1] += (char)H_count[index+1];
      }
   }
   else if (ap->query_H_count != 0)
   {
      if (ap->query_H_count > NONE+1  &&
          H_count[index+1] == 0       &&  expect_smarts)
      {
         strcat(buffer,";");
         for (i=ap->query_H_count-1;; i++)
         {
            strcat(buffer,"*H0,");
            buffer[strlen(buffer)-2] += i;
            if (i>=4) break;        /* trick to make it loop at least once */
         }
         buffer[strlen(buffer)-1] = '\0';
      }
      else if (!expect_smarts)
      {
         if (ap->query_H_count == 2)
         {
            strcat(buffer,"H");
         }
         else if (ap->query_H_count == 1)
         {
            /* NOP */
         }
         else
         {
            strcat(buffer,"H0");
            buffer[strlen(buffer)-1] += ap->query_H_count-1;
         }
      }
   }

   if (ap->charge > 0)
   {
      if (strchr(buffer,';')) strcat(buffer, ";*"); /* for complicated SMARTS */
      if (ap->charge == 1)
         strcat(buffer,"+");
      else if (ap->charge <= 9)
      {
         strcat(buffer,"+0");
         buffer[strlen(buffer)-1] += (char)ap->charge;
      }
      else
      {
         strcat(buffer,"+00");
         buffer[strlen(buffer)-2] += (char)((ap->charge/10)%10);
         buffer[strlen(buffer)-1] += (char)(ap->charge%10);
      }
   }
   else if (ap->charge < 0)
   {
      if (strchr(buffer,';')) strcat(buffer, ";*"); /* for complicated SMARTS */
      if (ap->charge == -1)
      {
         strcat(buffer,"-");
      }
      else if (ap->charge >= -9)
      {
         strcat(buffer,"-0");
         buffer[strlen(buffer)-1] -= (char)ap->charge;
      }
      else
      {
         strcat(buffer,"-00");
         buffer[strlen(buffer)-2] += (char)((-ap->charge/10)%10);
         buffer[strlen(buffer)-1] += (char)((-ap->charge)%10);
      }
   }

   if (ap->mapping > 0)         /* add reaction mapping info */
   {
      strcat(buffer,":");
      sprintf(tmp,"%d", ap->mapping-1);
      strcat(buffer, tmp);
   }

   if (ap->sub_desc != NONE  &&  expect_smarts)
   {
      if (ap->sub_desc > 0)                /* required substitution */
      {
         sprintf(tmp,";D%d", ap->sub_desc);
         strcat(buffer, tmp);
      }
      else if (ap->sub_desc == (-1))        /* no substitution */
      {
         strcat(buffer, ";D0");
      }
      else if (ap->sub_desc == (-2))        /* substitution as is */
      {
         sprintf(tmp,";D%d", nbp->n_ligands);
         strcat(buffer, tmp);
      }
   }

   strcat(buffer, "]");

   result = TypeAlloc(strlen(buffer)+1, char);
   strcpy(result, buffer);
   return (result);
}

char *BndAppend(char *result, char *cph, int btype, int from_first,
                int topography, int bcolor, int expect_smarts)
/*
 * Appends the string corresponding to the bond type btype followwed
 * by the one pointed to by cph to the string pointed to by result
 * and returns a pointer to the resulting string.
 * The string pointed to by cph is deallocated as a side-effect.
 *
 * The bond color is used to determine if stereo UP and DOWN
 * bonds are to be written instead of single bonds.
 */
{
   switch (btype)
   {
      case NONE:
            result = StrAppend(result, cph, ".");
         break;
      case SINGLE:      /* default => do nothing */
         if (bcolor == 0)        /* ordinary single bond */
         {
            if (topography == RING  &&  expect_smarts)
               result = StrAppend(result, cph, "@;-");
            else if (topography == CHAIN  &&  expect_smarts)
               result = StrAppend(result, cph, "!@;-");
            else if (!result) result = cph;
         }
         else if (from_first   &&  (bcolor & DY_ABOVE_FIRST))        /* OK */
            result = StrAppend(result, cph, "/");
         else if (from_first   &&  (bcolor & DY_BELOW_FIRST))        /* OK */
            result = StrAppend(result, cph, "\\");
         else if (from_first  &&  (bcolor & DY_ABOVE_SECOND))        /* OK */
            result = StrAppend(result, cph, "\\");
         else if (from_first  &&  (bcolor & DY_BELOW_SECOND))        /* OK */
            result = StrAppend(result, cph, "/");
         else if (!from_first   &&  (bcolor & DY_ABOVE_FIRST))        /* OK */
            result = StrAppend(result, cph, "\\");
         else if (!from_first   &&  (bcolor & DY_BELOW_FIRST))        /* OK */
            result = StrAppend(result, cph, "/");
         else if (!from_first  &&  (bcolor & DY_ABOVE_SECOND))        /* OK */
            result = StrAppend(result, cph, "/");
         else if (!from_first  &&  (bcolor & DY_BELOW_SECOND))        /* OK */
            result = StrAppend(result, cph, "\\");
         else
         {
            if (topography == RING  &&  expect_smarts)
               result = StrAppend(result, cph, "@;-");
            else if (topography == CHAIN  &&  expect_smarts)
               result = StrAppend(result, cph, "!@;-");
            else if (!result) result = cph;
         }
         break;
      case DOUBLE:
         if (topography == RING  &&  expect_smarts)
            result = StrAppend(result, cph, "@;=");
         else if (topography == CHAIN  &&  expect_smarts)
            result = StrAppend(result, cph, "!@;=");
         else
            result = StrAppend(result, cph, "=");
         break;
      case TRIPLE:
         result = StrAppend(result, cph, "#");
         break;
      case AROMATIC:
         if (!result)
         {
            result = cph;
         }
         else
         {
            if (!cph ||
                (cph[0] != 'c'  &&  cph[0] != 'n'))
               result = StrAppend(result, cph, "");        /* short-cut */
            else
               if (expect_smarts)
                  result = StrAppend(result, cph, ":");        /* explict aro */
               else
                  result = StrAppend(result, cph, ":");        /* explict aro */
         }
         break;
      case SINGLE_DOUBLE:
         if (expect_smarts)
         {
            if (topography == RING  &&  expect_smarts)
               result = StrAppend(result, cph, "@;=,@;-");
            else if (topography == CHAIN  &&  expect_smarts)
               result = StrAppend(result, cph, "!@;=,!@;-");
            else
               result = StrAppend(result, cph, "-,=");
         }
         else
            result = StrAppend(result, cph, "");        /* ignore query */
         break;
      case SINGLE_AROMATIC:
         if (expect_smarts)
         {
            if (topography == RING  &&  expect_smarts)
               result = StrAppend(result, cph, "@;-,@;:");
            else if (topography == CHAIN  &&  expect_smarts)
               result = StrAppend(result, cph, "!@;-,!@;:");
            else
               result = StrAppend(result, cph, "-,:");
         }
         else
            result = StrAppend(result, cph, "");        /* ignore query */
         break;
      case DOUBLE_AROMATIC:
         if (expect_smarts)
         {
            if (topography == RING  &&  expect_smarts)
               result = StrAppend(result, cph, "@;=,@;:");
            else if (topography == CHAIN  &&  expect_smarts)
               result = StrAppend(result, cph, "!@;=,!@;:");
            else
               result = StrAppend(result, cph, "=,:");
         }
         else
            result = StrAppend(result, cph, "");        /* ignore query */
         break;
      case ANY_BOND:
         if (expect_smarts)
         {
            if (topography == RING  &&  expect_smarts)
               result = StrAppend(result, cph, "@;~");
            else if (topography == CHAIN  &&  expect_smarts)
               result = StrAppend(result, cph, "!@;~");
            else
               result = StrAppend(result, cph, "~");
         }
         else
            result = StrAppend(result, cph, "");        /* ignore query */
         break;
      default:
         if (expect_smarts)
            result = StrAppend(result, cph, "?");       /* query bond */
         else
            result = StrAppend(result, cph, "");        /* ignore query */
         break;
   }

   return (result);
}

char *SmilesBranch(int istart, int iprev,
                   int open_branches[MAXOPEN],
                   int open_neighbours[MAXOPEN],
                   struct reaccs_molecule_t *mp,
                   neighbourhood_t *nbp,
                   int H_count[],
                   int iso,
                   int *nsmi_atoms,
                   int expect_smarts,
                   int aromatic_smiles,
                   int use_z)
/*
 * Traces the not-yet-visited branch starting at index istart of molecule
 * *mp. It returns a SMILES string describing this branch.
 * The array open_branches[0..MAXOPEN-1] lists the current set of not yet
 * finished branches, i.e. potential candidates for ring closures.
 * open_neighbours[i] stores the other atom of the bond from
 * open_branches[i].
 * iprev is the index of the atom from which the branch sprouted or (-1).
 * The function assumes that hydrogen counts have been perceived
 * and stored in H_count[1..mp->n_atoms] (note: index origin == 1!!!)
 *
 * The stereochemistry of the structure is included in the SMILES if
 * iso was set to TRUE.
 *
 * *nsmi_atoms is incremented by 1 whenever a new atom is written. This
 * should allow tracking of SMILES vs. MOL numbering.
 *
 * use_z indicates if 3D z-coordinates shall be used to infer stereochemistry.
 */
{
   char *result, *ligands, *cp, *cph, tmp[10];
   char *symbol, *closures;
   int nbranch, jbranch;
   int i, j;
   int lig_ind[MAXNEIGHBOURS], ilig;    // used to compute stereochemical parity from ligand indices

   result = NULL;

short_cut:
   ligands = NULL;
   closures = NULL;
   ilig = 0;
   if (iprev >= 0)                   /* atom from which we got here */
   {
      lig_ind[ilig] = iprev;
// fprintf(stderr,"~%d ",lig_ind[ilig]+1);
      ilig++;
   }
   if (nbp[istart].n_ligands == 3  &&      /* implicit hydrogen possibly needed for stereochemical parity */
       AtomSymbolMatch(mp->atom_array[istart].atom_symbol, "C,Si,Ge"))
   {
      lig_ind[ilig] = istart;
// fprintf(stderr,"!%d ",lig_ind[ilig]+1);
      ilig++;
   }
       
   (*nsmi_atoms)++;
   mp->atom_array[istart].color = (*nsmi_atoms);
// fprintf(stderr,"mp-atom #%d -> smiles #%d\n", istart+1, (*nsmi_atoms));

   nbranch = 0;                                      /* Count branching */
   for (j=0; j<nbp[istart].n_ligands; j++)
      if (nbp[istart].atoms[j] != iprev) /* don't trace back */
      {                                  /* where you came from */
         if (mp->atom_array[nbp[istart].atoms[j]].color != NONE) /* loop */
         {
            // first check if loop was already recorded
            for (i=0; i<MAXOPEN; i++)
               if (open_branches[i] == istart  &&  open_neighbours[i] == nbp[istart].atoms[j])    // free branching slot
                  break;
            if (i<MAXOPEN)      // found an already recorded loop
            {
               if (i >= 9)
               {
// if (i>=40) fprintf(stderr, "3: ring closure number '%d' needed\n", i+1);
                  if (i < 99)
                     sprintf(tmp, "%%%d", i+1);
                  else
                     sprintf(tmp, "%%{%d}", i+1);
               }
               else
                  sprintf(tmp, "%d", i+1);
               if (ilig < 4) /* save first kind of closure neighbours */
               {
                  lig_ind[ilig] = nbp[istart].atoms[j];
// fprintf(stderr,"1#%d ",lig_ind[ilig]+1);
                  ilig++;
               }
               closures = BndAppend(closures, NULL,
                             mp->bond_array[nbp[istart].bonds[j]].bond_type,
                             mp->bond_array[nbp[istart].bonds[j]].atoms[0] == istart+1,
                             mp->bond_array[nbp[istart].bonds[j]].topography,
                             mp->bond_array[nbp[istart].bonds[j]].color,
                             expect_smarts);
               closures = StrAppend(closures, NULL, tmp); /* output loop indicator */
//                open_branches[i] = -2;
               open_branches[i] = -2;
               open_neighbours[i] = -1;
               continue;
            }
            // clear branch recording if no further one is pending
            for (i=0; i<MAXOPEN; i++)
               if (open_branches[i] >= 0  ||  open_neighbours[i] >= 0) break;
            if (i == MAXOPEN) // no pending closure found => clear
            {
// fprintf(stderr, "1: clearing ring closure buffers\n");
               for (i=0; i<MAXOPEN; i++)
               {
                  open_branches[i] = -1;
                  open_neighbours[i] = -1;
               }
            }
            // not yet recorded => then record it
            for (i=0; i<MAXOPEN; i++)
               if (open_branches[i] == (-1))    // free branching slot
               {
                  if (ilig < 4) /* save first kind of closure neighbours */
                  {
                     lig_ind[ilig] = nbp[istart].atoms[j];
// fprintf(stderr,"2#%d ",lig_ind[ilig]+1);
                     ilig++;
                  }
                  if (i >= 9)
                  {
// if (i>=40) fprintf(stderr, "4: ring closure number '%d' needed\n", i+1);
                     if (i < 99)
                        sprintf(tmp, "%%%d", i+1);
                     else
                        sprintf(tmp, "%%{%d}", i+1);
                  }
                  else
                     sprintf(tmp, "%d", i+1);
                  closures = BndAppend(closures, NULL,
                                mp->bond_array[nbp[istart].bonds[j]].bond_type,
                                mp->bond_array[nbp[istart].bonds[j]].atoms[0] == istart+1,
                                mp->bond_array[nbp[istart].bonds[j]].topography,
                                mp->bond_array[nbp[istart].bonds[j]].color,
                                expect_smarts);
                  closures = StrAppend(closures, NULL, tmp); /* output loop indicator */
// fprintf(stderr, "2: opening branch '%s'\n", tmp);
                  open_branches[i] = nbp[istart].atoms[j];   /* save label */
                  open_neighbours[i] = istart;               /* save from atom */
                  break;
               }
         }
         else                                   /* normal branch */
         {
            jbranch = j;
            nbranch++;
         }
      }

   if (nbranch == 1)   /* simple chain => use short-cut */
   {
      lig_ind[ilig] = nbp[istart].atoms[jbranch];
// fprintf(stderr,"!%d ",lig_ind[ilig]+1);
      ilig++;
      symbol = SmilesAtom(mp, istart, lig_ind, ilig, nbp+istart,
                          H_count, expect_smarts, aromatic_smiles, use_z);
      ligands = BndAppend(ligands, NULL,
                          mp->bond_array[nbp[istart].bonds[jbranch]].bond_type,
                          mp->bond_array[nbp[istart].bonds[jbranch]].atoms[0] == istart+1,
                          mp->bond_array[nbp[istart].bonds[jbranch]].topography,
                          mp->bond_array[nbp[istart].bonds[jbranch]].color,
                          expect_smarts);
      result = StrAppend(result, symbol, "");
      result = StrAppend(result, closures, "");
      result = StrAppend(result, ligands, "");
      symbol = closures = ligands = NULL;
      iprev = istart;
      istart = nbp[istart].atoms[jbranch];
      goto short_cut;
   }
   else            /* regular processing */
   {
      cp = NULL; ligands = NULL;
      for (j=0; j<nbp[istart].n_ligands; j++)
      {
         if (nbp[istart].atoms[j] == iprev) continue;

         if (mp->atom_array[nbp[istart].atoms[j]].color == NONE)
         {                      /* loops have been dealt with before */
            cph = SmilesBranch(nbp[istart].atoms[j], istart,
                               open_branches, open_neighbours,
                               mp, nbp, H_count, iso,
                               nsmi_atoms,
                               expect_smarts, aromatic_smiles, use_z);
            cph = BndAppend(NULL, cph,
                            mp->bond_array[nbp[istart].bonds[j]].bond_type,
                            mp->bond_array[nbp[istart].bonds[j]].atoms[0] == istart+1,
                            mp->bond_array[nbp[istart].bonds[j]].topography,
                            mp->bond_array[nbp[istart].bonds[j]].color,
                            expect_smarts);
            if (cp)                             /* branch with parentesis */
            {
               ligands = StrAppend(ligands, NULL, "(");
               ligands = StrAppend(ligands, cp, "");
               ligands = StrAppend(ligands, NULL, ")");
            }
            cp = cph;
         }
      }
      if (cp) ligands = StrAppend(ligands, cp, "");    /* branch w/o parenthesis */
      for (i=0; i<MAXOPEN; i++)
         if (open_branches[i] == istart)
         {
            if (ilig < 4)
            {
               lig_ind[ilig] = open_neighbours[i];
// fprintf(stderr,"$%d ",lig_ind[ilig]+1);
               ilig++;
            }
            if (i >= 9)
            {
// FORTIFY if (i>=50) fprintf(stderr, "1: ring closure number '%d' needed\n", i+1);
               if (i < 99)
                  sprintf(tmp, "%%%d", i+1);
               else
                  sprintf(tmp, "%%{%d}", i+1);
            }
            else
               sprintf(tmp, "%d", i+1);
            closures = StrAppend(closures, NULL, tmp);
// fprintf(stderr, "closing branch '%s'\n", tmp);
//            open_branches[i] = -2;
            open_branches[i] = -2;
            open_neighbours[i] = -1;
         }
      // clear branch recording if no further one is pending
      for (i=0; i<MAXOPEN; i++)
         if (open_branches[i] >= 0  ||  open_neighbours[i] >= 0) break;
      if (i == MAXOPEN) // no pending closure found => clear
      {
fprintf(stderr, "2: clearing ring closure buffers\n");
         for (i=0; i<MAXOPEN; i++)
         {
            open_branches[i] = -1;
            open_neighbours[i] = -1;
         }
      }
      for (i=0; i<nbp[istart].n_ligands; i++) /* collect regular ligands */
      {
         for (j=0; j<ilig; j++)
            if (nbp[istart].atoms[i] == lig_ind[j]) break;
         if (j == ilig) /* new ligand atom */
            if (ilig < 4)
            {
               lig_ind[ilig] = nbp[istart].atoms[i];
// fprintf(stderr,"@%d ",lig_ind[ilig]+1);
               ilig++;
            }
      }
      symbol = SmilesAtom(mp, istart, lig_ind, ilig, nbp+istart,
                          H_count, expect_smarts, aromatic_smiles, use_z);
      result = StrAppend(result, symbol, NULL);
      result = StrAppend(result, closures, NULL);
      result = StrAppend(result, ligands, NULL);
      symbol = closures = ligands = NULL;
   }

   return (result);
}

#define DB_OPEN          16
#define DB_CLOSED        32

/**
 * Auxillary function used to generate a depth-first list of bonds
 * from *mp to be used to generate SMILES output.
 */
int DepthFirstVisit(int root_index, struct reaccs_molecule_t *mp,
                    neighbourhood_t nbp[], int *result, int ibond)
{
   int i;

   for (i=0; i<nbp[root_index].n_ligands; i++)
      if (mp->bond_array[nbp[root_index].bonds[i]].color == 0)
      {
         mp->atom_array[nbp[root_index].atoms[i]].color = VISITED;
         mp->bond_array[nbp[root_index].bonds[i]].color = VISITED;
         result[ibond++] = nbp[root_index].bonds[i];
         ibond = DepthFirstVisit(nbp[root_index].atoms[i],
                                 mp, nbp, result, ibond);
      }
   return (ibond);
}

int *DepthFirstBondList(struct reaccs_molecule_t *mp,
                        neighbourhood_t *nbp,
                        int numbering[])
/*
 * Returns a freshly allocated array of indices into mp->bond_array
 * such that the list of bonds is visited in a depth first order.
 * This is a function used to make sure that double bond stereo
 * assignments are not conflicting in the resulting SMILES.
 */
{
   int *result;
   int i, ibond, ibranch, best_number;

   result = (int *)calloc(mp->n_bonds, sizeof(int));
   for (i=0; i<mp->n_atoms; i++) mp->atom_array[i].color = 0;
   for (i=0; i<mp->n_bonds; i++) mp->bond_array[i].color = 0;

   ibond = 0;
   for (;;)
   {
      if (numbering == (int *)NULL)
      {
         for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].color == 0) break;
         ibranch = i;
      }
      else      /* Use the numbering to select next branch */
      {

         ibranch = -1; best_number = mp->n_atoms+1;
         for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].color == 0  &&
                numbering[i] < best_number)
            {
               ibranch = i; best_number = numbering[i];
// fprintf(stderr, "numbering[%d] = %d\n", ibranch, numbering[ibranch]);
            }
         if (ibranch == -1)
            ibranch = mp->n_atoms;
// else
// fprintf(stderr, "numbering[%d] = %d\n", ibranch, numbering[ibranch]);
      }
      if (ibranch >= mp->n_atoms) break;        /* no more open nodes */
      mp->atom_array[ibranch].color = VISITED;
      ibond = DepthFirstVisit(ibranch, mp, nbp, result, ibond);
   }

   for (i=0; i<mp->n_atoms; i++) mp->atom_array[i].color = 0;
   for (i=0; i<mp->n_bonds; i++) mp->bond_array[i].color = 0;

   return (result);
}

typedef unsigned pair_t[2];

static void PerceiveSmallRingBonds(struct reaccs_molecule_t *mp)
/*
 * Sets the topography of the bonds of *mp that are members of small rings.
 */
{
   unsigned int i;
   bond_set_node *ring_list, *plist;
   pair_t *graph;

   graph = TypeAlloc(mp->n_bonds+1, pair_t);
   for (i=0; i<mp->n_bonds; i++)
   {
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
   }
   ring_list = RingList(graph,mp->n_bonds);
   ring_list = CombineRings(ring_list);

   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (i=0; i<mp->n_bonds; i++)
         if (IsMember(plist->bond_set,i))
         {
            if (Cardinality(plist->bond_set) < 8)
               mp->bond_array[i].topography = RING;
         }

   DisposeBondSetList(ring_list);

   MyFree((char *) graph);
}

static void FixBisaryl(struct reaccs_molecule_t *mp)
/*
 * Change any non-ring aromatic bonds between ring atoms to s/a to catch bis-aryl cases.
 */
{
   unsigned int i;
   bond_set_node *ring_list, *plist;
   int *is_ring_bond;
   pair_t *graph;

   is_ring_bond = TypeAlloc(mp->n_bonds, int);
   for (i=0; i<mp->n_bonds; i++) is_ring_bond[i] = FALSE;
   graph = TypeAlloc(mp->n_bonds+1, pair_t);
   for (i=0; i<mp->n_bonds; i++)
   {
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
   }
   ring_list = RingList(graph,mp->n_bonds);

   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (i=0; i<mp->n_bonds; i++)
         if (IsMember(plist->bond_set,i))
         {
            is_ring_bond[i] = TRUE;
         }

   DisposeBondSetList(ring_list);
   for (i=0; i<mp->n_bonds; i++)
      if (mp->bond_array[i].bond_type == AROMATIC  &&  !is_ring_bond[i])
         mp->bond_array[i].bond_type = SINGLE_AROMATIC;

   MyFree((char *) graph);
   MyFree((char *) is_ring_bond);
}

void PerceiveDBStereo(struct reaccs_molecule_t *mp,
                      neighbourhood_t *nbp,
                      int *numbering)
/*
 * Perceives the DY_ABOVE_* and DY_BELOW_* properties of the bonds attached
 * to cis/trans defined double bonds. DY_ABOVE_* is set if the bond
 * is below the double bond directed left to right, and DY_ABOVE_*
 * is set otherwise.
 *
 * The *_FIRST flags are used when the double bond is at the first
 * atom of the bond and *_SECOND is used otherwise.
 *
 * After preliminary settings are defined, an attempt is made to work
 * around conflickting settings of the same single bond by flipping
 * the UP<=>DOWN settings for one of the offending double bonds.
 *
 * These settings can be directly translated into
 * UP ('/') and DOWN ('\') bonds while trancing the SMILES depending
 * on whether the ligand or the double bond is cited first on the path.
 *
 * The perceived bond direction is ORed into the bonds color field.
 *
 * If not NULL, numbering[] is the canonical numbering of mp. This numbering
 * is used to select between otherwise equivalent ABOVE/BELOW selections.
 */
{
   struct reaccs_atom_t *ap1, *ap2;
   struct reaccs_bond_t *bp;
   int i, j, ai1, ai2;
   int nmulti;
   int add_color;
   double dx1, dy1, dx2, dy2;
   int *topography_sav;

   int *df_bond_indices;

   int collision_found;

   if (mp->n_bonds <= 0) return;

   topography_sav = TypeAlloc(mp->n_bonds, int);
   for (i=0; i<mp->n_bonds; i++)
      topography_sav[i] = mp->bond_array[i].topography;

   df_bond_indices = DepthFirstBondList(mp, nbp, numbering);

   PerceiveSmallRingBonds(mp);

   for (i=0; i<mp->n_bonds; i++)
   {
      bp = &mp->bond_array[df_bond_indices[i]];

      if (bp->bond_type     != DOUBLE)           continue;
      if (bp->topography    == RING)             continue;
      if (bp->stereo_symbol == CIS_TRANS_EITHER) continue;
// if (0 && numbering)
// fprintf(stderr, "perceiving %d-%d\n",
// numbering[bp->atoms[0]-1],
// numbering[bp->atoms[1]-1]);

      /* check for trivial non-stereogenic bonds */
      ai1 = bp->atoms[0]-1; ai2 = bp->atoms[1]-1;
      if (nbp[ai1].n_ligands < 2  ||  nbp[ai1].n_ligands > 3) continue;
      if (nbp[ai2].n_ligands < 2  ||  nbp[ai2].n_ligands > 3) continue;
      ap1 = &mp->atom_array[ai1]; ap2 = &mp->atom_array[ai2];
      if (!AtomSymbolMatch(ap1->atom_symbol, "C,N")) continue;
      if (!AtomSymbolMatch(ap2->atom_symbol, "C,N")) continue;

      /* no allenes */
      nmulti = 0;
      for (j=0; j<nbp[ai1].n_ligands; j++)
         if (mp->bond_array[nbp[ai1].bonds[j]].bond_type != SINGLE) nmulti++;
      if (nmulti > 1) continue;
      nmulti = 0;
      for (j=0; j<nbp[ai2].n_ligands; j++)
         if (mp->bond_array[nbp[ai2].bonds[j]].bond_type != SINGLE) nmulti++;
      if (nmulti > 1) continue;

      collision_found = FALSE;
      /* set bond colors */
      dx1 = mp->atom_array[ai2].x - mp->atom_array[ai1].x;
      dy1 = mp->atom_array[ai2].y - mp->atom_array[ai1].y;
      for (j=0; j<nbp[ai1].n_ligands; j++)
      {
         if (nbp[ai1].atoms[j]+1 == bp->atoms[1]) continue;
         dx2 = mp->atom_array[nbp[ai1].atoms[j]].x - mp->atom_array[ai1].x;
         dy2 = mp->atom_array[nbp[ai1].atoms[j]].y - mp->atom_array[ai1].y;
         if (dx1*dy2-dx2*dy1 > 0)
         {
            if (ai1+1 == mp->bond_array[nbp[ai1].bonds[j]].atoms[0])
               add_color = DY_ABOVE_FIRST;
            else
               add_color = DY_ABOVE_SECOND;
            if (numbering != NULL  &&
                numbering[ai1] > numbering[ai2])
               add_color = SwapFirstSecond(add_color);
            mp->bond_array[nbp[ai1].bonds[j]].color |= add_color;
// fprintf(stderr,"bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai1].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai1].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai1].bonds[j]].color);
         }
         else
         {
            if (ai1+1 == mp->bond_array[nbp[ai1].bonds[j]].atoms[0])
               add_color = DY_BELOW_FIRST;
            else
               add_color = DY_BELOW_SECOND;
            if (numbering != NULL  &&
                numbering[ai1] > numbering[ai2])
               add_color = SwapFirstSecond(add_color);
            mp->bond_array[nbp[ai1].bonds[j]].color |= add_color;
// fprintf(stderr,"bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai1].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai1].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai1].bonds[j]].color);
         }
         if (mp->bond_array[nbp[ai1].bonds[j]].color == DY_ABOVE_BOTH)
            collision_found = TRUE;
         if (mp->bond_array[nbp[ai1].bonds[j]].color == DY_BELOW_BOTH)
            collision_found = TRUE;
      }
      for (j=0; j<nbp[ai2].n_ligands; j++)
      {
         if (nbp[ai2].atoms[j]+1 == bp->atoms[0]) continue;
         dx2 = mp->atom_array[nbp[ai2].atoms[j]].x - mp->atom_array[ai2].x;
         dy2 = mp->atom_array[nbp[ai2].atoms[j]].y - mp->atom_array[ai2].y;
         if (dx1*dy2-dx2*dy1 > 0)
         {
            if (ai2+1 == mp->bond_array[nbp[ai2].bonds[j]].atoms[0])
               add_color = DY_ABOVE_FIRST;
            else
               add_color = DY_ABOVE_SECOND;
            if (numbering != NULL  &&
                numbering[ai1] > numbering[ai2])
               add_color = SwapFirstSecond(add_color);
            mp->bond_array[nbp[ai2].bonds[j]].color |= add_color;
// fprintf(stderr,"bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai2].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai2].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai2].bonds[j]].color);
         }
         else
         {
            if (ai2+1 == mp->bond_array[nbp[ai2].bonds[j]].atoms[0])
               add_color = DY_BELOW_FIRST;
            else
               add_color = DY_BELOW_SECOND;
            if (numbering != NULL  &&
                numbering[ai1] > numbering[ai2])
               add_color = SwapFirstSecond(add_color);
            mp->bond_array[nbp[ai2].bonds[j]].color |= add_color;
// fprintf(stderr,"bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai2].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai2].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai2].bonds[j]].color);
         }
         if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_ABOVE_BOTH)
            collision_found = TRUE;
         if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_BELOW_BOTH)
            collision_found = TRUE;
      }
      if (collision_found)        /* retry with opposite db orientation */
      {
         collision_found = FALSE;
         dx1 = -dx1; dy1 = -dy1;
         for (j=0; j<nbp[ai1].n_ligands; j++)
         {
            if (nbp[ai1].atoms[j]+1 == bp->atoms[1]) continue;
            dx2 = mp->atom_array[nbp[ai1].atoms[j]].x - mp->atom_array[ai1].x;
            dy2 = mp->atom_array[nbp[ai1].atoms[j]].y - mp->atom_array[ai1].y;
            if (dx1*dy2-dx2*dy1 > 0)
            {
               if (ai1+1 == mp->bond_array[nbp[ai1].bonds[j]].atoms[0])
               {
                  mp->bond_array[nbp[ai1].bonds[j]].color &= ~DY_BELOW_FIRST;
                  mp->bond_array[nbp[ai1].bonds[j]].color |= DY_ABOVE_FIRST;
               }
               else
               {
                  mp->bond_array[nbp[ai1].bonds[j]].color &= ~DY_BELOW_SECOND;
                  mp->bond_array[nbp[ai1].bonds[j]].color |= DY_ABOVE_SECOND;
               }
// fprintf(stderr,"collision: bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai1].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai1].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai1].bonds[j]].color);
            }
            else
            {
               if (ai1+1 == mp->bond_array[nbp[ai1].bonds[j]].atoms[0])
               {
                  mp->bond_array[nbp[ai1].bonds[j]].color &= ~DY_ABOVE_FIRST;
                  mp->bond_array[nbp[ai1].bonds[j]].color |= DY_BELOW_FIRST;
               }
               else
               {
                  mp->bond_array[nbp[ai1].bonds[j]].color &= ~DY_ABOVE_SECOND;
                  mp->bond_array[nbp[ai1].bonds[j]].color |= DY_BELOW_SECOND;
               }
// fprintf(stderr,"collision: bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai1].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai1].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai1].bonds[j]].color);
            }
            if (mp->bond_array[nbp[ai1].bonds[j]].color == DY_ABOVE_BOTH)
               collision_found = TRUE;
            if (mp->bond_array[nbp[ai1].bonds[j]].color == DY_BELOW_BOTH)
               collision_found = TRUE;
         }
         for (j=0; j<nbp[ai2].n_ligands; j++)
         {
            if (nbp[ai2].atoms[j]+1 == bp->atoms[0]) continue;
            dx2 = mp->atom_array[nbp[ai2].atoms[j]].x - mp->atom_array[ai2].x;
            dy2 = mp->atom_array[nbp[ai2].atoms[j]].y - mp->atom_array[ai2].y;
            if (dx1*dy2-dx2*dy1 > 0)
            {
               if (ai2+1 == mp->bond_array[nbp[ai2].bonds[j]].atoms[0])
               {
                  mp->bond_array[nbp[ai2].bonds[j]].color &= ~DY_BELOW_FIRST;
                  mp->bond_array[nbp[ai2].bonds[j]].color |= DY_ABOVE_FIRST;
               }
               else
               {
                  mp->bond_array[nbp[ai2].bonds[j]].color &= ~DY_BELOW_SECOND;
                  mp->bond_array[nbp[ai2].bonds[j]].color |= DY_ABOVE_SECOND;
               }
// fprintf(stderr,"collision: bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai2].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai2].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai2].bonds[j]].color);
            }
            else
            {
               if (ai2+1 == mp->bond_array[nbp[ai2].bonds[j]].atoms[0])
               {
                  mp->bond_array[nbp[ai2].bonds[j]].color &= ~DY_ABOVE_FIRST;
                  mp->bond_array[nbp[ai2].bonds[j]].color |= DY_BELOW_FIRST;
               }
               else
               {
                  mp->bond_array[nbp[ai2].bonds[j]].color &= ~DY_ABOVE_SECOND;
                  mp->bond_array[nbp[ai2].bonds[j]].color |= DY_BELOW_SECOND;
               }
// fprintf(stderr,"collision: bond %d-%d colored to %d\n",
// mp->bond_array[nbp[ai2].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai2].bonds[j]].atoms[1],
// mp->bond_array[nbp[ai2].bonds[j]].color);
            }
            if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_ABOVE_BOTH)
               collision_found = TRUE;
            if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_BELOW_BOTH)
               collision_found = TRUE;
         }
      }
      if (collision_found)
         ShowMessage("Could not fix DB designation collisions", "PerceiveDBStereo");
      if (bp->stereo_symbol == CIS_TRANS_SWAPPED)
      {
         for (j=0; j<nbp[ai2].n_ligands; j++)
         {
            if (nbp[ai2].atoms[j]+1 == bp->atoms[0]) continue;
            if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_ABOVE_FIRST)
                mp->bond_array[nbp[ai2].bonds[j]].color = DY_BELOW_FIRST;
            else if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_BELOW_FIRST)
                mp->bond_array[nbp[ai2].bonds[j]].color = DY_ABOVE_FIRST;
            else if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_ABOVE_SECOND)
                mp->bond_array[nbp[ai2].bonds[j]].color = DY_BELOW_SECOND;
            else if (mp->bond_array[nbp[ai2].bonds[j]].color == DY_BELOW_SECOND)
                mp->bond_array[nbp[ai2].bonds[j]].color = DY_ABOVE_SECOND;
// fprintf(stderr, "swapped bond symbol %d of bond %d-%d\n",
// mp->bond_array[nbp[ai2].bonds[j]].color,
// mp->bond_array[nbp[ai2].bonds[j]].atoms[0],
// mp->bond_array[nbp[ai2].bonds[j]].atoms[1]);
         }
      }
   }
   /* restore topography */
   for (i=0; i<mp->n_bonds; i++)
      mp->bond_array[i].topography = topography_sav[i];
   MyFree((char*)topography_sav);
   MyFree((char *)df_bond_indices);
}

void SortNeighboursByShortcut(struct reaccs_molecule_t *mp, neighbourhood_t *nbp)
/*
 * Sort the neighbourhoods such that shortcut atoms get traversed first, which means they are last in the nb list but sorted ascending among themselves.
 */
{
    int ai, i, j, h;
    int is_shortcut1, is_shortcut2;

    for (ai=0; ai<mp->n_atoms; ai++, nbp++)
    {
        for (i=1; i<nbp->n_ligands; i++)
        {
            for (j=i-1; j>=0; j--)
            {
                is_shortcut1 = (0 == strcmp(mp->atom_array[nbp->atoms[j]].atom_symbol, "R"))    && (mp->atom_array[nbp->atoms[j]].atext[0] != '\0');
                is_shortcut2 = (0 == strcmp(mp->atom_array[nbp->atoms[j+1]].atom_symbol, "R"))  && (mp->atom_array[nbp->atoms[j+1]].atext[0] != '\0');
                if ((is_shortcut2 && !is_shortcut1)  ||  // shortcuts come last
                    (is_shortcut1 && is_shortcut2 &&     // NO_BONDs sort first
                     mp->bond_array[nbp->bonds[j]].bond_type < mp->bond_array[nbp->bonds[j+1]].bond_type)  ||
                    (is_shortcut1 && is_shortcut2 &&
                     mp->bond_array[nbp->bonds[j]].bond_type == mp->bond_array[nbp->bonds[j+1]].bond_type  &&
                     nbp->atoms[j] > nbp->atoms[j+1]))      // shortcuts with lower index go first
                {
                    h=nbp->bonds[j]; nbp->bonds[j]=nbp->bonds[j+1]; nbp->bonds[j+1]=h;
                    h=nbp->atoms[j]; nbp->atoms[j]=nbp->atoms[j+1]; nbp->atoms[j+1]=h;
                }
                else
                    break;
            }
        }
    }
}

char *MOLToSMIExt(struct reaccs_molecule_t *mp, int flags,
                  int numbering[],
                  char **coordpp)
/*
 * Returns a pointer to the SMILES string representing the connectivity
 * of the molecule *mp. The memory must be deallocated by the caller.
 *
 * flags can be ISOMERIC_SMILES which retains isomeric SMILES information
 * and/or SMARTS_PERCEPTION which performs a SMARTS_PERCEPTION before
 * *mp is translated into SMILES.
 *
 * The numbering[] array is used to break ties when generating the out SMILES.
 * It may be NULL which means an arbitrary order. Index origin of numbering[]
 * is 0.
 *
 * The coordinates of the atoms in the molecule are written as comma
 * separated values to a string and a pointer to this string is assigned
 * to (*coordpp).
 */
{
   char coordbuffer[30];
   char *result, *cp, *closures, *coordp;
   char num_buf[10];
   int i;
   int j;
   neighbourhood_t *nbp;
   int open_branches[MAXOPEN];
   int open_neighbours[MAXOPEN];
   int *H_count;
   int nsmi_atoms;  /* used to track number of SMILES atom already written */
   int iso = FALSE;
   int ibranch, first_shortcut, shortcut_bond_type, best_number;
   int best_degree;

   if (flags & ISOMERIC_SMILES) iso = TRUE;

   if (mp->n_atoms <= 0)
   {
      if (!IsNULL(coordpp)) (*coordpp) = (char *)NULL;
      return ((char *)NULL);
   }

   if (flags & SMARTS_PERCEPTION)
   {
      MakeHydrogensImplicit(mp);
   }
   else if (flags & DY_AROMATICITY)
   {
      MakeHydrogensImplicit(mp);
      for (i=0; i<mp->n_atoms; i++)
         mp->atom_array[i].query_H_count = NONE;
   }

   H_count = TypeAlloc(mp->n_atoms+1, int);
   if (flags & SMARTS_PERCEPTION)
   {
      for (i=0; i<mp->n_atoms; i++)
         if (mp->atom_array[i].query_H_count <= NONE)
            H_count[i+1] = 0;
         else
            H_count[i+1] = -100;
   }
   ComputeImplicitH(mp, H_count);
   if (flags & SMARTS_PERCEPTION)
   {
      for (i=0; i<mp->n_atoms; i++)
         if (H_count[i+1] < 0) H_count[i+1] = 0;
   }

   nsmi_atoms = 0;
   ResetColors(mp);
   nbp   = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   if (flags & SMARTS_PERCEPTION)
   {
      if (flags & CLASSIC_AROMATICITY)
         PerceiveAromaticBonds(mp);
      else
         PerceiveDYAromaticity(mp, nbp);
      PerceiveMarkush(mp, nbp);
   }
   else
   {
      if (flags & DY_AROMATICITY)
         PerceiveDYAromaticity(mp, nbp);
   }

   result = NULL;

   /* Use the numbering to sort neighbourhood */
   if (numbering != (int *)NULL)
   {
      SortBondsWithNumbering(mp, numbering);
      SetupNeighbourhood(mp,nbp,mp->n_atoms);
      SortNeighboursByNumbering(mp, nbp, numbering);
   }
   else
       SortNeighboursByShortcut(mp, nbp);

   PerceiveDBStereo(mp, nbp, numbering);

   for (i=0; i<MAXOPEN; i++) open_branches[i] = (-1);
   // first spit out strings of shortcuts while collecting open_branches[] and open_neighbours[]
   // for any branch points
   if (numbering == (int *)NULL)        // only do this if not canonicalized
   {
      for (;;)
      {
         // search for first shortcut 'atom' if any
         for (i=0; i<mp->n_atoms; i++)
            if (0 == strcmp(mp->atom_array[i].atom_symbol, "R")  &&
                mp->atom_array[i].atext[0] != '\0'  &&
                mp->atom_array[i].color == 0) break;
         if (i == mp->n_atoms)  // no further shortcut found => done
             break;
         ibranch = i;
// fprintf(stderr,"branch %d:\t", ibranch);
         // trace current sequence of shortcuts
         if (result) result = StrAppend(result, NULL, ".");
         result = StrAppend(result, NULL, "{");
         result = StrAppend(result, NULL, mp->atom_array[ibranch].atext);
         result = StrAppend(result, NULL, "}");
// fprintf(stderr, "{{%s:%d}}", mp->atom_array[ibranch].atext, ibranch+1);
         mp->atom_array[ibranch].color = ++nsmi_atoms;
         // now we trace the shortcut string
         for (;;)
         {
             first_shortcut = -1; shortcut_bond_type = NONE;
             closures = (char *)NULL;  // collect closures here
             for (j=0; j<nbp[ibranch].n_ligands; j++)
             {
                 if (mp->atom_array[nbp[ibranch].atoms[j]].color != NONE) continue; // already visited
                 // if (mp->bond_array[nbp[ibranch].bonds[j]].bond_type == NONE) continue;
                 if (first_shortcut < 0  &&
                     0 == strcmp(mp->atom_array[nbp[ibranch].atoms[j]].atom_symbol, "R") &&     // shortcut
                     mp->atom_array[nbp[ibranch].atoms[j]].atext[0] != '\0')
                 {
                     if (first_shortcut < 0)    // redundant
                     {
                        first_shortcut = nbp[ibranch].atoms[j];
                        shortcut_bond_type = mp->bond_array[nbp[ibranch].bonds[j]].bond_type;
// fprintf(stderr, "recording first shortcut (%d->%d) '%d'\n", ibranch, j, first_shortcut);
                     }
                 }
                 else // deal with SMILES section closures
                 {
                     // find a slot in open_branches
                     for (i=0; i<MAXOPEN; i++)
                         if (open_branches[i] == -1) break;
                     if (i >= 9)
                     {
// if (i>=40) fprintf(stderr, "2: ring closure number '%d' needed\n", i+1);
                        if (i < 99)
                           sprintf(num_buf, "%%%d", i+1);
                        else
                           sprintf(num_buf, "%%{%d}", i+1);
                     }
                     else
                        sprintf(num_buf, "%d", i+1);
                     /*
                     closures = BndAppend(closures, NULL,
                                mp->bond_array[nbp[ibranch].bonds[j]].bond_type,
                                mp->bond_array[nbp[ibranch].bonds[j]].atoms[0] == ibranch+1,
                                mp->bond_array[nbp[ibranch].bonds[j]].topography,
                                mp->bond_array[nbp[ibranch].bonds[j]].color,
                                FALSE);
                     */
                     closures = StrAppend(closures, NULL, num_buf); /* output loop indicator */
                     open_branches[i] = nbp[ibranch].atoms[j];   /* save label */
// fprintf(stderr, "1: opening branch '%s/%d' to '%s/%s' atom %d\n",
// num_buf, ibranch,
// mp->atom_array[open_branches[i]].atom_symbol,
// mp->atom_array[open_branches[i]].atext,
// 1+open_branches[i]);
                     open_neighbours[i] = ibranch;               /* save from atom */
                 }
             }
             result = StrAppend(result, closures, "");
             if (first_shortcut < 0) break;
             if (shortcut_bond_type == NONE) result = StrAppend(result, NULL, ".");
             ibranch = first_shortcut;
             result = StrAppend(result, NULL, "{");
             result = StrAppend(result, NULL, mp->atom_array[ibranch].atext);
             result = StrAppend(result, NULL, "}");
// fprintf(stderr, "{{%s:%d}}", mp->atom_array[ibranch].atext, ibranch+1);
            for (j=0; j<MAXOPEN; j++)
            {
                if (open_branches[j] == -1) continue;
                if (ibranch == open_branches[j])    // add ring digit
                {
                   sprintf(num_buf, "%d", j+1);
                   result = StrAppend(result, NULL, num_buf);
                   open_branches[j] = -1;
                   open_neighbours[j] = -1;
                }
             }
             mp->atom_array[ibranch].color = ++nsmi_atoms;
         }
// fprintf(stderr, "\n");
      }
// fprintf(stderr, "\n");
// fprintf(stderr, "Pending Closures\n");
// for (i=0; i<MAXOPEN; i++)
// {
//     if (open_branches[i] == -1) continue;
//     fprintf(stderr, "digit=%d, atom=%d, back-atom=%d\n", i+1, open_branches[i]+1, open_neighbours[i]+1);
//  }
//  fprintf(stderr, "\n");
   }
   // now, do regular atom processing
   for (;;)
   {
      if (numbering == (int *)NULL)
      {
         // search for pseudo-peripheral vertex
         ibranch = mp->n_atoms;
         // prefer terminal atoms
         best_degree = 100;
         for (i=0; i<mp->n_atoms; i++)
         {
            if (mp->atom_array[i].color == 0)
            {
               if (nbp[i].n_ligands < best_degree)
               {
                  ibranch = i;
                  best_degree = nbp[i].n_ligands;
               }
            }
         }
      }
      else      /* Use the numbering to select next branch */
      {

         ibranch = -1;
         best_number = mp->n_atoms+1;
         /* Scan R-atoms first */
         for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].color == 0  &&
                numbering[i] < best_number  &&
                mp->atom_array[i].atom_symbol[0] == 'R')
            {
               ibranch = i;
               best_number = numbering[i];
            }
         /* Don't prefer R-atoms if monomer */
         if (flags & TO_MONOMER) best_number = mp->n_atoms+1;
         /* Now, consider non R-atoms */
         for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].color == 0  &&
                numbering[i] < best_number  &&
                mp->atom_array[i].atom_symbol[0]!='R')
            {
               ibranch = i;
               best_number = numbering[i];
// fprintf(stderr, "numbering[%d] = %d\n", ibranch, numbering[ibranch]);
// fprintf(stderr, "atom_symbol = '%s'\n", mp->atom_array[i].atom_symbol);
            }
         if (ibranch == -1)
            ibranch = mp->n_atoms;
// else
// fprintf(stderr, "numbering[%d] = %d\n", ibranch, numbering[ibranch]);
      }

      if (ibranch == mp->n_atoms) break;     /* all atoms visited */
      /* make sure that all fragments start with fresh branching */
      /* This is needed to have disconnected canonicalized fragment */
      /* yield identical substrings, while making sure that shortcut sequences */
      /* are restarted and connected correctly. */
      for (i=0; i<MAXOPEN; i++)
          if (open_branches[i] < 0) open_branches[i] = (-1);
// fprintf(stderr, "Pending Closures at atom %d\n", ibranch+1);
// for (i=0; i<MAXOPEN; i++)
// {
//      if (open_branches[i] == -1) continue;
//      fprintf(stderr, "digit=%d, atom=%d, back-atom=%d\n", i+1, open_branches[i]+1, open_neighbours[i]+1);
// }
// fprintf(stderr, "\n");
      cp = SmilesBranch(ibranch, (-1),
                        open_branches, open_neighbours,
                        mp, nbp, H_count, iso, &nsmi_atoms,
                        flags & SMARTS_PERCEPTION,
                        flags & (AROMATIC_SMILES|DY_AROMATICITY),
                        flags & USE_Z);
      if (result)
         result = StrAppend(result, cp, ".");
      else
         result = StrAppend(result, cp, "");
   }

   if (!IsNULL(coordpp))
   {
      coordp = NULL;

      for (j=1; j<=mp->n_atoms; j++)
         for (i=0; i<mp->n_atoms; i++)
            if (j==mp->atom_array[i].color)
            {
               sprintf(coordbuffer,"%.2f,%.2f,",
                       mp->atom_array[i].x,mp->atom_array[i].y);
               coordp = StrAppend(coordp, (char *)NULL, coordbuffer);
            }
      /* Get rid of trailing ',' */
      if (coordp != NULL  &&
          strlen(coordp) > 0)
         coordp[strlen(coordp)-1] = '\0';
      (*coordpp) = coordp;
   }

   MyFree((char *)nbp);
   MyFree((char *)H_count);

   return (result);
}

char *MOLToSMI(struct reaccs_molecule_t *mp, int flags)
/*
 * Returns a pointer to the SMILES string representing the connectivity
 * of the molecule *mp. The memory must be deallocated by the caller.
 *
 * flags can be ISOMERIC_SMILES which retains isomeric SMILES information
 * and/or SMARTS_PERCEPTION which performs a SMARTS perception before
 * *mp is translated into SMILES.
 */
{
   char *coordp;
   char *result;

   result = MOLToSMIExt(mp, flags, (int *)NULL, &coordp);
   if (coordp) MyFree((char *)coordp);

   return (result);
}

#define MAXRINGS    70
#define MAXBRANCH   MAXATOMS
const char * SmilesToMDLAtom(const char *smiles,
                             char *symbol,
                             int *charge,
                             int *radical,
                             int *isotope,
                             int *hydrogen,
                             int *stereo,
                             int *mapping,
                             int *aromatic)
/*
 * Reads the SMILES atom symbol codes enclosed in square brackets.
 * The function assumes to be positioned at the '[' character.
 * It returns a pointer to the terminal square bracket and fills
 * *symbol, *charge, *radical, *hydrogen, and *stereo with the respective
 * values in MDL conventions. *aromatic is set if a lowercase atom symbol
 * is detected indicating an aromatic atom.
 *
 * If there is an error, symbol[0] will be set to '\0'.
 */
{
   int i, len;

   symbol[0] = '\0';
   if (smiles[0] == '[') smiles++;
   else                  return (smiles);
// fprintf(stderr, "SmilesToMDLAtom(0), smiles = '%s'\n", smiles);

   (*charge) = 0; (*radical) = 0;
   (*isotope) = 0; (*hydrogen) = 0;
   (*stereo) = 0; (*mapping) = 0;
   (*aromatic) = FALSE;

   while (isdigit(*smiles))  /* isotope label */
   {
      (*isotope) *= 10;
      (*isotope) += (*smiles)-'0';
      smiles++;
   }
// fprintf(stderr, "SmilesToMDLAtom(1)\n");

   /* atom symbol */
   if (('A' <= smiles[0]  &&  smiles[0] <= 'Z')  ||  '*' == smiles[0])
      symbol[0] = smiles[0];
   else if ('a' <= smiles[0]  &&  smiles[0] <= 'z')
   {
      (*aromatic) = TRUE;
      symbol[0] = (char)(smiles[0]+'A'-'a');
   }
   else
   {
// fprintf(stderr, "SmilesToMDLAtom(2), smiles = '%s'\n", smiles);
      return (smiles);
   }
   smiles++;
// fprintf(stderr, "SmilesToMDLAtom(2), symbol = '%s'\n", symbol);

   len = 1; symbol[len] = '\0';
   while ('a' <= smiles[0]  &&  smiles[0] <= 'z')
   {
      symbol[len] = smiles[0];
      len++; smiles++;
      symbol[len] = '\0';
   }

   /* stereo class */
   while (smiles[0] == '@')
   {
      (*stereo)++; smiles++;
   }
   if ((*stereo) == 1  && '0' <= smiles[0]  &&  smiles[0] <= '9')
   {
      (*stereo) = 0;
      while ('0' <= smiles[0]  &&  smiles[0] <= '9')
      {
         (*stereo) *= 10;
         (*stereo) += (smiles[0] - '0');
         smiles++;
      }
   }
// fprintf(stderr, "SmilesToMDLAtom(3), symbol = '%s'\n", symbol);

   if ((*isotope) > 0  &&  'A' <= symbol[0]  &&  symbol[0] <= 'Z')  /* convert isotope to MDL form */
   {
      for (i=0; ptable[i].symbol; i++)
      if (0 == strcmp(ptable[i].symbol, symbol)) break;
      if (ptable[i].symbol == NULL)   /* element not found */
      {
         symbol[0] = '\0';
         return (smiles);
      }
      if (0 == strcmp(symbol,"H"))
      {
         if ((*isotope) == 2)
            symbol[0] = 'D';
         else if ((*isotope) == 3)
            symbol[0] = 'T';
         (*isotope) = 0;
      }
      else
         (*isotope) = (*isotope) - (int)(ptable[i].mass + 0.5);
   }
// fprintf(stderr, "SmilesToMDLAtom(4)\n");

   if (smiles[0] == 'H')      /* hydrogen count specification */
   {
      smiles++;
      if ('0' <= smiles[0]  &&  smiles[0] <= '9')
      {
         (*hydrogen) = (smiles[0] - '0');
         smiles++;
      }
      else
         (*hydrogen) = 1;
   }
   (*hydrogen)++;          /* hydrogen is always explicit!! */

   while (smiles[0] == '+'  ||  smiles[0] == '-')
   {
      if      ((*charge) >= 0  &&  smiles[0] == '+') (*charge)++;
      else if ((*charge) <= 0  &&  smiles[0] == '-') (*charge)--;
      else                                         break;
      smiles++;
   }
   if ('1' <= smiles[0]  &&  smiles[0] <= '9'  &&
       ((*charge) == (-1)  ||  (*charge) == 1))
   {
      (*charge) *= smiles[0]-'0';
      smiles++;
   }
   if ('0' <= smiles[0]  &&  smiles[0] <= '9'  && (*charge != 0))
   {
      if (*charge < 0)
        (*charge) = (*charge)*10 - (smiles[0]-'0');
      else
         (*charge) = (*charge)*10 + (smiles[0]-'0');
      smiles++;
   }
// fprintf(stderr, "SmilesToMDLAtom(5)\n");

   /* reaction mapping */
   if (smiles[0] == ':')
   {
      smiles++;
      while ('0' <= smiles[0]  &&  smiles[0] <= '9')
      {
         (*mapping) *= 10;
         (*mapping) += smiles[0]-'0';
         smiles++;
      }
      (*mapping)++;     /* MDL mappings start with 1 */
   }

   if (smiles[0] != ']') symbol[0] = '\0';       /* parsing error */

   if (0 == strcmp(symbol,"*")  ||  0 == strcmp(symbol,"R"))
   {
      if ((*isotope) > 0  &&  (*isotope) <= 99)     // treat the same as charged '*' atoms
      {
         symbol[0] = 'R';
         symbol[1] = '#';
         symbol[2] = '\0';
         (*charge) = (*isotope); (*isotope) = 0;
      }
      else if ((*charge) > 0  &&  (*charge) <= 99)
      {
         symbol[0] = 'R';
         symbol[1] = '#';
         symbol[2] = '\0';
      }
      else if ((*charge) == (-1))
      {
         strcpy(symbol, "Frm");
         (*charge) = 0;
      }
      else if ((*charge) == (-2))
      {
         strcpy(symbol, "To");
         (*charge) = 0;
      }
      else
      {
         symbol[0] = 'R';
         symbol[1] = '\0';
      }
   }
// fprintf(stderr, "SmilesToMDLAtom(6)\n");

   return (smiles);
}

void ComputeSmilesValenceAndHCount(struct reaccs_molecule_t *mp)
/*
 * Puts the valence and hydrogen count information consistent
 * with SMILES conventions into the dummy1 and dummy2 fields
 * of the atoms of *mp.
 * This is a preprocessing step for denormalization.
 */
{
   struct reaccs_atom_t *ap;
   int i, j;
   neighbourhood_t *nbp;
   int valence, charge, mdl_valence;
   int nsingle, ndouble, ntriple, naromatic;

   nbp   = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      for (j=0; ptable[j].symbol; j++)        /* what is MDL's valence */
      if (0 == strcmp(ptable[j].symbol, ap->atom_symbol)) break;
      if (ptable[j].symbol == NULL)  /* element not found */
         mdl_valence = 0;
      else
         mdl_valence = ptable[j].valence;

      nsingle = ndouble = ntriple = naromatic = 0;
      for (j=0; j<nbp[i].n_ligands; j++)
      switch (mp->bond_array[nbp[i].bonds[j]].bond_type)
      {
          case SINGLE: nsingle++; break;
          case DOUBLE: ndouble++; break;
          case TRIPLE: ntriple++; break;
          default:     naromatic++; break;
      }
      charge = ap->charge;
      valence = (nsingle*2 + ndouble*4 + ntriple*6 + naromatic*3)/2;
      valence += ABS(charge);
      if (ap->query_H_count > 0)
         valence += ap->query_H_count-1;
      if (valence == mdl_valence)
      {
         if (ap->query_H_count > 0)             /* [OH]CC */
            ap->dummy2 = ap->query_H_count-1;
         else
            ap->dummy2 = 0;                     /* O(CC)CC */
         ap->dummy1 = mdl_valence;
      }
      else if (valence < mdl_valence)
      {
         if (ap->query_H_count > 0)             /* [OH]CC */
         {
            ap->dummy2 = ap->query_H_count-1;
            if (naromatic == 2  &&  nsingle == 0  &&  ndouble == 0  &&
                mdl_valence == 4)
               ap->dummy1 = 2;
            else
               ap->dummy1 = valence;
         }
         else
         {
            ap->dummy2 = mdl_valence-valence;      /* OCC */
            ap->dummy1 = mdl_valence;
         }
      }
      else if (naromatic > 0)
      {
         if (ap->query_H_count > 0)             /* [O-][nH0+]1ccccc1 */
         {
            ap->dummy2 = ap->query_H_count-1;
            if (valence > mdl_valence)
               ap->dummy1 = mdl_valence + 2*((valence-mdl_valence)/2);
            else
               ap->dummy1 = valence;
         }
         else if (valence-mdl_valence  ==  3)    /* O=s1cccc1 */
         {
            ap->dummy2 = 0;
            ap->dummy1 = mdl_valence+2;
         }
         else if (valence-mdl_valence  >  1)    /* [O-][n+]1ccccc1 */
         {
            ap->dummy2 = 0;
            ap->dummy1 = valence;
         }
         else                                   /* s1cccc1 */
         {
            ap->dummy2 = 0;
            ap->dummy1 = mdl_valence;
         }
      }
      else /* valence > mdl_valence */
      {
         if (ap->query_H_count > 0)     /* probably an error */
         {
            ap->dummy2 = ap->query_H_count-1;
            ap->dummy1 = valence;
         }
         else                           /* P(=O)(-O)(-O) */
         {
            ap->dummy2 = (valence - mdl_valence)%2;
            ap->dummy1 = valence + ap->dummy2;
         }
      }
   }

   if (nbp) MyFree((char *)nbp);
}

#define SMILES_UC          1
#define SMILES_LC          2
#define SMILES_CLOSE       4
#define SMILES_OPEN        8
#define SMILES_DIGIT      16
#define SMILES_BOND       32
#define SMILES_AROMATIC   64
#define SMILES_ARO_BOND  128
#define SMILES_OTHER    1024
#define SMILES_BEFORE_RING (SMILES_UC | SMILES_AROMATIC | SMILES_ARO_BOND | SMILES_BOND | SMILES_CLOSE | SMILES_DIGIT)

static int smiles_class[256];
static int smiles_class_init = FALSE;
void InitSMILESClasses()
{
   int i;
   for (i=0; i<256; i++) smiles_class[i] = 0;
   for (i='A'; i<='Z'; i++) smiles_class[i] |= SMILES_UC;
   for (i='a'; i<='z'; i++) smiles_class[i] |= SMILES_LC;
   smiles_class[')'] |= SMILES_CLOSE; smiles_class[']'] |= SMILES_CLOSE; smiles_class['}'] |= SMILES_CLOSE;
   smiles_class['('] |= SMILES_OPEN; smiles_class['['] |= SMILES_OPEN; smiles_class['{'] |= SMILES_OPEN;
   for (i='0'; i<='9'; i++) smiles_class[i] |= SMILES_DIGIT;
   smiles_class['-']  |= SMILES_BOND; smiles_class['=']  |= SMILES_BOND; smiles_class['#']  |= SMILES_BOND;
   smiles_class['\\'] |= SMILES_BOND; smiles_class['/']  |= SMILES_BOND;
   smiles_class['.']  |= SMILES_BOND;
   smiles_class[':']  |= SMILES_ARO_BOND;
   smiles_class['c']  |= SMILES_AROMATIC; smiles_class['n']  |= SMILES_AROMATIC;
   smiles_class['o']  |= SMILES_AROMATIC; smiles_class['s']  |= SMILES_AROMATIC; smiles_class['p']  |= SMILES_AROMATIC;
   smiles_class['*']  |= SMILES_AROMATIC;
   smiles_class['.']  |= SMILES_OTHER;
   smiles_class_init = TRUE;
}

static int debug_smiles = TRUE;

/**
 * Perform a quick check if this string can be a legal SMILES.
 * It basically amounts to checking matching '[]' and '()' pairs, and testing if ring digits match up
 * and follow a legal character.
 */
int QuickCheck(const char *smiles)
{
    int i, n;
    const char *cp;
    int bra, ket;
    int closure_count[100];
    int in_atom, in_shortcut;

    if (!smiles_class_init) InitSMILESClasses();
    /* empty SMILES */
    if (!smiles  ||  !smiles[0]) return FALSE;
    bra = ket = 0;
    // check for matching '[]' pairs
    for (cp=smiles; (*cp) && (*cp) != '|'; cp++)
    {
        if ((*cp) == '[') bra++;
        else if ((*cp) == ']') ket++;
        if (bra < ket)
        {
           if (debug_smiles) fprintf(stderr, "'[]' mismatch in SMILES '%s'\n", smiles);
           return FALSE;
        }
    }
    if (bra != ket)
    {
       if (debug_smiles) fprintf(stderr, "'[]' mismatch in SMILES '%s'\n", smiles);
       return FALSE;
    }
    bra = ket = 0;
    // check for matching '()' pairs
    for (cp=smiles; (*cp) && (*cp) != '|'; cp++)
    {
        if ((*cp) == '(') bra++;
        else if ((*cp) == ')') ket++;
        if (bra < ket)
        {
           if (debug_smiles) fprintf(stderr, "'()' mismatch in SMILES '%s'\n", smiles);
           return FALSE;
        }
    }
    if (bra != ket)
    {
       if (debug_smiles) fprintf(stderr, "'()' mismatch in SMILES '%s'\n", smiles);
       return FALSE;
    }
    bra = ket = 0;
    // check for matching '{}' pairs
    for (cp=smiles; (*cp) && (*cp) != '|'; cp++)
    {
        if ((*cp) == '{') bra++;
        else if ((*cp) == '}') ket++;
        if (bra < ket)
        {
           if (debug_smiles) fprintf(stderr, "'{}' mismatch in SMILES '%s'\n", smiles);
           return FALSE;
        }
    }
    if (bra != ket)
    {
       if (debug_smiles) fprintf(stderr, "'{}' mismatch in SMILES '%s'\n", smiles);
       return FALSE;
    }
    // SMILES must not start with a digit
    if (isdigit(*smiles)) return FALSE;
    if ((*smiles) == '%') return FALSE;
    // check for ring digit (or %99 patterns)
    for (i=0; i<100; i++)
       closure_count[i] = 0;
    // return FALSE if ring closures don't match up
    in_atom = 0;
    in_shortcut = 0;
    for (cp=smiles; (*cp) && (*cp) != '|'; cp++)
    {
       if (*cp == '{') in_shortcut++;
       else if (*cp == '}') in_shortcut--;
       if (*cp == '[') in_atom++;
       else if (*cp == ']') in_atom--;
       if (in_atom > 0  ||  in_shortcut > 0) continue;       // ring digits can only be outside of atoms
       if (smiles_class[*cp] & SMILES_DIGIT)
       {
          if ((smiles_class[*(cp-1)] & SMILES_BEFORE_RING) == 0)
          {
             if (debug_smiles) fprintf(stderr, "illegal character '%c' before ring in SMILES '%s'\n", *(cp-1), smiles);
             return FALSE;
          }
          closure_count[(*cp)-'0']++;
       }
       else if (cp[0] == '%')
       {
          if ((smiles_class[*(cp-1)] & SMILES_BEFORE_RING) == 0)
          {
             if (debug_smiles) fprintf(stderr, "illegal character '%c' before two-digit ring closure in SMILES '%s'\n", *(cp-1), smiles);
             return FALSE;
          }
          if (!(smiles_class[cp[1]]&SMILES_DIGIT)  ||  !(smiles_class[cp[2]]&SMILES_DIGIT))
          {
             if (debug_smiles) fprintf(stderr, "illegal two-digit ring closure in SMILES '%s'\n", smiles);
             return FALSE;
          }
          closure_count[10*(cp[1]-'0')+(cp[2]-'0')]++;
          cp+=2;
       }
    }
    for (i=0; i<100; i++)
       if (closure_count[i] % 2 != 0)
       {
          if (debug_smiles) fprintf(stderr, "unmatched ring closure %d in SMILES '%s'\n", i, smiles);
          return FALSE;
       }
    return TRUE;
}

struct reaccs_molecule_t *SMIToMOL(const char *smiles, int flags)
/*
 * Converts the molecule described by smiles[] to the corresponding
 * MDL data structure and returns a pointer to it.
 *
 * The color field of the atom array contains the ordinal number of the
 * explicit atoms in SMILES order.
 *
 * flags is used to modify the behaviour of the conversion. Currently,
 * it is used to request a layout calculation and to remove non-essential
 * hydrogens.
 *
 * In addition to the regular SMILES syntax, smiles[] may contain a
 * coordinate section that follows the regular smiles after a '|' character.
 * The coordinate string contains a comma separated list of
 * 2D coordinates, i.e. 2 values per atom. They are used to place the
 * the atoms of the SMILES that are not implicit hydrogen atoms.
 * The implicit ones are then positioned with the layout algorithm.
 * If coordinatestring == NULL, all atoms will be layed out.
 * Atoms with missing coordinates are also layed out. Only relative
 * positioning is considered and only connected atoms are layed out
 * as one unit.
 * E.g.: "1.54,1.54,,,3.08,3.08"
 * Note: The smiles[] will be modified if there is a coordinate section.
 * I can, therefore, not be a read-only buffer like argv[*].
 */
{
   struct reaccs_molecule_t *result, *mp;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp, *bph;
   char symbol[4];
   char atext[MDL_MAXLINE+1];
   int open_rings[MAXRINGS], open_aromatic[MAXRINGS], open_bonds[MAXRINGS];
   int iring, idigit;
   int branch_stack[MAXBRANCH], aromatic_stack[MAXBRANCH], top_of_stack;
   int i, j, k, h;
   int prev_atom, bond_type, prev_aromatic;
   int charge, radical, isotope, hydrogen, stereo, mapping, aromatic;
   int *old_H_count, *new_H_count;
   int atno, *atnos;

   int bond_dir;
   int isomeric_smiles;
   int chiral;
   int ligands[MAXBRANCH], n_ligands, parity;
   int at1, at2, lig1[2], nlig1, lig2[2], nlig2;
   char *coordp, *cp;
   int novalue;
   double x, y;
   double value1, value2;
   
   if (!QuickCheck(smiles)) return ((struct reaccs_molecule_t *)NULL);
// fprintf(stderr, "QuickCheck of '%s' succeeded\n", smiles);

   coordp = strstr(smiles, "|");
   if (coordp)
   {
      coordp[0] = '\0'; coordp++;
   }
// fprintf(stderr,"1\n"); // Fortify_ListAllMemory();
   result = NewMolecule(strlen(smiles), strlen(smiles)*2); /* bonds are temp. stored twice */
// fprintf(stderr,"1a\n"); CheckFortifyMemory();
   ap = result->atom_array;
   bp = result->bond_array;
   strcpy(result->version, "V2000");
   atno = 1;
   for (i=0; i<MAXRINGS; i++) open_rings[i] = 0;
   for (i=0; i<MAXRINGS; i++) open_aromatic[i] = 0;
   for (i=0; i<MAXRINGS; i++) open_bonds[i] = 0;
   top_of_stack = 0;
   prev_atom = 0; prev_aromatic = FALSE;
   bond_type = NONE; bond_dir = NONE;
   if (flags & TRUE_DB_STEREO)
      isomeric_smiles = TRUE;
   else
      isomeric_smiles = FALSE;
   chiral = FALSE;
   while (*smiles)
   {
      atext[0] = '\0';
      switch (*smiles)
      {
         case 'C': case 'B': case 'N': case 'O':
         case 'P': case 'S': case 'F': case 'I':
         case '*': case '{':
            if (smiles[0] == '{')       // atom text
            {
                smiles++;
                for (i=0; i<MDL_MAXLINE; i++)
                    if (smiles[0] == '}'  ||  smiles[0] == '\0') break;
                    else
                    {
                        atext[i] = smiles[0]; smiles++;
                    }
                if (smiles[0] != '}' ||  i == 0) goto error_cleanup;
                atext[i] = '\0';
                strcpy(symbol, "R");
            }
            else if (smiles[0] == 'C'  &&  smiles[1] == 'l')
            {
               strcpy(symbol,"Cl");
               smiles++;
            }
            else if (smiles[0] == 'B' && smiles[1] == 'r')
            {
               strcpy(symbol,"Br");
               smiles++;
            }
            else
            {
               symbol[0] = *smiles;
               symbol[1] = '\0';
            }
            if (symbol[0] == '*')        /* user R-atom for dummies */
               symbol[0] = 'R';
            strcpy(ap->atom_symbol, symbol);
            if (atext[0] != '\0') strcpy(ap->atext, atext);
            ap->color = atno; atno++;
            ap++;
            if (prev_atom > 0)
            {
               if (bond_type == NONE)
                  bp->bond_type = SINGLE;
               else
                  bp->bond_type = bond_type;
               if (bond_dir == NONE)    bp->value =  0.0;
               else if (bond_dir == UP) bp->value =  1.0;
               else                     bp->value = -1.0;
               bp->stereo_symbol = NONE;
               bp->atoms[0] = prev_atom;
               bp->atoms[1] = (int)(ap-result->atom_array);
// fprintf(stderr, "case(1): "); PrintREACCSBond(stderr, bp);
               bp++;
               bp[0] = bp[-1];
               bp[0].value = -bp[-1].value;
               bp[0].atoms[0] = bp[-1].atoms[1];
               bp[0].atoms[1] = bp[-1].atoms[0];
// fprintf(stderr, "case(2): "); PrintREACCSBond(stderr, bp);
               bp++;
            }
            bond_type = NONE; bond_dir = NONE;
            prev_atom = (int)(ap-result->atom_array);
            prev_aromatic = FALSE;
            break;
         case 'c': case 'b': case 'n':
         case 'o': case 'p': case 's':
            symbol[0] = (char)((*smiles) + 'A'-'a');
            symbol[1] = '\0';
            strcpy(ap->atom_symbol, symbol);
            ap->color = atno; atno++;
            ap++;
            if (prev_atom > 0)
            {
               if (prev_aromatic  &&  bond_type == NONE)
                  bp->bond_type = AROMATIC;
               else if (bond_type == NONE)
                  bp->bond_type = SINGLE;
               else
                  bp->bond_type = bond_type;
               if (bond_dir == NONE)    bp->value =  0.0;
               else if (bond_dir == UP) bp->value =  1.0;
               else                     bp->value = -1.0;
               bp->stereo_symbol = NONE;
               bp->atoms[0] = prev_atom;
               bp->atoms[1] = (int)(ap-result->atom_array);
// fprintf(stderr, "case(3): "); PrintREACCSBond(stderr, bp);
               bp++;
               bp[0] = bp[-1];
               bp[0].value = -bp[-1].value;
               bp[0].atoms[0] = bp[-1].atoms[1];
               bp[0].atoms[1] = bp[-1].atoms[0];
// fprintf(stderr, "case(4): "); PrintREACCSBond(stderr, bp);
               bp++;
            }
            bond_type = NONE; bond_dir = NONE;
            prev_atom = (int)(ap-result->atom_array);
            prev_aromatic = TRUE;
            break;
         case '-' :
            bond_type = SINGLE; bond_dir = NONE; break;
         case '/':
            isomeric_smiles = TRUE;
            bond_type = SINGLE; bond_dir = UP; break;
         case '\\':
            isomeric_smiles = TRUE;
            bond_type = SINGLE; bond_dir = DOWN; break;
         case '=' :
            bond_type = DOUBLE; bond_dir = NONE; break;
         case '#' :
            bond_type = TRIPLE; bond_dir = NONE; break;
         case ':' :
            bond_type = AROMATIC; bond_dir = NONE; break;
         case '%':
         case '1': case '2': case '3':
         case '4': case '5': case '6':
         case '7': case '8': case '9':
            if ((*smiles) == '%')
            {
               smiles++;
               iring = 0;
               idigit = 0;
               while (idigit < 2 && '0' <= (*smiles)  &&  (*smiles) <= '9')
               {
                  iring = 10*iring + (*smiles)-'0';
                  smiles++;
                  idigit++;
               }
               smiles--;
            }
            else
               iring = (*smiles)-'0';
            if (iring >= MAXRINGS)
            {
               ShowMessage("Too many rings","SMIToMOL");
// fprintf(stderr,"SMIToMOL(5)\n");
               goto error_cleanup;
            }
            if (open_rings[iring] > 0)  /* closing ring */
            {
               bph = result->bond_array+open_bonds[iring]; /* fixing bond */
               /* largest bond type gets precedence */
               if (open_aromatic[iring]  &&
                   prev_aromatic         &&
                   bond_type == NONE)
                  bph->bond_type = AROMATIC;
               else if (bond_type == NONE)
               {
                  if (bph->bond_type <= SINGLE) bph->bond_type = SINGLE;
               }
               else
               {
                  if (bph->bond_type <= bond_type) bph->bond_type = bond_type;
               }
                 /* bond dir might be defined on either side of ring closure */
               bph->atoms[0] = open_rings[iring];
               bph->atoms[1] = prev_atom;
               open_rings[iring] = 0;
               bp[0] = bph[0];
               bp[0].atoms[0] = bph[0].atoms[1];
               bp[0].atoms[1] = bph[0].atoms[0];
               if (bph->value == 0.0)   /* if closure defines direction */
               {                        /* then opening bond gets opposite */
                  if      (bond_dir == UP)   bph->value = -1.0;
                  else if (bond_dir == DOWN) bph->value =  1.0;
               }
               bp->value = -bph->value;
// fprintf(stderr, "case(6): "); PrintREACCSBond(stderr, bp);
// fprintf(stderr, "case(6a): "); PrintREACCSBond(stderr, bph);
               bp++;
            }
            else                                        /* opening ring */
            {
               open_rings[iring] = prev_atom;           /* saving state */
               open_aromatic[iring] = prev_aromatic;
               open_bonds[iring] = (int)(bp - result->bond_array);

               bp->bond_type = bond_type;               /* creating bond */
               bp->stereo_symbol = NONE;
               if (bond_dir == NONE)    bp->value =  0.0;
               else if (bond_dir == UP) bp->value =  1.0;
               else                     bp->value = -1.0;
// fprintf(stderr, "case(7): "); PrintREACCSBond(stderr, bp);
               bp++;
            }
            bond_type = NONE; bond_dir = NONE;
            break;
         case '[':
            smiles = SmilesToMDLAtom(smiles,
                                     symbol,
                                    &charge,
                                    &radical,
                                    &isotope,
                                    &hydrogen,
                                    &stereo,
                                    &mapping,
                                    &aromatic);
            if ((*symbol) == '\0')
            {
// fprintf(stderr,"SMIToMOL(2)\n");
               goto error_cleanup;
            }
            strcpy(ap->atom_symbol, symbol);
            if (0 == strcmp(symbol, "R#"))      /* R-Group */
            {
               AddNumProperty(result,
                              "M  RGP",
                             (int)(ap-result->atom_array)+1,
                             charge);
               ap->charge  = 0;
            }
            else
               ap->charge  = charge;
            ap->radical = radical;
            ap->mass_difference = isotope;
            if (isotope != 0) isomeric_smiles = TRUE;
            ap->query_H_count = hydrogen;
            ap->mapping = mapping;
            ap->stereo_parity = stereo;
            if (stereo != NONE) isomeric_smiles = TRUE;
            if (stereo != NONE) chiral = TRUE;
            ap->color = atno; atno++;
            ap++;
            if (prev_atom > 0)
            {
               if (aromatic  &&  prev_aromatic  && bond_type == NONE)
                  bp->bond_type = AROMATIC;
               else if (bond_type == NONE)
                  bp->bond_type = SINGLE;
               else
                  bp->bond_type = bond_type;
               if (bond_dir == NONE)    bp->value =  0.0;
               else if (bond_dir == UP) bp->value =  1.0;
               else                     bp->value = -1.0;
               bp->stereo_symbol = NONE;
               bp->atoms[0] = prev_atom;
               bp->atoms[1] = (int)(ap-result->atom_array);
// fprintf(stderr, "case(8): "); PrintREACCSBond(stderr, bp);
               bp++;
               bp[0] = bp[-1];
               bp[0].value = -bp[-1].value;
               bp[0].atoms[0] = bp[-1].atoms[1];
               bp[0].atoms[1] = bp[-1].atoms[0];
// fprintf(stderr, "case(9): "); PrintREACCSBond(stderr, bp);
               bp++;
            }
            bond_type = NONE; bond_dir = NONE;
            prev_atom = (int)(ap-result->atom_array);
            prev_aromatic = aromatic;
            if (stereo != NONE  &&  hydrogen == AT_LEAST_1) /* add hydrogen */
            {
               bp->bond_type = SINGLE;
               bp->value = 0.0;
               bp->stereo_symbol = NONE;
               bp->atoms[0] = prev_atom;
               bp->atoms[1] = prev_atom+1;
// fprintf(stderr, "case(10): "); PrintREACCSBond(stderr, bp);
               bp++;
               bp[0] = bp[-1];
               bp[0].value = -bp[-1].value;
               bp[0].atoms[0] = bp[-1].atoms[1];
               bp[0].atoms[1] = bp[-1].atoms[0];
// fprintf(stderr, "case(11): "); PrintREACCSBond(stderr, bp);
               bp++;
               strcpy(ap->atom_symbol, "H");
               (ap-1)->query_H_count = NONE;
               ap++;
            }
            break;
         case '(':
            branch_stack[top_of_stack] = prev_atom;
            aromatic_stack[top_of_stack] = prev_aromatic;
            top_of_stack++;
            break;
         case ')':
            top_of_stack--;
            prev_atom = branch_stack[top_of_stack];
            prev_aromatic = aromatic_stack[top_of_stack];
            break;
         case '.':
            if (isdigit(smiles[1]) ||  smiles[1] == '%')
            {
               bond_type = ANY_BOND; bond_dir = NONE;
// fprintf(stderr, "setting bond_type to AROMATIC at smiles = '%s'\n", smiles);
            }
            else
            {
               prev_atom = 0;
               prev_aromatic = FALSE;
               bond_type = NONE; bond_dir = NONE;
            }
            break;
         default :              /* unexpected character */
// fprintf(stderr,"SMIToMOL(3)\n");
            goto error_cleanup;
      }
      smiles++;
   }
// fprintf(stderr,"2\n"); Fortify_ListAllMemory();
   if (chiral) result->chiral_flag = TRUE;
   else        result->chiral_flag = FALSE;
   result->n_atoms = (int)(ap-result->atom_array);
   result->n_bonds = (int)(bp-result->bond_array);

   /* here, we have bonds duplicated to point in both directions.  */
   /* This is necessary to convert stereochmistry conventions from */
   /* SMILES (which is bond based) to MOL-Files which is set to be */
   /* numbering based.                                             */

   if (result->n_atoms == 0)
   {
// fprintf(stderr,"SMIToMOL(4)\n");
      goto error_cleanup;
   }
   for (i=0; i<MAXRINGS; i++)
      if (open_rings[i] != 0)
      {
// fprintf(stderr,"SMIToMOL(5)\n");
         goto error_cleanup;
      }
// fprintf(stderr,"SMIToMOL(5a)\n");

         /* convert SMILES atom stereo parities to atom number parities */
   for (i=0, ap=result->atom_array; i<result->n_atoms; i++, ap++)
   {
      if (ap->stereo_parity == NONE) continue;
      n_ligands = 0;
      for (j=0, bp=result->bond_array; j<result->n_bonds; j++, bp++)
         if (bp->atoms[0] == i+1)
         {
            ligands[n_ligands] = bp->atoms[1];
            n_ligands++;
            if (n_ligands > 4) break;
         }
      if (j < result->n_bonds  ||  n_ligands < 4)/* wrong number of ligands */
      {
         ap->stereo_parity = NONE;
         continue;
      }
      parity = 1;
      for (j=1; j<n_ligands; j++)
         for (k=j-1; k>=0; k--)
            if (ligands[k] > ligands[k+1])
            {
               h = ligands[k]; ligands[k] = ligands[k+1]; ligands[k+1] = h;
               parity *= (-1);
            }
      if (parity < 0)      /* odd permutation -> switch atom stereo parity */
         ap->stereo_parity = 3 - ap->stereo_parity;
   }
// fprintf(stderr,"3\n"); // Fortify_ListAllMemory();
   /* convert SMILES bond stereo parities to atom numbering base  */
   /* stereo information is coded in CIS_MASK or TRANS_MASK ORed  */
   /* into the stereo_symbol field of the double bond. This mask  */
   /* is assumed to be removed when the layout is made.           */
   for (i=0, bp=result->bond_array; i<result->n_bonds; i++, bp++)
      if (bp->bond_type == DOUBLE)
      {
         at1 = bp->atoms[0]; at2 = bp->atoms[1];
         if (at1 > at2) continue;
         nlig1 = nlig2 = 0;
         /* collect bonds sprouting off the double bond */
         for (j=0, bph=result->bond_array; j<result->n_bonds; j++, bph++)
         {
            if (bph->atoms[0] == at1  &&  bph->atoms[1] != at2)
            {
               if (bph->bond_type == AROMATIC)  /* make sure the algorithm does not invent stereochemistry during denormalization */
               {
                  bp->stereo_symbol = CIS_TRANS_EITHER;
                  break;
               }
               if (bph->bond_type != SINGLE)    /* ignoring allenes */
               {
                  bp->stereo_symbol = NONE;
                  break;
               }
               if (nlig1 >= 2) break;
               lig1[nlig1] = j; nlig1++;
            }
            if (bph->atoms[0] == at2  &&  bph->atoms[1] != at1)
            {
               if (bph->bond_type == AROMATIC)  /* make sure the algorithm does not invent stereochemistry during denormalization */
               {
                  bp->stereo_symbol = CIS_TRANS_EITHER;
                  break;
               }
               if (bph->bond_type != SINGLE)    /* ignoring allenes */
               {
                  bp->stereo_symbol = NONE;
                  break;
               }
               if (nlig2 >= 2) break;
               lig2[nlig2] = j; nlig2++;
            }
         }
         if (j < result->n_bonds)  /* ignore non-fatal problems */
            continue;

         if (nlig1 == 0  ||  nlig2 == 0)        /* terminal double bond */
            continue;

         /* Check for EITHER double bonds */
         if (nlig1 == 1  &&  result->bond_array[lig1[0]].value == 0.0)
         {
// fprintf(stderr,"CIS_TRANS_EITHER(1)\n");
            bp->stereo_symbol = CIS_TRANS_EITHER;
            continue;
         }
         if (nlig2 == 1  &&  result->bond_array[lig2[0]].value == 0.0)
         {
// fprintf(stderr,"CIS_TRANS_EITHER(2)\n");
            bp->stereo_symbol = CIS_TRANS_EITHER;
            continue;
         }
         if (nlig1 == 2  &&
             result->bond_array[lig1[0]].value == 0.0  &&
             result->bond_array[lig1[1]].value == 0.0)
         {
// fprintf(stderr,"CIS_TRANS_EITHER(3)\n");
            bp->stereo_symbol = CIS_TRANS_EITHER;
            continue;
         }
         if (nlig2 == 2  &&
             result->bond_array[lig2[0]].value == 0.0  &&
             result->bond_array[lig2[1]].value == 0.0)
         {
// fprintf(stderr,"CIS_TRANS_EITHER(4)\n");
            bp->stereo_symbol = CIS_TRANS_EITHER;
            continue;
         }

         /* make sure that lower numbered ligand comes first */
         if (nlig1 == 2  &&
             result->bond_array[lig1[0]].atoms[0] +        /* This assumes,      */
             result->bond_array[lig1[0]].atoms[1] >        /* correctly, that    */
             result->bond_array[lig1[1]].atoms[0] +        /* one atom is shared */
             result->bond_array[lig1[1]].atoms[1])

         {
            h = lig1[0]; lig1[0] = lig1[1]; lig1[1] = h;
         }
         if (nlig2 == 2  &&
             result->bond_array[lig2[0]].atoms[0] +        /* This assumes,      */
             result->bond_array[lig2[0]].atoms[1] >        /* correctly, that    */
             result->bond_array[lig2[1]].atoms[0] +        /* one atom is shared */
             result->bond_array[lig2[1]].atoms[1])

         {
            h = lig2[0]; lig2[0] = lig2[1]; lig2[1] = h;
         }

         /* Check for consistency of bond values */
         if (nlig1 == 2  &&
             result->bond_array[lig1[0]].value ==
               result->bond_array[lig1[1]].value)
         {
// fprintf(stderr,"CIS_TRANS_EITHER(5)\n");
            bp->stereo_symbol = CIS_TRANS_EITHER;
            continue;
         }
         if (nlig2 == 2  &&
             result->bond_array[lig2[0]].value ==
               result->bond_array[lig2[1]].value)
         {
// fprintf(stderr,"CIS_TRANS_EITHER(6)\n");
            bp->stereo_symbol = CIS_TRANS_EITHER;
            continue;
         }

         value1 = result->bond_array[lig1[0]].value;
         if (nlig1 == 2) value1 -= result->bond_array[lig1[1]].value;

         value2 = result->bond_array[lig2[0]].value;
         if (nlig2 == 2) value2 -= result->bond_array[lig2[1]].value;

         if (value1*value2 < 0)
            bp->stereo_symbol |= TRANS_MASK;
         else
            bp->stereo_symbol |= CIS_MASK;
      }

   for (i=0, j=0; i<result->n_bonds; i++)  /* remove dupl. bonds */
      if (result->bond_array[i].atoms[0] < result->bond_array[i].atoms[1])
      {
         result->bond_array[j] = result->bond_array[i]; j++;
      }
   result->n_bonds = j;

             /* ignore stereo double bonds if not isomeric SMILES */
   if (!isomeric_smiles)
      for (i=0, bp=result->bond_array; i<result->n_bonds; i++, bp++)
         if (bp->bond_type == DOUBLE) bp->stereo_symbol = NONE;

   atnos = TypeAlloc(result->n_atoms, int);
   for (i=0; i<result->n_atoms; i++)
      atnos[i] = result->atom_array[i].color;

// FORTIFY CheckFortifyMemory();
// fprintf(stderr,"4\n"); // Fortify_ListAllMemory();
   if (!CheckNeighbourhood(result))
   {
// fprintf(stderr,"4a\n");
// CheckFortifyMemory();
      FreeMolecule(result);
// fprintf(stderr,"4b\n");
// CheckFortifyMemory();
      return ((struct reaccs_molecule_t *)NULL);
   }
// if (result->n_atoms > 50) PrintREACCSMolecule(stderr, result, "after parsing");

// fprintf(stderr,"4c\n"); CheckFortifyMemory();
   if (!(flags & EXPECT_SMARTS)) ComputeSmilesValenceAndHCount(result);
// fprintf(stderr,"5\n"); // Fortify_ListAllMemory();
   if (flags & EXPECT_SMARTS) FixBisaryl(result);
   if (!(flags & EXPECT_SMARTS)) CCTDenormalize(result);
   for (i=0; i<result->n_atoms; i++)
      result->atom_array[i].dummy1 = result->atom_array[i].dummy2 = 0;
// fprintf(stderr,"6\n"); // Fortify_ListAllMemory();

   /* fix hydrogen count display */
   old_H_count = TypeAlloc(result->n_atoms+1, int);
   new_H_count = TypeAlloc(result->n_atoms+1, int);
   for (i=0; i<result->n_atoms; i++)
   {
      old_H_count[i+1] = result->atom_array[i].query_H_count;
      result->atom_array[i].query_H_count = 0;
   }
   ComputeImplicitH(result, new_H_count);
// fprintf(stderr,"7\n"); // Fortify_ListAllMemory();
   for (i=0; i<result->n_atoms; i++)
      if (old_H_count[i+1]   != 0  &&
          new_H_count[i+1]+1 != old_H_count[i+1])
      {
// fprintf(stderr, "%s atom %d has old_H_count=%d, new_H_count=%d, and query_H_count=%d\n", result->atom_array[i].atom_symbol, i+1, old_H_count[i+1], new_H_count[i+1], result->atom_array[i].query_H_count);
         if (AtomSymbolMatch(result->atom_array[i].atom_symbol, "C,N,O")  &&
             new_H_count[i+1]==old_H_count[i+1])
         {      // one hydrogen less than default => doublet radical
             result->atom_array[i].radical = DOUBLET;
         }
         else if (AtomSymbolMatch(result->atom_array[i].atom_symbol, "C,N")  &&
             new_H_count[i+1]==old_H_count[i+1]+1)
         {      // one hydrogen less than default => doublet radical
             result->atom_array[i].radical = TRIPLET;
         }
         else   // H count mismatch of unknown origin
             result->atom_array[i].query_H_count = old_H_count[i+1];
      }
   if (old_H_count) MyFree((char *)old_H_count);
   if (new_H_count) MyFree((char *)new_H_count);

// fprintf(stderr,"8\n"); Fortify_ListAllMemory();
   if (flags & DO_LAYOUT)
   {
      RecolorMolecule(result);
      mp = result;

      if (coordp)
      {
         cp = coordp;
         atno = 0;
         /* parse coordinate string */
         for (;;)
         {
            atno++; novalue = FALSE;
            if (1 != sscanf(cp, "%lf", &x))        /* no value found */
            {
               x = 0.0; novalue = TRUE;
               /* skip field */
            }
            while ((*cp) && (*cp) != ',') cp++; if (*cp) cp++;
            if (1 != sscanf(cp, "%lf", &y))        /* no value found */
            {
               y = 0.0; novalue = TRUE;
               /* skip field */
            }
            while ((*cp) && (*cp) != ',') cp++; if (*cp) cp++;

            if (!novalue)
               for (i=0; i<mp->n_atoms; i++)
                  if (atnos[i] == atno)
                  {
                      mp->atom_array[i].x = x;
                      mp->atom_array[i].y = y;
                      mp->atom_array[i].color = KEEP_POSITION;
                      break;
                  }

            if (*cp == '\0') break;        /* end of coordinate string */
            if (atno > mp->n_atoms) break;        /* safeguard agains too long input */
         }

         /* fuse connected fragments with coordinates */
         /* TO BE IMPLEMENTED */

         for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
            if (mp->atom_array[bp->atoms[0]-1].color & KEEP_POSITION  &&
                mp->atom_array[bp->atoms[1]-1].color & KEEP_POSITION)
               bp->bond_type |= DONT_FLIP_BOND;
      }
      result = LayoutMolecule(mp);
      if (0  &&  (flags & DROP_TRIVIAL_HYDROGENS)  &&  coordp)
      {
         MakeSomeHydrogensImplicit(result,
                                   NON_STEREO_HYDROGENS |
                                   ANCILLARY_STEREO |
                                   NO_QUERY_H_COUNT);
      }
      FreeMolecule(mp);
   }

   for (i=0; i<result->n_atoms; i++)
      result->atom_array[i].color = atnos[i];
   MyFree((char *)atnos);

// fprintf(stderr,"9\n"); Fortify_ListAllMemory();
   strncpy(result->program_name, "NSMI2MOL", NPROGNAME);
   result->program_name[NPROGNAME] = '\0';
// if (result->n_atoms > 50) PrintREACCSMolecule(stderr, result, "EXIT");
   return (result);

error_cleanup:
   FreeMolecule(result);
   return (NULL);
#undef SLASH
#undef BACK_SLASH
}

struct reaccs_reaction_t *SMIToRXN(char *smiles)
/*
 * Converts a reactions smiles to a REACCS reaction data structure.
 *
 * It uses the now abandoned '+' notation to reconstruct the
 * different molecules on reactant and product side of the >..> symbols.
 * This translation should transparently also work with reactions that
 * don't use the '+' of couse.
 */
{
   struct reaccs_molecule_t *tmpmol;
   struct reaccs_molecule_t *reac;
   struct reaccs_molecule_t *prod;
   struct reaccs_reaction_t *rxn;
   char *reac_stop, *prod_start, *buffer;
   char *cp, *cph;
   int in_symbol;
   int nmol;
   int i;

   if (!QuickCheck(smiles))
       return ((struct reaccs_reaction_t *)NULL);

   reac_stop  = strchr(smiles, '>');        /* end of reactant SMILES */
   prod_start = strrchr(smiles, '>');       /* start of product SMILES */

   if (!reac_stop || !prod_start || prod_start <= reac_stop)
   {
      return ((struct reaccs_reaction_t *)NULL);
   }

   buffer = (char *)calloc(strlen(smiles)+1, sizeof(char));
   cp = smiles;
   /* First, we collect the reactants. */
   reac = (struct reaccs_molecule_t *)NULL;
   nmol = 0;
   for(;;)
   {
      /* Fill buffer[] with the SMILES of the reactant. */
      cph = buffer;
      for (cph = buffer, in_symbol = FALSE;
           (*cp) != '\0'  &&  (*cp) != '>';
           (*cph++) = (*cp++))
      {
         if (cp[0] == '[')        /* entering atom symbol */
            in_symbol = TRUE;
         else if (cp[0] == ']')
            in_symbol = FALSE;
         else if (cp[0] == '+'  &&  !in_symbol)
            break;
      }
      cph[0] = '\0';
      tmpmol = SMIToMOL(buffer, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
      if (!tmpmol  ||  !cp[0])
      {
         free(buffer);
         FreeMoleculeList(tmpmol);
         FreeMoleculeList(reac);
         return ((struct reaccs_reaction_t *)NULL);
      }
      /* patch for MOL file output to ISISDraw */
      for (i=0; i<tmpmol->n_atoms; i++)
      {
         tmpmol->atom_array[i].x *= 0.5;
         tmpmol->atom_array[i].y *= 0.5;
      }
      sprintf(tmpmol->name, "reactant%d", ++nmol);
      strcpy(tmpmol->program_name, "SMIToRXN");
      strncpy(tmpmol->comment, buffer, MDL_MAXLINE); tmpmol->comment[MDL_MAXLINE] = '\0';
      tmpmol->next = reac; reac = tmpmol;

      if ((*cp++) == '>') break;
   }

   /* Now, we skip any reagents. */
   while ((*cp) && (*cp) != '>')
      cp++;
   if (!(*cp))
   {
      free(buffer);
      FreeMoleculeList(reac);
      return ((struct reaccs_reaction_t *)NULL);
   }
   cp++;

   /* Finally, we collect the reactants. */
   prod = (struct reaccs_molecule_t *)NULL;
   nmol = 0;
   for(;;)
   {
      /* Fill buffer[] with the SMILES of the reactant. */
      cph = buffer;
      for (cph = buffer, in_symbol = FALSE;
           (*cp) != '\0';
           (*cph++) = (*cp++))
      {
         if (cp[0] == '[')        /* entering atom symbol */
            in_symbol = TRUE;
         else if (cp[0] == ']')
            in_symbol = FALSE;
         else if (cp[0] == '+'  &&  !in_symbol)
            break;
      }
      cph[0] = '\0';
      tmpmol = SMIToMOL(buffer, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);
      if (!tmpmol  ||  cp[0] == '>')
      {
         free(buffer);
         FreeMoleculeList(tmpmol);
         FreeMoleculeList(reac);
         FreeMoleculeList(prod);
         return ((struct reaccs_reaction_t *)NULL);
      }
      /* patch for MOL file output to ISISDraw */
      for (i=0; i<tmpmol->n_atoms; i++)
      {
         tmpmol->atom_array[i].x *= 0.5;
         tmpmol->atom_array[i].y *= 0.5;
      }
      sprintf(tmpmol->name, "product%d", ++nmol);
      strcpy(tmpmol->program_name, "SMIToRXN");
      strncpy(tmpmol->comment, buffer, MDL_MAXLINE); tmpmol->comment[MDL_MAXLINE] = '\0';
      tmpmol->next = prod; prod = tmpmol;
      if (cp[0] == '\0') break;
      cp++;
   }

   rxn = NewReaction();
   rxn->n_reactants = 0; rxn->n_products  = 0;
   /* rxn->reactants = reac; */
   for (rxn->n_reactants = 0; reac != NULL; rxn->n_reactants++)
   {
      tmpmol = reac->next; reac->next = rxn->reactants;
      rxn->reactants = reac; reac = tmpmol;
   }
   /* rxn->products  = prod; */
   for (rxn->n_products  = 0; prod != NULL; rxn->n_products++)
   {
      tmpmol = prod->next; prod->next = rxn->products;
      rxn->products = prod; prod = tmpmol;
   }

   PerceiveReactionCenter(rxn);
   strncpy(rxn->comment, smiles, MDL_MAXLINE); rxn->comment[MDL_MAXLINE] = '\0';

   return (rxn);
}
