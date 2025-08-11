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
/*      File:           utilities.c                                     */
/*                                                                      */
/*      Purpose:        This file implements some utility functions     */
/*                      dealing with REACCS data structures. They don't */
/*                      perform any I/O except for debugging purposes.  */
/*                      I/O-functions are in reaccsio.c.                */
/*                                                                      */
/*      History:        04-Jan-93       Completed ANSI prototypes.      */
/*                                                                      */
/************************************************************************/

#include "utilities.h"

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "graph.h"
#include "local.h"
#include "reaccs.h"
#include "set.h"
#include "symbol_lists.h"

int CountSTextLines(struct stext_line_t *stext_lines)
/*
 * Counts the length of the list stext_lines.
 */
{
   int i;

   for (i=0; !IsNULL(stext_lines); i++)
   {
      stext_lines = stext_lines->next;
   }

   return (i);
}

void FreeSTextLines(struct stext_line_t *stext_lines)
/*
 * Frees the storage allocated to the STEXT lines stext_lines.
 */
{
   struct stext_line_t *hp;

   while (stext_lines)
   {
      hp = stext_lines->next; MyFree((char *)stext_lines); stext_lines = hp;
   }
}

struct stext_line_t *CopySTextLines(struct stext_line_t *stext_lines)
/*
 * Copies the STEXT lines pointed to by stext_lines.
 */
{
   struct stext_line_t *hp, *result;

   result = (struct stext_line_t *)NULL; /* create list in reverse order */
   while (stext_lines)
   {
      hp = TypeAlloc(1,struct stext_line_t);
      *hp = *stext_lines; hp->next = result; result = hp;
      stext_lines = stext_lines->next;
   }

   /* restore order of list */
   stext_lines = result; result = (struct stext_line_t *)NULL;
   while (stext_lines)
   {
      hp = stext_lines; stext_lines = hp->next;
      hp->next = result; result = hp;
   }

   return (result);
}

int CountPropLines(struct prop_line_t *prop_lines)
/*
 * Counts the length of the list prop_lines.
 */
{
   int i;

   for (i=0; !IsNULL(prop_lines); i++)
      prop_lines = prop_lines->next;

   return (i);
}

void FreePropLines(struct prop_line_t *prop_lines)
/*
 * Frees the storage allocated to the property lines prop_lines.
 */
{
   struct prop_line_t *hp;

   while (prop_lines)
   {
      hp = prop_lines->next; MyFree((char *)prop_lines); prop_lines = hp;
   }
}

struct prop_line_t *CopyPropLines(struct prop_line_t *prop_lines)
/*
 * Copies the property lines pointed to by prop_lines.
 */
{
   struct prop_line_t *hp, *result;

   result = (struct prop_line_t *)NULL; /* create list in reverse order */
   while (prop_lines)
   {
      hp = TypeAlloc(1,struct prop_line_t);
      *hp = *prop_lines; hp->next = result; result = hp;
      prop_lines = prop_lines->next;
   }

   /* restore order of list */
   prop_lines = result; result = (struct prop_line_t *)NULL;
   while (prop_lines)
   {
      hp = prop_lines; prop_lines = hp->next;
      hp->next = result; result = hp;
   }

   return (result);
}

void AddNumProperty(struct reaccs_molecule_t *mp,
                    char *tag,
                    int atno,
                    int prop)
/*
 * Adds the atom property prop for atom atno to the list of property
 * lines in mp using the prefix tag[].
 */
{
   struct prop_line_t *hp;

   hp = TypeAlloc(1,struct prop_line_t);
   sprintf(hp->text, "%s  1 %3d %3d", tag, atno, prop);
   hp->next = mp->prop_lines;
   mp->prop_lines = hp;
   mp->n_props++;
}

int GetNumProperty(struct prop_line_t *prop_lines,
                   char *tag,
                   int atno,
                   int *prop)
/*
 * Fetches the property identified by tag[] for atom atno from the
 * property strings linked at prop_lines. Returns TRUE if a property
 * entry was found for that atom and FALSE otherwise. *prop is set to
 * the value found.
 */
{
   int nentries, n;
   int tmp_ats[8];
   int tmp_vals[8];

   while (prop_lines)
   {
      if (strncmp(prop_lines->text, tag, strlen(tag)) == 0)
      {
	 sscanf(prop_lines->text+strlen(tag), "%3d", &nentries);
	 n = sscanf(prop_lines->text+strlen(tag)+3,
		   " %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
		   &tmp_ats[0],  &tmp_vals[0],  &tmp_ats[1],  &tmp_vals[1],
		   &tmp_ats[2],  &tmp_vals[2],  &tmp_ats[3],  &tmp_vals[3],
		   &tmp_ats[4],  &tmp_vals[4],  &tmp_ats[5],  &tmp_vals[5],
		   &tmp_ats[6],  &tmp_vals[6],  &tmp_ats[7],  &tmp_vals[7]);
	 if (n != nentries*2)
	 {
	    ShowMessageI("n = %d","GetNumProperty",n);
	    ShowMessageI("nentries = %d","GetNumProperty",nentries);
	    ShowMessageS("buffer = '%s'\n","GetNumProperty",prop_lines->text);
	 }
	 for (n=0; n<nentries; n++)
	    if(tmp_ats[n] == atno)
	    {
	       (*prop) = tmp_vals[n];
	       return (TRUE);
	    }
      }
      prop_lines = prop_lines->next;
   }

   return (FALSE);
}

void FreeMoleculeChildObjects(struct reaccs_molecule_t *mp)
/*
 * Deallocates the storage for the REACCS molecule *mp leaving the struct itself allocated.
 */
{
   if ((char *)mp != NULL)
   {
      FreeSymbolLists(mp->symbol_lists);
      mp->symbol_lists = (struct symbol_list_t *) NULL;
      FreeSTextLines(mp->stext_lines);
       mp->stext_lines = (struct stext_line_t *)NULL;
      FreePropLines(mp->prop_lines);
       mp->prop_lines = (struct prop_line_t *)NULL;
      if ((char *)mp->atom_array != NULL) MyFree((char *)mp->atom_array);
      if ((char *)mp->bond_array != NULL) MyFree((char *)mp->bond_array);
   }
}

void FreeMolecule(struct reaccs_molecule_t *mp)
/*
 * Deallocates the storage for the REACCS molecule *mp.
 */
{
   if ((char *)mp != NULL)
   {
      FreeMoleculeChildObjects(mp);
      MyFree((char *)mp);
   }
}

struct reaccs_molecule_t *NewMolecule(int n_atoms, int n_bonds)
/*
 * Allocates memory for a molecule with n_atoms atoms and n_bonds bonds.
 */
{
   struct reaccs_molecule_t *mp;

   mp = TypeAlloc(1,struct reaccs_molecule_t);
   mp->name[0] = mp->user_initials[0] = mp->program_name[0] = '\0';
   mp->date[0] = mp->time[0] = '\0';
   strcpy(mp->dimensionality, "2D");    /* default is 2D file */
   mp->scale1 = 1; mp->scale2 = 1.0;
   mp->registry_number = 0;
   mp->comment[0] = '\0';
   mp->n_atom_lists = 0;
   mp->symbol_lists = (struct symbol_list_t *)NULL;
   mp->stext_lines = (struct stext_line_t *)NULL;
   mp->prop_lines = (struct prop_line_t *)NULL;
   mp->dummy1 = 0; mp->chiral_flag = ACHIRAL;
   mp->n_props = 0; mp->version[0] = '\0';
   strcpy(mp->version, "V2000");
   mp->atom_array = TypeAlloc(n_atoms,struct reaccs_atom_t);
   mp->n_atoms = n_atoms;
   mp->bond_array = TypeAlloc(n_bonds,struct reaccs_bond_t);
   mp->n_bonds = n_bonds;
   mp->next = (struct reaccs_molecule_t *)NULL;
   mp->color = 0;

   return (mp);
}

void FreeMoleculeList(struct reaccs_molecule_t *mp)
/*
 * Deallocates the storage occupied by the molecules listed
 * in mp, mp->next, ...
 */
{
   struct reaccs_molecule_t *mph;

   while (!IsNULL(mp))
   {
      mph = mp->next; FreeMolecule(mp); mp = mph;
   }
}

struct reaccs_reaction_t *NewReaction()
/*
 * Allocates memory for an empty reaction.
 */
{
   struct reaccs_reaction_t *result;

   result = TypeAlloc(1,struct reaccs_reaction_t);
   result->name[0] = '\0'; result->user_initials[0] = '\0';
   result->program_name[0] = '\0';
   result->date[0] = '\0'; result->time[0] = '\0';
   result->registry_number = 0; result->comment[0] = '\0';
   result->molecularity = 0; result->next = (struct reaccs_reaction_t *)NULL;

   return (result);
}

void FreeReaction(struct reaccs_reaction_t *rp)
/*
 * Deallocates the storage for the REACCS reaction *rp.
 */
{
   FreeMoleculeList(rp->reactants);
   FreeMoleculeList(rp->products);

   MyFree((char *)rp);
}

struct reaccs_molecule_t *CopyMolecule(struct reaccs_molecule_t *mp)
/*
 * Allocates space for a new molecule and copies the contents of *mp
 * to that space. The function returns a pointer to the new molecule.
 */
{
   unsigned int i;
   struct reaccs_molecule_t *result;

   result = TypeAlloc(1, struct reaccs_molecule_t);/* allocate space for */
   if ((char *)result == (char *)NULL)          /* molecule data        */
      return(result);
   (*result) = (*mp);       /* Copy direct data */
   result->stext_lines  = CopySTextLines(mp->stext_lines);
   result->prop_lines   = CopyPropLines(mp->prop_lines);
                                      /* allocating space for atom entries */
   result->atom_array = TypeAlloc(mp->n_atoms,struct reaccs_atom_t);
   for (i=0; i<mp->n_atoms; i++)
      result->atom_array[i] = mp->atom_array[i];
                                      /* allocating space for bond entries */
   result->bond_array = TypeAlloc(mp->n_bonds,struct reaccs_bond_t);
   for (i=0; i<mp->n_bonds; i++)
      result->bond_array[i] = mp->bond_array[i];

   result->symbol_lists = CopySymbolLists(mp->symbol_lists);
   result->next = (struct reaccs_molecule_t *)NULL;

   return(result);
}

int AddAtomToMolecule(struct reaccs_molecule_t *mp,
                      double x, double y, double z,
                      char *symbol)
/*
 * Adds a new atom to the structure pointed to by *mp. This atom is
 * not yet connected. It is placed at the given coordinates x, y, z,
 * and is of type *symbol. The function returns the number of the new
 * atom.
 */
{
   unsigned int i;
   struct reaccs_atom_t *new_atoms;

   new_atoms = TypeAlloc(mp->n_atoms+1,struct reaccs_atom_t);
   for (i=0; i<mp->n_atoms; i++)
      new_atoms[i] = mp->atom_array[i];
   new_atoms[mp->n_atoms].x = x;
   new_atoms[mp->n_atoms].y = y;
   new_atoms[mp->n_atoms].z = z;
   strncpy(new_atoms[mp->n_atoms].atom_symbol, symbol, 3);

   MyFree((char *)mp->atom_array); mp->atom_array = new_atoms;
   mp->n_atoms++;

   return(mp->n_atoms);
}

void AddBondToMolecule(struct reaccs_molecule_t *mp,
                       int from_atom, int to_atom,
                       int bond_type, int stereo_symbol)
/*
 * Adds a new bond from from_atom to to_atom to the structure pointed
 * to by *mp. The bond get the type bond_type and the stereosymbol
 * stereo_symbol.
 */
{
   unsigned int i;
   struct reaccs_bond_t *new_bonds;

   new_bonds = TypeAlloc(mp->n_bonds+1,struct reaccs_bond_t);
   for (i=0; i<mp->n_bonds; i++)
      new_bonds[i] = mp->bond_array[i];
   new_bonds[mp->n_bonds].atoms[0] = from_atom;
   new_bonds[mp->n_bonds].atoms[1] = to_atom;
   new_bonds[mp->n_bonds].bond_type = bond_type;
   new_bonds[mp->n_bonds].stereo_symbol = stereo_symbol;

   MyFree((char *)mp->bond_array); mp->bond_array = new_bonds;
   mp->n_bonds++;
}

struct reaccs_reaction_t *CopyReaction(struct reaccs_reaction_t *rp)
/*
 * Allocates space for a new reaction and copies the contents of *mp
 * to that space. The function returns a pointer to the new reaction.
*/
{
   struct reaccs_reaction_t *result;
   struct reaccs_molecule_t *mp, *mph;

   result = TypeAlloc(1,struct reaccs_reaction_t);
   if ((char *)result == (char *)NULL)
      return(result);
   (*result) = (*rp);
   result->reactants = (struct reaccs_molecule_t *)NULL;
   result->products  = (struct reaccs_molecule_t *)NULL;
   result->next = (struct reaccs_reaction_t *)NULL;

   for (mph=rp->reactants; !IsNULL(mph); mph=mph->next)
   {
      mp = CopyMolecule(mph);
      if (IsNULL(mp))
      {
         FreeReaction(result);
         return((struct reaccs_reaction_t *)NULL);
      }
      mp->next = result->reactants; result->reactants = mp;
   }
   mp = (struct reaccs_molecule_t *)NULL; /* fix order of reactants */
   while (!IsNULL(result->reactants))
   {
      mph = result->reactants; result->reactants = mph->next;
      mph->next = mp; mp = mph;
   }
   result->reactants = mp;

   for (mph=rp->products; !IsNULL(mph); mph=mph->next)
   {
      mp = CopyMolecule(mph);
      if (IsNULL(mp))
      {
         FreeReaction(result);
         return((struct reaccs_reaction_t *)NULL);
      }
      mp->next = result->products; result->products = mp;
   }
   mp = (struct reaccs_molecule_t *)NULL; /* fix order of products */
   while (!IsNULL(result->products))
   {
      mph = result->products; result->products = mph->next;
      mph->next = mp; mp = mph;
   }
   result->products = mp;

   return(result);
}

int NewParity(unsigned int         atn,
              int                  old_parity,
              unsigned int         n_bonds,
              struct reaccs_bond_t bonds[],
              int                  renumbering[])
/*
 * Computes the new parity for atom number atn when the structure
 * is renumbered using renumbering[].
 */
{
   unsigned int i;
   int j, n, h, orient[4], index;

   for (index=0; index<4; index++) orient[index]=0; /* clear tuple */

   for (i=0, index=0; i<n_bonds && index<4; i++)    /* get neighbours */
   {
      if (bonds[i].atoms[0] == atn && renumbering[bonds[i].atoms[1]] != 0)
         orient[index++] = bonds[i].atoms[1];
      if (bonds[i].atoms[1] == atn && renumbering[bonds[i].atoms[0]] != 0)
         orient[index++] = bonds[i].atoms[0];
   }

   if (index < 3) return(NONE); /* too few ligands for stereocenter */

   for (i=1; i<index; i++)              /* sort neighbours in old numbering */
      for (j=i; j>0; j--)
         if (orient[j-1] < orient[j])
         {
            h = orient[j-1]; orient[j-1] = orient[j]; orient[j] = h;
         }

   for (i=0; i<index; i++)              /* renumber neighbours */
      orient[i] = renumbering[orient[i]];

   for (i=1,n=0; i<index; i++)          /* compute parity of renumbering */
      for (j=i; j>0; j--)
         if (orient[j-1] < orient[j])
         {
            n++;
            h = orient[j-1]; orient[j-1] = orient[j]; orient[j] = h;
         }

   if (n%2 == 0 || (old_parity!=ODD && old_parity!=EVEN))
      return (old_parity);
   else
      if (old_parity == ODD) return (EVEN);
      else                   return (ODD);

}

struct reaccs_molecule_t *ThreadQueryAtoms(struct reaccs_molecule_t *query)
/*
 * Re-sorts the atoms in query such that (if possible) for all
 * n < query->n_atoms, atom n+1 is connected to one from 1..n.
 *
 * This is needed to improve the recursive phase of substruture matching and
 * shall only be done for atoms that don't have complications like symbol
 * lists and the like. Flag-based stereochemistry is ignored!
 *
 * Note: any neighbourhood arrays may become invalid when this method is called.
 */
{
   int i, j, n;
   struct reaccs_atom_t atom;
   struct reaccs_bond_t *bp;

   if (!IsNULL(query->symbol_lists)) return query;
   // PrintREACCSMolecule(stderr, query, "");
   for (n=1; n<query->n_atoms; n++)
   {
      // test if this atom is connected to an atom with a lower index
      for (j=0, bp=query->bond_array; j<query->n_bonds; j++, bp++)
      {
         if (bp->atoms[0] == n+1  &&  bp->atoms[1] <= n) break;
         if (bp->atoms[1] == n+1  &&  bp->atoms[0] <= n) break;
      }
      if (j < query->n_bonds) continue; // is connected => next atom
      // search for atom beyond n+1 that is connected and swap it with n+1.
      for (i=n+1; i<query->n_atoms; i++)
      {
         for (j=0, bp=query->bond_array; j<query->n_bonds; j++, bp++)
         {
            if (bp->atoms[0] == i+1  &&  bp->atoms[1] <= n) break;
            if (bp->atoms[1] == i+1  &&  bp->atoms[0] <= n) break;
         }
         if (j < query->n_bonds)  // atom i+1 is connected => swap i+1 with n+1 and goto next atom
         {
// fprintf(stderr, "swapping atom %d with atom %d\n", n+1, i+1);
            // swap atoms
            atom = query->atom_array[i];
            query->atom_array[i] = query->atom_array[n];
            query->atom_array[n] = atom;
            // swap atom numbers of bonds
            for (j=0, bp=query->bond_array; j<query->n_bonds; j++, bp++)
            {
               if (bp->atoms[0] == i+1) bp->atoms[0] = n+1;
               else if (bp->atoms[0] == n+1) bp->atoms[0] = i+1;
               if (bp->atoms[1] == i+1) bp->atoms[1] = n+1;
               else if (bp->atoms[1] == n+1) bp->atoms[1] = i+1;
            }
            break;
         }
      }
   }
   // PrintREACCSMolecule(stderr, query, "");
   // exit (EXIT_FAILURE);
   return query;
}

void RemoveAtomFromMolecule(struct reaccs_molecule_t *mp, int atom)
/*
 * Removes atom from molecule *mp.
 *
 * The atoms are numbered starting from 1!
 * This function has only been tested with isolated atoms although it
 * should also work with multiply connected ones.
 */
{
   int* good_atoms;
   int* good_bonds;
   int i;

   if (mp == (struct reaccs_molecule_t *)NULL) return;
   if (mp->n_atoms <= 1) return;        /* do not remove last atom */
   /* set up arrays of atom and bond numbers to be deleted */
   good_atoms = TypeAlloc(mp->n_atoms+1, int);
   good_bonds = TypeAlloc(mp->n_bonds, int);
   for (i=1; i<mp->n_atoms+1; i++)
      good_atoms[i] = i != atom;
   for (i=0; i<mp->n_bonds; i++)
   {
      if (mp->bond_array[i].atoms[0] == atom  ||
          mp->bond_array[i].atoms[1] == atom)
         good_bonds[i] = FALSE;
      else
         good_bonds[i] = TRUE;
   }
   StripMolecule(mp, good_atoms, good_bonds);
   MyFree((char *) good_bonds);
   MyFree((char *) good_atoms);
}

void StripMolecule(struct reaccs_molecule_t *mp,
                   int good_atoms[], int good_bonds[])
/*
 * Deletes from molecule *mp all atoms i and bonds j for which
 * good_atoms[i] == FALSE, and good_bonds[j] == FALSE, respectively.
 * The affected stereosymbols are adjusted.
 * sizeof(good_atoms) = n_atoms+1 !!
 */
{
   unsigned int i, ii;
   int *new_number,n,h;
   struct symbol_list_t *slp, *slph;
   int had_parity;
   struct prop_line_t *known_prop_lines;
   struct prop_line_t *prop_lines, *hp;
   char *cp;
   int nrgp, ratom, rval;

   if (NULL == strstr(mp->name, "stripped"))
       strncat(mp->name," stripped",MAXNAME);
   strncpy(mp->user_initials,"BR",2);
   strncpy(mp->program_name,"StripMol",NPROGNAME);
   strncpy(mp->comment,"This file has been stripped",MDL_MAXLINE);

   new_number = TypeAlloc(mp->n_atoms+1,int);
   for (i=0; i<=mp->n_atoms; i++) new_number[i] = 0;
   had_parity = FALSE;
   for (i=0; i<mp->n_atoms; i++)
      if (mp->atom_array[i].stereo_parity != NONE)
         had_parity = TRUE;
   for (i=1,n=1; i<=mp->n_atoms; i++)           /* new_numbers maps old */
      if (good_atoms[i]) new_number[i] = n++;   /* to new numbers       */

   for (i=0; i<mp->n_bonds; i++)        /* delete stereo information */
      if (!good_bonds[i])               /* related to deleted bonds  */
      {
         mp->atom_array[mp->bond_array[i].atoms[0]-1].stereo_parity = NONE;
         mp->atom_array[mp->bond_array[i].atoms[1]-1].stereo_parity = NONE;
      }

   for (i=0,n=0; i<mp->n_bonds; i++)    /* update bond array */
      if (good_bonds[i])
         mp->bond_array[n++] = mp->bond_array[i];
   mp->n_bonds = n;

   for (i=0; i<mp->n_atoms; i++)        /* adjust stereo parities */
      if (mp->atom_array[i].stereo_parity != NONE)
         mp->atom_array[i].stereo_parity =
            NewParity(i+1,
                      mp->atom_array[i].stereo_parity,
                      mp->n_bonds, mp->bond_array,
                      new_number);

   /* Only do this if original moleule had parities perceived */
   if (had_parity)
   {
      h = FALSE;
      for (i=0; i<mp->n_atoms; i++)        /* adjust chiral flag */
         if (mp->atom_array[i].stereo_parity != NONE)
            h = TRUE;
      if (!h) mp->chiral_flag = FALSE;     /* no stereo center -> not chiral */

      for (i=0; i<mp->n_bonds; i++)                /* fix stereo bonds */
      {                                            /* stereo bond symbol */
         h = mp->bond_array[i].stereo_symbol;      /* makes sense only   */
         if (h == UP || h == DOWN)                 /* if at least one    */
         {                                         /* atom is a center   */
            h = mp->atom_array[mp->bond_array[i].atoms[0]-1].stereo_parity;
            if (h == EVEN || h == ODD)
               continue;

            h = mp->atom_array[mp->bond_array[i].atoms[1]-1].stereo_parity;
            if (h == EVEN || h == ODD)
               continue;

            mp->bond_array[i].stereo_symbol = NONE;
         }
      }
   }

   for (i=0,n=0; i<mp->n_atoms; i++)            /* squeeze atoms */
      if (good_atoms[i+1] != 0)
         mp->atom_array[n++] = mp->atom_array[i];
   mp->n_atoms = n;

   for (i=0; i<mp->n_bonds; i++)        /* renumber bond atom pairs */
   {
      mp->bond_array[i].atoms[0] = new_number[mp->bond_array[i].atoms[0]];
      mp->bond_array[i].atoms[1] = new_number[mp->bond_array[i].atoms[1]];
   }
   for (i=0,n=0; i<mp->n_bonds; i++)
      if (mp->bond_array[i].atoms[0] != 0 &&
          mp->bond_array[i].atoms[1] != 0)
         mp->bond_array[n++] = mp->bond_array[i];
   mp->n_bonds = n;

   if (mp->symbol_lists)
   {
      mp->n_atom_lists = 0;
      slp = (struct symbol_list_t *)NULL;
      while (mp->symbol_lists)
      {
         if (new_number[mp->symbol_lists->atom] == 0)  /* deleted atom */
         {
            slph = mp->symbol_lists->next;
            mp->symbol_lists->next = (struct symbol_list_t *)NULL;
            FreeSymbolLists(mp->symbol_lists);
            mp->symbol_lists = slph;
         }
         else                                   /* renumbered atom */
         {
            mp->symbol_lists->atom = new_number[mp->symbol_lists->atom];
            slph = mp->symbol_lists->next;
            mp->symbol_lists->next = slp;
            slp = mp->symbol_lists;
            mp->symbol_lists = slph;
            mp->n_atom_lists++;
         }
      }
      mp->symbol_lists = slp;
   }

   /* [TODO] Need to implement this more completely for toolkit use */
   if (mp->prop_lines)
   {
      known_prop_lines = (struct prop_line_t *)NULL;
      mp->n_props = 0;
      for (prop_lines = mp->prop_lines; !IsNULL(prop_lines); /* within loop */)
      {
         // renumber r-group properties
         if (0 == strncmp(prop_lines->text, "M  RGP", strlen("M  RGP")))
         {
            sscanf(prop_lines->text+strlen("M  RGP"), "%d", &nrgp);
            cp = prop_lines->text+strlen("M  RGP")+3;
            if (nrgp > 4)       // safeguard against read errors
            {
               ShowMessageS("Could not handle property line '%s'","StripMolecule", prop_lines->text);
               hp = prop_lines; prop_lines = prop_lines->next; hp->next = (struct prop_line_t *)NULL;
               FreePropLines(hp);
               continue;
            }
            for (ii=0; ii<nrgp; ii++)
            {
               sscanf(cp, "%d %d", &ratom, &rval);
               sprintf(cp, " %3d %3d", new_number[ratom], rval);
               cp += 2*(1+3);
            }
            hp = prop_lines; prop_lines = prop_lines->next;
            hp->next = known_prop_lines; known_prop_lines = hp;
            mp->n_props++;
         }
         else   // unknown properties
         {
            ShowMessageS("Could not handle property line '%s'","StripMolecule", prop_lines->text);
            hp = prop_lines; prop_lines = prop_lines->next; hp->next = (struct prop_line_t *)NULL;
            FreePropLines(hp);
         }
      }
      mp->prop_lines = known_prop_lines;
   }

   MyFree((char *)new_number);
}

struct reaccs_molecule_t *SplitBond(struct reaccs_molecule_t *mp,
                                    int ibond,
                                    char *atsym1, int charge1, int mdiff1,
                                    char *atsym2, int charge2, int mdiff2)
/*
 * Splits the bond with index ibond attaching new atoms with the given
 * characteristics to the residual. The coordinates of the new atoms will be
 * on top of the old ones. It may therefore be needed to convert to SMILES or use
 * SplitMolecule to separate and visualize the result.
 *
 * The change is affected in-place and the modified molecule is returned.
 */
{
    int i, ai1, ai2;
    struct reaccs_bond_t *new_bonds, *bp;
    struct reaccs_atom_t *new_atoms, *ap;

    // PrintREACCSMolecule(stderr,mp,"ORIGINAL");
    /* need room for one more bond */
    new_bonds = TypeAlloc(mp->n_bonds+1, struct reaccs_bond_t);
    for (i=0; i<mp->n_bonds; i++)
        new_bonds[i] = mp->bond_array[i];
    MyFree((char *)mp->bond_array);
    mp->bond_array = new_bonds;
    mp->n_bonds++;
    /* need room for two more atoms */
    new_atoms = TypeAlloc(mp->n_atoms+2, struct reaccs_atom_t);
    for (i=0; i<mp->n_atoms; i++)
        new_atoms[i] = mp->atom_array[i];
    MyFree((char *)mp->atom_array);
    mp->atom_array = new_atoms;
    mp->n_atoms+=2;
    /* Start modification */
    ai1 = mp->bond_array[ibond].atoms[0]-1;
    ai2 = mp->bond_array[ibond].atoms[1]-1;
    /* copy and update atom properties */
    ap = mp->atom_array+(mp->n_atoms-2);
    (*ap) = mp->atom_array[ai1];
    strncpy(ap->atom_symbol, atsym1, 3);
    ap->charge = charge1;
    ap->mass_difference = mdiff1;
    ap = mp->atom_array+(mp->n_atoms-1);
    (*ap) = mp->atom_array[ai2];
    strncpy(ap->atom_symbol, atsym2, 3);
    ap->charge = charge2;
    ap->mass_difference = mdiff2;
    /* copy bond properties */
    mp->bond_array[mp->n_bonds-1] = mp->bond_array[ibond];
    mp->bond_array[ibond].atoms[0] = mp->n_atoms;
    mp->bond_array[mp->n_bonds-1].atoms[1] = mp->n_atoms-1;
    mp->bond_array[mp->n_bonds-1].atoms[0] = ai1+1;
    // PrintREACCSMolecule(stderr,mp,"RESULT");
    return mp;
}

struct reaccs_molecule_t *SplitMolecule(struct reaccs_molecule_t *mp)
/*
 * Checks if mp contains more than one fragment, modifies *mp to be
 * one of them and returns a pointer to a structure containing the others
 */
{
   struct reaccs_molecule_t *result;
   int *marked_atoms, *marked_bonds;
   int i;
   int changed;

   marked_atoms = TypeAlloc(mp->n_atoms+1, int);
   marked_bonds = TypeAlloc(mp->n_bonds+1, int);

   for (i=0; i<=mp->n_atoms; i++) marked_atoms[i] = FALSE;
   for (i=0; i<mp->n_bonds; i++)  marked_bonds[i] = FALSE;
                                /* start graph search at atom one */
   marked_atoms[1] = TRUE;
   do
   {
      changed = FALSE;
      for (i=0; i<mp->n_bonds; i++)
         if (marked_atoms[mp->bond_array[i].atoms[0]] !=
             marked_atoms[mp->bond_array[i].atoms[1]])
         {
            marked_atoms[mp->bond_array[i].atoms[0]] = TRUE;
            marked_atoms[mp->bond_array[i].atoms[1]] = TRUE;
            changed = TRUE;
         }
   } while (changed);

   for (i=1; i<=mp->n_atoms; i++)
      if (!marked_atoms[i]) break;
   if (i == mp->n_atoms+1)              /* only one fragment -> return */
   {
      MyFree((char *)marked_atoms); MyFree((char *)marked_bonds);
      return((struct reaccs_molecule_t *) NULL);
   }

   for (i=0; i<mp->n_bonds; i++)        /* mark bonds of component */
      if (marked_atoms[mp->bond_array[i].atoms[0]])
         marked_bonds[i] = TRUE;
   result = CopyMolecule(mp);
   StripMolecule(result,marked_atoms,marked_bonds);
   for (i=1; i<=mp->n_atoms; i++) marked_atoms[i] = !marked_atoms[i];
   for (i=0; i<mp->n_bonds; i++)  marked_bonds[i] = !marked_bonds[i];
   StripMolecule(mp,marked_atoms,marked_bonds);

   MyFree((char *)marked_atoms); MyFree((char *)marked_bonds);

   return(result);
}

struct reaccs_molecule_t *SplitMoleculeList(struct reaccs_molecule_t *mp)
/*
 * Scans the molecule list pointed to by mp and returns a list of
 * of all the fragments found therein.
 * The original molecules are destroyed.
 */
{
   struct reaccs_molecule_t *mph, *mphh, *result;

   result = (struct reaccs_molecule_t *)NULL;
   while (!IsNULL(mp))
   {
      mph = mp->next; mp->next = (struct reaccs_molecule_t *)NULL;
      for (mphh = SplitMolecule(mp); !IsNULL(mphh); mphh = SplitMolecule(mp))
      {
         mphh->next = result; result = mphh;
      }
      mp->next = result; result = mp; mp = mph;
   }
   return(result);
}

struct reaccs_molecule_t *RemoveEmptyMolecules(struct reaccs_molecule_t *mp)
/*
 * Scans the molecule list pointed to by mp, removes the empty molecules
 * and returns the resulting list. The original list is destroyed.
 */
{
   struct reaccs_molecule_t *mph, *result;

   result = (struct reaccs_molecule_t *)NULL;
   while (!IsNULL(mp))
   {
      if (mp->n_atoms == 0)
      {
         mph = mp; mp = mp->next;
         FreeMolecule(mph);
      }
      else
      {
         mph = mp; mp = mp->next;
         mph->next = result; result = mph;
      }
   }

   return(result);
}

void SplitUnconnectedMolecules(struct reaccs_reaction_t *rp)
/*
 * Splits the molecules listed as reactants and products of *rp
 * into their connected components.
 */
{
   struct reaccs_molecule_t *mp;

   rp->reactants = RemoveEmptyMolecules(rp->reactants);
   rp->reactants = SplitMoleculeList(rp->reactants);
   rp->n_reactants = 0;
   for (mp=rp->reactants; !IsNULL(mp); mp=mp->next) rp->n_reactants++;

   rp->products  = RemoveEmptyMolecules(rp->products);
   rp->products  = SplitMoleculeList(rp->products);
   rp->n_products = 0;
   for (mp=rp->products; !IsNULL(mp); mp=mp->next) rp->n_products++;
}

struct reaccs_reaction_t *
ConnectedReactionComponents(struct reaccs_reaction_t *rp)
/*
 * Transforms the list of reactions rp into a list of all connected
 * subrections of those reactions. *rp is distroyed during that process.
 */
{
   int i,j,k,l,n,h;
   struct reaccs_reaction_t *result, *rhp;
   struct reaccs_molecule_t *mp, *mph;
   int new_label, old_label;
   int changed;
   int found, label;

   int *reactant_label;
   int *product_label;

   result = (struct reaccs_reaction_t *)NULL;

   while (!IsNULL(rp))
   {
      SplitUnconnectedMolecules(rp);
                                        /* allocate temporary storage */
      reactant_label = TypeAlloc(rp->n_reactants, int);
      product_label  = TypeAlloc(rp->n_products, int);

        /* at first all reactands and product form their own reaction */
      for (i=0; i<rp->n_reactants; i++)
         reactant_label[i] = i;
      for (i=0; i<rp->n_products;  i++)
         product_label[i]  = i+rp->n_reactants;

      do
      {
         changed = FALSE;
         for (i=0, mp=rp->reactants; !IsNULL(mp); i++, mp=mp->next)
            for (j=0, mph=rp->products; !IsNULL(mph); j++, mph=mph->next)
               if (reactant_label[i] != product_label[j])
               {
                  for (k=0; k<mp->n_atoms; k++) /* look for linking mapping */
                     for (l=0; l<mph->n_atoms; l++)
                        if (mp->atom_array[k].mapping ==
                            mph->atom_array[l].mapping &&
                            mp->atom_array[k].mapping != 0)
                        {               /* new link found -> merge labels */
                           new_label = reactant_label[i];
                           old_label = product_label[j];
                           if (new_label > old_label) {h=new_label; new_label=old_label; old_label=h;}
                           for (n=0; n<rp->n_reactants; n++)
                              if (reactant_label[n] == old_label)
                                 reactant_label[n] = new_label;
                           for (n=0; n<rp->n_products; n++)
                              if (product_label[n] == old_label)
                                 product_label[n] = new_label;
                           changed = TRUE;
                        }
               }
      } while (changed);

      for (label=0; label<rp->n_reactants+rp->n_products; label++)
      { /* check if a reactant or product is marked with "label" */
         found = FALSE;
         for (i=0; i<rp->n_reactants; i++)
            if (reactant_label[i] == label) found = TRUE;
         for (i=0; i<rp->n_products; i++)
            if (product_label[i] == label) found = TRUE;
         if (!found) continue;

         rhp = TypeAlloc(1,struct reaccs_reaction_t);
         if ((char *)rhp == (char *)NULL)
            return(result);
         *rhp = *rp;

         rhp->reactants = (struct reaccs_molecule_t *)NULL;
         for (i=0, n=0, mp=rp->reactants; i<rp->n_reactants; i++, mp=mp->next)
            if (reactant_label[i] == label)
            {
               n++;
               mph = CopyMolecule(mp);
               mph->next = rhp->reactants;
               rhp->reactants = mph;
            }
         rhp->n_reactants = n;

         rhp->products = (struct reaccs_molecule_t *)NULL;
         for (i=0, n=0, mp=rp->products; i<rp->n_products; i++, mp=mp->next)
            if (product_label[i] == label)
            {
               n++;
               mph = CopyMolecule(mp);
               mph->next = rhp->products;
               rhp->products = mph;
            }
         rhp->n_products = n;

         rhp->next = result; result = rhp;
      }

      rhp = rp->next; FreeReaction(rp); rp = rhp;

      if (reactant_label) MyFree((char *)reactant_label);
      if (product_label)  MyFree((char *)product_label);
   }

   return(result);
}

void UNUSED_AddClosureBonds(struct reaccs_molecule_t *mp,
                     int RT_atom[], /* 1..n */
                     int RT_bond[])
/*
 * Adds to the reaction type defined by the boolean arrays RT_atom[] and
 * RT_bond[] those bonds that connect atoms in RT_atom.
 */
{
   int i;
   struct reaccs_bond_t *bp;

   for (i=0,bp=mp->bond_array; i<mp->n_bonds; i++,bp++)
      if (RT_atom[bp->atoms[0]] && RT_atom[bp->atoms[1]])
         RT_bond[i] = TRUE;
}

struct reaccs_bond_t *UNUSED_SearchBondPtr(int at1, int at2,
                                    struct reaccs_molecule_t *mp)
/*
 * purpose:     searches the pointer to the bond in the molecule *mp
 *              that connects at1 and at2.
 * returns:     pointer to that bond or NULL if no such bond exists.
 */
{
   struct reaccs_bond_t *bp;
   int i;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if ((bp->atoms[0]==at1 && bp->atoms[1]==at2) ||
          (bp->atoms[0]==at2 && bp->atoms[1]==at1))
         return (bp);
   return ((struct reaccs_bond_t *)NULL);
}

void UNUSED_SearchNeighbours(int atnr,
                      struct reaccs_molecule_t *mp,
                      int neigh[MAXNEIGHBOURS],
                      int *nneighp)
/*
 * purpose:     searches the neighbours of atom atnr in *mp and stores
 *              them in neigh[0..*nneigh-1].
 */
{
   struct reaccs_bond_t *bp;
   int i;

   *nneighp = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->atoms[0]==atnr) neigh[(*nneighp)++] = bp->atoms[1];
      if (bp->atoms[1]==atnr) neigh[(*nneighp)++] = bp->atoms[0];
   }
}

// hydrogen that can be dropped
#define TRIVIAL   1
// stereo bond reassigned to allow drop of hydrogen
#define ANCILLARY 2
// hydrogen should not be dropped
#define ESSENTIAL 3

void PrepareCenterForHDrop(struct reaccs_molecule_t *mp,
                           neighbourhood_t nbp[],
                           int iatom,
                           int atom_status[])
/*
 * Prepares the hash/wedge bonds to hydrogen atoms and their neigbours
 * at atom index iatom for the removal of the bond to the hydrogen.
 *
 * nbp[iatom] is the neighbourhood of the central atom in question.
 */
{
   int nup = 0;
   int ndown = 0;
   int nhydrogens = 0;
   int i, j;
   int ipartner;
   int last_up, last_down, last_h;
   int stereo_symbol;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp, *bph;
   double xh, yh, x1, y1, x2, y2, xi, yi, norm, ref_cosa, cosa;
   int problem;
   int best_index;
   double best_cosa;

   /* return quickly for trivial cases */
   if (0 == strcmp("H", mp->atom_array[iatom].atom_symbol)) return;
   if (nbp[iatom].n_ligands < 3) return;

   /* analyze up/down configuration */
   last_up = last_down = -1;
   stereo_symbol = NONE; // stereo symbol to attached hydrogen
   bph = NULL;
   for (i=0; i<nbp[iatom].n_ligands; i++)
   {
      ap = &mp->atom_array[nbp[iatom].atoms[i]];
      bp = &mp->bond_array[nbp[iatom].bonds[i]];
      if (bp->stereo_symbol == UP  &&  bp->atoms[0] == iatom+1)
      {
         nup++;
         last_up = nbp[iatom].atoms[i];
      }
      if (bp->stereo_symbol == DOWN  &&  bp->atoms[0] == iatom+1)
      {
         ndown++;
         last_down = nbp[iatom].atoms[i];
      }
      if (0 == strcmp("H", ap->atom_symbol)  &&
          ap->mass_difference == NONE        &&
          ap->charge          == NONE        &&
          ap->radical         == NONE)
      {
         nhydrogens++;
         last_h = nbp[iatom].atoms[i];
         stereo_symbol = bp->stereo_symbol;
         bph=bp;
      }
   }

   if (nhydrogens == 0)         /* no hydrogen => NOP */
   {
      return;
   }
   else if (nhydrogens > 0  &&  nup+ndown == 0) /* no stereo bonds =>trivial */
   {
      for (i=0; i<nbp[iatom].n_ligands; i++)
      {
         ap = &mp->atom_array[nbp[iatom].atoms[i]];
         if (0 == strcmp("H", ap->atom_symbol))
         {
            if (ap->mass_difference == NONE &&
                ap->charge          == NONE &&
                ap->radical         == NONE)
               ap->color = TRIVIAL;
            else
            {
               ap->color = ESSENTIAL;
// fprintf(stderr, "1: hydrogen atom %d is essential\n", iatom+1);
            }
         }
      }
      return;
   }
   else if (nhydrogens >= 2  &&  nup+ndown != 0)        /* superflous stereo */
   {
      for (i=0; i<nbp[iatom].n_ligands; i++)
      {
         ap = &mp->atom_array[nbp[iatom].atoms[i]];
         bp = &mp->bond_array[nbp[iatom].bonds[i]];
         if (bp->atoms[0] == iatom+1) bp->stereo_symbol = NONE;
         ap->color = TRIVIAL;
      }
      return;
   }
   /* nhydrogens must be == 1 from here on */
   else if (stereo_symbol != UP  &&  stereo_symbol != DOWN)
   {                                  /* undesignated hydrogen => not so trivial */
// fprintf(stderr, "PrepareCenterForHDrop: clearing non-stereo (%d) bond %d-%d\n", stereo_symbol, last_h+1, iatom+1);
      /* Check if all remaining angles would be <<180 degrees */
      best_cosa = 2;
      for (i=0; i<nbp[iatom].n_ligands; i++)
      {
         if (nbp[iatom].atoms[i] == last_h) continue;
         ap = &mp->atom_array[nbp[iatom].atoms[i]];
         x1 = ap->x - mp->atom_array[iatom].x;
         y1 = ap->y - mp->atom_array[iatom].y;
         for (j=i+1; j<nbp[iatom].n_ligands; j++)
         {
            if (nbp[iatom].atoms[j] == last_h) continue;
            ap = &mp->atom_array[nbp[iatom].atoms[j]];
            xh = ap->x - mp->atom_array[iatom].x;
            yh = ap->y - mp->atom_array[iatom].y;
            cosa = xh*x1+yh*y1;
            norm = sqrt(xh*xh+yh*yh)*sqrt(x1*x1+y1*y1);
            if (norm > 0.0001) cosa /= norm;
            else              continue;
            if (cosa < best_cosa) best_cosa = cosa;
         }
      }
      if (nup+ndown == 0)
      {
         mp->atom_array[last_h].color = TRIVIAL;
      }
      else if (best_cosa > -0.8)
         mp->atom_array[last_h].color = TRIVIAL;
      else
      {
// fprintf(stderr, "2: hydrogen atom %d is ancillary\n", last_h+1);
         mp->atom_array[last_h].color = ANCILLARY;
      }
   }
   /* stereo_symbol is UP or DOWN from here on */
   else if (nup+ndown == 1)     /* only hydrogen has symbol */
   {    /* find closest neighbour bond and apply opposite wedge to it */

      xh = mp->atom_array[last_h].x - mp->atom_array[iatom].x;
      yh = mp->atom_array[last_h].y - mp->atom_array[iatom].y;
      best_index = 0;
      best_cosa = -2;
      for (i=0; i<nbp[iatom].n_ligands; i++)
      {
         if (nbp[iatom].atoms[i] == last_h) continue;
         ap = &mp->atom_array[nbp[iatom].atoms[i]];
         bp = &mp->bond_array[nbp[iatom].bonds[i]];
         if (bp->stereo_symbol != NONE) continue;  /* avoid conflicts */
         if (bp->bond_type != SINGLE) continue;    /* avoid conflicts */
         x1 = ap->x - mp->atom_array[iatom].x;
         y1 = ap->y - mp->atom_array[iatom].y;
         cosa = xh*x1+yh*y1;
         norm = sqrt(xh*xh+yh*yh)*sqrt(x1*x1+y1*y1);
         if (norm > 0.001) cosa /= norm;
         else              continue;
         /* prefer wedge to hetero atoms */
         if (0 != strcmp("C", ap->atom_symbol)) cosa *= 1.0001;
         /* prefer non_ring wedges */
         if (atom_status[nbp[iatom].atoms[i]] == 0)
            cosa *= 1.0002;
         if (cosa > best_cosa)
         {
            best_cosa = cosa;
            best_index = i;
         }
      }

      if (best_cosa < 0.05)
      {
// fprintf(stderr, "PrepareCenterForHDrop: Could not find normal bond close to hydrogen %d\n", last_h+1);
         mp->atom_array[last_h].color = ESSENTIAL;
// fprintf(stderr, "2: hydrogen atom %d is essential\n", last_h+1);
      }
      else
      {
         /* Check if all remaining angles would be <<180 degrees */
         best_cosa = 2;
         for (i=0; i<nbp[iatom].n_ligands; i++)
         {
            if (nbp[iatom].atoms[i] == last_h) continue;
            ap = &mp->atom_array[nbp[iatom].atoms[i]];
            x1 = ap->x - mp->atom_array[iatom].x;
            y1 = ap->y - mp->atom_array[iatom].y;
            for (j=i+1; j<nbp[iatom].n_ligands; j++)
            {
               if (nbp[iatom].atoms[j] == last_h) continue;
               ap = &mp->atom_array[nbp[iatom].atoms[j]];
               xh = ap->x - mp->atom_array[iatom].x;
               yh = ap->y - mp->atom_array[iatom].y;
               cosa = xh*x1+yh*y1;
               norm = sqrt(xh*xh+yh*yh)*sqrt(x1*x1+y1*y1);
               if (norm > 0.0001) cosa /= norm;
               else              continue;
               if (cosa < best_cosa) best_cosa = cosa;
            }
         }
         bp = &mp->bond_array[nbp[iatom].bonds[best_index]];
         if (stereo_symbol == UP) bp->stereo_symbol = DOWN;
         else                     bp->stereo_symbol = UP;
         if (bp->atoms[1] == iatom+1)   /* bond is wrong way =>swap it */
         {
            i = bp->atoms[0]; bp->atoms[0] = bp->atoms[1]; bp->atoms[1] = i;
         }
         if (bph) bph->stereo_symbol = NONE;
         if (best_cosa > -0.8)
         {
            mp->atom_array[last_h].color = TRIVIAL;
         }
         else
         {
// fprintf(stderr, "1: hydrogen atom %d is ancillary, best_cosa = %g\n", last_h+1, best_cosa);
            mp->atom_array[last_h].color = ANCILLARY;
         }
      }
   }
   else if (nup+ndown == 2  &&  nup*ndown==1)   /* two different wedges */
   {
      if (stereo_symbol == UP)   ipartner = last_down;
      else                       ipartner = last_up;
      xh = mp->atom_array[last_h].x - mp->atom_array[iatom].x;
      yh = mp->atom_array[last_h].y - mp->atom_array[iatom].y;
      x1 = mp->atom_array[ipartner].x - mp->atom_array[iatom].x;
      y1 = mp->atom_array[ipartner].y - mp->atom_array[iatom].y;
      /* Now, checking if there is plain bond pointing between those wedges */
      xi = (xh+x1)/2; yi = (yh+y1)/2;
      ref_cosa = xi*x1+yi*y1;
      norm = sqrt(xi*xi+yi*yi)*sqrt(x1*x1+y1*y1);
      if (norm > 0.001) ref_cosa /= norm;
      problem = FALSE;
      for (i=0; i<nbp[iatom].n_ligands; i++)
      {
         if (nbp[iatom].atoms[i] == last_h) continue;
         if (nbp[iatom].atoms[i] == ipartner) continue;
         ap = &mp->atom_array[nbp[iatom].atoms[i]];
         bp = &mp->bond_array[nbp[iatom].bonds[i]];
         if (bp->bond_type != SINGLE) continue;  /* avoid conflicts */
         x2 = ap->x - mp->atom_array[iatom].x;
         y2 = ap->y - mp->atom_array[iatom].y;
         cosa = xi*x2+yi*y2;
         norm = sqrt(xi*xi+yi*yi)*sqrt(x2*x2+y2*y2);
         if (norm > 0.001) cosa /= norm;
         if (cosa > ref_cosa) problem = TRUE;
      }
      if (problem)
      {
// fprintf(stderr, "PrepareCenterForHDrop: normal bond between wedges\n");
         mp->atom_array[last_h].color = ESSENTIAL;
// fprintf(stderr, "3: hydrogen atom %d is essential\n", last_h+1);
      }
      else
      {
// fprintf(stderr, "0: hydrogen atom %d is ancillary\n", last_h+1);
         mp->atom_array[last_h].color = ANCILLARY;
         if (bph) bph->stereo_symbol = NONE;
      }
   }
   else if (nup+ndown == 2)                    /* partner with same wedge */
   {
// fprintf(stderr, "PrepareCenterForHDrop: same wedge => cannot remove\n");
         mp->atom_array[last_h].color = ESSENTIAL;
// fprintf(stderr, "4: hydrogen atom %d is essential\n", last_h+1);
   }
   else                         /* more than one additional stereo bond */
   {
// fprintf(stderr, "PrepareCenterForHDrop: > 2 wedges => cannot remove\n");
         mp->atom_array[last_h].color = ESSENTIAL;
// fprintf(stderr, "5: hydrogen atom %d is essential\n", last_h+1);
   }
}

void MakeSomeHydrogensImplicit(struct reaccs_molecule_t *mp,
			       int which_hydrogens)
/*
 * purpose:     removes explicit hydrogens from molecule *mp;
 *              stereochemistry is preserved as far as possible.
 *		The bit flags in 'which_hydrogens' define which hydrogens
 *		will be candidates for removal.
 */
{
   int i, n;
   int hnr;
   int nlig;
   int partner, bond_index;
   struct symbol_list_t *slp;
   struct prop_line_t *hp;
   int *no_center;
   int *degree;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   neighbourhood_t *nbp;
   int *atom_status, *bond_status;
   int nrgp, ratom, rval;

   char *cp;
   int ii;
   int essential_hydrogens = 0;

// fprintf(stderr, "MakeSomeHydrogensImplicit: which_hydrogens = %d\n", which_hydrogens);
   for (i=0; i<mp->n_atoms; i++)
      mp->atom_array[i].color = NONE;
   nbp   = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);
   for (i=0; i<mp->n_atoms; i++)
      PrepareCenterForHDrop(mp, nbp, i, atom_status);

                                        /* find possible stereocenters */
   no_center = TypeAlloc(mp->n_atoms+1,int);
   for (i=1; i<=mp->n_atoms; i++) no_center[i] = FALSE;
   degree    = TypeAlloc(mp->n_atoms+1,int);
   for (i=1; i<=mp->n_atoms; i++) degree[i] = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->bond_type != SINGLE)
      {
         no_center[bp->atoms[0]] = TRUE;
         no_center[bp->atoms[1]] = TRUE;
      }
      degree[bp->atoms[0]]++; degree[bp->atoms[1]]++;
   }
   for (i=1; i<=mp->n_atoms; i++)
      if (degree[i] < 3) no_center[i] = TRUE;

   /*
    * loop over all atoms;
    *    increment hnr if not hydrogen or if unremovable hydrogen;
    *    remove atom and decrement mp->n_atoms if removable hydrogen.
    */
   for (hnr=1; hnr<=mp->n_atoms; /* void */ )
   {
      ap = mp->atom_array+hnr-1;
      if (0 != strcmp(ap->atom_symbol,"H"))	/* Not a hydrogen => keep it */
      {
         hnr++; continue;
      }

      /* Check if isotopically labelled hydrogens are to be removed. */
      if (!(which_hydrogens & ISOTOPE_HYDROGENS)  &&
          ap->mass_difference != 0)
      {
// fprintf(stderr, "1: hnr = %d\n", hnr);
         essential_hydrogens++;
         hnr++; continue;
      }

      nlig=0;                           /* search bonding partner */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (bp->atoms[0] == hnr)
         {
            partner = bp->atoms[1];
            bond_index = i;
            nlig++;
         }
         if (bp->atoms[1] == hnr)
         {
            partner = bp->atoms[0];
            bond_index = i;
            nlig++;
         }
      }

                        /* treat only mono-valent hydrogens */
      if (nlig >  1 || nlig == 0)
      {
         ShowMessageI("nlig == %d","MakeHydrogensImplicit",nlig);
// fprintf(stderr, "2: hnr = %d\n", hnr);
         essential_hydrogens++;
         hnr++; continue;
      }

      /* Don't remove explicit hydrogens from R atoms */
      if (mp->atom_array[partner-1].atom_symbol[0] == 'R')
      {
// fprintf(stderr, "0: essential bond (%d-%d) to hydrogen\n", bp->atoms[0], bp->atoms[1]);
// fprintf(stderr, "3: hnr = %d\n", hnr);
         essential_hydrogens++;
         hnr++; continue;
      }

      bp=mp->bond_array+bond_index;

      /* Check if any stereo hydrogens are to be removed. */
      if (!(which_hydrogens & (ANCILLARY_STEREO | ESSENTIAL_STEREO))  	&&
          bp->stereo_symbol != EITHER 					&&
          bp->stereo_symbol != NONE)
      {
// fprintf(stderr, "1: essential bond (%d-%d) to hydrogen\n", bp->atoms[0], bp->atoms[1]);
// fprintf(stderr, "4: hnr = %d\n", hnr);
         essential_hydrogens++;
         hnr++; continue;
      }

      /* Keep essential stereo atoms */
      if (!(which_hydrogens & ESSENTIAL_STEREO)  	&&
          ap->color == ESSENTIAL)
      {
// fprintf(stderr, "2: essential bond (%d-%d) to hydrogen\n", bp->atoms[0], bp->atoms[1]);
// fprintf(stderr, "5: hnr = %d\n", hnr);
         essential_hydrogens++;
         hnr++; continue;
      }

      /* Keep ancillary stereo atoms */
      if (!(which_hydrogens & ANCILLARY_STEREO)  	&&
          ap->color == ANCILLARY)
      {
// fprintf(stderr, "3: essential bond (%d-%d) to hydrogen, which_hydrogen = %d\n", bp->atoms[0], bp->atoms[1], which_hydrogens);
// fprintf(stderr, "6: hnr = %d\n", hnr);
         essential_hydrogens++;
         hnr++; continue;
      }

      if (bp->stereo_symbol != NONE)
      {
         ShowMessage("Ignoring stereo bond","MakeSomeHydrogensImplicit");
      }

      if (!(which_hydrogens & NO_QUERY_H_COUNT))
      {
	 if (mp->atom_array[partner-1].query_H_count == 0)
	    mp->atom_array[partner-1].query_H_count = 2;
	 else
	    mp->atom_array[partner-1].query_H_count++;
      }

      for (i=0, n=0; i<mp->n_atoms; i++)        /* remove H-atom */
         if (i != hnr-1) mp->atom_array[n++] = mp->atom_array[i];
      for (i=1, n=1; i<=mp->n_atoms; i++)       /* update no_center flags */
         if (i != hnr) no_center[n++] = no_center[i];
      mp->n_atoms--;
      for (i=0, n=0; i<mp->n_bonds; i++)        /* remove corresp. bond */
         if (i != bond_index) mp->bond_array[n++] = mp->bond_array[i];
      mp->n_bonds--;
                                                /* renumber bonds */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (bp->atoms[0] > hnr) bp->atoms[0]--;
         if (bp->atoms[1] > hnr) bp->atoms[1]--;
      }

      /*
       * Here, we need to loop over atom property lists to fix their
       * atom numbering.
       */
      for (slp = mp->symbol_lists; !IsNULL(slp); slp=slp->next)
      {
          if (slp->atom > hnr) slp->atom--;
      }
      /*
       * Here, we should fix the 'M  ALS' lines.
       */
      for (hp = mp->prop_lines; !IsNULL(hp); hp=hp->next)
      {
         /* Atom symbol list properties have been treated above. */
         /* Other properties should be handled as well though!   */
         if (0 == strncmp(hp->text, "M  RGP", strlen("M  RGP")))
         {
            sscanf(hp->text+strlen("M  RGP"), "%d", &nrgp);
            cp = hp->text+strlen("M  RGP")+3;
            for (ii=0; ii<nrgp; ii++)
            {
               sscanf(cp, "%d %d", &ratom, &rval);
               if (ratom > hnr)
               {
                   ratom--;
                   sprintf(cp, " %3d %3d", ratom, rval);
               }
               cp += 2*(1+3);
            }
         }
         else if (0 != strncmp(hp->text, "M  ALS", strlen("M  ALS")))
         {
            fprintf(stderr, "need to fix line '%s'\n", hp->text);
         }
      }
   }

   MyFree((char *)no_center);
   MyFree((char *)degree);
   if (nbp) MyFree((char *)nbp);
   MyFree((char *)atom_status);
   MyFree((char *)bond_status);
   for (i=0; i<mp->n_atoms; i++)
      mp->atom_array[i].color = NONE;
}

void MakeHydrogensImplicit(struct reaccs_molecule_t *mp)
/*
 * purpose:     removes all explicit hydrogens from molecule *mp;
 *              stereochemistry is preserved as far as possible
 */
{
   MakeSomeHydrogensImplicit(mp, ALL_HYDROGENS);
   // MakeSomeHydrogensImplicit(mp, (ALL_HYDROGENS & (~ANCILLARY_STEREO)));
   // fprintf(stderr, "MakeHydrogensImplicit: ALL_HYDROGENS = %d, ANCILLARY_STEREO = %d, which_hydrogens = %d\n", ALL_HYDROGENS, ANCILLARY_STEREO, (ALL_HYDROGENS & (~ANCILLARY_STEREO)));
}

void RemoveQueryOptions(struct reaccs_molecule_t *mp)
/*
 * Removes all query options found in the molecule *mp.
 * Currently, it only deals with HCounts.
 */
{
   int i;
   struct reaccs_atom_t *ap;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      ap->query_H_count = NONE;
}

struct reaccs_molecule_t *FuseMoleculeList(struct reaccs_molecule_t *mlp)
/*
 * Merges all molecules on list *mlp into one and returns the result.
 * Global structure components are taken from the first list element.
 */
{
   struct reaccs_molecule_t *resultp, *mph;
   int i, j, total_atoms;

   if (IsNULL(mlp)) return (mlp);

   resultp = TypeAlloc(1,struct reaccs_molecule_t);

   (*resultp) = (*mlp);     /* copy non-dynamic structure parts */
   if (resultp->stext_lines)
   {
      ShowMessage("Could not handle STEXT lines","FuseMoleculeList");
      ShowMessage("STEXT lines have been removed","FuseMoleculeList");
      resultp->n_atom_lists = 0;
      resultp->stext_lines = (struct stext_line_t *)NULL;
   }
   if (resultp->prop_lines)
   {
      ShowMessage("Could not handle property lines","FuseMoleculeList");
      ShowMessage("property lines have been removed","FuseMoleculeList");
      resultp->n_atom_lists = 0;
      resultp->prop_lines = (struct prop_line_t *)NULL;
   }

                                                /* new atoms and bonds */
   for (mph=mlp, resultp->n_atoms=0; !IsNULL(mph); mph=mph->next)
      resultp->n_atoms += mph->n_atoms;
   resultp->atom_array =
      TypeAlloc(resultp->n_atoms, struct reaccs_atom_t);
   for (mph=mlp, resultp->n_bonds=0; !IsNULL(mph); mph=mph->next)
      resultp->n_bonds += mph->n_bonds;
   resultp->bond_array =
      TypeAlloc(resultp->n_bonds, struct reaccs_bond_t);

   for (mph=mlp,j=0; !IsNULL(mph); mph=mph->next)       /* copy atoms */
      for (i=0; i<mph->n_atoms; i++, j++)
         resultp->atom_array[j] = mph->atom_array[i];

   for (mph=mlp,j=0,total_atoms=0; !IsNULL(mph); mph=mph->next) /* copy bonds */
   {
      for (i=0; i<mph->n_bonds; i++, j++)
      {
         resultp->bond_array[j] = mph->bond_array[i];
         resultp->bond_array[j].atoms[0] += total_atoms;
         resultp->bond_array[j].atoms[1] += total_atoms;
      }
      total_atoms += mph->n_atoms;
   }

   if (resultp->symbol_lists)
   {
      ShowMessage("Could not handle atom lists","FuseMoleculeList");
      ShowMessage("atom list has been removed","FuseMoleculeList");
      resultp->n_atom_lists = 0;
      resultp->symbol_lists = (struct symbol_list_t *)NULL;
   }

   resultp->next = (struct reaccs_molecule_t *)NULL;

   return (resultp);
}

struct reaccs_molecule_t *LinkMolecules(struct reaccs_molecule_t *mp1,
                                        unsigned int              at1,
                                        int                       bd,
                                        unsigned int              at2,
                                        struct reaccs_molecule_t *mp2)
/*
 * Creates a new molecule which combines *mp1 and *mp2. The two
 * fragments are linked at the former atoms at1 and at2 with a
 * bond of type bd. The return value is a pointer to the constructed
 * molecule. Global structure components are taken from mp1.
 */
{
   struct reaccs_molecule_t *resultp;
   int i;
   struct reaccs_bond_t *bp;

   resultp = TypeAlloc(1,struct reaccs_molecule_t);

   *resultp = *mp1;     /* copy non-dynamic structure parts */
   resultp->stext_lines = (struct stext_line_t *)NULL;
   resultp->prop_lines = (struct prop_line_t *)NULL;
   resultp->symbol_lists = (struct symbol_list_t *)NULL;
   resultp->next = (struct reaccs_molecule_t *)NULL;

   resultp->n_atoms = mp1->n_atoms+mp2->n_atoms; /* allocate space for  */
   resultp->atom_array =                         /* new atoms and bonds */
      TypeAlloc(resultp->n_atoms, struct reaccs_atom_t);
   resultp->n_bonds = mp1->n_bonds+mp2->n_bonds+1;
   resultp->bond_array =
      TypeAlloc(resultp->n_bonds, struct reaccs_bond_t);

   for (i=0; i<mp1->n_atoms; i++)       /* copy atoms */
      resultp->atom_array[i] = mp1->atom_array[i];
   for (i=0; i<mp2->n_atoms; i++)
      resultp->atom_array[i+mp1->n_atoms] = mp2->atom_array[i];

   for (i=0; i<mp1->n_bonds; i++)       /* copy bonds */
      resultp->bond_array[i] = mp1->bond_array[i];
   for (i=0; i<mp2->n_bonds; i++)
   {
      bp = &resultp->bond_array[i+mp1->n_bonds];
      *bp = mp2->bond_array[i];
      bp->atoms[0] += mp1->n_atoms;
      bp->atoms[1] += mp1->n_atoms;
   }

   bp = &resultp->bond_array[resultp->n_bonds-1]; /* make new bond */
   bp->atoms[0] = at1; bp->atoms[1] = at2+mp1->n_atoms;
   bp->bond_type = bd;
   bp->stereo_symbol = NONE; bp->dummy = 0;
   bp->topography = NONE; bp->reaction_mark = MAKE_BREAK;

   return (resultp);
}

// Standard PTABLE data structure as used by MDL/Accelrys/Symyx ...
// see utilities.h for definition of entry struct
struct ptable_entry ptable[] =
   {
      { "C"  ,  4 ,   12.01115 ,  2.55 , 11.26 },  /* copied for speed */
      { "N"  ,  3 ,   14.00670 ,  3.04 , 14.53 },  /* copied for speed */
      { "O"  ,  2 ,   15.99940 ,  3.44 , 13.62 },  /* copied for speed */
      { "H"  ,  1 ,    1.00797 ,  2.20 , 13.60 },
      { "He" ,  0 ,    4.00260 ,  0.00 , 24.59 },
      { "Li" ,  1 ,    6.93900 ,  0.98 ,  5.39 },
      { "Be" ,  2 ,    9.01220 ,  1.57 ,  9.32 },
      { "B"  ,  3 ,   10.81100 ,  2.04 ,  8.30 },
      { "C"  ,  4 ,   12.01115 ,  2.55 , 11.26 },
      { "N"  ,  3 ,   14.00670 ,  3.04 , 14.53 },
      { "O"  ,  2 ,   15.99940 ,  3.44 , 13.62 },
      { "F"  ,  1 ,   18.99840 ,  3.98 , 17.42 },
      { "Ne" ,  0 ,   20.18300 ,  0.00 , 21.56 },
      { "Na" ,  1 ,   22.98980 ,  0.93 ,  5.14 },
      { "Mg" ,  2 ,   24.31200 ,  1.31 ,  7.65 },
      { "Al" ,  3 ,   26.98150 ,  1.61 ,  5.99 },
      { "Si" ,  4 ,   28.08600 ,  1.90 ,  8.15 },
      { "P"  ,  5 ,   30.97380 ,  2.19 , 10.49 },
      { "S"  ,  2 ,   32.06400 ,  2.58 , 10.36 },
      { "Cl" ,  1 ,   35.45300 ,  3.16 , 12.97 },
      { "Ar" ,  0 ,   39.94800 ,  0.00 , 15.76 },
      { "K"  ,  1 ,   39.10200 ,  0.82 ,  4.34 },
      { "Ca" ,  2 ,   40.08000 ,  1.00 ,  6.11 },
      { "Sc" ,  3 ,   44.95600 ,  1.36 ,  6.54 },
      { "Ti" ,  3 ,   47.90000 ,  1.54 ,  6.82 },
      { "V"  ,  3 ,   50.94200 ,  1.63 ,  6.74 },
      { "Cr" ,  3 ,   51.99600 ,  1.66 ,  6.77 },
      { "Mn" ,  4 ,   54.93800 ,  1.55 ,  7.44 },
      { "Fe" ,  2 ,   55.84700 ,  1.83 ,  7.87 },
      { "Co" ,  2 ,   58.93320 ,  1.88 ,  7.86 },
      { "Ni" ,  2 ,   58.71000 ,  1.91 ,  7.64 },
      { "Cu" ,  1 ,   63.54600 ,  1.90 ,  7.73 },
      { "Zn" ,  2 ,   65.37000 ,  1.65 ,  9.39 },
      { "Ga" ,  2 ,   69.72000 ,  1.81 ,  6.00 },
      { "Ge" ,  4 ,   72.59000 ,  2.01 ,  7.90 },
      { "As" ,  3 ,   74.92160 ,  2.18 ,  9.81 },
      { "Se" ,  4 ,   78.96000 ,  2.55 ,  9.75 },
      { "Br" ,  1 ,   79.90400 ,  2.96 , 11.81 },
      { "Kr" ,  0 ,   83.80000 ,  0.00 , 14.00 },
      { "Rb" ,  1 ,   85.47000 ,  0.82 ,  4.18 },
      { "Sr" ,  2 ,   87.62000 ,  0.95 ,  5.70 },
      { "Y"  ,  3 ,   88.90500 ,  1.22 ,  6.38 },
      { "Zr" ,  4 ,   91.22000 ,  1.33 ,  6.84 },
      { "Nb" ,  3 ,   92.90600 ,  1.60 ,  6.88 },
      { "Mo" ,  4 ,   95.94000 ,  2.16 ,  7.10 },
      { "Tc" ,  6 ,   98.90620 ,  1.90 ,  7.28 },
      { "Ru" ,  4 ,  101.07000 ,  2.20 ,  7.37 },
      { "Rh" ,  3 ,  102.90500 ,  2.28 ,  7.46 },
      { "Pd" ,  4 ,  106.40000 ,  2.20 ,  8.34 },
      { "Ag" ,  1 ,  107.86800 ,  1.93 ,  7.58 },
      { "Cd" ,  2 ,  112.40000 ,  1.69 ,  8.99 },
      { "In" ,  3 ,  114.82000 ,  1.78 ,  5.79 },
      { "Sn" ,  2 ,  118.69000 ,  1.96 ,  7.34 },
      { "Sb" ,  3 ,  121.75000 ,  2.05 ,  8.64 },
      { "Te" ,  4 ,  127.60000 ,  2.10 ,  9.01 },
      { "I"  ,  1 ,  126.90440 ,  2.66 , 10.45 },
      { "Xe" ,  0 ,  131.30000 ,  0.00 , 12.13 },
      { "Cs" ,  1 ,  132.90500 ,  0.79 ,  3.89 },
      { "Ba" ,  2 ,  137.33000 ,  0.89 ,  5.21 },
      { "La" ,  3 ,  138.91000 ,  1.10 ,  5.58 },
      { "Ce" ,  3 ,  140.12000 ,  1.12 ,  5.47 },
      { "Pr" ,  3 ,  140.90700 ,  1.13 ,  5.42 },
      { "Nd" ,  3 ,  144.24000 ,  1.14 ,  5.49 },
      { "Pm" ,  3 ,  145.00000 ,  1.20 ,  5.55 },
      { "Sm" ,  2 ,  150.35000 ,  1.17 ,  5.63 },
      { "Eu" ,  2 ,  151.96000 ,  1.20 ,  5.67 },
      { "Gd" ,  3 ,  157.25000 ,  1.20 ,  6.14 },
      { "Tb" ,  3 ,  158.92400 ,  1.20 ,  5.85 },
      { "Dy" ,  3 ,  162.50000 ,  1.22 ,  5.93 },
      { "Ho" ,  3 ,  164.93000 ,  1.23 ,  6.02 },
      { "Er" ,  3 ,  167.26000 ,  1.24 ,  6.10 },
      { "Tm" ,  3 ,  168.93400 ,  1.25 ,  6.18 },
      { "Yb" ,  2 ,  173.04000 ,  1.10 ,  6.25 },
      { "Lu" ,  3 ,  174.97000 ,  1.27 ,  5.43 },
      { "Hf" ,  4 ,  178.49000 ,  1.30 ,  7.00 },
      { "Ta" ,  5 ,  180.94800 ,  1.50 ,  7.89 },
      { "W"  ,  6 ,  183.85000 ,  2.36 ,  7.98 },
      { "Re" ,  0 ,  186.20000 ,  1.90 ,  7.88 },
      { "Os" ,  3 ,  190.20000 ,  2.20 ,  8.70 },
      { "Ir" ,  3 ,  192.20000 ,  2.20 ,  9.10 },
      { "Pt" ,  2 ,  195.09000 ,  2.28 ,  9.00 },
      { "Au" ,  1 ,  196.96700 ,  2.54 ,  9.23 },
      { "Hg" ,  2 ,  200.59000 ,  2.00 , 10.44 },
      { "Tl" ,  1 ,  204.37000 ,  2.04 ,  6.11 },
      { "Pb" ,  2 ,  207.19000 ,  2.33 ,  7.42 },
      { "Bi" ,  3 ,  208.98000 ,  2.02 ,  7.29 },
      { "Po" ,  0 ,  209.00000 ,  2.00 ,  8.42 },
      { "At" ,  3 ,  210.00000 ,  2.20 ,  0.00 },
      { "Rn" ,  0 ,  222.00000 ,  0.00 , 10.75 },
      { "Fr" ,  1 ,  223.00000 ,  0.70 ,  0.00 },
      { "Ra" ,  2 ,  226.03000 ,  0.90 ,  5.28 },
      { "Ac" ,  0 ,  227.00000 ,  1.10 ,  6.90 },
      { "Th" ,  4 ,  232.03800 ,  1.30 ,  0.00 },
      { "Pa" ,  0 ,  231.04000 ,  1.50 ,  0.00 },
      { "U"  ,  6 ,  238.03000 ,  1.38 ,  0.00 },
      { "Np" ,  5 ,  237.05000 ,  1.36 ,  0.00 },
      { "Pu" ,  4 ,  244.00000 ,  1.28 ,  5.80 },
      { "Am" ,  4 ,  243.00000 ,  1.30 ,  6.00 },
      { "Cm" ,  3 ,  247.00000 ,  1.30 ,  0.00 },
      { "Bk" ,  3 ,  247.00000 ,  1.30 ,  0.00 },
      { "Cf" ,  0 ,  251.00000 ,  1.30 ,  0.00 },
      { "Es" ,  0 ,  254.00000 ,  1.30 ,  0.00 },
      { "Fm" ,  0 ,  257.00000 ,  1.30 ,  0.00 },
      { "Md" ,  0 ,  258.00000 ,  0.00 ,  0.00 },
      { "No" ,  0 ,  259.00000 ,  0.00 ,  0.00 },
      { "Lr" ,  0 ,  260.00000 ,  0.00 ,  0.00 },
      { "D"  ,  1 ,    2.01400 ,  2.20 , 13.60 },
      { "R"  ,  1 ,    0.00000 ,  0.00 , 00.00 },
      { (char *)NULL ,  0 ,  0.00000 ,  0.00 ,  0.00 }  // sentinel
   };

/**
 * Converts the mass value for an atomic symbol into a MOL file mass difference
 * value.
 */
int LookupMassDifference(double mass, char *symbol)
{
    int i;
    for (i=0; i<sizeof(ptable); i++)
        if (0 == strcmp(ptable[i].symbol, symbol))
        {
            double diff = mass-ptable[i].mass;
            if (diff > 0) return (int)(diff+0.49);
            else          return (int)(diff-0.49);
        }
    // just return rounded value if not found in ptable
    if (mass > 0) return (int)(mass+0.49);
    else          return (int)(mass-0.49);
}

string_int_table periodic_table[] =
   {
      {"H",             1},
      {"D",             1},
      {"T",             1},
      {"He",            2},
      {"Li",            3},
      {"Be",            4},
      {"B",             5},
      {"C",             6},
      {"N",             7},
      {"O",             8},
      {"F",             9},
      {"Ne",           10},
      {"Na",           11},
      {"Mg",           12},
      {"Al",           13},
      {"Si",           14},
      {"P",            15},
      {"S",            16},
      {"Cl",           17},
      {"Ar",           18},
      {"K",            19},
      {"Ca",           20},
      {"Sc",           21},
      {"Ti",           22},
      {"V",            23},
      {"Cr",           24},
      {"Mn",           25},
      {"Fe",           26},
      {"Co",           27},
      {"Ni",           28},
      {"Cu",           29},
      {"Zn",           30},
      {"Ga",           31},
      {"Ge",           32},
      {"As",           33},
      {"Se",           34},
      {"Br",           35},
      {"Kr",           36},
      {"Rb",           37},
      {"Sr",           38},
      {"Y",            39},
      {"Zr",           40},
      {"Nb",           41},
      {"Mo",           42},
      {"Tc",           43},
      {"Ru",           44},
      {"Rh",           45},
      {"Pd",           46},
      {"Ag",           47},
      {"Cd",           48},
      {"In",           49},
      {"Sn",           50},
      {"Sb",           51},
      {"Te",           52},
      {"I",            53},
      {"Xe",           54},
      {"Cs",           55},
      {"Ba",           56},
      {"La",           57},
      {"Ce",           58},
      {"Pr",           59},
      {"Nd",           60},
      {"Pm",           61},
      {"Sm",           62},
      {"Eu",           63},
      {"Gd",           64},
      {"Tb",           65},
      {"Dy",           66},
      {"Ho",           67},
      {"Er",           68},
      {"Tm",           69},
      {"Yb",           70},
      {"Lu",           71},
      {"Hf",           72},
      {"Ta",           73},
      {"W",            74},
      {"Re",           75},
      {"Os",           76},
      {"Ir",           77},
      {"Pt",           78},
      {"Au",           79},
      {"Hg",           80},
      {"Tl",           81},
      {"Pb",           82},
      {"Bi",           83},
      {"Po",           84},
      {"At",           85},
      {"Rn",           86},
      {"Fr",           87},
      {"Ra",           88},
      {"Ac",           89},
      {"Th",           90},
      {"Pa",           91},
      {"U",            92},
      {"Np",           93},
      {"Pu",           94},
      {"X",           120},
      {"Q",           121},
      {"M",           122},
      {"R",           123},
      {"A",           124},
      {(char *)NULL,    0}
   };

/*
 * predefined generic atom type sets for use in STRUCHK
 */
char *HC_table[] =                      /* pseudosymbol "G" */
   {"H", "C", (char *)NULL};

char *non_metal_hetero_elements[] =     /* pseudosymbol "Q" */
   {
      "He",
      "B", "N", "O", "F", "Ne",
      "Si", "P", "S", "Cl", "Ar",
      "As", "Se", "Br", "Kr",
      "Sb", "Te", "I", "Xe",
      "At",                     /* "Rn", This element must be removed */
      (char *)NULL,             /* because of a trick in utils.c */
   };

char *metals[] =                /* pseudosymbol "M" */
   {
      "Li", "Be",
      "Na", "Mg", "Al",
      "K", "Ca", "Sc",
        "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
           "Ga",
      "Rb", "Sr", "Y",
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
           "In", "Sn",
      "Cs", "Ba", "La",
       "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
       "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
        "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
           "Tl", "Pb", "Bi", "Po",
      "Fr", "Ra", "Ac",
       "Th", "Pa", "U", "Np", "Pu",
      (char *)NULL,
  };

char *non_metal_small_solution[] =      /* pseudosymbol "Qs" */
   {
      "H",
      "B",  "C", "N", "O", "F",
           "Si", "P", "S", "Cl",
                     "Se", "Br",
                            "I",
      (char *)NULL,
   };

char *alkali_metals[] =                /* pseudosymbol "alk" */
   {
      "Li", "Na", "K", "Rb", "Cs", "Fr",
      (char *)NULL,
  };

char *gr2[] =                              /* pseudosymbol "gr2" */
   {
      "Be", "Mg", "Ca", "Sr", "Ba", "Ra",
      (char *)NULL,
  };

char *gr3[] =                               /* pseudosymbol "gr3" */
   {
      "B", "Al", "Ga", "In", "Tl",
      (char *)NULL,
  };

char *gr4[] =                              /* pseudosymbol "gr4" */
   {
      "C", "Si", "Ge", "Sn", "Pb",
      (char *)NULL,
  };

char *ONS_table[] =                     /* pseudosymbol "ONS" or "ons" */
   {"O", "N", "S", (char *)NULL};

char *on2[] =                            /* pseudosymbol "on2" */
   {
      "O", "N", "S", "P", "Se", "Te", "Po",
      (char *)NULL,
  };

char *halogenes[] =                     /* pseudosymbol "X" or "hal" */
   {"F", "Cl", "Br", "I", "At", (char *)NULL};

char *ha2[] =                     /* pseudosymbol "ha2" */
   {"Cl", "Br", "I", "At", (char *)NULL};

char *transition_metals[] =                /* pseudosymbol "trn" */
   {
      "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
      "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
      "La", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
      (char *)NULL,
  };

char *tra[] =                /* pseudosymbol "tra" */
   {
      "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni",
      "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
      "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt",
      (char *)NULL,
  };

char *trb[] =                /* pseudosymbol "trb" */
   {
      "Cu", "Zn", "Ag", "Cd", "Au", "Hg",
      (char *)NULL,
  };

char *tm1[] =                /* pseudosymbol "tm1" */
   {
      "Cu", "Ag", "Au",
      (char *)NULL,
  };

char *tm2[] =                /* pseudosymbol "tm2" */
   {
      "Zn", "Cd", "Hg",
      (char *)NULL,
  };

char *tm3[] =                /* pseudosymbol "tm3" */
   {
      "Sc", "Y", "La",
      (char *)NULL,
  };

char *tm4[] =                /* pseudosymbol "tm4" */
   {
      "Ti", "Zr", "Hf",
      (char *)NULL,
  };

char *tm5[] =                /* pseudosymbol "tm5" */
   {
      "V", "Nb", "Ta",
      (char *)NULL,
  };

char *tm6[] =                /* pseudosymbol "tm6" */
   {
      "Cr", "Mo", "W",
      (char *)NULL,
  };

char *tm7[] =                /* pseudosymbol "tm7" */
   {
      "Mn", "Tc", "Re",
      (char *)NULL,
  };

char *tm8[] =                /* pseudosymbol "tm8" */
   {
      "Fe", "Co", "Ni",
      "Ru", "Rh", "Pd",
      "Os", "Ir", "Pt",
      (char *)NULL,
  };

char *lanthanoids[] =                /* pseudosymbol "lan" */
   {
       "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
       "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
      (char *)NULL,
  };

char *amino_acids[] =                   /* pseudosymbol "Ami" or "ami"*/
   {
      "Ala", "Arg", "Asn", "Asp", "Cys",
      "Gln", "Glu", "Gly", "His", "Ile",
      "Leu", "Lys", "Met", "Phe", "Pro",
      "Ser", "Thr", "Trp", "Tyr", "Val",
      (char *)NULL,
   };

static int IsInStringTable(char *symbol, char *table[])
/*
 * Checks if the string symbol is listed in table[] and returns
 * TRUE if it is and FALSE otherwise.
 */
{
   char **stringp;

   for (stringp = table;
        !IsNULL(*stringp);
        stringp++)
      if (0 == strcmp(*stringp,symbol)) return (TRUE);
   return (FALSE);
}

int AtomSymbolMatch(char *atsym, char *pattern)
/*
 * Returns TRUE if atsym is in the comma delimited list of atom symbols
 * stored in pattern and FALSE otherwise.
 * There are also a number of standard atom type lists like "alk" for alkali metals or
 * "Q" for non-C/non-H defined above as arrays of strings.
 */
{
   /* static */ char pat_buf[400];
   char *tokp;

   strcpy(pat_buf,pattern);
   for (tokp = strtok(pat_buf,",");
        !IsNULL(tokp);
        tokp = strtok((char *)NULL,","))
   {
      if (islower(*tokp))
      {
         if (0 == strcmp("alk",tokp))
         {
            if (IsInStringTable(atsym,alkali_metals)) return (TRUE);
         }
         else if (0 == strcmp("gr2",tokp))
         {
            if (IsInStringTable(atsym,gr2)) return (TRUE);
         }
         else if (0 == strcmp("gr3",tokp))
         {
            if (IsInStringTable(atsym,gr3)) return (TRUE);
         }
         else if (0 == strcmp("gr4",tokp))
         {
            if (IsInStringTable(atsym,gr4)) return (TRUE);
         }
         else if (0 == strcmp("ons",tokp))
         {
            if (IsInStringTable(atsym,ONS_table)) return (TRUE);
         }
         else if (0 == strcmp("on2",tokp))
         {
            if (IsInStringTable(atsym,on2)) return (TRUE);
         }
         else if (0 == strcmp("hal",tokp))
         {
            if (IsInStringTable(atsym,halogenes)) return (TRUE);
         }
         else if (0 == strcmp("ha2",tokp))
         {
            if (IsInStringTable(atsym,ha2)) return (TRUE);
         }
         else if (0 == strcmp("trn",tokp))
         {
            if (IsInStringTable(atsym,transition_metals)) return (TRUE);
         }
         else if (0 == strcmp("tra",tokp))
         {
            if (IsInStringTable(atsym,tra)) return (TRUE);
         }
         else if (0 == strcmp("trb",tokp))
         {
            if (IsInStringTable(atsym,trb)) return (TRUE);
         }
         else if (0 == strcmp("tm1",tokp))
         {
            if (IsInStringTable(atsym,tm1)) return (TRUE);
         }
         else if (0 == strcmp("tm2",tokp))
         {
            if (IsInStringTable(atsym,tm2)) return (TRUE);
         }
         else if (0 == strcmp("tm3",tokp))
         {
            if (IsInStringTable(atsym,tm3)) return (TRUE);
         }
         else if (0 == strcmp("tm4",tokp))
         {
            if (IsInStringTable(atsym,tm4)) return (TRUE);
         }
         else if (0 == strcmp("tm5",tokp))
         {
            if (IsInStringTable(atsym,tm5)) return (TRUE);
         }
         else if (0 == strcmp("tm6",tokp))
         {
            if (IsInStringTable(atsym,tm6)) return (TRUE);
         }
         else if (0 == strcmp("tm7",tokp))
         {
            if (IsInStringTable(atsym,tm7)) return (TRUE);
         }
         else if (0 == strcmp("tm8",tokp))
         {
            if (IsInStringTable(atsym,tm8)) return (TRUE);
         }
         else if (0 == strcmp("lan",tokp))
         {
            if (IsInStringTable(atsym,lanthanoids)) return (TRUE);
         }
         else if (0 == strcmp("ami",tokp))
         {
            if (IsInStringTable(atsym,amino_acids)) return (TRUE);
         }
      }
      if (0 == strcmp(atsym,tokp)) return (TRUE);
   }

   if (0 == strcmp("A",pattern))
      return (0 != strcmp("H",atsym));
   else if (0 == strcmp("Qs",pattern))
      return (IsInStringTable(atsym,non_metal_small_solution));
   else if (0 == strcmp("G",pattern))
      return (IsInStringTable(atsym,HC_table));
   else if (0 == strcmp("ONS",pattern))
      return (IsInStringTable(atsym,ONS_table));
   else if (0 == strcmp("X",pattern))
      return (IsInStringTable(atsym,halogenes));
   else if (0 == strcmp("Q",pattern))
      return (IsInStringTable(atsym,non_metal_hetero_elements));
   else if (0 == strcmp("M",pattern))
      return (IsInStringTable(atsym,metals));
   else if (0 == strcmp("Ami",pattern))
      return (IsInStringTable(atsym,amino_acids));
   return (FALSE);
}

int MapToQueryAtomType(char *symbol,
                       char *source_table[],
                       char *target_symbol)
/*
 * If the string symbol is listed in source_table[], then it is
 * replaced by the string target_symbol.
 * This function is used to map sets of atom types to a query atom
 * type symbol.
 */
{
   int i;

   for (i=0; !IsNULL(source_table[i]); i++)
      if (0 == strcmp(symbol,source_table[i]))
      {
         strcpy(symbol,target_symbol);
         return (TRUE);
      }
   return (FALSE);
}

void ResetColors(struct reaccs_molecule_t *mp)
/*
 *  Purpose:    Sets the color fields of all bond and atom structures in
 *              the molecule *mp to 0.
 */
{
   int i;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      ap->color = NO_COLOR;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      bp->color = NO_COLOR;
}

void ResetValues(struct reaccs_molecule_t *mp)
/*
 *  Purpose:    Sets the value fields of all bond and atom structures in
 *              the molecule *mp to 0.0.
 */
{
   int i;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      ap->value = 0.0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      bp->value = 0.0;
}

void ApplyToAllMolecules(struct reaccs_reaction_t *rp,
                         void (*func)(struct reaccs_molecule_t *))
/*
 *  Purpose:    Executes func(mp) for all molecules, i.e. reactants and
 *              products, in *rp.
 */
{
   struct reaccs_molecule_t *mp;

   for (mp=rp->reactants; !IsNULL(mp); mp=mp->next)
      (*func)(mp);
   for (mp=rp->products; !IsNULL(mp); mp=mp->next)
      (*func)(mp);
}

void ApplyToAllProducts(struct reaccs_reaction_t *rp,
                        void (*func)(struct reaccs_molecule_t *))
/*
 *  Purpose:    Executes func(mp) for all products in *rp.
 */
{
   struct reaccs_molecule_t *mp;

   for (mp=rp->products; !IsNULL(mp); mp=mp->next)
      (*func)(mp);
}

int CheckNeighbourhood(struct reaccs_molecule_t *mp)
/*
 * Checks if the number of neighbours of the atoms in *mp is within the allowed limits.
 */
{
   int i;
   int at0, at1;
   int result;
   int* nbcounts;

   result = TRUE;
   if (IsNULL(mp)) return FALSE;


   nbcounts = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_bonds; i++)
   {
      at0 = mp->bond_array[i].atoms[0]-1;
      at1 = mp->bond_array[i].atoms[1]-1;
      if (at0 <0  ||  at0 >= mp->n_atoms  ||
          at1 < 0  ||  at1 >= mp->n_atoms)
      {
         ShowMessageI("illegal atom numbers in bond %d ", "CheckNeighbourhood", i);
         // result = FALSE;
         // break;
      }
      nbcounts[at0]++;
      nbcounts[at1]++;
   }
   for (i=0; i<mp->n_atoms; i++)
      if (nbcounts[i] > MAXNEIGHBOURS)
      {
          ShowMessageI("atom %d has too man ligands", "CheckNeighbourhood", i+1);
          result = FALSE;
      }
//  FORTIFY CheckFortifyMemory();
   MyFree((char *)nbcounts);
   return result;
}

int SetupNeighbourhood(struct reaccs_molecule_t *mp,
                       neighbourhood_t *neighbour_array,
                       int nlimit)
/*
 * Computes the array of atom and bond neighbourhoods of the
 * atoms in *mp which have an atom number less than nlimit.
 *
 * Returns FALSE in case of failure.
 */
{
   int i;
   int at0, at1;

   for (i=0; i<mp->n_atoms; i++)        /* setup neighbourhood */
      neighbour_array[i].n_ligands = 0;
   for (i=0; i<mp->n_bonds; i++)
   {
      at0 = mp->bond_array[i].atoms[0]-1;
      at1 = mp->bond_array[i].atoms[1]-1;
      if (at0 >= nlimit  ||  at1 >= nlimit) continue;
      neighbour_array[at0].atoms[neighbour_array[at0].n_ligands] = at1;
      neighbour_array[at0].bonds[neighbour_array[at0].n_ligands] = i;
      neighbour_array[at0].n_ligands++;
      if (neighbour_array[at0].n_ligands > MAXNEIGHBOURS)
      {
         fprintf(stderr,"Too many ligands at atom %d\n", at0);
         ShowMessageI("Too many neighbours at atom %d\n",
                      "SetupNeighbourhood",
                      at0+1);
         return FALSE;
      }
      neighbour_array[at1].atoms[neighbour_array[at1].n_ligands] = at0;
      neighbour_array[at1].bonds[neighbour_array[at1].n_ligands] = i;
      neighbour_array[at1].n_ligands++;
      if (neighbour_array[at1].n_ligands > MAXNEIGHBOURS)
      {
         fprintf(stderr,"Too many ligands at atom %d\n", at1);
         ShowMessageI("Too many neighbours at atom %d\n",
                      "SetupNeighbourhood",
                      at1+1);
         return FALSE;
      }
   }
   return TRUE;
}

// bond type strings used in STRUCHK rules
symbol_entry_t bond_to_string[] =
   {
      {NONE,            "?"},
      {SINGLE,          "-"},
      {DOUBLE,          "="},
      {TRIPLE,          "#"},
      {AROMATIC,        "~"},
      {SINGLE_DOUBLE,   "-="},
      {SINGLE_AROMATIC, "-~"},
      {DOUBLE_AROMATIC, "=~"},
      {ANY_BOND,        "*"},
      {(-1), (char *)NULL}
   };

// charge strings used in STRUCHK rules
symbol_entry_t charge_to_string[] =
   {
      {NONE,   ""},
      {7,      "+7"},
      {6,      "+6"},
      {5,      "+5"},
      {4,      "+4"},
      {3,      "+3"},
      {2,      "+2"},
      {1,      "+1"},
      {-1,     "-1"},
      {-2,     "-2"},
      {-3,     "-3"},
      {-4,     "-4"},
      {-5,     "-5"},
      {-6,     "-6"},
      {-7,     "-7"},
      {ANY_CHARGE, "+-"},
      {0, (char *)NULL}
   };

// radical definition strings used in STRUCHK rules
symbol_entry_t radical_to_string[] =
   {
      {NONE,       ""},
      {SINGLET,    "|"},
      {DOUBLET,    "."},
      {TRIPLET,    ":"},
      {0, (char *)NULL}
   };

void SortNeighbourhood(neighbourhood_t *np, struct reaccs_molecule_t *mp)
/*
 * Sorts the neighbourhood *np with respect to the bond orders and
 * atom symbols defined in *mp.
 */
{
   int i,j,h;
   int btj, btj1;
   int hcj, hcj1;
   char *cpj, *cpj1;

   for (i=1; i<np->n_ligands; i++)
      for (j=i-1; j>=0; j--)
      {
         btj  = mp->bond_array[np->bonds[j]].bond_type;
         btj1 = mp->bond_array[np->bonds[j+1]].bond_type;
         cpj  = mp->atom_array[np->atoms[j]].atom_symbol;
         cpj1 = mp->atom_array[np->atoms[j+1]].atom_symbol;
         hcj  = mp->atom_array[np->atoms[j]].query_H_count;
         hcj1 = mp->atom_array[np->atoms[j+1]].query_H_count;
         if (btj < btj1 ||
             (btj == btj1  &&  0 <  strcmp(cpj,cpj1))  ||
             (btj == btj1  &&  0 == strcmp(cpj,cpj1)  &&  hcj < hcj1))
         {
            h=np->bonds[j]; np->bonds[j]=np->bonds[j+1]; np->bonds[j+1]=h;
            h=np->atoms[j]; np->atoms[j]=np->atoms[j+1]; np->atoms[j+1]=h;
         }
         else
             break;
      }
}

void SortBondsWithNumbering(struct reaccs_molecule_t *mp,
                            int numbering[])
/*
 * Sorts the bond_array member of *mp using numbering.
 * Note: the sorting preserves invariants with stereo-bonds.
 */
{
   int i, itmp, j;
   struct reaccs_bond_t bond, *bp;

   if (!numbering) return;
   /* Sort the bond atoms into the correct order if invariant */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      /* check if really swappable */
      if (bp->bond_type == SINGLE  &&  bp->stereo_symbol != NONE) continue;
      if (numbering[bp->atoms[0]-1] > numbering[bp->atoms[1]-1])
      {
         itmp = bp->atoms[0]; bp->atoms[0] = bp->atoms[1]; bp->atoms[1] = itmp;
      }
   }
   /* Now, we sort the bonds in ascending order */
   for (i=1; i<mp->n_bonds; i++)
      for (j=i-1; j>=0; j--)
      {
         if (numbering[mp->bond_array[j].atoms[0]-1] > 
               numbering[mp->bond_array[j+1].atoms[0]-1]  ||
             (numbering[mp->bond_array[j].atoms[0]-1] == 
                numbering[mp->bond_array[j+1].atoms[0]-1]  &&
              numbering[mp->bond_array[j].atoms[1]-1] > 
                numbering[mp->bond_array[j+1].atoms[1]-1]))
         {
            bond = mp->bond_array[j];
            mp->bond_array[j] = mp->bond_array[j+1];
            mp->bond_array[j+1] = bond;
         }
         else
            break;
      }
}

void SortNeighboursByNumbering(struct reaccs_molecule_t *mp,
                               neighbourhood_t *nbp,
                               int numbering[])
/*
 * Sorts the atom neighbourhood_t objects of nbp[0..mp->n_atoms] using the
 * atom ranking in numbering[].
 */
{
   int i, j, jj, tmp;

   for (i=0; i<mp->n_atoms; i++)
   {
      for (j=1; j<nbp[i].n_ligands; j++)
         for (jj=j-1; jj>=0; jj--)
            if (numbering[nbp[i].atoms[jj]] > numbering[nbp[i].atoms[jj+1]])
            {
               tmp = nbp[i].atoms[jj];
               nbp[i].atoms[jj] = nbp[i].atoms[jj+1];
               nbp[i].atoms[jj+1] = tmp;
               tmp = nbp[i].bonds[jj];
               nbp[i].bonds[jj] = nbp[i].bonds[jj+1];
               nbp[i].bonds[jj+1] = tmp;
            }
            else
               break;
   }
}

struct valence_table_entry      /* This is the table of possible valences */
   {                            /* of the chemical elements. Elements not */
      char *atom_type;          /* in the table are assumed to have no    */
      int   from_valence,       /* valence constraints at all.            */
            to_valence,
            step_valence,
            lone_pairs,
	    pair_deficit;
   } valence_table[] =  {
                        {"C",  4, 4, 1},        // make common elements hit first
                        {"H",  -1, 1, 2},
                        {"N",  3, 5, 2, 1},
                        {"O",  2, 2, 1, 2},
                        {"Cl", 1, 7, 2},
                        {"P",  3, 5, 2, 1},
                        {"S",  2, 6, 2, 2},
                        {"F",  1, 1, 1, 1},

                        {"H",  -1, 1, 2},

                        {"Li", -1, 1, 2},
                        {"Na", -1, 1, 2},
                        {"K",  -1, 1, 2},
                        {"Rb", -1, 1, 2},
                        {"Cs", -1, 1, 2},

                        {"Be", -2, 2, 2},
                        {"Mg", -2, 2, 2},

                        {"B",   3, 3, 1, 0, 1},

                        {"Al", -3, 3, 2, 0, 1},
                        {"Ga", -3, 3, 2, 0, 1},
                        {"In", -3, 3, 2, 0, 1},
                        {"Tl", -3, 3, 2},

                        // {"C",  4, 4, 1},
                        {"Si", 4, 4, 1},

                        // {"N",  3, 5, 2, 1},
                        // {"P",  3, 5, 2, 1},
                        {"As", 3, 5, 2},
                        {"Sb", 3, 5, 2},
                        {"Bi", 3, 5, 2},

                        // {"O",  2, 2, 1, 2},
                        // {"S",  2, 6, 2, 2},
                        {"Se", 2, 6, 2},
                        {"Te", 2, 6, 2},

                        {"La", -3, 3, 2, 0, 1},
                        {"Ce", -3, 3, 2, 0, 1},
                        {"Pr", -3, 3, 2, 0, 1},
                        {"Nd", -3, 3, 2, 0, 1},
                        {"Pm", -3, 3, 2, 0, 1},
                        {"Sm", -3, 3, 2, 0, 1},
                        {"Eu", -3, 3, 2, 0, 1},
                        {"Gd", -3, 3, 2, 0, 1},
                        {"Tb", -3, 3, 2, 0, 1},
                        {"Dy", -3, 3, 2, 0, 1},
                        {"Ho", -3, 3, 2, 0, 1},
                        {"Er", -3, 3, 2, 0, 1},
                        {"Tm", -3, 3, 2, 0, 1},
                        {"Yb", -3, 3, 2, 0, 1},
                        {"Lu", -3, 3, 2, 0, 1},

                        // {"F",  1, 1, 1, 1},
                        // {"Cl", 1, 7, 2},
                        {"Br", 1, 7, 2},
                        {"I",  1, 7, 2},

                        {(char *)NULL, 0}};

int
ImplicitHydrogens(char *symbol,
                  int   nsingle,
                  int   naromatic,
                  int   ndouble,
                  int   ntriple,
                  int   radical,
                  int   charge)
/*
 * Computes the number of implicit hydrogens attached to an atom of type
 * symbol with nsingle single bonds, naromatic aromatic bonds, ndouble
 * double bonds, ntriple triple bonds, and radical and charge state
 * radical and charge, resp.
 */
{
   int i,val,h;
   int bond_electrons;

   bond_electrons = nsingle+2*ndouble+3*ntriple;
   if (radical) bond_electrons++;
   switch (naromatic)
   {
      case 0: break;
      case 1: /* one aromatic bond could be an error */
              bond_electrons+=2;
              break;
      case 2: bond_electrons+=3;
              break;
      case 3: bond_electrons+=4;
              break;
      default:ShowMessageI("atom with %d aromatic bonds",
                           "ImplicitHydrogens",
                           naromatic);
              bond_electrons += naromatic+1;
              break;
   }

   for (i=0; valence_table[i].atom_type!=(char *)NULL; i++)
      if (0 == strcmp(valence_table[i].atom_type,symbol))
      {
         if (charge == 0)       /* Easy case */
         {
            for (val=valence_table[i].from_valence;
                 val<=valence_table[i].to_valence;
                 val+=valence_table[i].step_valence)
               if (0 <= (h=val-bond_electrons)) return(h);
         }
         else if (charge > 0)
         {
            for (val=valence_table[i].from_valence;
                 val<=valence_table[i].to_valence;
                 val+=valence_table[i].step_valence)
            {
               h=val-bond_electrons+charge;
               if (h < 0) continue;
               if (valence_table[i].lone_pairs > 0)
                  return (h);
               else
                  return (0);
            }
         }
         else /* charge < 0 */
         {
            for (val=valence_table[i].from_valence;
                 val<=valence_table[i].to_valence;
                 val+=valence_table[i].step_valence)
            {
               h=val-bond_electrons-charge;
               if (h < 0) continue;
               if (valence_table[i].pair_deficit<h)
                      h=valence_table[i].pair_deficit;
               /* hydrid ions are not stable if there are lone-pairs. */
               if (valence_table[i].lone_pairs > 0)
                  return (0);
               else
                  return (h);
            }
         }
      }

   return(0);
}

void
ComputeImplicitH(struct reaccs_molecule_t *mp,
                 int H_count[])
/*
 * Computes the implicit hydrogen counts for the atoms in *mp. The
 * result is returned in H_count[] starting at H_count[1] for atom
 * number 1. H_count[0] is not changed.
 *
 * Note: The elements of H_count are changed only if the are zero!!
 *       They need to be initialized by the caller!
 */
{
   int i;
   /* static */
   int *single_bond,      /* Counts of attached bonds of a type, */
       *aromatic_bond,    /* <array>[0] is unused */
       *double_bond,
       *triple_bond,
       *radical,	  /* total count of property */
       *charge;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;

   single_bond   = TypeAlloc(mp->n_atoms+1, int);
   aromatic_bond = TypeAlloc(mp->n_atoms+1, int);
   double_bond   = TypeAlloc(mp->n_atoms+1, int);
   triple_bond   = TypeAlloc(mp->n_atoms+1, int);
   radical       = TypeAlloc(mp->n_atoms+1, int);
   charge        = TypeAlloc(mp->n_atoms+1, int);

   for (i=0; i<mp->n_atoms; i++)
   {
      radical[i+1] = mp->atom_array[i].radical == DOUBLET;
      charge[i+1]  = mp->atom_array[i].charge;
   }

   for (i=0; i<=mp->n_atoms; i++)
      single_bond[i] = aromatic_bond[i] = double_bond[i] = triple_bond[i] = 0;

   for (i=0,bp=mp->bond_array; i<mp->n_bonds; i++,bp++)
   {
      switch (bp->bond_type)
      {
         case SINGLE: single_bond[bp->atoms[0]]++;
                      single_bond[bp->atoms[1]]++;
                      break;
         case DOUBLE: double_bond[bp->atoms[0]]++;
                      double_bond[bp->atoms[1]]++;
                      break;
         case TRIPLE: triple_bond[bp->atoms[0]]++;
                      triple_bond[bp->atoms[1]]++;
                      break;
         case AROMATIC: aromatic_bond[bp->atoms[0]]++;
                        aromatic_bond[bp->atoms[1]]++;
                        break;
         default :
            single_bond[bp->atoms[0]]++;
            single_bond[bp->atoms[1]]++;
            break;
      }
   }

   for (i=0,ap=mp->atom_array; i<mp->n_atoms; i++,ap++)
      if (H_count[i+1] == 0)
      {
	 H_count[i+1] = ImplicitHydrogens(ap->atom_symbol,
					  single_bond[i+1],
					  aromatic_bond[i+1],
					  double_bond[i+1],
					  triple_bond[i+1],
					  radical[i+1],
					  charge[i+1]);
	 if (H_count[i+1] < 0) H_count[i+1] = 0;
      }

   MyFree((char *)single_bond); MyFree((char *)aromatic_bond);
   MyFree((char *)double_bond); MyFree((char *)triple_bond);
   MyFree((char *)radical); MyFree((char *)charge);
}

int RemoveIsolatedAtoms(struct reaccs_molecule_t *mp)
/*
 * Scans *mp for atoms with no bonds attached to them and then deletes
 * those atoms.
 * Renumbers the bonds accordingly.
 * It returns the number of atoms removed.
 */
{
   int i, n;
   int nremove;
   struct reaccs_bond_t *bp;

   ResetColors(mp);

   for (i=0, bp=mp->bond_array;                 /* find connected atoms */
        i<mp->n_bonds;
        i++, bp++)
   {
      mp->atom_array[bp->atoms[0]-1].color++;
      mp->atom_array[bp->atoms[1]-1].color++;
   }

   nremove = 0;
   for (i=0, n=1; i<mp->n_atoms; i++)                 /* find new numbers */
      if (mp->atom_array[i].color != NO_COLOR)
         mp->atom_array[i].color = n++;
      else
         nremove++;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++) /* renumber bonds */
   {
      bp->atoms[0] = mp->atom_array[bp->atoms[0]-1].color;
      bp->atoms[1] = mp->atom_array[bp->atoms[1]-1].color;
   }

   for (i=0, n=0; i<mp->n_atoms; i++)                   /* squeeze atom list */
      if (mp->atom_array[i].color != NO_COLOR)
         mp->atom_array[n++] = mp->atom_array[i];

   if (n == 0  &&  mp->n_atoms != 0)            /* keep at least one atom */
   {
      mp->n_atoms = 1;
      return (nremove-1);
   }
   else
   {
      mp->n_atoms = n;
      return (nremove);
   }
}

typedef unsigned atom_pair[2];
void RingState(struct reaccs_molecule_t *mp,
               int atom_status[],
               int bond_status[])
/*
 * Computes how many basis rings each bond shares and how many
 * ring bonds are attached to an atom. The results are stored in
 * atom_status[] and bond_status[] respectively.
 */
{
   bond_set_node *rph, *ring_list;
   atom_pair *bonds;
   struct reaccs_bond_t *bp;
   int i;

   for (i=0; i<mp->n_atoms; i++)
      atom_status[i] = 0;
   for (i=0; i<mp->n_bonds; i++)
      bond_status[i] = 0;

   if (mp->n_bonds == 0) return;

   bonds = TypeAlloc(mp->n_bonds, atom_pair); /* get basis rings */
   for (i=0; i<mp->n_bonds; i++)
   {
      bonds[i][0] = mp->bond_array[i].atoms[0];
      bonds[i][1] = mp->bond_array[i].atoms[1];
   }
   ring_list = RingList(bonds,mp->n_bonds);
   ring_list = CombineRings(ring_list);
   MyFree((char *)bonds);

   for (rph=ring_list; rph; rph=rph->next)
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (IsMember(rph->bond_set,i))
          bond_status[i]++;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bond_status[i] > 0)
      {
         atom_status[bp->atoms[0]-1]++;
         atom_status[bp->atoms[1]-1]++;
      }

   DisposeBondSetList(ring_list);
}

void MakeRingsClockwise(struct reaccs_molecule_t *mp)
/*
 * Orients the bonds in *mp such that the small basis rings are
 * traced in a clockwise direction.
 *
 */
{
   struct reaccs_bond_t *bp;
   bit_set_t *hset;
   bond_set_node *ring_list, *plist;
   /* static */ unsigned graph[MAXBONDS][2];
   /* static */ unsigned bonds[MAXBONDS][2];
   /* static */ unsigned ring[MAXBONDS];
   double angle;
   int at1, at2, at3;
   int i, j, nsize;
   int changed;

   if (mp->n_bonds > MAXBONDS  ||  mp->n_atoms > MAXATOMS) return;

   for (i=0; i<mp->n_bonds; i++)        /* find all small basis rings */
   {
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
      mp->bond_array[i].topography = CHAIN;
      mp->bond_array[i].dummy      = 0;
   }
   ring_list = RingList(graph,mp->n_bonds);
   ring_list = CombineRings(ring_list);

   if (ring_list  &&  ring_list->next)  /* at least two rings */
   {
      for (plist=ring_list; plist; plist=plist->next)
      { /* misuse cardinality to store # of double bonds in ring */
         plist->cardinality = 0;
         for (i=0; i<mp->n_bonds; i++)
            if (IsMember(plist->bond_set,i)  &&
                mp->bond_array[i].bond_type != SINGLE)
               plist->cardinality++;
      }

      do        /* Bubble sort by increasing "cardinality" */
      {
         changed = FALSE;
	 for (plist=ring_list; plist->next; plist = plist->next)
            if (plist->cardinality > plist->next->cardinality)
            {
               hset = plist->bond_set;
               plist->bond_set = plist->next->bond_set;
               plist->next->bond_set = hset;
               nsize = plist->cardinality;
               plist->cardinality = plist->next->cardinality;
					plist->next->cardinality = nsize;
               changed = TRUE;
            }
      } while (changed);
   }

   for (plist=ring_list;        /* for all rings make bond chains */
        plist!=(bond_set_node *)NULL;
        plist=plist->next)
   {
      nsize = 0;                /* fetch bonds of ring */
      for (i=0; i<mp->n_bonds; i++)
         if (IsMember(plist->bond_set,i))
         {
            mp->bond_array[i].topography = RING;
            bonds[nsize][0] = mp->bond_array[i].atoms[0];
            bonds[nsize][1] = mp->bond_array[i].atoms[1];
            nsize++;
         }
      /* line up ring atoms */
      ring[0] = bonds[0][0]; ring[1] = bonds[0][1];
      for (i=2; i<nsize; i++)
         for (j=1; j<nsize; j++)
         {
            if (bonds[j][0] == ring[i-1]  &&  bonds[j][1] != ring[i-2])
            {
               ring[i] = bonds[j][1]; break;
            }
            if (bonds[j][1] == ring[i-1]  &&  bonds[j][0] != ring[i-2])
            {
               ring[i] = bonds[j][0]; break;
            }
         }
      angle = 0.0;      /* Sum up angles to test if clockwise */
      for (i=0; i<nsize; i++)
      {
         at1 = ring[i]; at2 = ring[(i+1)%nsize]; at3 = ring[(i+2)%nsize];
         angle += Angle(mp->atom_array[at1-1].x - mp->atom_array[at2-1].x,
                        mp->atom_array[at1-1].y - mp->atom_array[at2-1].y,
                        mp->atom_array[at3-1].x - mp->atom_array[at2-1].x,
                        mp->atom_array[at3-1].y - mp->atom_array[at2-1].y);
      }
      if (angle > (nsize-1.5)*3.1415926)
      {                 /* counter clockwise -> swap direction */
         for (i=0; i<nsize/2; i++)
         {
            j = ring[i]; ring[i] = ring[nsize-1-i]; ring[nsize-1-i] = j;
         }
      }
      for (i=0; i<nsize; i++)   /* fix bonds */
      {
         at1 = ring[i]; at2 = ring[(i+1)%nsize];
         for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
            if (bp->atoms[1]      == at1  &&
                bp->atoms[0]      == at2  &&
                (bp->stereo_symbol == NONE  || /* don't spoil stereo bonds!!! */
		 bp->stereo_symbol == CIS_TRANS_EITHER))
            {
               bp->atoms[0] = at1; bp->atoms[1] = at2; break;
            }
      }
   }

   DisposeBondSetList(ring_list);
}

double MolecularWeight(struct reaccs_molecule_t *mp)
/*
 * Computes and returns the molecular weight of the molecule *mp.
 */
{
   double result;
   int i;
   struct reaccs_atom_t *ap;
   int *H_counts;
   struct ptable_entry *ptp;

   result = 0;
   H_counts = TypeAlloc(mp->n_atoms+1, int);
   for (i=0; i<=mp->n_atoms; i++)
         H_counts[i] = 0;
   ComputeImplicitH(mp, H_counts);
   /* Take care of already fixed hydrogen counts */
   /* Note: Index origin of H_counts is 1 */
   for (i=0; i<mp->n_atoms; i++)
      if (mp->atom_array[i].query_H_count != NONE)
         H_counts[i+1] = mp->atom_array[i].query_H_count-ZERO_COUNT;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      for (ptp = ptable; ptp->symbol != NULL; ptp++)
	 if (strcmp(ptp->symbol, ap->atom_symbol) == 0) break;
      if (ptp->symbol != NULL)
	 result += ptp->mass;
      result += H_counts[i+1]*1.0080;
   }

   MyFree((char *)H_counts);
   return (result);
}

int MolecularFormula(struct reaccs_molecule_t *mp,
                     char formula[],
                     int bufsize)
/*
 * Collects the molecular formula of *mp and writes it to formula[].
 * Disconnected fragments are separated by a "." character.
 * It returns the number of fragments found.
 * bufsize is the usable length of the formula buffer including the terminal '\0'.
 */
{
   int i, j;
   int nfrag;
   int iseed, seed, first, color, changed;
   struct reaccs_bond_t *bp;
   struct reaccs_atom_t *ap;
   int *H_counts;
   int count, charge;

   static struct
   {
      char symbol[4];
      int  count;
   } element_counts[20] = {{"C",0}};
   int maxelem;
   char tmp_sym[4];
   int h;

   H_counts = TypeAlloc(mp->n_atoms+1, int);

   for (i=0; i<mp->n_atoms; i++) mp->atom_array[i].color = NO_COLOR;
   for (i=0; i<=mp->n_atoms; i++) H_counts[i] = 0;
   ComputeImplicitH(mp, H_counts);
   /* Take care of already fixed hydrogen counts */
   /* Note: Index origin of H_counts is 1 */
   for (i=0; i<mp->n_atoms; i++)
      if (mp->atom_array[i].query_H_count != NONE)
         H_counts[i+1] = mp->atom_array[i].query_H_count-ZERO_COUNT;
   nfrag = 0;
   formula[0] = '\0'; bufsize--;

   first = TRUE; color = NO_COLOR;
   for (iseed=(-1); iseed<(int)mp->n_atoms; iseed++)
   {
      if (iseed == (-1))  /* try to seed at carbon atom for first fragment */
      {
	 for (i=0; i<mp->n_atoms; i++)
	    if (0 == strcmp(mp->atom_array[i].atom_symbol, "C"))
	       break;
	 if (i < mp->n_atoms)
	 {
	    seed = i;
	 }
	 else
	    seed = 0;
      }
      else
	 seed = iseed;
      if (seed >= mp->n_atoms) break;

      if (mp->atom_array[seed].color != NO_COLOR)
         continue;
      else
      {
         nfrag++;
         mp->atom_array[seed].color = ++color;              /* next fragment */
         do
         {
            for (i=0, bp=mp->bond_array, changed=FALSE;     /* color fragment */
                 i<mp->n_bonds;
                 i++, bp++)
               if (mp->atom_array[bp->atoms[0]-1].color == NO_COLOR &&
                   mp->atom_array[bp->atoms[1]-1].color == color       ||
                   mp->atom_array[bp->atoms[1]-1].color == NO_COLOR &&
                   mp->atom_array[bp->atoms[0]-1].color == color)
               {
                  changed = TRUE;
                  mp->atom_array[bp->atoms[0]-1].color = color;
                  mp->atom_array[bp->atoms[1]-1].color = color;
               }
         } while (changed);

	 /* count explicit hydrogens */
         for (count=0, i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
            if (ap->color == color  &&
	        0 == strcmp("H", ap->atom_symbol))
	       count++;

	 /* count implicit hydrogens */
         for (charge=0, i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].color == color)
            {
               count += H_counts[i+1];
               charge += mp->atom_array[i].charge;
            }

         maxelem = 1; element_counts[0].count = 0;
         for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
            if (ap->color == color  &&
	        0 != strcmp("H", ap->atom_symbol)) /* Hs are special */
            {
               for (j=0; j<maxelem; j++)
                  if (0 == strcmp(element_counts[j].symbol,ap->atom_symbol))
                  {
                     element_counts[j].count++;
                     break;
                  }
               if (j == maxelem)
               {
                  strcpy(element_counts[j].symbol,ap->atom_symbol);
                  element_counts[j].count = 1;
                  maxelem++;
               }
            }
	 for (i=2; i<maxelem; i++)
	    for (j=i-1; j>=1; j--)
	       if (strcmp(element_counts[j].symbol,
	                  element_counts[j+1].symbol) > 0)
	       {
		  strcpy(tmp_sym, element_counts[j].symbol);
		  strcpy(element_counts[j].symbol, element_counts[j+1].symbol);
		  strcpy(element_counts[j+1].symbol, tmp_sym);
		  h = element_counts[j].count;
		  element_counts[j].count = element_counts[j+1].count;
		  element_counts[j+1].count = h;
	       }
	       else
		  break;

         if (first)                            /* . if more than one fragment */
            first = FALSE;
         else if (bufsize > 0)
	 {
	    sprintf(formula,".");
            bufsize -= strlen(formula);
            formula += strlen(formula);
	 }

         if (bufsize <= 0) break;

         if (element_counts[0].count != 0)     /* 'C' atoms */
            if (element_counts[0].count == 1)
	    {
               if (strlen(element_counts[0].symbol) < bufsize)
               {
                   sprintf(formula,"%s",  element_counts[0].symbol);
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }
            else
	    {
               if (element_counts[0].count < 1000  &&  strlen(element_counts[0].symbol)+3 < bufsize)
               {
                   sprintf(formula,"%s%d",element_counts[0].symbol,
                                          element_counts[0].count);
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }

         if (count != 0)                     /* 'H' atoms */
            if (count == 1)
	    {
               if (bufsize > 0)
               {
                   sprintf(formula,"H");
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }
            else
	    {
                if (bufsize > 1+3  &&  count < 3)
                {
                   sprintf(formula,"H%d",count);
                   bufsize -= strlen(formula);
	           formula += strlen(formula);
                }
	    }

         for (i=1; i<maxelem; i++)          /* all other atoms */
            if (element_counts[i].count == 1)
	    {
               if (bufsize > strlen(element_counts[i].symbol))
               {
                   if (0 == strcmp(element_counts[i].symbol, "R#"))
                      sprintf(formula,"%s",  "R");
                   else
                      sprintf(formula,"%s",  element_counts[i].symbol);
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }
            else
	    {
               if (bufsize > strlen(element_counts[i].symbol)+3  &&  element_counts[i].count < 1000)
               {
                   if (0 == strcmp(element_counts[i].symbol, "R#"))
                      sprintf(formula,"%s%d","R",
                                             element_counts[i].count);
                   else
                      sprintf(formula,"%s%d",element_counts[i].symbol,
                                             element_counts[i].count);
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }

         if (charge != 0)
            if (charge == 1)
	    {
               if (bufsize > 3)
               {
                   sprintf(formula,"(+)");
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }
            else if (charge == (-1))
	    {
               if (bufsize > 3)
               {
                   sprintf(formula,"(-)");
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }
            else if (charge > 1)
	    {
               if (bufsize > 3+2)
               {
                   sprintf(formula,"(%d+)",charge);
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }
            else
	    {
               if (bufsize > 3+2)
               {
                   sprintf(formula,"(%d-)",-charge);
                   bufsize -= strlen(formula);
                   formula += strlen(formula);
               }
	    }
      }
   }
   MyFree((char *)H_counts);
   return (nfrag);
}

int TotalCharge(struct reaccs_molecule_t *mp)
/*
 * Returns the total charge of the atoms in *mp.
 */
{
   struct reaccs_atom_t *ap;
   int i, result;

   for (i=result=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      result += ap->charge;

   return (result);
}

int CombineChargeRadical(int charge, int radical)
/*
 * Returns MDL's combined charge/radical code.
 */
{
   if (radical == NONE  &&
       (-3) <= charge   &&
       charge <= 3   &&
       charge != 0)
      return (RADICAL-charge);
   else if (charge == 0  &&  radical == NONE)    return (NONE);
   else if (charge == 0  &&  radical == DOUBLET) return (RADICAL);
   else
      return (ANY_CHARGE);
}

void SplitChargeRadical(int charge_radical, int *chargep, int *radicalp)
/*
 * Splits MDL's combined charge/radical code into it's components.
 */
{
   if (charge_radical == NONE)
   {
      *chargep = NONE; *radicalp = NONE;
   }
   else if (charge_radical == RADICAL)
   {
      *chargep = NONE; *radicalp = DOUBLET;
   }
   else
   {
      *chargep = RADICAL-charge_radical;
      *radicalp = NONE;
   }
}

int IsFieldHeader(char *header, char *field)
/*
 * Returns TRUE if header is a valid SD file data section header
 * for the field *field.
 */
{
   char *cp;

   if (*header != '>') return (FALSE);

   for (cp=header; *cp; cp++)
      if (*cp == '<') break;
   cp++;

   if (!*cp) return (FALSE);

   while (*field)
      if (toupper(*cp) == toupper(*field))
      {
         cp++; field++;
      }
      else
         break;

   return (!*field);
}

void FlipStereoSymbols(struct reaccs_molecule_t *mp, int color)
/*
 * Flips the stereo symboles of the parts of the molecule *mp that are
 * colored with color.
 */
{
   int i;
   struct reaccs_bond_t *bp;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (color == NONE  ||
          (mp->atom_array[bp->atoms[0]-1].color == color &&
           mp->atom_array[bp->atoms[1]-1].color == color))
         if (bp->stereo_symbol == UP)
       bp->stereo_symbol = DOWN;
    else if (bp->stereo_symbol == DOWN)
        bp->stereo_symbol = UP;
}

void FlipMolecule(struct reaccs_molecule_t *mp, int color)
/*
 * Flips the parts of the molecule *mp that are colored with color
 * with respect to the vertical (Y) axis. If color is 'NONE',
 * the whole molecule is flipped.
 */
{
   int i;
   int natoms;
   double xcenter;

   xcenter = 0.0; natoms = 0;
   for (i=0; i<mp->n_atoms; i++)
      if (color == NONE || mp->atom_array[i].color == color)
      {
    xcenter += mp->atom_array[i].x;
         natoms++;
      }

   if (natoms == 0) return;

   xcenter /= natoms;

   for (i=0; i<mp->n_atoms; i++)
      if (color == NONE || mp->atom_array[i].color == color)
    mp->atom_array[i].x = xcenter-mp->atom_array[i].x;

   FlipStereoSymbols(mp, color);
}

void FloodFillFragments(struct reaccs_molecule_t *mp)
/*
 * The color of colored atoms is smeared out to the whole fragment
 * containig that atom. If there are more than one distict color,
 * the winning one on each fragment is arbitrary.
 */
{
   int changed;
   int i, j;
   int col1, col2;
   struct reaccs_bond_t *bp;

   do
   {
      changed = FALSE;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->bond_type & 0xF)
         {
            col1 = mp->atom_array[bp->atoms[0]-1].color;
            col2 = mp->atom_array[bp->atoms[1]-1].color;
            if (col1 == col2) continue;
            if (col1 != NONE)
            {
               if (col2 == NONE)        /* not yet colored */
                  mp->atom_array[bp->atoms[1]-1].color = col1;
               else                     /* different colors */
                  for (j=0; j<mp->n_atoms; j++)
                     if (mp->atom_array[j].color == col2)
                        mp->atom_array[j].color = col1;
               changed = TRUE;
            }
         else if (col2 != NONE)
         {
            mp->atom_array[bp->atoms[0]-1].color = col2;
            changed = TRUE;
         }
      }
   } while (changed);
}

int FloodWithColor(struct reaccs_molecule_t *mp,
                   neighbourhood_t *nbp,
	               int aindex, int color)
/*
 * Recursively colors uncolored parts of *mp starting with the atom with
 * index aindex with color. It uses a flood fill algorithm.
 * It returns the number of atoms recolored during the process.
 */
{
   int result;
   int i;

   mp->atom_array[aindex].color = color; result = 1;
   for (i=0; i<nbp[aindex].n_ligands; i++)
      if (mp->atom_array[nbp[aindex].atoms[i]].color == NO_COLOR)
	 result += FloodWithColor(mp, nbp, nbp[aindex].atoms[i], color);

   return (result);
}

void StripColoredPart(struct reaccs_molecule_t *mp, int color)
/*
 * Removes all atoms and adjacent bonds from *mp that are colored with color.
 */
{
    neighbourhood_t *nbp;
    int *good_atoms, *good_bonds;
    struct reaccs_atom_t *ap;
    struct reaccs_bond_t *bp;
    int i;

    good_atoms = TypeAlloc(mp->n_atoms+1, int);
    good_bonds = TypeAlloc(mp->n_bonds, int);
    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
        good_atoms[i+1] = ap->color != color;
    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
        good_bonds[i] = mp->atom_array[bp->atoms[0]-1].color != color  &&
                        mp->atom_array[bp->atoms[1]-1].color != color;
    StripMolecule(mp, good_atoms, good_bonds);
    MyFree((char *)good_atoms);
    MyFree((char *)good_bonds);
}

#ifdef __TURBOC__
#include <setjmp.h>
void ShowStack(char header[])
{
#ifndef __WIN32__
	jmp_buf jmpb;
   setjmp(jmpb);
   printf("%s:  stack pointer = %ld\n",header,(long)(jmpb->j_sp));
#endif
	return;
}
#else
void ShowStack(char header[])
{
   return;
}
#endif
