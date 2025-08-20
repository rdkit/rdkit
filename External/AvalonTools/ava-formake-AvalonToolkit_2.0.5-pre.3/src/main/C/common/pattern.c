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
/*    File:           pattern.c                                         */
/*                                                                      */
/*    Purpose:        This file implements the functions to translate   */
/*                    one augmented atom pattern to another. This       */
/*                    service is mainly used for structure              */
/*                    standardization e.g. by struchk.exe.              */
/*                                                                      */
/*    History:        23-Dec-1992     Switched to ANSI prototypes.      */
/*                                                                      */
/************************************************************************/
 
#include "pattern.h"

#include <string.h>
#include <ctype.h>

#include "local.h"
#include "reaccs.h"
#include "symboltable.h"
#include "utilities.h"

int RecMatch(struct reaccs_molecule_t *mp,
             unsigned int match[],
             unsigned int level,
             augmented_atom_t *aap,
             neighbourhood_t nbp[])
/*
 * Recursively searches *mp for the next atom to be appended to
 * match[0..level]. Returns TRUE if all ligands in aap are mapped
 * successfully. nbp[i] describes the neighbour bonds and atoms
 * of atom i+1.
 */
{
   register unsigned short i, j;
   int is_new;
   neighbourhood_t *nbph;

   if (level == aap->n_ligands) return(TRUE);

   nbph = &nbp[match[0]];
   for (i=0; i<nbph->n_ligands; i++)
      if ((mp->atom_array[nbph->atoms[i]].charge == aap->ligands[level].charge   ||
           aap->ligands[level].charge == ANY_CHARGE)                                &&
          (mp->atom_array[nbph->atoms[i]].radical == aap->ligands[level].radical ||
           aap->ligands[level].radical == ANY_RADICAL)                              &&
          mp->bond_array[nbph->bonds[i]].bond_type == aap->ligands[level].bond_type &&
          (aap->ligands[level].s_count == NONE                                   ||
           nbp[nbph->atoms[i]].n_ligands == aap->ligands[level].s_count)            &&
         AtomSymbolMatch(mp->atom_array[nbph->atoms[i]].atom_symbol, aap->ligands[level].atom_symbol))
      {
         is_new = TRUE; match[level+1] = nbph->atoms[i];
         for (j=0; j<=level; j++)
           if (nbph->atoms[i] == match[j]) is_new = FALSE;
         if (is_new && RecMatch(mp,match,level+1,aap,nbp))
            return(TRUE);
      }

   return(FALSE);
}

int AAMatch(struct reaccs_molecule_t *mp,
            unsigned int i,
            unsigned int match[],
            augmented_atom_t *aap,
            int atom_ring_status[],
            neighbourhood_t nbp[])
/*
 * Tests if atom i in *mp matches the augmented atom description
 * *aap. nbp[] is used to speed up access to neighbour atoms and
 * bonds. The first matching atom mapping is placed into match[1..].
 * i is stored in match[0].
 */
{
   if (nbp[i].n_ligands == aap[0].n_ligands                                            &&
       (aap[0].charge == ANY_CHARGE  || mp->atom_array[i].charge == aap[0].charge)     &&
       (aap[0].radical == ANY_RADICAL  || mp->atom_array[i].radical == aap[0].radical) &&
       AtomSymbolMatch(mp->atom_array[i].atom_symbol,aap[0].atom_symbol))
   {
      if (!IsNULL(atom_ring_status)  &&  aap->topography == RING   &&  atom_ring_status[i] == 0) return (FALSE);
      if (!IsNULL(atom_ring_status)  &&  aap->topography == CHAIN  &&  atom_ring_status[i] != 0) return (FALSE);
      match[0] = i;
      return (RecMatch(mp,match,0,aap,nbp));
   }
   else return (FALSE);
}

void TransformAA(struct reaccs_molecule_t *mp,
                 int i,
                 unsigned int match[MAXNEIGHBOURS+1],
	         augmented_atom_t *aapin,
	         augmented_atom_t *aapout,
                 int atom_ring_state[],
	         neighbourhood_t nbp[])
/*
 * Transforms the augmented atom with index i in molecule *mp. match[]
 * contains the mapping of the ligands in the augmented atom to the
 * corresponding ligands in *mp. *aapin is the source augmented atom
 * and aapout is the target. Atom and bond attributes are changed only
 * if they differ in aapin and aapout, thus allowing for 'generic' attributes
 * in the matching augmented atom descriptions. nbp gives the atom and bond
 * neighbouring to speed-up access.
 */
{
   int j, k;
                                        /* fix central atom */
   if (aapin->charge != aapout->charge)
      mp->atom_array[i].charge = aapout->charge;
   if (aapin->radical != aapout->radical)
      mp->atom_array[i].radical = aapout->radical;
   if (0 != strcmp(aapin->atom_symbol,aapout->atom_symbol))
      strcpy(mp->atom_array[i].atom_symbol,aapout->atom_symbol);
                                        /* fix ligand atoms/bonds */
   for (j=0; j<nbp[i].n_ligands; j++)
      for (k=0; k<aapout->n_ligands; k++)
         if (nbp[i].atoms[j] == match[k+1])
         {
            if (aapin->ligands[k].charge != aapout->ligands[k].charge)
               mp->atom_array[nbp[i].atoms[j]].charge =
            aapout->ligands[k].charge;
            if (aapin->ligands[k].radical != aapout->ligands[k].radical)
               mp->atom_array[nbp[i].atoms[j]].radical =
                  aapout->ligands[k].radical;
            if (0 != strcmp(aapin->ligands[k].atom_symbol,
                            aapout->ligands[k].atom_symbol))
               strcpy(mp->atom_array[nbp[i].atoms[j]].atom_symbol,
                      aapout->ligands[k].atom_symbol);
            if (aapin->ligands[k].bond_type !=
                aapout->ligands[k].bond_type)
               mp->bond_array[nbp[i].bonds[j]].bond_type =
                  aapout->ligands[k].bond_type;
         }
}

int AAFix(struct reaccs_molecule_t *mp,
          augmented_atom_t aap[2],
          neighbourhood_t nbp[])
/*
 * Searches for all occurences of augmented atom aap[0] in *mp and
 * replaces them with the augmented atom aap[1]. The function 
 * returns TRUE if there was a match for aap[0]. nbp[] is used
 * for rapid access of neighboour bonds and atoms.
 * Atom symbol strings are only changed if corresponding atom symbols in
 * aap[0] and aap[1] don't match. This allows atom type lists with-out
 * introducing them into the fixed molecule.
 */
{
   unsigned int i;
   int matched; 
   unsigned int match[MAXNEIGHBOURS+1];
   char buffer[1024];
   int *atom_status, *bond_status;

   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);

   matched = FALSE;
   buffer[0] = '\0';
   for (i=0; i<mp->n_atoms; i++)
   {
      if (AAMatch(mp,i,match,aap,atom_status,nbp))
      {
         matched = TRUE;

         TransformAA(mp,i,match,&aap[0],&aap[1],atom_status,nbp);
         snprintf(buffer, 1023, "%s %d : '%s' -> '%s'",
                 mp->name, i+1,
                 aap[0].short_name,
                 aap[1].short_name);
         AddMsgToList(buffer);
      }
   }
   MyFree((char *)bond_status);
   MyFree((char *)atom_status);

   return(matched);
}

int TransformAugmentedAtoms(struct reaccs_molecule_t *mp,
                            augmented_atom_t table[][2],
                            int nentries)
/*
 * Scans the molecule *mp for occurrences of the augmented atoms in
 * table[0..nentries-1][0] and replaces the matches with the corresponding
 * atom in table[0..nentries-1][1].
 */
{
   int i, n;
   neighbourhood_t *neighbour_array;
   int result;

   if (mp->n_atoms <= 0) return (FALSE);

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

   result = FALSE;
   for (i=0; i<nentries; i++)
   {
      if (AAFix(mp,table[i],neighbour_array)) result = TRUE;
   }

   for (i=0,n=0; i<mp->n_bonds; i++)     /* remove deleted bonds */
      if (mp->bond_array[i].bond_type != NONE)
         mp->bond_array[n++] = mp->bond_array[i];
   mp->n_bonds = n;

   MyFree((char *) neighbour_array);

   return (result);
}

char *AddSymbol(char *cp, const char *symbol, int *nbufp)
/*
 * Copies the string symbol[] to position cp, decrements *nbufp with the
 * length of symbol[], and returns the position in cp[] after the copy.
 * If nbufp would become less than 0, an error is reported and cp is returned.
 */
{
   int len;

   if (!symbol) return (cp);
   len = strlen(symbol);
   if (*nbufp <= len)
   {
      fprintf(stderr,"AddSymbol: buffer overflow by symbol '%s'\n",symbol);
      return (cp);
   }
   strcpy(cp,symbol); return (cp+len);
}

#define MAX_LOCAL_BUFFER      1280

char *H_counttext[] =
   {
      "",        /* ANY_COUNT  */
      "H0",      /* ZERO_COUNT */
      "H",       /* AT_LEAST_1 */
      "H",       /* AT_LEAST_2 */
      "H",       /* AT_LEAST_3 */
                 /*  pseudo counts */
      "H",       /* AT_LEAST_4 */
      "H",       /* AT_LEAST_5 */
      "H",       /* AT_LEAST_6 */
   };

#define  NONE              0
#define  SINGLE            1
#define  DOUBLE            2
#define  TRIPLE            3
#define  AROMATIC          4
#define  SINGLE_DOUBLE     5
#define  SINGLE_AROMATIC   6
#define  DOUBLE_AROMATIC   7

 /* only used in ConvertMoleculeToString */
static
char *bondtext[9] = {"?",       /* NONE   */
                     "-",       /* SINGLE */
                     "=",       /* DOUBLE */
                     "^",       /* TRIPLE */
                     "%",       /* AROMATIC */
                     ":",       /* SINGLE_DOUBLE */
                     "?",       /* SINGLE_AROMATIC */
                     "#",       /* DOUBLE_AROMATIC */
                     "~"};      /* ANY */

void TraceBranch(char buffer[MAX_LOCAL_BUFFER],
                 int root,
                 int *nmappedp,
                 unsigned short branch_bond,
                 unsigned short branch_atom,
                 struct reaccs_molecule_t *mp,
                 neighbourhood_t nbarr[])
/*
 * Traces the branch starting at atom root through bond branch_bond
 * and branch_atom. It writes the corresponding string to buffer.
 * every atom visited during this (possibly recursive) process is
 * colored by its sequence number starting with *nmappedp+1.
 * *nvisited is incremented accordingly.
 */
{
   neighbourhood_t *nbp;
   int nbuf;
   int i;
   char local_buffer[MAX_LOCAL_BUFFER];

            /* add starting bond and atom to path */
   nbuf = MAX_LOCAL_BUFFER;
   buffer = AddSymbol(buffer,
                      bondtext[mp->bond_array[branch_bond].bond_type],
                      &nbuf);
   buffer = AddSymbol(buffer,
                     mp->atom_array[branch_atom].atom_symbol,
                     &nbuf);
   buffer = AddSymbol(buffer,
                     H_counttext[mp->atom_array[branch_atom].query_H_count],
                     &nbuf);
   buffer = AddSymbol(buffer,
                     SAFEGUARD(IdToString(charge_to_string, mp->atom_array[branch_atom].charge)),
                     &nbuf);
   buffer = AddSymbol(buffer,
                     SAFEGUARD(IdToString(radical_to_string, mp->atom_array[branch_atom].radical)),
                     &nbuf);

            /* mark new atom as already mapped */
   mp->atom_array[branch_atom].color = ++(*nmappedp);
      
                /* look at neighbours of new atom */
   nbp = &nbarr[branch_atom];
   local_buffer[0] = '\0';
   for (i=0; i<nbp->n_ligands; i++)
      if (nbp->atoms[i] != root &&
          mp->atom_array[nbp->atoms[i]].color < mp->atom_array[branch_atom].color)
      {
         if (local_buffer[0] != '\0')   /* append previous branch */
         {
            buffer = AddSymbol(buffer, "(", &nbuf);
            buffer = AddSymbol(buffer, local_buffer, &nbuf);
            buffer = AddSymbol(buffer, ")", &nbuf);
            local_buffer[0] = '\0';
         }

         if (mp->atom_array[nbp->atoms[i]].color != NO_COLOR)
         {                              /* ring closure */
            sprintf(local_buffer,"%s@%d",
                    bondtext[mp->bond_array[nbp->bonds[i]].bond_type],
                    mp->atom_array[nbp->atoms[i]].color);
         }
         else if (nbp->atoms[i] != root)
         {
            TraceBranch(local_buffer,
                        (int) branch_atom,
                        nmappedp,
                        nbp->bonds[i],
                        nbp->atoms[i],
                        mp,
                        nbarr);
         }
      }
   if (local_buffer[0] != '\0')
      buffer = AddSymbol(buffer, local_buffer, &nbuf);
}

int ConvertMoleculeToString(struct reaccs_molecule_t *mp,
                            char string[], int nbuf)
/*
 * Computes a string representation of molecule *mp. The result is
 * put into string[0..nbuf-1]. The function returns TRUE if it
 * successfully translated the structure and FALSE otherwise.
 */
{
   neighbourhood_t *neighbour_array, *nbp;
   char *continuep;
   int root, minlig;
   int i;
   int nmapped;         /* number of atoms already mapped in current string */
   char local_buffer[MAX_LOCAL_BUFFER];

   if (IsNULL(mp))
   {
      strcpy(string,"/");
      return (TRUE);
   }

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);

   SetupNeighbourhood(mp, neighbour_array,mp->n_atoms);
   ResetColors(mp);

   continuep = string;
   nmapped = 0;
   for (;;)
   {
      minlig = MAXNEIGHBOURS;
      for (i=0, root=(-1); i<mp->n_atoms; i++)
         if (mp->atom_array[i].color == NO_COLOR &&
             neighbour_array[i].n_ligands < minlig)
         {
             minlig = neighbour_array[i].n_ligands;
             root = i;
         }
      if (root == (-1)) break;    /* no new root */
      if (string != continuep)
         continuep = AddSymbol(continuep, " ", &nbuf);

      mp->atom_array[root].color = ++nmapped;
      continuep = AddSymbol(continuep,
                            mp->atom_array[root].atom_symbol,
                            &nbuf);
      continuep = AddSymbol(continuep,
                            H_counttext[mp->atom_array[root].query_H_count],
                            &nbuf);
      continuep = AddSymbol(continuep,
                            SAFEGUARD(IdToString(charge_to_string,
                            mp->atom_array[root].charge)),
                            &nbuf);
      continuep = AddSymbol(continuep,
                            SAFEGUARD(IdToString(radical_to_string,
                            mp->atom_array[root].radical)),
                            &nbuf);
      nbp = &neighbour_array[root];
      local_buffer[0] = '\0';
      for (i=0; i<nbp->n_ligands; i++)
         if (mp->atom_array[nbp->atoms[i]].color == NO_COLOR)
         {
            if (local_buffer[0] != '\0')
            {
               continuep = AddSymbol(continuep, "(", &nbuf);
               continuep = AddSymbol(continuep, local_buffer, &nbuf);
               continuep = AddSymbol(continuep, ")", &nbuf);
               local_buffer[0] = '\0';
            }
            TraceBranch(local_buffer,
                        root,
                        &nmapped,
                        nbp->bonds[i],
                        nbp->atoms[i],
                        mp,
                        neighbour_array);
         }
      if (local_buffer[0] != '\0') continuep = AddSymbol(continuep, local_buffer, &nbuf);
   }
   continuep = AddSymbol(continuep, "/", &nbuf);
   *continuep = '\0';
   
   MyFree((char *)neighbour_array);

   return (TRUE);
}

int AAPrint(FILE *fp,
            struct reaccs_molecule_t *mp,
            int i,
            neighbourhood_t *np)
/*
 * Prints a string representation of the neighbourhood of atom
 * i in molecule *mp to the file *fp.
 * It uses the neighbourhood information stored in *np.
 */
{
   int j;
   char buffer[100];
   int len;

   SortNeighbourhood(np,mp);

   sprintf(buffer,"%s%s%s",
    mp->atom_array[i].atom_symbol,
          SAFEGUARD(IdToString(charge_to_string, mp->atom_array[i].charge)),
      SAFEGUARD(IdToString(radical_to_string, mp->atom_array[i].radical)));
   fprintf(fp,"%s",buffer);
   len = strlen(buffer);
   
   for (j=0; j<np->n_ligands; j++)
   {
      sprintf(buffer,"(%s%s%s%s)",
              SAFEGUARD(IdToString(bond_to_string, mp->bond_array[np->bonds[j]].bond_type)),
              mp->atom_array[np->atoms[j]].atom_symbol,
              SAFEGUARD(IdToString(charge_to_string, mp->atom_array[np->atoms[j]].charge)),
              SAFEGUARD(IdToString(radical_to_string, mp->atom_array[np->atoms[j]].radical)));
      fprintf(fp,"%s",buffer);
      len += strlen(buffer);
   }
    
   return (len);
}

char *AAToString(char buffer[],
                 struct reaccs_molecule_t *mp,
                 int i,
                 neighbourhood_t *np)
/*
 * Prints a string representation of the neighbourhood of atom
 * i in molecule *mp to the string buffer.
 * It uses the neighbourhood information stored in *np.
 */
{
   int j;
   char *old_buffer;

   old_buffer = buffer;
   SortNeighbourhood(np,mp);

   sprintf(buffer,"%s%s%s",
           mp->atom_array[i].atom_symbol,
           SAFEGUARD(IdToString(charge_to_string, mp->atom_array[i].charge)),
           SAFEGUARD(IdToString(radical_to_string, mp->atom_array[i].radical)));
   buffer += strlen(buffer);
   
   for (j=0; j<np->n_ligands; j++)
   {
      sprintf(buffer,"(%s%s%s%s)",
              SAFEGUARD(IdToString(bond_to_string, mp->bond_array[np->bonds[j]].bond_type)),
              mp->atom_array[np->atoms[j]].atom_symbol,
              SAFEGUARD(IdToString(charge_to_string, mp->atom_array[np->atoms[j]].charge)),
              SAFEGUARD(IdToString(radical_to_string, mp->atom_array[np->atoms[j]].radical)));
      buffer += strlen(buffer);
   }
    
   return (old_buffer);
}

char *AtsymCopy(char *symbol)
/*
 * Returns a pointer to a copy of *symbol. This function is intended
 * to save space by returning pointers to static strings for
 * common atom types.
 */
{
   char *result;
   static char *Carbon = "C";
   static char *Nitrogen = "N";
   static char *Oxigen = "O";
   static char *Phosphorus = "P";
   static char *Boron = "B";

   /*
   if (symbol[0] == 'C'  &&  symbol[1] == '\0') return (Carbon);
   if (symbol[0] == 'N'  &&  symbol[1] == '\0') return (Nitrogen);
   if (symbol[0] == 'O'  &&  symbol[1] == '\0') return (Oxigen);
   if (symbol[0] == 'P'  &&  symbol[1] == '\0') return (Phosphorus);
   if (symbol[0] == 'B'  &&  symbol[1] == '\0') return (Boron);
   */

   result = TypeAlloc(1+strlen(symbol),char);
   strcpy(result,symbol);

   return (result);
}

int StringToAugmentedAtom(char *string,
                          augmented_atom_t *aap)
/*
 * Reads string[] and constructs the corresponding augmented atom *aap.
 * Returns TRUE if scan was successful and FALSE otherwise.
 * The syntax of a augmented atom string is as follows:
 *
 *   <AAString> ::= ['@'|'!@'] <Atom Symbol> [<Charge>]
 *                 {'(' <Bond Symbol> <Atom Symbol> [<Charge>] ')'}.
 * Bond symbols and charge descriptors are defined in the symbol tables
 * 'bond_to_string' and 'charge_to_string'.
 * '@' means central atom is in ring, '!@' means central atom is not in ring, omitting means any topography could match.
 */
{
   int i;
   char *oldstring;
   char atsym[300];
   char ligsymbols[MAXNEIGHBOURS][200];
   symbol_entry_t *stp;

   oldstring = string;

   aap->topography = NONE;
   if (string[0] == '@')
   {
      aap->topography = RING;
      string++;
   }
   else if (string[0] == '!'  &&  string[1] == '@')
   {
      aap->topography = CHAIN;
      string++; string++;
   }
   aap->n_ligands = 0;                       /* initially no ligands */
   if (!isalpha(string[0])) return (FALSE);  /* atom symbol is required */

   i = 0;                                    /* fetch atom symbol */
   while (isalnum(*string)  ||  *string == ','  ||  *string == '#')
      atsym[i++] = *string++;
   atsym[i] = '\0';

   if (string == strpbrk(string,"+-"))            /* charge definition */
   {
      for (stp = charge_to_string;
           !IsNULL(stp->symbol_string);
           stp++)
         if (stp->symbol_string[0] != '\0'  &&
             0 == strncmp(string,
                          stp->symbol_string,
                          strlen(stp->symbol_string)))
         {
            aap->charge = stp->symbol_id;
            string += strlen(stp->symbol_string);
            break;
         }
         if (IsNULL(stp->symbol_string))    /* no match */
             return (FALSE);
   }
   else aap->charge = NONE;

   if (string == strpbrk(string,".|:")  &&  /* radical definition */
       !isdigit(string[1]))                 /* don't confuse with degree specification */
   {
      for (stp = radical_to_string; !IsNULL(stp->symbol_string); stp++)
         if (stp->symbol_string[0] != '\0'  &&
             0 == strncmp(string,
                          stp->symbol_string,
                          strlen(stp->symbol_string)))
         {
            aap->radical = stp->symbol_id;
            string += strlen(stp->symbol_string);
            break;
         }
      if (IsNULL(stp->symbol_string))    /* no match */
         return (FALSE);
   }
   else aap->radical = NONE;

   while (*string == '(')                 /* read ligand descriptions */
   {
      string++;                               /* skip '(' */

      for (stp = bond_to_string;          /* read bond type */
          !IsNULL(stp->symbol_string);
          stp++)
         if (0 == strncmp(string,
                          stp->symbol_string,
                          strlen(stp->symbol_string)))
         {
            aap->ligands[aap->n_ligands].bond_type = stp->symbol_id;
            string += strlen(stp->symbol_string);
            break;
         }
      if (IsNULL(stp->symbol_string))    /* no match */
         return (FALSE);

      if (!isalpha(string[0]))           /* there must be an atom symbol */
         return (FALSE);
      i = 0;                             /* fetch atom symbol */
      while (isalnum(*string)  ||  *string == ','  ||  *string == '#')
         ligsymbols[aap->n_ligands][i++] = *string++;
      ligsymbols[aap->n_ligands][i] = '\0';

      if (string == strpbrk(string,"+-"))            /* charge definition */
      {
         for (stp = charge_to_string; !IsNULL(stp->symbol_string); stp++)
            if (stp->symbol_string[0] != '\0'  &&
                0 == strncmp(string,
                             stp->symbol_string,
                             strlen(stp->symbol_string)))
            {
               aap->ligands[aap->n_ligands].charge = stp->symbol_id;
               string += strlen(stp->symbol_string);
               break;
            }
         if (IsNULL(stp->symbol_string))    /* no match */
            return (FALSE);
      }
      else aap->ligands[aap->n_ligands].charge = NONE;

      if (string == strpbrk(string,".|:")  &&    /* radical definition */
          !isdigit(string[1]))   /* don't confuse with degree specification */
      {
         for (stp = radical_to_string; !IsNULL(stp->symbol_string); stp++)
            if (stp->symbol_string[0] != '\0'  &&
                0 == strncmp(string,
                             stp->symbol_string,
                             strlen(stp->symbol_string)))
            {
               aap->ligands[aap->n_ligands].radical = stp->symbol_id;
               string += strlen(stp->symbol_string);
               break;
            }
         if (IsNULL(stp->symbol_string))    /* no match */
            return (FALSE);
      }
      else aap->ligands[aap->n_ligands].radical = NONE;

      /* substitution count descriptor */
      if (string[0] == ':'  &&  isdigit(string[1]))
      {
         aap->ligands[aap->n_ligands].s_count = string[1]-'0';
         string += 2;
      }
      else
         aap->ligands[aap->n_ligands].s_count = 0;

      if (*string != ')')           /* look for ')' */
         return (FALSE);
      else
      {
         string++;
         aap->n_ligands++;
      }
   }

   if (*string != '\0')                    /* additional chars. after scan? */
      return (FALSE);
   else                                 /* correct AA description */
   {                                        /*  -> allocate memory    */
      aap->atom_symbol = AtsymCopy(atsym);
      aap->short_name  = TypeAlloc(1+strlen(oldstring),char);
      strcpy(aap->short_name,oldstring);

      for (i=0; i<aap->n_ligands; i++)
        aap->ligands[i].atom_symbol = AtsymCopy(ligsymbols[i]);
   }

   return (TRUE);
}

char aa_check_version[10] = {'\0'};

augmented_atom_t *ReadAugmentedAtoms(FILE *fp, int *natomp)
/*
 * Reads file fp and constructs an array of augmented atom
 * descriptions.
 * The function expects the number of augmented atom description
 * on the first line and the strings corresponding to the data structures
 * on the *natomp following lines.
 * Only the portion of the strings between the two '"' characters is used,
 * other characters are regarded as comments.
 */
{
   char atom_string[500];
   augmented_atom_t *atoms;
   int i;
   char buffer[80];
   char *cp;
   int vers_len;
   int success = TRUE;

   fscanf(fp,"%d",natomp);

   fgets(buffer,80,fp);
   cp = strchr(buffer,'_');
   if (cp && strchr(cp+1,'_'))
   {
      vers_len = (int)(strchr(cp+1,'_')-cp)-1;
      strncpy(aa_check_version, cp+1, vers_len);
      aa_check_version[vers_len] = '\0';
   }

   if (log_file)
      fprintf(log_file, "augmented atom check version = %s\n", aa_check_version);

   atoms = TypeAlloc(*natomp,augmented_atom_t);

   for (i=0; i<*natomp; i++)
   {
      SearchChar(fp,'"');

      fscanf(fp,"%499[^\"]",atom_string); SearchChar(fp,'\n');
      if (!StringToAugmentedAtom(atom_string,&atoms[i]))
      {
         fprintf(stderr, "ReadAugmentedAtoms: unsuccessful translation of %s\n", atom_string);
         success = FALSE;
      }
   }
   if (success) return (atoms);

   MyFree((char *)atoms); (*natomp) = 0;
   return ((augmented_atom_t *)NULL);
}

char aa_trans_version[10] = {'\0'};

aa_pair *ReadAAPairs(FILE *fp, int *ntransp)
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
{
   char trans_string[1000];
   aa_pair *trans_pairs;
   int i;
   char buffer[80];
   char *cp;
   int vers_len;
   int success = TRUE;

   fscanf(fp,"%d",ntransp);

   fgets(buffer,80,fp);
   cp = strchr(buffer,'_');
   if (cp && strchr(cp+1,'_'))
   {
      vers_len = (int)(strchr(cp+1,'_')-cp)-1;
      strncpy(aa_trans_version, cp+1, vers_len);
      aa_trans_version[vers_len] = '\0';
   }

   if (log_file)
      fprintf(log_file, "augmented atom transformation version = %s\n", aa_trans_version);

   trans_pairs = TypeAlloc(*ntransp,aa_pair);

   for (i=0; i<*ntransp; i++)
   {
      SearchChar(fp,'"'); fscanf(fp,"%999[^\"]",trans_string); SearchChar(fp,'"');
      if (!StringToAugmentedAtom(trans_string,&trans_pairs[i][0]))
      {
         fprintf(stderr, "ReadAAPairs: unsuccessful translation of from-pattern[%d]\n", i);
         success = FALSE;
      }
      
      SearchChar(fp,'"'); fscanf(fp,"%999[^\"]",trans_string); SearchChar(fp,'\n');
      if (!StringToAugmentedAtom(trans_string,&trans_pairs[i][1]))
      {
         fprintf(stderr, "ReadAAPairs: unsuccessful translation of to-pattern[%d]\n", i);
         success = FALSE;
      }
   }
   if (success) return (trans_pairs);

   MyFree((char *)trans_pairs); (*ntransp) = 0;
   return ((aa_pair *)NULL);
}

int IsUseAtom(int n, struct reaccs_molecule_t *mp)
/*
 * Tests if atom number n has a reaction center bond attached to it in
 * molecule *mp.
 */
{
   int i;
   
   for (i=0; i<mp->n_bonds; i++)
      if (mp->bond_array[i].reaction_mark != NONE       &&
          mp->bond_array[i].reaction_mark != NOT_CENTER &&
          mp->bond_array[i].reaction_mark != UNCHANGED  &&
          (mp->bond_array[i].atoms[0]==n || mp->bond_array[i].atoms[1]==n))
         return (TRUE);

   return (FALSE);
}

#define TRANSFORMED     1

struct reaccs_molecule_t *
AlternativeAATransformation(struct reaccs_molecule_t *mp,
                            int                       next_atom,
                            aa_pair                   tfm_table[],
                            int                       ntfm,
                            struct reaccs_molecule_t *old_list,
                            int                       use_old)
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
 * This function was designed to be used to genereate input to the CASP
 * program.
 */
{
   int i, j, nmatch, jmatch;
   unsigned int match[MAXNEIGHBOURS+1];
   struct reaccs_molecule_t *result, *mph;
   neighbourhood_t          *neighbour_array;
   int *atom_status, *bond_status;

   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);

   SetupNeighbourhood(mp, neighbour_array,mp->n_atoms);
tail_recursion_short_cut:
                                /* Search for atom to be transformed */
   for (i=next_atom; i<mp->n_atoms; i++)
   {
      nmatch = 0;
      for (j=0; j<ntfm; j++)
         if (IsUseAtom(i+1,mp) &&
             AAMatch(mp,i,match,&tfm_table[j][0],atom_status,neighbour_array))
         {
            jmatch = j;
            nmatch++;
         }
      if (nmatch > 0) break;
   }
   
   if (use_old)
   {
      if (i == mp->n_atoms)
      {
         mp->next = old_list;
	 MyFree((char *)neighbour_array);
         MyFree((char *)bond_status);
         MyFree((char *)atom_status);
         return (mp);
      }
      else
      {
         result = old_list;
         for (j=0; j<ntfm; j++)
            if (IsUseAtom(i+1,mp) &&
                AAMatch(mp,i,match,&tfm_table[j][0],atom_status,neighbour_array))
            {
               mph = CopyMolecule(mp);
               TransformAA(mph,i,match,
                          &tfm_table[j][0],
                          &tfm_table[j][1],
                           atom_status,
                           neighbour_array);
               mph->color = TRANSFORMED;
               if (!IsNULL(log_file))
                  fprintf(log_file,"atom #%d: '%s' -> '%s'\n",
                          i+1,
                          tfm_table[j][0].short_name,
                          tfm_table[j][1].short_name);
               result =
                  AlternativeAATransformation(mph,i+1, tfm_table, ntfm,
                                              result, use_old);
            }
         result = AlternativeAATransformation(mp, i+1, tfm_table, ntfm,
                                              result, use_old);
         MyFree((char *)bond_status);
         MyFree((char *)atom_status);
	 MyFree((char *)neighbour_array);
         return (result);
      }
   }
   else if (i == mp->n_atoms)        /* no additional transformations */
   {
      if (next_atom == 0)
      {
         FreeMolecule(mp);
         MyFree((char *)neighbour_array);
         MyFree((char *)bond_status);
         MyFree((char *)atom_status);
         return (old_list);
      }
      else
         mp->next = old_list;

      MyFree((char *)bond_status);
      MyFree((char *)atom_status);
      MyFree((char *)neighbour_array);
      return (mp);
   }
   else if (nmatch == 1) /* only one transformation -> tail recursion */
   {
      AAMatch(mp,i,match,&tfm_table[jmatch][0],atom_status,neighbour_array);
      TransformAA(mp,i,match,
                 &tfm_table[jmatch][0],
                 &tfm_table[jmatch][1],
                  atom_status,
                  neighbour_array);
      mp->color = TRANSFORMED;
      if (!IsNULL(log_file))
         fprintf(log_file,"atom #%d: '%s' -> '%s'\n",
                 i+1,
                 tfm_table[jmatch][0].short_name,
                 tfm_table[jmatch][1].short_name);
      next_atom = i+1;
      goto tail_recursion_short_cut;
   }
   else
   {
      result = old_list;
      for (j=0; j<ntfm; j++)
         if (IsUseAtom(i+1,mp) &&
             AAMatch(mp,i,match,&tfm_table[j][0],atom_status,neighbour_array))
         {
            mph = CopyMolecule(mp);
            TransformAA(mph,i,match,
                       &tfm_table[j][0],
                       &tfm_table[j][1],
                        atom_status,
                        neighbour_array);
            mph->color = TRANSFORMED;
            if (!IsNULL(log_file))
               fprintf(log_file,"atom #%d: '%s' -> '%s'\n",
                       i+1,
                       tfm_table[j][0].short_name,
                       tfm_table[j][1].short_name);
            result =
               AlternativeAATransformation(mph,i+1, tfm_table, ntfm,
                                           result, use_old);
         }
      FreeMolecule(mp);
      MyFree((char *)bond_status);
      MyFree((char *)atom_status);
      MyFree((char *)neighbour_array);
      return (result);
   }
}
