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
/*      File:           ssmatch.c                                       */
/*                                                                      */
/*      Purpose:        Implements substructure matching for reaccs     */
/*                      molecules. It also implements fingerprinting.   */
/*                                                                      */
//      History:        creation:           23-Feb-1990                 */
//                      added new hashing   14-May-2011                 */
//                                                                      */
/************************************************************************/

#include <string.h>

#include "hashcode.h"
#include "local.h"
#include "reaccs.h"
#include "smi2mol.h"
#include "string_counter.h"
#include "symbol_lists.h"
#include "utilities.h"

#include "perceive.h"

#include "ssmatch.h"

int CountMatches(ssmatch_t *smp)
/*
 * Counts the number of substructure matches listed in smp.
 */
{
   int result;

   for (result=0; smp; smp=smp->next) result++;

   return (result);
}

static ssmatch_t *free_ss_matches = (ssmatch_t *)NULL;

ssmatch_t *NewSSMatch(void)
/*
 * Creates an empty new substructure mapping structure. It either allocates
 * it or takes it from the free list.
 */
{
   register ssmatch_t *result;

   if (IsNULL(free_ss_matches))
   {
      result = TypeAlloc(1,ssmatch_t);
      if (IsNULL(result))
      {
         fprintf(stderr,"Out of memory when allocating ssmatch_t\n");
         exit (EXIT_FAILURE);
      }
   }
   else
   {
      result = free_ss_matches;
      free_ss_matches = result->next;
      result->next = (ssmatch_t *)NULL;
   }

   return (result);
}

void FreeSSMatch(ssmatch_t *matchp)
/*
 * Places the argument into the list of free substructure match records.
 */
{
   if (IsNULL(matchp))
      fprintf(stderr,"FreeSSMatch: NULL argument\n");
   else
   {
      matchp->next = free_ss_matches;
      free_ss_matches = matchp;
   }
}

void FreeSSMatchHeap(void)
/*
 * Frees all the ssmatch_t structures in the list free_ss_matches.
 */
{
   ssmatch_t *ssp;

   while (free_ss_matches)
   {
      ssp = free_ss_matches->next;
      MyFree((char *)free_ss_matches);
      free_ss_matches = ssp;
   }
}


static int AtomMatch(struct reaccs_atom_t *map,
                     struct reaccs_atom_t *ssap,
                     struct symbol_list_t *symbol_lists,
                     int                   atom)
/*
 * Tests if the substructure atom *ssap matches the structure atom map.
 * 'symbol_list' contains the list of query atom type sets for the substructure
 * and atom is the current substructure atom (this is necessary to pick the
 * correct list entry).
 */
{
   int result;

   if ((map->rsize_flags & ssap->rsize_flags)  !=  ssap->rsize_flags)
      return FALSE;
// if (ssap->query_H_count != NONE)
// fprintf(stderr, "atom %d has hydrogen count %d\n", atom, ssap->query_H_count);
// if (map->radical != NONE  &&  ssap->radical != NONE)
// fprintf(stderr, "trying to match mradical=%X with sradical=%X\n",
// map->radical, ssap->radical);
   if (0 == strcmp(ssap->atom_symbol,"L"))
   {
      while (!IsNULL(symbol_lists))
      {
         if (symbol_lists->atom == atom)
         {
            result =
               (ssap->charge == NONE  || map->charge == ssap->charge)     &&
               (ssap->radical == NONE  ||  map->radical == ssap->radical) &&
               (ssap->query_H_count == NONE   ||
                (ssap->query_H_count == ZERO_COUNT &&
                 map->query_H_count == ZERO_COUNT) ||
                ssap->query_H_count <= map->query_H_count);
            if (!result) return (result);
            result = AtomSymbolMatch(map->atom_symbol,symbol_lists->string);
            if (!symbol_lists->logic) result = !result;
            return (result);
         }
         symbol_lists = symbol_lists->next;
      }
      return (FALSE);
   }
   else
      return ((ssap->charge == NONE  ||
               map->charge == ssap->charge) &&
              (ssap->radical == NONE  ||
               map->radical == ssap->radical) &&
              (ssap->query_H_count == NONE   ||
              (ssap->query_H_count == ZERO_COUNT &&
               map->query_H_count == ZERO_COUNT) ||
               ssap->query_H_count <= map->query_H_count) &&
              AtomSymbolMatch(map->atom_symbol,ssap->atom_symbol));
}

// #define BOND_QUICK_TEST(mbt1, sbt2) (sbt2>AROMATIC  || (mbt1) == (sbt2))
#define BOND_QUICK_TEST(mbt1, sbt2) (sbt2>=AROMATIC  || (mbt1) == (sbt2))
#define BOND_QUICK_TEST_(bt1, bt2) (FALSE)

static int BondMatch(struct reaccs_bond_t *mbp,
                     struct reaccs_bond_t *ssbp)
/*
 * Tests if the substructure bond *ssbp matches the structure bond *mbp.
 */
{
   if (ssbp->topography != NONE  &&
       ssbp->topography != mbp->topography) return FALSE;
   if ((mbp->rsize_flags & ssbp->rsize_flags) != ssbp->rsize_flags)
      return FALSE;
   switch (ssbp->bond_type)
   {
      case  SINGLE:
      case  DOUBLE:
      case  TRIPLE:
// fprintf(stderr, "testing SDT qry %d bond (%d-%d) flags %X against mol %d bond (%d-%d) with flags %X. Result = %X\n",
// ssbp->bond_type, ssbp->atoms[0], ssbp->atoms[1], ssbp->bond_type_flags,
// mbp->bond_type, mbp->atoms[0], mbp->atoms[1], mbp->bond_type_flags, 
// (ssbp->bond_type_flags & mbp->bond_type_flags));
         return (ssbp->bond_type == mbp->bond_type);
      case  AROMATIC:
         if (ssbp->bond_type == mbp->bond_type) return TRUE;
// fprintf(stderr, "testing AROMATIC qry bond (%d-%d) flags %X against mol bond (%d-%d) with flags %X. Result = %X\n",
// ssbp->atoms[0], ssbp->atoms[1], ssbp->bond_type_flags,
// mbp->atoms[0], mbp->atoms[1], mbp->bond_type_flags, 
// (ssbp->bond_type_flags & mbp->bond_type_flags));
         return ((ssbp->bond_type_flags & mbp->bond_type_flags) != 0);
      case  SINGLE_DOUBLE:
         return (mbp->bond_type == SINGLE  ||
                 mbp->bond_type == DOUBLE);
      case  DOUBLE_AROMATIC:
// fprintf(stderr, "testing DOUBLE_AROMATIC qry bond (%d-%d) flags %X against mol bond (%d-%d) with flags %X. Result = %X\n",
// ssbp->atoms[0], ssbp->atoms[1], ssbp->bond_type_flags,
// mbp->atoms[0], mbp->atoms[1], mbp->bond_type_flags, 
// (ssbp->bond_type_flags & mbp->bond_type_flags));
         return (mbp->bond_type == AROMATIC  ||
                 mbp->bond_type == DOUBLE);
      case  ANY_BOND:
         return (TRUE);
      case  SINGLE_AROMATIC:
         return (mbp->bond_type == AROMATIC  ||
                 mbp->bond_type == SINGLE);
      default:
         fprintf(stderr,"BondMatch: Illegal bond type %d\n",ssbp->bond_type);
         return (FALSE);
   }
}

static
ssmatch_t *RecSSMatch(struct reaccs_molecule_t *mp,
                      struct reaccs_molecule_t *ssp,
                      int atom_match[],
                      int bond_match[],
                      int n, int nlimit,
                      int single_match,
                      neighbourhood_t nbp[],
                      neighbourhood_t ssnbp[],
                      ssmatch_t *old_matches,
                      uint64_t *match_signature)
/*
 * Recursively constructs the matchings of *ssp in *mp. It assumes that
 * the fields atom_match[0..n-1] contain already matched atoms and that
 * bond_match[0..ssp->n_bonds-1] contains matching bonds. Unassigned fields
 * in atom_match and bond_match are set to NOT_ASSIGNED (-1).
 * The next atom to be assigned a matching atom in *mp is atom number n+1
 * (atom number == 1+atom array index!!). nbp and ssnbp are used to speed up
 * access to neighbouring and the resulting matches are prepended to
 * old_matches.
 * Only 'nlimit' substructure atoms are considered. This feature is used
 * to test for CASP functional groups.
 * If 'single_match' is set TRUE, the recursion ends after the first match
 * has been found.
 *
 * if not NULL, match_sgnature[i] contains a bit pattern that is created by hashing all compatible ssp atoms into the range [0..30].
 * This can be used to quickly check if an ssp atom can match to the mp atom while doing recursive ABAS.
 */
{
   int *memory;
   int *new_atom_match, *new_bond_match;
   int *atom_candidate;
   unsigned int i, j, k;
   int false_match;
   ssmatch_t *result;
   int lig;
   int ncand, *candidates;
   struct reaccs_atom_t *ap1;
   int is_ok;
   struct reaccs_bond_t *sbp, *mbp;

   if (n == nlimit)     /* match found -> create matching */
   {
      // check if bond match is compatible with either pure standard or pure Daylight aromaticity
      int dy_mismatch = FALSE;
      int dy_required_match = FALSE;
      for (j=0; j<ssp->n_bonds; j++)
         if (bond_match[j] >= 0)
         {
            sbp = ssp->bond_array+j;
            mbp = mp->bond_array+bond_match[j];
// fprintf(stderr, "%d-(%o|%o)-%d => %d-(%o|%o)-%d\n",
// sbp->atoms[0], sbp->bond_type, sbp->bond_type_flags, sbp->atoms[1],
// mbp->atoms[0], mbp->bond_type, mbp->bond_type_flags, mbp->atoms[1]);
            // regular bond types take precedence
            if (sbp->bond_type < AROMATIC  &&  sbp->bond_type == mbp->bond_type) continue;
            if (sbp->bond_type_flags & BTF_AROMATIC)
            {
                if (mbp->bond_type_flags & BTF_DYAROMATIC)
		{
// fprintf(stderr,"dy_required_match = TRUE\n");
		    dy_required_match = TRUE;
		}
            }
            else
            {
                if (mbp->bond_type_flags & BTF_DYAROMATIC)
		{
// fprintf(stderr,"dy_mismatch = TRUE\n");
		    dy_mismatch = TRUE;
		}
            }
         }
// fprintf(stderr,"dy_required_match = %d, dy_mismatch = %d\n", dy_required_match, dy_mismatch);
      if (dy_required_match  && dy_mismatch) return old_matches;
      result = NewSSMatch();
      for (j=0; j<n; j++) result->match_atoms[j] = atom_match[j];
      result->n_match = n;
      result->next = old_matches;
      return (result);
   }

// ncand = 0;
// if (match_signature != NULL)
// {
//   /* take advantage of signature if it's defined */
//   for (i=0; i<mp->n_atoms; i++)
//      if (0 != (match_signature[i]&(1L<<(n%31))))
//         ncand++;
// }
// for (i=0; i<n; i++) fprintf(stderr, "%d->%d ", i+1, atom_match[i]+1); fprintf(stderr, "\tncand=%d\n", ncand);

   memory = TypeAlloc(ssp->n_atoms+ssp->n_bonds+mp->n_atoms+mp->n_atoms, int);
   new_atom_match = memory;
   new_bond_match = memory+ssp->n_atoms;
   atom_candidate = new_bond_match+ssp->n_bonds;
   candidates     = atom_candidate+mp->n_atoms;

   if (match_signature != NULL)
   {
      /* take advantage of signature if it's defined */
      for (i=0; i<mp->n_atoms; i++)
         if (0 != (match_signature[i]&(1L<<(n%31))))
	    atom_candidate[i] = TRUE;
   }
   else
      for (i=0; i<mp->n_atoms; i++) atom_candidate[i] = TRUE;
   /* already mapped atoms cannot be candidates */
   for (i=0; i<n; i++) atom_candidate[atom_match[i]] = FALSE;

   /* If atom n+1 of ss is connected to the already mapped part of the  */
   /* ss, then the candidate for the mapping must be a neighbour of the */
   /* maps of the neighbours of ss atom n+1.                            */
   ncand = 0;
   for (i=0; i<ssnbp[n].n_ligands; i++)
      if (ssnbp[n].atoms[i] < n)   /* link to already mapped part of ss */
      {                            /* fs atom mapped to ss ligand i */
         lig = atom_match[ssnbp[n].atoms[i]];
         for (j=0; j<nbp[lig].n_ligands; j++)
            if (atom_candidate[nbp[lig].atoms[j]])  /* unmapped fs lig. */
               candidates[ncand++] = nbp[lig].atoms[j];
      }
   if (ncand > 0)
   {
      for (i=0; i<mp->n_atoms; i++) atom_candidate[i] = FALSE;
      for (i=0; i<ncand; i++) atom_candidate[candidates[i]] = TRUE;
   }
// for (i=0; i<n; i++) fprintf(stderr, "%d->%d ", i+1, atom_match[i]+1); fprintf(stderr, "\tncand=%d\n", ncand);

                        /* look for candidate matches of ss atom #n+1 */
   for (i=0; i<mp->n_atoms; i++)
      if (atom_candidate[i] &&
          nbp[i].n_ligands >= ssnbp[n].n_ligands &&
          AtomMatch(&mp->atom_array[i],
                    &ssp->atom_array[n], ssp->symbol_lists, n+1))
   {
// fprintf(stderr, "matching SS '%s' atom %d to atom %d\n", ssp->atom_array[n].atom_symbol, n+1, i+1);
      ap1 = &ssp->atom_array[n];     /* Check substitution description */
      if (ap1->sub_desc != NONE)
      {
         if (ap1->sub_desc == SUB_AS_IS)
            is_ok = (nbp[i].n_ligands == ssnbp[n].n_ligands);
         else if (ap1->sub_desc == SUB_ZERO)
            is_ok = (nbp[i].n_ligands == 0);
         else if (ap1->sub_desc == SUB_MORE)
            is_ok = (nbp[i].n_ligands >= SUB_MORE);
         else
            is_ok = (nbp[i].n_ligands == ap1->sub_desc);
      }
      else
         is_ok = TRUE;
      if (!is_ok) continue;
// for (j=0; j<n; j++) fprintf(stderr," ");
// fprintf(stderr, "matching atom %d\n", i+1);

      for (j=0; j<ssp->n_atoms; j++)                /* initialize matching */
         new_atom_match[j] = atom_match[j];
      new_atom_match[n] = i;
      for (j=0; j<ssp->n_bonds; j++)
         new_bond_match[j] = bond_match[j];
                                           /* look if neighbourhood matches */
      for (j=0; j<ssnbp[n].n_ligands; j++)   /* for all ss ligands */
      {
         false_match = FALSE;
         if (new_atom_match[ssnbp[n].atoms[j]] != NOT_ASSIGNED)
         {           /* ligand atom already matched -> check bond value */
            if (new_bond_match[ssnbp[n].bonds[j]] != NOT_ASSIGNED)
               fprintf(stderr,"unexp. bond match\n"); /* should be new */
            else
            {
               for (k=0; k<nbp[i].n_ligands; k++) /* for all fs ligands */
                  if (nbp[i].atoms[k] == new_atom_match[ssnbp[n].atoms[j]])
                     if (BondMatch(&mp->bond_array[nbp[i].bonds[k]],
                                   &ssp->bond_array[ssnbp[n].bonds[j]]))
                     {
                        new_bond_match[ssnbp[n].bonds[j]] = nbp[i].bonds[k];
                        break;
                     }
                     else            /* fs bond does not match ss bond */
                        false_match = TRUE;
               if (k == nbp[i].n_ligands)         /* no corresponding bond */
                  false_match = TRUE;
            }
         }
         if (false_match) break;
      }
      if (j == ssnbp[n].n_ligands)   /* no ligand violation */
      {
         old_matches =
            RecSSMatch(mp,ssp,new_atom_match,new_bond_match,
                       n+1,nlimit,single_match,
                       nbp,ssnbp,old_matches, match_signature);
      }
      if (!IsNULL(old_matches) && single_match)
         break;
   }

   MyFree((char *)memory);
// MyFree((char *)candidates); MyFree((char *)atom_candidate);
// MyFree((char *)new_bond_match); MyFree((char *)new_atom_match);
   return (old_matches);
}

static int GoodFGEnvironment(ssmatch_t *match,
                             struct reaccs_molecule_t *mp,
                             neighbourhood_t nbp[],
                             struct reaccs_molecule_t *ssp,
                             neighbourhood_t  ssnbp[],
                             int nlimit)
/*
 * Checks if *match defines a mapping of the first nlimit atoms of
 * the substructure *ssp onto *mp which has an environment (alpha atoms)
 * compatible with the ligands allowed by the additional atoms in *ssp.
 * nbp[] and ssnbp[] define the neighbourhood structures of *mp and *ssp,
 * respectively, to speed-up access.
 */
{
   int i, j;
   int mindex, sindex;  /* index of current structure and substructure node */
   int mligand, sligand; /* indices of ligands of struc. and substruc. node */

   for (sindex=0; sindex<nlimit; sindex++)
   {
      mindex = match->match_atoms[sindex];
      for (i=0; i<nbp[mindex].n_ligands; i++)
      {
         mligand = nbp[mindex].atoms[i];
         for (j=0; j<ssnbp[sindex].n_ligands; j++)
         {
            sligand = ssnbp[sindex].atoms[j];
            if (sligand < nlimit && match->match_atoms[sligand] == mligand)
               break;                             /* check if already mapped */
            if (sligand < nlimit) continue;       /* mapped to another atom */

            if (AtomMatch(&mp->atom_array[mligand],
                          &ssp->atom_array[sligand],
                           ssp->symbol_lists, sligand+1) &&
                BondMatch(&mp->bond_array[nbp[mindex].bonds[i]],
                          &ssp->bond_array[ssnbp[sindex].bonds[j]]))
               break;
         }
         if (j == ssnbp[sindex].n_ligands) /* no matching ligand found */
            return (FALSE);
      }
   }

   return (TRUE);
}

void PrintMatch(FILE *fp, ssmatch_t *match)
/*
 * Prints a substructure match onto file *fp.
 */
{
   int i;

   for (i=0; i<match->n_match; i++)
      fprintf(fp,"%d->%d ",i+1,match->match_atoms[i]+1);
   fprintf(fp,"\n");
}

static
ssmatch_t *PruneFGMatches(ssmatch_t *old_matches,
                          struct reaccs_molecule_t *mp,
                          neighbourhood_t nbp[],
                          struct reaccs_molecule_t *ssp,
                          neighbourhood_t ssnbp[],
                          int nlimit)
/*
 * Checks the matches from the list 'old_matches' and prunes away those
 * which do not fulfill the 'optional neighbourhood' constraint defined
 * by the (ssp->n_atoms-nlimit) additional atoms in ssp.
 */
{
   ssmatch_t *result, *h;

   result = (ssmatch_t *)NULL;
   while (!IsNULL(old_matches))
   {
      if (GoodFGEnvironment(old_matches,mp,nbp,ssp,ssnbp,nlimit))
      {
         h = old_matches->next; old_matches->next = result;
         result = old_matches; old_matches = h;
      }
      else
      {
         h = old_matches; old_matches = h->next; FreeSSMatch(h);
      }
   }

   return (result);
}

ssmatch_t *FGMatch(struct reaccs_molecule_t *mp,
                   struct reaccs_molecule_t *ssp,
                   int nlimit)
/*
 * Returns the list of matchings of the substructure *ssp in structure *mp.
 * 'nlimit' atoms form the match, the remaining atoms define the allowed
 * neighbourhood of the matching atoms (optional neighbours).
 */
{
   ssmatch_t *old_matches;
   neighbourhood_t nbp[MAXATOMS], ssnbp[MAXATOMS];
   int i;
   int atom_match[MAXATOMS], bond_match[MAXBONDS];

   old_matches = (ssmatch_t *)NULL;

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      return old_matches;
   }

   if (ssp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",ssp->name);
      return old_matches;
   }

   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   SetupNeighbourhood(ssp,ssnbp,nlimit);
   for (i=0; i<ssp->n_atoms; i++) atom_match[i] = NOT_ASSIGNED;
   for (i=0; i<ssp->n_bonds; i++) bond_match[i] = NOT_ASSIGNED;

   for (i=0; i<mp->n_atoms; i++)
      if (nbp[i].n_ligands >= ssnbp[0].n_ligands &&
          AtomMatch(&mp->atom_array[i],
                    &ssp->atom_array[0], ssp->symbol_lists, 1))
      {
         atom_match[0] = i;
         old_matches = RecSSMatch(mp,ssp,atom_match,bond_match,
                                  1,nlimit,MULTIPLE_MATCHES,
                                  nbp,ssnbp,old_matches, NULL);
      }

   SetupNeighbourhood(ssp,ssnbp,ssp->n_atoms);
   return (PruneFGMatches(old_matches, mp, nbp, ssp, ssnbp, nlimit));
}

ssmatch_t *RSSMatch(struct reaccs_molecule_t *mp,
                    struct reaccs_molecule_t *ssp)
/*
 * Returns the list of matchings of the substructure *ssp in structure *mp.
 * Atoms with a zero mapping field must not have any connections other than
 * than to also mapped atoms.
 */
{
   ssmatch_t *old_matches, *matchph, *result;
   neighbourhood_t nbp[MAXATOMS], ssnbp[MAXATOMS];
   int i;
   int atom_match[MAXATOMS], bond_match[MAXBONDS];

   old_matches = (ssmatch_t *)NULL;

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      return old_matches;
   }

   if (ssp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",ssp->name);
      return old_matches;
   }

   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   SetupNeighbourhood(ssp,ssnbp,ssp->n_atoms);
   for (i=0; i<ssp->n_atoms; i++) atom_match[i] = NOT_ASSIGNED;
   for (i=0; i<ssp->n_bonds; i++) bond_match[i] = NOT_ASSIGNED;

   for (i=0; i<mp->n_atoms; i++)
      if (nbp[i].n_ligands >= ssnbp[0].n_ligands &&
          AtomMatch(&mp->atom_array[i],
                    &ssp->atom_array[0], ssp->symbol_lists, 1))
      {
         atom_match[0] = i;
         old_matches = RecSSMatch(mp,ssp,atom_match,bond_match,
                                  1,ssp->n_atoms,MULTIPLE_MATCHES,
                                  nbp,ssnbp,old_matches, NULL);
      }

   result = (ssmatch_t *)NULL;
   while (!IsNULL(old_matches))
   {
      for (i=0; i<ssp->n_atoms; i++)
         if (ssp->atom_array[i].mapping == NONE  &&
             ssnbp[i].n_ligands != nbp[old_matches->match_atoms[i]].n_ligands)
            break;
      if (i == ssp->n_atoms)      /* all matches OK */
      {
         matchph = old_matches->next; old_matches->next = result;
         result = old_matches; old_matches = matchph;
      }
      else                  /* bad match for atom i+1 */
      {
         matchph = old_matches; old_matches = old_matches->next;
         FreeSSMatch(matchph);
      }
   }

   return (result);
}

int Refine(struct reaccs_molecule_t *mp,  neighbourhood_t nbp[], unsigned long aa_mask[],
           struct reaccs_molecule_t *ssp, neighbourhood_t ssnbp[], unsigned long saa_mask[],
           int is_candidate[],
           uint64_t *match_signature)
/*
 * Checks for a possible matching of *ssp into *mp by a refinement
 * process. is_candidate[0..mp->n_atoms] is filled with TRUE or FALSE
 * if the corresponding atom index of *mp can match to index 0 of *ssp.
 *
 * match_sgnature[i] will contain a bit pattern that is created by hashing all compatible ssp atoms into the range [0..30].
 * This can be used to quickly check if an ssp atom can match to the mp atom while doing recursive ABAS.
 * The array is filled by Refine() if match_signature is != NULL, but it should be provided by the caller.
 *
 * If bot arrays are not NULL, the bit(s) in sa_class will be requed to match the mapping bits in aa_mask[]. This is needed
 * to implement shortcut matching.
 *
 * The function returns the number of candidates for this matching.
 */
{
   int i, j, k;
   int result, ok;
   struct reaccs_atom_t *ap1, *ap2;
   struct reaccs_bond_t *bp1, *bp2;
//   /* static */
   int tmp_atom1[MAXATOMS], ai1;
//   /* static */
   int tmp_atom2[MAXATOMS], ai2;
   int changed;
   int nmatch;
   int refine_ok;

   int *atom_memory, *bond_memory;
//   /* static */
   int *atom_cand[MAXATOMS], ncand[MAXATOMS];
//   /* static */
   int *bond_cand[MAXBONDS];

   int sscounts[AROMATIC+1], mcounts[AROMATIC+1];

static int nrefine=0, npass1=0, npass2=0, npass3=0;

int iloop;
   int nbond_cand;
   int quick_path;

   /* Double checking. The calling function should have done it, too */
   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      return 0;
   }
   if (ssp->n_atoms <= 0  ||  ssp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many (or too few) atoms in query '%s'\n",ssp->name);
      return 0;
   }

   /* Check compatibility of bond type counts */
   for (i=0; i<=AROMATIC; i++)
      sscounts[i] = mcounts[i] = 0;
   for (i=0, bp1=ssp->bond_array; i<ssp->n_bonds; i++, bp1++)
      if (bp1->bond_type <= AROMATIC) sscounts[bp1->bond_type]++;
   for (j=0, bp2=mp->bond_array; j<mp->n_bonds; j++, bp2++)
      if (bp2->bond_type <= AROMATIC) mcounts[bp2->bond_type]++;
   for (i=0; i<AROMATIC; i++)	// only check non-aromatic bond type counts
      if (sscounts[i] > mcounts[i]) return (0);

   /* Allocate memory as one big chunk to improve performance */
   atom_memory = TypeAlloc(mp->n_atoms*ssp->n_atoms, int);
   for (i=0; i<ssp->n_atoms; i++)
   {
      atom_cand[i] = atom_memory+mp->n_atoms*i;
      ncand[i] = 0;
   }

   bond_memory = TypeAlloc(mp->n_bonds*ssp->n_bonds, int);
   for (i=0; i<ssp->n_bonds; i++)
      bond_cand[i] = bond_memory+mp->n_bonds*i;

   nrefine++;
   refine_ok = FALSE;
                          /* Set up possible atom matches */
   for (i=0, ap1=ssp->atom_array; i<ssp->n_atoms; i++, ap1++)
   {
      if (ap1->atom_symbol[1] == '\0'   &&  ap1->atom_symbol[0] != 'A'  &&
          ap1->atom_symbol[0] != 'L'    &&  ap1->atom_symbol[0] != 'Q'  &&
          ap1->atom_symbol[0] != 'R'    &&  ap1->sub_desc == NONE       &&
	  ap1->query_H_count == NONE)	// query_H_count is special
      {
	 for (j=0, ap2=mp->atom_array; j<mp->n_atoms; j++, ap2++)
         {         /* Quick atom symbol test possible */
            if (ap1->atom_symbol[0] == ap2->atom_symbol[0]        &&
                ap2->atom_symbol[1] == '\0'                  &&
                (ap2->rsize_flags & ap1->rsize_flags)  ==  ap1->rsize_flags  &&
                nbp[j].n_ligands >= ssnbp[i].n_ligands)
            {
               atom_cand[i][j] = TRUE;
               ncand[i]++;
            }
            else
               atom_cand[i][j] = FALSE;
         }
         if (!IsNULL(saa_mask)  &&  !IsNULL(aa_mask))
             for (j=0; j<mp->n_atoms; j++)
                if ((saa_mask[i]&aa_mask[j]) != saa_mask[i])
                    atom_cand[i][j] = FALSE;
      }
      else
      {
         for (j=0, ap2=mp->atom_array; j<mp->n_atoms; j++, ap2++)
         {      /* need full atom symbol comparison */
            if (nbp[j].n_ligands >= ssnbp[i].n_ligands  &&
                (ap2->rsize_flags & ap1->rsize_flags)  ==  ap1->rsize_flags  &&
                AtomMatch(ap2, ap1, ssp->symbol_lists, i+1))
            {
               if (ap1->sub_desc != NONE)
               {
                  if (ap1->sub_desc == SUB_AS_IS)
                     atom_cand[i][j] = (nbp[j].n_ligands == ssnbp[i].n_ligands);
                  else if (ap1->sub_desc == SUB_ZERO)
                     atom_cand[i][j] = (nbp[j].n_ligands == 0);
                  else if (ap1->sub_desc == SUB_MORE)
                     atom_cand[i][j] = (nbp[j].n_ligands >= SUB_MORE);
                  else
                     atom_cand[i][j] = (nbp[j].n_ligands == ap1->sub_desc);
               }
               else
                  atom_cand[i][j] = TRUE;
               if (!IsNULL(saa_mask)  &&  !IsNULL(aa_mask))
                  if (atom_cand[i][j]  &&  (saa_mask[i]&aa_mask[j]) != saa_mask[i])
                  {
                      atom_cand[i][j] = FALSE;
                  }
               if (atom_cand[i][j])
               {
                  ncand[i]++;
               }
            }
            else
               atom_cand[i][j] = FALSE;
         }
      }
      if (ncand[i] == 0)
      {
          goto violation;
      }
   }
   npass1++;

                    /* Set up possible bond matches */
   for (i=0, bp1=ssp->bond_array; i<ssp->n_bonds; i++, bp1++)
   {
      nbond_cand = 0;
      if (bp1->bond_type < AROMATIC  &&		// aromatics are special
          bp1->rsize_flags  == 0     &&
          bp1->topography == NONE)
      {
         for (j=0, bp2=mp->bond_array; j<mp->n_bonds; j++, bp2++)
         {
            if (bp2->bond_type == bp1->bond_type)
            {
               bond_cand[i][j] = TRUE;
               nbond_cand++;
            }
            else
            {
               bond_cand[i][j] = FALSE;
            }
         }
      }
      else
      {
         for (j=0, bp2=mp->bond_array; j<mp->n_bonds; j++, bp2++)
            if (BOND_QUICK_TEST(bp2->bond_type, bp1->bond_type)  &&
                (bp2->rsize_flags & bp1->rsize_flags) == bp1->rsize_flags  &&
                BondMatch(bp2, bp1))
            {
               bond_cand[i][j] = TRUE;
               nbond_cand++;
            }
            else
	    {
               bond_cand[i][j] = FALSE;
            }
      }
      if (nbond_cand <= 0)
      {
         ai1 = bp1->atoms[0]-1; ai2 = bp1->atoms[1]-1;
         for (j=0; j<mp->n_atoms; j++)  /* clear away proven non-matches */
         {
            atom_cand[ai1][j] = FALSE;
            atom_cand[ai2][j] = FALSE;
         }
         goto violation;
      }
   }
   npass2++;

iloop=0;

loop_back:
   do             /* do refinement */
   {
iloop++;

// fprintf(stderr, "atom map %d\n", iloop);
// for (i=0; i<ssp->n_atoms; i++)
// {
//    for (j=0; j<mp->n_atoms; j++)
//       if (atom_cand[i][j]) fprintf(stderr, "1"); else fprintf(stderr, ".");
//    fprintf(stderr, "\n");
// }
// fprintf(stderr, "\nbond map\n");
// for (i=0; i<ssp->n_bonds; i++)
// {
//    for (j=0; j<mp->n_bonds; j++)
//       if (bond_cand[i][j]) fprintf(stderr, "1"); else fprintf(stderr, ".");
//    fprintf(stderr, "\n");
// }
// fprintf(stderr, "\n");

      changed = FALSE;
      for (i=0, bp1=ssp->bond_array; i<ssp->n_bonds; i++, bp1++)
      {
         for (j=0; j<mp->n_atoms; j++) tmp_atom1[j] = FALSE;
         for (j=0; j<mp->n_atoms; j++) tmp_atom2[j] = FALSE;
         ai1 = bp1->atoms[0]-1; ai2 = bp1->atoms[1]-1;
         for (j=0, bp2=mp->bond_array; j<mp->n_bonds; j++, bp2++)
            if (bond_cand[i][j])
            {
               ok = FALSE;
               if (atom_cand[ai1][bp2->atoms[0]-1]  && /* normal match */
                   atom_cand[ai2][bp2->atoms[1]-1])
               {
                  tmp_atom1[bp2->atoms[0]-1] = TRUE;
                  tmp_atom2[bp2->atoms[1]-1] = TRUE;
                  ok = TRUE;
               }
               if (atom_cand[ai1][bp2->atoms[1]-1]  && /* switched match */
                   atom_cand[ai2][bp2->atoms[0]-1])
               {
                  tmp_atom1[bp2->atoms[1]-1] = TRUE;
                  tmp_atom2[bp2->atoms[0]-1] = TRUE;
                  ok = TRUE;
               }
               if (!ok)                 /* bond candidate does not match */
               {
                  changed = TRUE;
                  bond_cand[i][j] = FALSE;
               }
            }
         for (j=0; j<mp->n_atoms; j++)  /* clear away proven non-matches */
            if (atom_cand[ai1][j] && !tmp_atom1[j])
            {
               atom_cand[ai1][j] = FALSE;
               ncand[ai1]--;
               if (ncand[ai1] == 0) goto violation;
               changed = TRUE;
            }
         for (j=0; j<mp->n_atoms; j++)
            if (atom_cand[ai2][j] && !tmp_atom2[j])
            {
               atom_cand[ai2][j] = FALSE;
               ncand[ai2]--;
               if (ncand[ai2] == 0) goto violation;
               changed = TRUE;
            }
      }
      for (i=0; i<ssp->n_atoms; i++)
         if (ncand[i] == 0)
         {
            changed = FALSE;
            goto violation;
         }

   } while (changed);
   npass3++;
// fprintf(stderr, "nrefine=%d, npass1=%d, npass2=%d, npass3=%d\n", nrefine, npass1, npass2, npass3);

   for (i=0; i<ssp->n_atoms; i++)
      if (ncand[i] == 1)          /* Exploit completely refined atoms */
      {
         for (j=0; j<mp->n_atoms; j++)
            if (atom_cand[i][j]) break;
         nmatch = j;              /* matching superstructure atom */
         for (j=0; j<ssp->n_atoms; j++)
            if (j != i  &&  atom_cand[j][nmatch])
            {
               ncand[j]--;
               if (ncand[j] == 0) goto violation;
               atom_cand[j][nmatch] = FALSE;
               changed = TRUE;
            }
      }
   if (changed) goto loop_back;
   refine_ok = TRUE;

violation:
   result = 0;
   for (i=0; i<mp->n_atoms; i++)
   {
      if (match_signature != NULL)
      {
         match_signature[i] = 0;
         for (j=0; j<ssp->n_atoms; j++)
            if (refine_ok  &&  atom_cand[j][i])
	       match_signature[i] |= ((uint64_t)1L<<(j%31));
      }
      if (refine_ok && atom_cand[0][i])
      {
         result++;
         is_candidate[i] = TRUE;
      }
      else
         is_candidate[i] = FALSE;
// if (match_signature != NULL) fprintf(stderr, "match_signature[%d] = \t%09lX\n", i, match_signature[i]);
   }

   for (i=0; i<ssp->n_atoms; i++)
      if (ncand[i] == 0)
      {
         result = 0;
         break;
      }

   if (result != 0)
      for (i=0; i<ssp->n_atoms; i++)   /* count neighbour candidates */
      {
         for (j=0; j<mp->n_atoms; j++)
            tmp_atom1[j] = FALSE;
         for (j=0; j<ssnbp[i].n_ligands; j++)
            for (k=0; k<mp->n_atoms; k++)
               if (atom_cand[ssnbp[i].atoms[j]][k])
                  tmp_atom1[k] = TRUE;
         for (j=0, nmatch=0; j<mp->n_atoms; j++)
            if (tmp_atom1[j]) nmatch++;
               if (nmatch < ssnbp[i].n_ligands)       /* cardinality violation */
               {
                  result = 0;
                  break;
               }
      }

//   if (FALSE && result)
//      for (i=0; i<ssp->n_bonds; i++)
//      {
//	 fprintf(stderr, "(%d-%d): ", ssp->bond_array[i].atoms[0], ssp->bond_array[i].atoms[1]);
//	 for (j=0; j<mp->n_bonds; j++)
//	    if (bond_cand[i][j])
//	       fprintf(stderr, "(%d-%d) ", mp->bond_array[j].atoms[0], mp->bond_array[j].atoms[1]);
//	 fprintf(stderr, "\n");
//      }
//
   MyFree((char *)bond_memory);
   MyFree((char *)atom_memory);

   return (result);
}

void PerceiveAromaticity(struct reaccs_molecule_t *mp, neighbourhood_t *nbp)
/*
 * Perceive standard 4n+2 alternating cycle aromaticity and store it with the
 * bond_type_flag member for matching.
 */
{
   int i;
   struct reaccs_bond_t *bp;
   int *h_counts;
   int *bond_types;

   PerceiveAromaticBonds(mp);
   h_counts = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++) h_counts[i] = mp->atom_array[i].query_H_count;
   bond_types = TypeAlloc(mp->n_bonds, int);
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      bond_types[i] = bp->bond_type;
      if (bp->bond_type == SINGLE) bp->bond_type_flags |= BTF_SINGLE;
      if (bp->bond_type == DOUBLE) bp->bond_type_flags |= BTF_DOUBLE;
      if (bp->bond_type == TRIPLE) bp->bond_type_flags |= BTF_TRIPLE;
      if (bp->bond_type == AROMATIC) bp->bond_type_flags |= BTF_AROMATIC;
      if (bp->bond_type == SINGLE_DOUBLE)
      {
	     bp->bond_type_flags |= BTF_SINGLE;
	     bp->bond_type_flags |= BTF_DOUBLE;
      }
      if (bp->bond_type == SINGLE_AROMATIC)
      {
	     bp->bond_type_flags |= BTF_SINGLE;
	     bp->bond_type_flags |= BTF_AROMATIC;
      }
      if (bp->bond_type == DOUBLE_AROMATIC)
      {
	     bp->bond_type_flags |= BTF_AROMATIC;
	     bp->bond_type_flags |= BTF_DOUBLE;
      }
   }
   for (i=0; i<mp->n_atoms; i++) mp->atom_array[i].query_H_count = h_counts[i];
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      bp->bond_type = bond_types[i];
   }
   PerceiveDYAromaticity(mp, nbp);
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->bond_type == AROMATIC  &&
          (bp->bond_type_flags & BTF_AROMATIC) == 0)
          bp->bond_type_flags |= BTF_DYAROMATIC;
      if (bp->bond_type == AROMATIC) bp->bond_type_flags |= BTF_AROMATIC;
      bp->bond_type = bond_types[i];
// fprintf(stderr, "(%d-%d) has flags %X\n", bp->atoms[0], bp->atoms[1], bp->bond_type_flags);
   }
   for (i=0; i<mp->n_atoms; i++) mp->atom_array[i].query_H_count = h_counts[i];
   MyFree((char *)h_counts);
   MyFree((char *)bond_types);
}

int refine_success = FALSE;

ssmatch_t *SSMatchExt(struct reaccs_molecule_t *mp,
                      unsigned long aa_mask[],      // perceived set if augnmented atom classes for shortcut matching or NULL
                      struct reaccs_molecule_t *ssp,
                      unsigned long saa_mask[],     // set of required augmented atom classes for the SS atom or NULL
                      int single_match)
/*
 * Returns the list of matchings of the substructure *ssp in structure *mp.
 * 'single_match' nedds to be TRUE if only the existence of a substructure
 * is to be checked.
 */
{
   ssmatch_t *old_matches;
   int i,j;
   int ncand;

   neighbourhood_t *nbp, *ssnbp;
   int *atom_status, *bond_status;
   struct reaccs_bond_t *bp;
   // storage used to fix fields modified by perception
   int hcounts[MAXATOMS];
   int bond_orders[MAXBONDS];
   int topographies[MAXBONDS];

// /* static */
   int atom_match[MAXATOMS], bond_match[MAXBONDS];
// /* static */
   int is_candidate[MAXATOMS];

   uint64_t *match_signature;
   match_signature = TypeAlloc(mp->n_atoms+1, uint64_t);

   old_matches = (ssmatch_t *)NULL;

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      return old_matches;
   }
   if (mp->n_bonds > MAXBONDS)
   {
      fprintf(stderr,"Too many bonds in molecule '%s'\n",mp->name);
      return old_matches;
   }

   if (ssp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",ssp->name);
      return old_matches;
   }

refine_success = FALSE;
   nbp   = TypeAlloc(mp->n_atoms, neighbourhood_t);
   ssnbp = TypeAlloc(ssp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   SetupNeighbourhood(ssp,ssnbp,ssp->n_atoms);

   for (j=0; j<mp->n_atoms; j++)
   {
       hcounts[j] = mp->atom_array[j].query_H_count;
   }
   for (j=0; j<mp->n_bonds; j++)
   {
      bond_orders[j] = mp->bond_array[j].bond_type;
      topographies[j] = mp->bond_array[j].topography;
   }
   PerceiveAromaticity(mp,nbp);
   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bond_status[i] > 0) bp->topography = RING;
      else                    bp->topography = CHAIN;
   MyFree((char *)atom_status);
   MyFree((char *)bond_status);

   PerceiveAromaticity(ssp,ssnbp);
   atom_status = TypeAlloc(ssp->n_atoms, int);
   bond_status = TypeAlloc(ssp->n_bonds, int);
   RingState(ssp, atom_status, bond_status);
   for (i=0, bp=ssp->bond_array; i<ssp->n_bonds; i++, bp++)
      if (bond_status[i] > 0) bp->topography = RING;
   MyFree((char *)atom_status);
   MyFree((char *)bond_status);

// if (mp->atom_array[0].x == 0) mp = LayoutMolecule(mp);
// PrintREACCSMolecule(stderr, mp, "MOL");
// PrintREACCSMolecule(stderr, ssp, "QRY");

   ncand = Refine(mp, nbp, aa_mask, ssp, ssnbp, saa_mask, is_candidate, match_signature);
// fprintf(stderr, "ncand = %d\n", ncand);

   if (ncand == 0)
   {
      MyFree((char *)ssnbp);
      MyFree((char *)nbp);
      MyFree((char *)match_signature);
      for (j=0; j<mp->n_bonds; j++)
      {
         mp->bond_array[j].bond_type = bond_orders[j];
         mp->bond_array[j].topography = topographies[j];
      }
      for (j=0; j<mp->n_atoms; j++)
         // mp->atom_array[j].query_H_count = NONE;
         mp->atom_array[j].query_H_count = hcounts[j];
      return (old_matches);
   }

refine_success = TRUE;

   for (i=0; i<ssp->n_atoms; i++) atom_match[i] = NOT_ASSIGNED;
   for (i=0; i<ssp->n_bonds; i++) bond_match[i] = NOT_ASSIGNED;

   for (i=0; i<mp->n_atoms; i++)
      if (is_candidate[i] &&
          nbp[i].n_ligands >= ssnbp[0].n_ligands &&
          AtomMatch(&mp->atom_array[i],
                    &ssp->atom_array[0], ssp->symbol_lists, 1))
      {
// fprintf(stderr, "matching SS '%s' atom 1 to atom %d\n", ssp->atom_array[0].atom_symbol, i+1);
         atom_match[0] = i;
         old_matches = RecSSMatch(mp,ssp,atom_match,bond_match,
                                  1,ssp->n_atoms,single_match,
                                  nbp,ssnbp,old_matches,
                                  match_signature);
         if (!IsNULL(old_matches) && single_match)
         {
            MyFree((char *)ssnbp);
            MyFree((char *)nbp);
            MyFree((char *)match_signature);
            for (j=0; j<mp->n_bonds; j++)
            {
               mp->bond_array[j].bond_type = bond_orders[j];
               mp->bond_array[j].topography = topographies[j];
            }
            for (j=0; j<mp->n_atoms; j++)
               // mp->atom_array[j].query_H_count = NONE;
               mp->atom_array[j].query_H_count = hcounts[j];
            return (old_matches);
         }
      }

   MyFree((char *)ssnbp);
   MyFree((char *)nbp);
   MyFree((char *)match_signature);
   for (i=0; i<mp->n_bonds; i++)
   {
      mp->bond_array[i].bond_type = bond_orders[i];
      mp->bond_array[i].topography = topographies[i];
   }
   for (j=0; j<mp->n_atoms; j++)
      // mp->atom_array[j].query_H_count = NONE;
      mp->atom_array[j].query_H_count = hcounts[j];
//   /* if (!old_matches) fprintf(stderr, "#"); */
   return (old_matches);
}

ssmatch_t *SSMatch(struct reaccs_molecule_t *mp,
                   struct reaccs_molecule_t *ssp,
                   int single_match)
{
    return SSMatchExt(mp, (unsigned long *)NULL, ssp, (unsigned long *)NULL, single_match);
}

void GuessHCountsFromSubstitution(struct reaccs_molecule_t *mp,
                                  neighbourhood_t nbp[])
/*
 * Uses substitution count query options to guess the number of
 * required hydrogens on oxigen and nitrogen atoms. This can speed
 * up screening for substructure matches.
 */
{
   int i, j, nsingle, ndouble, ntriple, naromatic, nother, nexplicit;
   struct reaccs_atom_t *ap;

   if (!mp) return;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      /* ignore if already defined */
      if (ap->query_H_count != NONE) continue;
      /* Only easy cases */
      if (ap->charge != 0) continue;
      if (ap->radical != 0) continue;
      nsingle = ndouble = naromatic= ntriple = nother = 0;
      nexplicit = 0;
      for (j=0; j<nbp[i].n_ligands; j++)
         if (AtomSymbolMatch(mp->atom_array[nbp[i].atoms[j]].atom_symbol, "H,D,T"))
            nexplicit++;
      /* Collect bond orders of ligands */
      for (j=0; j<nbp[i].n_ligands; j++)
         switch (mp->bond_array[nbp[i].bonds[j]].bond_type)
         {
            case SINGLE: nsingle++; break;
            case DOUBLE: ndouble++; break;
            case AROMATIC: naromatic++; break;
            case TRIPLE: ntriple++; break;
            default: nother++;
         }
      /* Ignore if there are query bonds */
      if (nother > 0) continue;
      /* for now, only handle AS_IS case */
      if (0 == strcmp("C", ap->atom_symbol))
      {
         if (naromatic == 2)
         {
            nsingle++; ndouble++; naromatic=0;
         }
         if (naromatic > 0) continue;
         if (ap->sub_desc == NONE)      /* degree from explicit hydrogens */
         {
            if (nsingle+2*ndouble+3*ntriple == 4)
            {
// fprintf(stderr, "%d: s=%d, d=%d, t=%d, sum=%d\n",
// i+1, nsingle, ndouble, ntriple, nsingle+2*ndouble+3*ntriple);
               ap->sub_desc = SUB_AS_IS;
            }
         }
         else if (ap->sub_desc == SUB_AS_IS)    /* H_count from degree */
         {
            ap->query_H_count =
               nexplicit+ZERO_COUNT+4-nsingle-2*ndouble-3*ntriple;
         }
      }
      else if (0 == strcmp("O", ap->atom_symbol))
      {
         if (ntriple+naromatic != 0) continue;
         if (ap->sub_desc == SUB_AS_IS)    /* H_count from degree */
            ap->query_H_count = nexplicit+ZERO_COUNT+2-nsingle-2*ndouble;
      }
      else if (0 == strcmp("N", ap->atom_symbol))
      {
         if (naromatic != 0) continue;
         if (ap->sub_desc == SUB_AS_IS)    /* H_count from degree */
            ap->query_H_count = nexplicit+ZERO_COUNT+3-nsingle-2*ndouble;
      }
// fprintf(stderr,"UTILITIES:\t'%s' atom %d has h_count=%d\n", ap->atom_symbol, i+1, ap->query_H_count);
   }
}

#define RING_PATTERN_SEED      11
#define RING_PATH_SEED         13
#define ATOM_SYMBOL_PATH_SEED  17
#define ATOM_CLASS_PATH_SEED   23
#define ATOM_COUNT_SEED        31
#define AUGMENTED_ATOM_SEED    37
#define HCOUNT_PATH_SEED       41
#define HCOUNT_CLASS_PATH_SEED 43
#define HCOUNT_PAIR_SEED       47
#define BOND_PATH_SEED         53
#define AUGMENTED_BOND_SEED    61
#define RING_SIZE_SEED         67
#define DEGREE_PATH_SEED       71
#define CLASS_SPIDER_SEED      79
#define RING_CLOSURE_SEED      101

#define NON_SSS_SEED          179

// /* macro to convert the current seed value into the 'incremented' one */
// #define OLD_NEXT_SEED(seed, increment) \
//         (((seed)&0x000FFFFL)*(unsigned)(increment)+(((seed)&0xFFF0000L)>>16)+7) 
// #define OLD_ADD_BIT(counts, ncounts, seed) (counts[(seed)%ncounts]++)

/* new macro to convert the current seed value into the 'incremented' one */
#define NEXT_SEED(seed, increment) next_hash(seed, increment)
#define ADD_BIT(counts, ncounts, seed) (counts[hash_position(seed,ncounts)]++)
#define ADD_BIT_COUNT(counts, ncounts, seed, count) (counts[hash_position(seed,ncounts)]+=count)
#define SET_BIT(bytes, nbytes, seed) \
        (bytes[((seed)/8)%nbytes] |= 0xFF&(1<<((seed)%8)))

/* Flags to be used to control recursive processing */
#define PROCESS_RING_CLOSURES 0x0001
#define PROCESS_CHAINS        0x0002
#define FORCED_HETERO_END     0x0004
#define IGNORE_PATH_SYMBOL    0x0008
#define IGNORE_TERM_SYMBOL    0x0010
#define FORCED_RING_PATH      0x0020
#define STOP_AT_HEAVY_ATOM    0x0040
#define DEBUG_PATH            0x0100

#define ANY_COLOR       113

#define CSP3      19
#define HETERO    23
#define GENERIC   (-1)

#define SPECIAL_RING (0xFC&~(1<<6))

static
void SetPathLengthFlags(struct reaccs_molecule_t *mp,
                        int touched_indices[],
                        int start_index, int path_length, int current_index,
                        int max_size,
                        int **length_matrix,
                        neighbourhood_t nbp[],
			int exclude_atom,
			char *prefix)
/*
 * Recursively traces the neighbouring of an atom (start_index+1)
 * collecting the path_lengths as bit flags in length_matrix[][].
 */
{
   int i, ai;

   for (i=0; i<nbp[current_index].n_ligands; i++)
   {
      if (path_length+1 > max_size)             /* don't go too far */
         continue;
      ai = nbp[current_index].atoms[i];
      if (ai+1 == exclude_atom) continue;
      if (touched_indices[ai]) continue;   /* don't walk backwards */
      if (mp->atom_array[ai].color == 0) continue;
      touched_indices[ai] = 1;  /* updating */
      length_matrix[start_index][ai] |= 1<<(path_length+1);
      SetPathLengthFlags(mp, touched_indices,
                         start_index, path_length+1, ai,
                         max_size, length_matrix, nbp,
			 exclude_atom,
			 prefix);
      touched_indices[ai] = 0;  /* down-dating */
   }
}

static
void SpecialNeighboursRec(struct reaccs_molecule_t *mp,
                          int touched_indices[],
                          int path_length,
                          int current_index,
                          int max_size,
                          /* count of sp3 carbons with >= 3 C neighbours */
                          int csp3[],
                          /* count of hetero atoms */
                          int hetero[],
                          neighbourhood_t nbp[],
			  int exclude_atom,
			  char *prefix)
/*
 * Recursively traces the neighbouring of an atom collecting
 * counts of special atoms at certain graph distances.
 *
 * exclude_atom terminates neighbourhood search paths.
 */
{
   int i, ai;

   for (i=0; i<nbp[current_index].n_ligands; i++)
   {
      if (path_length+1 > max_size)             /* don't go too far */
         continue;
      ai = nbp[current_index].atoms[i];
      if (ai+1 == exclude_atom) continue;
      if (touched_indices[ai]) continue;   /* don't walk backwards */
      if (mp->atom_array[ai].color <= 0) continue;
      if (mp->atom_array[ai].color == CSP3) csp3[path_length]++;
      else if (mp->atom_array[ai].color == HETERO) hetero[path_length]++;
      /* only carbon-connected SPIDERS (beta atoms are always included) */
      if (mp->atom_array[ai].color == HETERO  &&  path_length > 1)
         continue;
      touched_indices[ai] = 1;  /* updating */
      SpecialNeighboursRec(mp, touched_indices,
                           path_length+1, ai,
                           max_size, csp3, hetero, nbp,
			   exclude_atom,
			   prefix);
      touched_indices[ai] = 0;  /* down-dating */
   }
}

void printTouched(struct reaccs_molecule_t *mp,
		  int *touched_indices,
		  int ai,
		  long seed,
		  char *message)
/*
 * Debugging utility
 */
{
   int i;
   fprintf(stderr, "%s: ", message);
   for (i=0; i<mp->n_atoms; i++)
      if (touched_indices[i] != 0)
      {
	 fprintf(stderr, "%s|%d(%d)",
			 mp->atom_array[i].atom_symbol,
			 mp->atom_array[i].color,
			 i+1);
      }
   fprintf(stderr, " %ld\n", seed);
   for (i=0; i<mp->n_bonds; i++)
      if (touched_indices[mp->bond_array[i].atoms[0]-1] != 0  &&
          touched_indices[mp->bond_array[i].atoms[1]-1] != 0)
      {
	 fprintf(stderr, "(%d-%d-%d)",
			 mp->bond_array[i].atoms[0],
			 mp->bond_array[i].bond_type,
			 mp->bond_array[i].atoms[1]);
      }
   fprintf(stderr, " %ld\n", seed);
}

int SetPathBitsRec(struct reaccs_molecule_t *mp,
                   neighbourhood_t *nbp,
                   int *fp_counts, int ncounts,
                   uint64_t seed,
                   int *touched_indices,
                   int nbonds,
                   int minbonds, int maxbonds,
                   int sprout_index,
                   int first_index,
                   int last_index,
                   int flags,
		   int exclude_atom,
		   char *prefix)
/*
 * Recursively enumerates the paths through *mp. The next sprouting
 * step is done on the atom (sprout_index+1). seed represents the
 * hash_code processing value up to and including this atom. touched_indices[i]
 * is > 0 if atom (i+1) is already included in the parent path.
 *
 * Paths through exclude_atom are terminated.
 *
 * If prefix is != NULL, then it points to a string buffer to be used to construct the string representation of the
 * current paths. This is not yet implemented!
 */
{
   int result;
   int i, ai, bi, acolor, bcolor;
   uint64_t old_seed;
   struct reaccs_atom_t *ap;

   result = 0;
   old_seed = seed;
   if (nbonds > maxbonds) return (result);
   for (i=0; i<nbp[sprout_index].n_ligands; i++)
   {
      ai = nbp[sprout_index].atoms[i];
      if (ai == last_index) continue;
      if (ai+1 == exclude_atom) continue;
      bi = nbp[sprout_index].bonds[i];
      ap = mp->atom_array+ai;
      if ((flags&FORCED_RING_PATH) && mp->bond_array[bi].rsize_flags == 0)
	 continue;
      if (ap->color >= 18  &&  ap->color != ANY_COLOR  &&  0 != (flags & STOP_AT_HEAVY_ATOM)) continue;
      if (touched_indices[ai] > 0)      /* ring closure */
      {
         if (0 == (flags & PROCESS_RING_CLOSURES)) continue;
         if (touched_indices[ai] > 1) continue; // not just a plain ring
         bcolor = mp->bond_array[bi].color;
         if (bcolor == 0) continue;
         acolor = mp->atom_array[ai].color;
         if (acolor ==  0  &&  0 == (flags & IGNORE_PATH_SYMBOL)) continue;
         old_seed = seed;
         seed = NEXT_SEED(seed, bcolor*16);
         if (0 == (flags & IGNORE_PATH_SYMBOL))
	 {
	    seed = NEXT_SEED(seed, acolor);
	 }
         seed = NEXT_SEED(seed, touched_indices[ai]);  // codes how far this closure points back
         if (flags & DEBUG_PATH) printTouched(mp, touched_indices, ai, seed, "ring");
         ADD_BIT(fp_counts, ncounts, seed);
         seed = NEXT_SEED(seed, RING_CLOSURE_SEED*(nbonds-touched_indices[ai]));  // codes ring size
         ADD_BIT(fp_counts, ncounts, seed);
         result++;
         seed = old_seed;
      }
      else                              /* normal path */
      {
         bcolor = mp->bond_array[bi].color;
         if (bcolor == 0) continue;
         acolor = mp->atom_array[ai].color;
         if (acolor == 0  &&  0 == (flags & IGNORE_PATH_SYMBOL)) continue;
         /* saving and updating */
         touched_indices[ai] = nbonds+1;
         old_seed = seed;
         seed = NEXT_SEED(seed, bcolor*16);
         if (0 == (flags & IGNORE_PATH_SYMBOL))
	 {
            seed = NEXT_SEED(seed, acolor);
	 }
         if (nbonds >= minbonds  &&
             (acolor > 0  ||  0 != (flags & IGNORE_TERM_SYMBOL))  &&
             0 != (flags & PROCESS_CHAINS))
         {
            if (0 == (flags & FORCED_HETERO_END)  ||
                0 != (flags & IGNORE_TERM_SYMBOL) ||
                acolor != 6)
            {
               if (acolor > 1)         /* don't hit hydrogens! */
               {
// fprintf(stderr,"Setting chain bit for %d bonds, term color = %d\n", nbonds, acolor);
                  if (0 != (flags & IGNORE_TERM_SYMBOL))
		  {
                     ADD_BIT(fp_counts, ncounts, seed);
		  }
                  else
		  {
                     ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed,17*acolor));
		  }
                  result++;
               }
            }
         }
         if (nbonds+1 <= maxbonds)       /* continue recursion */
         {
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     nbonds+1,
                                     minbonds, maxbonds,
                                     ai,
                                     first_index, sprout_index,
				     flags, exclude_atom,
				     prefix);
         }

         /* restoring and down dating */
         touched_indices[ai] = 0;
         seed = old_seed;
      }
   }

   return result;
}

#define HETERO_FLAG     0x0100
#define RING_SUBST_FLAG 0x0200
#define QUART_FLAG      0x0400
#define CSP3_FLAG       0x0800
#define RS_SPECIAL_FLAG 0x1000
#define TYPE_MASK       0x00FF
#define C_FLAG          0x0001
#define O_FLAG          0x0002
#define N_FLAG          0x0003
#define S_FLAG          0x0004
#define P_FLAG          0x0005
#define X_FLAG          0x0006

int SetFeatureBits(struct reaccs_molecule_t *mp,
                   int *fp_counts,
                   int ncounts,
                   int start_flags,
                   int end_flags,
                   int path_min,
                   int path_max,
                   int use_counts,
                   int use_atom_types,
                   int **length_matrix,
                   uint64_t start_seed,
		   int exclude_atom,
		   char *prefix)
{
   int result = 0;
   int coli, colj;
   int i, j, k;
   uint64_t seed_i, seed;
   int *counts;

   counts = TypeAlloc(ncounts*4, int);	/* allocate tmp array for counts */
   for (i=0; i<mp->n_atoms; i++)
   {
      if (i+1 == exclude_atom) continue;
      coli = mp->atom_array[i].color;
      if (0 == (coli&start_flags)) continue;
      if (use_atom_types)
      {
	 if (0 == (coli&TYPE_MASK)) continue;	// ignore generic atoms
         seed_i = NEXT_SEED(start_seed, coli&TYPE_MASK);
      }
      else                seed_i = start_seed;
      for (j=0; j<mp->n_atoms; j++)
      {
	 if (j+1 == exclude_atom) continue;
         colj = mp->atom_array[j].color;
         if (0 == (colj&end_flags)) continue;
         if (use_atom_types)
         {
	    if (0 == (colj&TYPE_MASK)) continue;	// ignore generic atoms
	    seed = NEXT_SEED(seed_i, colj&TYPE_MASK);
         }
         else                seed = seed_i;
         for (k=path_min; k<=path_max; k++)
            if ((1<<k)&length_matrix[i][j])
            {
	       /* count the features */
	       counts[(k*19+seed)%(ncounts*4)]++;
               // ADD_BIT(fp_counts, ncounts, k*19+seed);
               result++;
// fprintf(stderr, "setting bit for (%d/%d)-%d-(%d/%d)\n",i+1,coli,k,j+1,colj);
            }
      }
   }
   /* Set the bits */
   for (i=0; i<ncounts*4; i++)
   {
      if (counts[i] > 0)
      {
	 ADD_BIT(fp_counts, ncounts, i);
	 result++;
      }
      if (use_counts  &&  counts[i] > 1)
      {
         ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,19), counts[i]);
	 result++;
      }
      if (use_counts  &&  counts[i] > 2)
      {
         ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,29), counts[i]);
	 result++;
      }
      if (use_counts  &&  counts[i] > 4)
      {
         ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,59), counts[i]);
	 result++;
      }
      if (use_counts  &&  counts[i] > 8)
      {
         ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i,79), counts[i]);
	 result++;
      }
   }
   MyFree((char *)counts);
   return (result);
}

int SetRingSizePairBits(struct reaccs_molecule_t *mp,
                        int *fp_counts,
                        int ncounts,
			int start_flags,
			int end_flags,
                        int path_min,
                        int path_max,
                        int **length_matrix,
                        uint64_t start_seed,
			int exclude_atom,
			char *prefix)
{
   int result = 0;
   int coli, colj;
   int rsizei, rsizej;
   int i, j, k;
   int ri, rj;
   uint64_t seed_i, seed;
   for (i=0; i<mp->n_atoms; i++)
   {
      if (i+1 == exclude_atom) continue;
      coli = mp->atom_array[i].color;
      if (0 == (coli&start_flags)) continue;
      rsizei = mp->atom_array[i].rsize_flags;
      if (0 == rsizei) continue;
      seed_i = start_seed;
      for (j=0; j<mp->n_atoms; j++)
      {
	 if (j+1 == exclude_atom) continue;
         colj = mp->atom_array[j].color;
         if (0 == (colj&end_flags)) continue;
         rsizej = mp->atom_array[j].rsize_flags;
         if (0 == rsizej) continue;
         seed = seed_i;
         for (k=path_min; k<=path_max; k++)
            if ((1<<k)&length_matrix[i][j])
            {
	       for (ri=3; ri<15; ri++)
	       {
		  if (0 == (rsizei&(1<<ri))) continue;
	          for (rj=3; rj<15; rj++)
	          {
		     if (0 == (rsizej&(1<<rj))) continue;
		     ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed, (k/2)*119), 37*ri*rj));
		     result++;
// fprintf(stderr, "setting bit for (%d/%d)-%d-(%d/%d)\n",i+1,ri,k,j+1,rj);
		  }
	       }
            }
      }
   }
   return (result);
}

/*
 * Sets the bits in fingerprint[0..nbytes-1] that correspond to the paths
 * with up to some defined size. which_bits is a set of flags that defines
 * which algorithms are triggered. as_query must be set to TRUE if the
 * fingerprint is to be computed for a query, which is mostly for hydrogen
 * counts.
 */
int SetFingerprintBits(struct reaccs_molecule_t *mp,
                       char *fingerprint, int nbytes,
                       int which_bits,
                       int as_query,
                       int fpflags)
{
   int* fp_counts;
   int result;
   int i;

   if (IsNULL(mp)) return 0;

// fprintf(stderr, "start SetFingerprintBits()\n");
   fp_counts = TypeAlloc(nbytes*8, int);
   result = SetFingerprintCountsWithFocus(mp, fp_counts, nbytes*8,
	                                  which_bits, as_query, fpflags, 0);
// fprintf(stderr, "middle SetFingerprintBits()\n");
   for (i=0; i<nbytes*8; i++)
      if (fp_counts[i] > 0) SET_BIT(fingerprint, nbytes, i);
   MyFree((char *)fp_counts);
   return result;
// fprintf(stderr, "done SetFingerprintBits()\n");
}

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
                                       int focus_atom)
{
   int* fp_counts, *fp_counts_excluded;
   int result;
   int i;

   if (IsNULL(mp)) return 0;

   fp_counts          = TypeAlloc(nbytes*8, int);
   fp_counts_excluded = TypeAlloc(nbytes*8, int);
   result = SetFingerprintCountsWithFocus(mp, fp_counts, nbytes*8,
	                                      which_bits, as_query, fpflags, 0);
   result -= SetFingerprintCountsWithFocus(mp, fp_counts_excluded, nbytes*8,
	                                       which_bits, as_query, fpflags, focus_atom);
   for (i=0; i<nbytes*8; i++)
      if (fp_counts[i]-fp_counts_excluded[i] > 0)
          SET_BIT(fingerprint, nbytes, i);
   MyFree((char *)fp_counts);
   MyFree((char *)fp_counts_excluded);
   return result;
}

/*
 * Sets the counts in fp_counts[0..nbytes-1] that correspond to the paths
 * with up to some defined size. which_bits is a set of flags that defines
 * which algorithms are triggered. as_query must be set to TRUE if the
 * fp_counts are to be computed for a query, which is mostly for hydrogen
 * counts.
 *
 * exclude_atom is the number of the atom that must not be touched by each
 * fragment. exclude_atom == 0 means don't exclude any atom.
 *
 * Derived from SetFingerprintCountsWithFocus to which it defaults if sc == NULL.
 */
int CountFingerprintPatterns(struct reaccs_molecule_t *mp,
                             int *fp_counts, int ncounts,
                             int which_bits,
                             int as_query,
                             int fpflags,
                             int exclude_atom,
                             struct string_counter_t *sc)
{
   neighbourhood_t *nbp;
   uint64_t seed, old_seed;
   int *touched_indices; /* this array is up- and down-dated during recursion */
   int *old_bond_types;
   struct reaccs_atom_t *ap, *ap1, *ap2, *ap3;
   struct reaccs_bond_t *bp;
   int i, j, j1, j2, k, tmp;
   int i1, i2, i3;
   uint64_t prod, sum, sumi, prodi, sumj, prodj;
   int ai, ai1, ai2;
   int qq_count;
   int qc_count;
   int nrbonds;
   int nrare_atoms;
   int is_rare;
   int nqtmp, nmulti;
   int hash;
   int *H_count;
   int nbits;
   int result;  /* the number of paths enumerated */
   struct symbol_list_t *symbol_lists;
   char prefix[100];	// string buffer used to print/construct text representation of FP feature
   char pat_buf[400];
   char *tokp;
#define NCOUNT_HASH 128
#define NCOUNT_SEED_HASH 128*128
   int atom_type_count_hash[NCOUNT_HASH];
   int atom_type_count_seed_hash[NCOUNT_SEED_HASH];
   int *atom_status, *bond_status, *degree, *cdegree, *unsaturated, *nspecial;
#define MAX_SPIDER 7
   int csp3[MAX_SPIDER+1], hetero[MAX_SPIDER+1];
   int tmp1, tmp2;
   int ndouble, naromatic;
   int rscounts[15][15];
   int nringch2, nfusionch, nspiro, nfusionb;
   int *length_tmp;
   int *extcon, *extcon2;
   int **length_matrix;
   int flags;
   int jj;
   int changed;
// FORTIFY    long total_bytes_allocated;
// FORTIFY    long bytes_allocated;
// FORTIFY    total_bytes_allocated = fortifySet();

   if (IsNULL(mp)) return 0;

   strcpy(prefix, "NONE:");
   result = 0;
   if (0 == (fpflags & ACCUMULATE_BITS))
      for (i=0; i<ncounts; i++) fp_counts[i] = 0;
   nbp     = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);
   touched_indices = TypeAlloc(mp->n_atoms, int);

   old_bond_types = TypeAlloc(mp->n_bonds, int);
   for (i=0; i<mp->n_bonds; i++)
      old_bond_types[i] = mp->bond_array[i].bond_type;

   /* collect hydrogen counts */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   if (as_query)
   {
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
            H_count[bp->atoms[1]]++;
         else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
            H_count[bp->atoms[0]]++;
      }
      GuessHCountsFromSubstitution(mp, nbp);
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->query_H_count != NONE)
         {
            H_count[i+1] = ap->query_H_count-ZERO_COUNT;
// fprintf(stderr,"SSMATCH:\t'%s' atom %d has h_count=%d\n", ap->atom_symbol, i+1, ap->query_H_count);
         }
   }
   else
   {
      ComputeImplicitH(mp, H_count);
      /* Add the explicit hydrogens to the implicit counts */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
            H_count[bp->atoms[1]]++;
         else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
            H_count[bp->atoms[0]]++;
      }
   }

   /* Do some perception */
   if (fpflags & USE_DY_AROMATICITY)
   {
      PerceiveDYAromaticity(mp, nbp);
   }
   else
   {
      PerceiveAromaticBonds(mp);
   }
   degree = TypeAlloc(mp->n_atoms, int);
   cdegree = TypeAlloc(mp->n_atoms, int);
   nspecial = TypeAlloc(mp->n_atoms, int);
   unsaturated = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   atom_status = TypeAlloc(mp->n_atoms, int);
   RingState(mp, atom_status, bond_status);
   SetRingSizeFlags(mp, 14, nbp);

   nrare_atoms = 0;
   /* Set the color property to represent all different atom types */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      unsaturated[i] = FALSE;
      ap->color = StringToInt(periodic_table, ap->atom_symbol);
      if (ap->color <= 1) ap->color = 0;        /* ignore hydrogens */
      /* mark special atom types */
      if (ap->color > 115) ap->color = -1;
      if (0 == strcmp("A", ap->atom_symbol)) ap->color = -1;
      if (0 == strcmp("L", ap->atom_symbol))
      {
         ap->color = -1;
         symbol_lists = mp->symbol_lists;
         while (!IsNULL(symbol_lists))
         {
            if (symbol_lists->atom == i+1)
            {
               is_rare = TRUE;
               if (!symbol_lists->logic) is_rare = FALSE;
               is_rare = is_rare && !AtomSymbolMatch("C",symbol_lists->string);
	           is_rare = is_rare && !AtomSymbolMatch("H",symbol_lists->string);
	           is_rare = is_rare && !AtomSymbolMatch("O",symbol_lists->string);
	           is_rare = is_rare && !AtomSymbolMatch("N",symbol_lists->string);
	           is_rare = is_rare && !AtomSymbolMatch("Cl",symbol_lists->string);
	           is_rare = is_rare && !AtomSymbolMatch("F",symbol_lists->string);
	        }
            symbol_lists = symbol_lists->next;
         }
      }
      else
	     is_rare = ap->color > 0  &&
		           !AtomSymbolMatch(ap->atom_symbol,"C,H,O,N,S,P,Cl,F");
      if (is_rare)
	     if (exclude_atom != i+1  ||  exclude_atom <= 0) nrare_atoms++;
      // color shortcut atoms if any
      if ((which_bits & USE_SHORTCUT_LABELS)  &&  strcmp(ap->atom_symbol,"R") == 0  &&  ap->atext[0] != '\0')
      {
         ap->color = (0xFFFF00&hash_string(ap->atext)) | 119;
      }
   }
// if (nrare_atoms > 0) fprintf(stderr, "nrare_atoms = %d\n", nrare_atoms);

   ndouble = 0;
   naromatic = 0;
   nfusionb = 0;
   /* Set the color property to represent the different bond type classes */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->bond_type == SINGLE)        bp->color = 1;
      else if (bp->bond_type == DOUBLE)
      {
         bp->color = 2;
         if (bp->atoms[0] != exclude_atom  &&  bp->atoms[1] != exclude_atom)
	    ndouble++;
      }
      else if (bp->bond_type == TRIPLE)   bp->color = 3;
      else if (bp->bond_type == AROMATIC)
      {
         bp->color = 4;
         if (bp->atoms[0] != exclude_atom  &&  bp->atoms[1] != exclude_atom)
	    naromatic++;
      }
      else                                bp->color = 0;
      if (bp->color > 1)
      {
         unsaturated[bp->atoms[0]-1] = TRUE;
         unsaturated[bp->atoms[1]-1] = TRUE;
      }

      /* Count non-hydrogen degree */
      if (mp->atom_array[bp->atoms[0]-1].color != 0  &&
          mp->atom_array[bp->atoms[1]-1].color != 0)
      {
         degree[bp->atoms[0]-1]++;
         degree[bp->atoms[1]-1]++;
      }
      /* Count carbon degree */
      if (bp->bond_type == DOUBLE)
      {
	     nspecial[bp->atoms[0]-1]++;
	     nspecial[bp->atoms[1]-1]++;
      }
      else if (bp->bond_type == TRIPLE)
      {
	     nspecial[bp->atoms[0]-1]+=2;
	     nspecial[bp->atoms[1]-1]+=2;
      }
      if (bp->bond_type != SINGLE) continue;
      if (mp->atom_array[bp->atoms[0]-1].color != 0  &&
          mp->atom_array[bp->atoms[1]-1].color == 6)
	     cdegree[bp->atoms[0]-1]++;
      if (mp->atom_array[bp->atoms[1]-1].color != 0  &&
          mp->atom_array[bp->atoms[0]-1].color == 6)
	     cdegree[bp->atoms[1]-1]++;
      
      if (exclude_atom <= 0  ||
	      (bp->atoms[0] != exclude_atom  &&  bp->atoms[1] != exclude_atom))
	     if (atom_status[bp->atoms[0]-1] > 2  &&
	         atom_status[bp->atoms[1]-1] > 2)	// ring fusion
	        nfusionb++;
   }
   /* ignore special atom types for further processing */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color < 0) ap->color = 0;

// fprintf(stderr, "USE_ATOM_COUNT\n");
   if (which_bits & USE_ATOM_COUNT)
   {
      strcpy(prefix, "AC:");
      /* Collect hashed counts of atom types with hydrogen counts */
      for (i=0; i<NCOUNT_HASH; i++)
         atom_type_count_hash[i] = 0;
      for (i=0; i<NCOUNT_SEED_HASH; i++)
         atom_type_count_seed_hash[i] = 0;
      nringch2=0; nfusionch=0; nspiro=0;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
	     if (i+1 == exclude_atom) continue;
         if (ap->color == 0) continue;
         if (ap->color == 6)
         {
            /* count CH2 in ring */
            if (atom_status[i] > 0  &&  H_count[i+1] >= 2) nringch2++;
            /* count ring fusion CH atoms */
            if (atom_status[i] > 2  &&  H_count[i+1] >= 1) nfusionch++;
            if (atom_status[i] > 3) nspiro++;
            continue;   // carbon has a special count processing
         }
         hash = 0;
         seed = NEXT_SEED(317*ATOM_COUNT_SEED,507);
         seed = NEXT_SEED(seed,ap->color+17);
         atom_type_count_seed_hash[seed%NCOUNT_SEED_HASH]++;
         /* normal hetero only with hydrogen */
         if ((ap->color == 7  ||  ap->color == 8)  &&  H_count[i+1] <= 0)
            continue;
         hash = hash*7 + ap->color+13;
         atom_type_count_hash[hash%NCOUNT_HASH]++;
         seed = NEXT_SEED(seed,ap->color+13);
         atom_type_count_seed_hash[seed%NCOUNT_SEED_HASH]++;
         /* one more bit for rare types */
         if (ap->color != 7  &&  ap->color != 8)
         {
            hash = hash*7 + ap->color+2*13;
            atom_type_count_hash[hash%NCOUNT_HASH]++;
            seed = NEXT_SEED(seed,ap->color+2*13);
            atom_type_count_seed_hash[seed%NCOUNT_SEED_HASH]++;
         }
      }
      nbits = 0;
      /* Now, we set the corresponding bits */

      if (ncounts <= 2048)      // old atom count fingerprints
          for (i=0; i<NCOUNT_HASH; i++)
          {
             if (atom_type_count_hash[i] > 0)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*19,3), atom_type_count_hash[i]);
                nbits++;
             }
             if (atom_type_count_hash[i] > 1)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*23,5), atom_type_count_hash[i]);
                nbits++;
             }
             if (atom_type_count_hash[i] > 2)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*29,7), atom_type_count_hash[i]);
                nbits++;
             }
             if (atom_type_count_hash[i] > 4)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*31,11), atom_type_count_hash[i]);
                nbits++;
             }
             if (atom_type_count_hash[i] > 8)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*37,13), atom_type_count_hash[i]);
                nbits++;
             }
             if (atom_type_count_hash[i] > 16)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*41,17), atom_type_count_hash[i]);
                nbits++;
             }
          }
      else      // new atom count fingerprints
          for (i=0; i<NCOUNT_SEED_HASH; i++)
          {
             if (atom_type_count_seed_hash[i] > 0)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*19,3), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 1)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*23,5), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 2)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*29,7), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 3)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*29,71), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 4)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*31,11), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 6)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*31,113), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 8)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*37,13), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 12)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*37,133), atom_type_count_seed_hash[i]);
                nbits++;
             }
             if (atom_type_count_seed_hash[i] > 16)
             {
                ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(i*41,17), atom_type_count_seed_hash[i]);
                nbits++;
             }
          }

      if (naromatic >= 10) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 341), naromatic); nbits++;}
      if (naromatic >= 14) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 441), naromatic); nbits++;}
      if (naromatic >= 18) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 541), naromatic); nbits++;}
      if (naromatic >= 22) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 641), naromatic); nbits++;}
      if (naromatic >= 26) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 741), naromatic); nbits++;}

      if (ndouble >= 3) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 371), naromatic); nbits++;}
      if (ndouble >= 5) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 471), naromatic); nbits++;}
      if (ndouble >= 8) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 571), naromatic); nbits++;}

      if (nringch2 >= 6)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 411), nringch2); nbits++;}
      if (nringch2 >= 12) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 511), nringch2); nbits++;}
      if (nringch2 >= 22) {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 611), nringch2); nbits++;}

      if (nfusionch >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 421), nfusionch); nbits++;}
      if (nfusionch >= 4)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 521), nfusionch); nbits++;}

      // Note spiro atoms
      if (nspiro >= 1)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 322), nspiro); nbits++;}
      if (nspiro >= 1)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 422), nspiro); nbits++;}
      if (nspiro >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 522), nspiro); nbits++;}
      if (nspiro >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 622), nspiro); nbits++;}

      if (nfusionb >= 1)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 425), nfusionb); nbits++;}
      if (nfusionb >= 2)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 525), nfusionb); nbits++;}
      if (nfusionb >= 3)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 625), nfusionb); nbits++;}
      if (nfusionb >= 5)  {ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(ATOM_COUNT_SEED, 825), nfusionb); nbits++;}
      if (nrare_atoms > 0)
      {
         seed = NEXT_SEED(317*ATOM_COUNT_SEED,5);
         ADD_BIT_COUNT(fp_counts, ncounts, seed, nrare_atoms); nbits++;
         ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,101), nrare_atoms); nbits++;
         ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,103), nrare_atoms); nbits++;
         if (nrare_atoms > 1)
         {
            ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,203), nrare_atoms);
            nbits++;
         }
         if (nrare_atoms > 3)
         {
            ADD_BIT_COUNT(fp_counts, ncounts, NEXT_SEED(seed,211), nrare_atoms);
            nbits++;
         }
      }
      result += nbits;
   }

// fprintf(stderr, "USE_ATOM_SYMBOL_PATH\n");
   if (which_bits & USE_ATOM_SYMBOL_PATH)
   {
      strcpy(prefix, "ASP:");
      seed = ATOM_SYMBOL_PATH_SEED;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
	 if (i+1 == exclude_atom) continue;
         if (ap->color <= 0) continue;
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
	 /* Ignore common atom types */
         if (ap->color != 6  &&  ap->color != 7  &&  ap->color != 8)
	 {
	    ADD_BIT(fp_counts, ncounts, seed);
	    result++;
	 }
         /* fingerprint two atom pairs if substitution count is defined */
if (1) // [TODO] 1
         if (ap->color !=6  &&
             (!as_query                 ||
              ap->sub_desc == SUB_AS_IS ||
              (ap->sub_desc != NONE     &&
               ap->sub_desc != SUB_MORE &&
               ap->sub_desc == degree[i]+SUB_ONE-1)))
         {
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed+12347, touched_indices,
                                     1,
                                     1, 1,  /* path length 1 to 1 */
                                     i, 0, -1,
                                     FORCED_HETERO_END |
                                     PROCESS_CHAINS,
                                     exclude_atom,
				     prefix);
         }
         /* 2 bond paths for not very common starts */
if (1) // [TODO] 2
         if (ap->color >= 10  ||
             (ap->color == 7  &&  degree[i] > 0)  ||
             (ap->color == 8  &&  degree[i] > 1)  ||
             (ap->color == 6  &&  degree[i] > 2  &&  atom_status[i] > 0))
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     1, 2,  /* path length 1 to 2 */
                                     i, 0, -1,
                                     STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
                                     FORCED_HETERO_END |
                                     PROCESS_CHAINS,
                                     exclude_atom,
				     prefix);
         /* Add more paths starting at special atoms */
if (1) // [TODO] 3
         if (ap->color > 6)
         {
	       result += SetPathBitsRec(mp, nbp,
				     fp_counts, ncounts,
				     217*ap->color+seed, touched_indices,
				     1,
				     3, 4,  /* path length 3 to 4 */
				     i, 0, -1,
				     IGNORE_PATH_SYMBOL |
				     PROCESS_RING_CLOSURES |
                                     STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
				     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);

               if (ap->color > 10  &&  ap->color <= 18)     // only third row of periodic table
                    result += SetPathBitsRec(mp, nbp,
					     fp_counts, ncounts,
					     17+seed, touched_indices,
					     1,
					     5, 7,
					     i, 0, -1,
					     FORCED_HETERO_END |
					     IGNORE_PATH_SYMBOL |
                                             STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
					     // PROCESS_RING_CLOSURES |      // rings containing metal cause too many bits
					     PROCESS_CHAINS,
					     exclude_atom,
					     prefix);

         }
if (1) // [TODO] 4
         if ((ap->color == 7 ||  ap->color == 8) &&      // only do this for common hetero elements
             cdegree[i]>2)
         {
            seed = NEXT_SEED(seed, ap->color*23);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     2, 6,  /* path length 1 to 4 */
                                     i, 0, -1,
                                     STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
                                     IGNORE_TERM_SYMBOL |
                                     // IGNORE_PATH_SYMBOL |
                                     PROCESS_CHAINS,
                                     exclude_atom,
                                     prefix);
         }

         seed = old_seed;
         /* Add another set of bits for paths starting at more crowded atoms */
if (1) // [TODO] 5
         if (degree[i] >= 4  &&
             ap->color != 5  &&         // don't do it for boron
             ap->color < 18)            // don't do it for transition metals
         {
            seed = NEXT_SEED(seed, ap->color);
            if (1)
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     NEXT_SEED(seed,3*107), touched_indices,
                                     1,
                                     2, 4,  /* path length 2 to 4 */
                                     i, 0, -1,
                                     STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
                                     IGNORE_TERM_SYMBOL |     /* gedeck.mol */
                                     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);
         }

         seed = old_seed;
         /* Add another set of bits for paths starting at spiro atoms */
if (1) // [TODO] 6
         if (atom_status[i] >= 4 &&
             ap->color != 5      &&         // don't do it for boron
             ap->color < 18)                // don't do it for transition metals
         {
// fprintf(stderr, "spiro atom %s(%d) found\n", ap->atom_symbol, i+1);
            seed = NEXT_SEED(seed, ap->color+55);
// if (1)
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     NEXT_SEED(seed,3*109), touched_indices,
                                     1,
                                     2, 4,  /* path length 1 to 4 */
                                     i, 0, -1,
                                     STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
                                     IGNORE_TERM_SYMBOL |
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_RING_CLOSURES |
                                     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);
         }

         seed = old_seed;
         /* Add bits for paths starting with rare bond orders */
         if (1)
         for (j=0; j<nbp[i].n_ligands; j++)
         {
            bp = &mp->bond_array[nbp[i].bonds[j]];
            ai = nbp[i].atoms[j];
            if (ai+1 == exclude_atom) continue;
            if (bp->color == 0) continue;
            if (bp->bond_type != DOUBLE  &&
                bp->bond_type != TRIPLE)
               continue;
            seed = old_seed;
            seed = NEXT_SEED(seed, ap->color);
            seed = NEXT_SEED(seed, bp->color*613);
            touched_indices[ai] = 1;    /* updating */
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     2,
                                     5, 5,
                                     ai, 0, i,
                                     STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
                                     IGNORE_PATH_SYMBOL |
                                     // IGNORE_TERM_SYMBOL |
                                     PROCESS_RING_CLOSURES |
                                     (0*PROCESS_CHAINS),
                                     exclude_atom,
                                     prefix);
            touched_indices[ai] = 0;    /* down-dating */
         }

         seed = old_seed;
         touched_indices[i] = 0; /* down-dating */
      }
   }

// fprintf(stderr, "USE_AUGMENTED_ATOM\n");
   if (which_bits & USE_AUGMENTED_ATOM)
   {
      strcpy(prefix, "AA:");
      /* Set bits for all triples of atoms connected to a common atom */
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
	 if (i+1 == exclude_atom) continue;

	 /* Set bit for atoms with more than one double or a triple bond */
	 if (nspecial[i] >= 2)
	 {
	    ADD_BIT(fp_counts, ncounts, NEXT_SEED(ap->color*AUGMENTED_ATOM_SEED,101));
	    // ADD_BIT(fp_counts, ncounts, NEXT_SEED(ap->color*AUGMENTED_ATOM_SEED,301));
	    result += 1;
	 }
         old_seed = seed;
         // Add some bits for hydrogen counted or hetero central atoms with
         // hetero neighbours
         if ((H_count[i+1] > 0 || ap->color != 6)  && degree[i] >= 2)
            for (i1=0; i1<nbp[i].n_ligands; i1++)
	    {
	       if (mp->atom_array[nbp[i].atoms[i1]].color == 0) continue;
	       if (nbp[i].atoms[i1]+1 == exclude_atom) continue;
	       if (mp->bond_array[nbp[i].bonds[i1]].color == 0) continue;
               for (i2=i1+1; i2<nbp[i].n_ligands; i2++)
               {
                  seed = NEXT_SEED(AUGMENTED_ATOM_SEED, 97);
                  seed = NEXT_SEED(seed, ap->color);
                  if (mp->atom_array[nbp[i].atoms[i2]].color == 0) continue;
		  if (nbp[i].atoms[i2]+1 == exclude_atom) continue;
                  if (mp->bond_array[nbp[i].bonds[i2]].color == 0) continue;
                  if (mp->atom_array[nbp[i].atoms[i1]].color == 6  &&
                      mp->atom_array[nbp[i].atoms[i2]].color == 6) continue;
                  sum = 0;
                  sum += mp->atom_array[nbp[i].atoms[i1]].color *
                         mp->bond_array[nbp[i].bonds[i1]].color;
                  sum += mp->atom_array[nbp[i].atoms[i2]].color *
                         mp->bond_array[nbp[i].bonds[i2]].color;
                  prod = 1;
                  prod *= mp->atom_array[nbp[i].atoms[i1]].color; prod &= 0xFFF;
                  prod *= mp->atom_array[nbp[i].atoms[i2]].color; prod &= 0xFFF;
                  seed = NEXT_SEED(seed, sum);
                  seed = NEXT_SEED(seed, prod);
                  ADD_BIT(fp_counts, ncounts, seed);
                  // seed = NEXT_SEED(seed, (sum*prod)&0xFFF);
                  // ADD_BIT(fp_counts, ncounts, seed);
		  result += 1;
		  nmulti = 0;
                  if (mp->bond_array[nbp[i].bonds[i1]].color >= 2) nmulti++;
                  if (mp->bond_array[nbp[i].bonds[i2]].color >= 2) nmulti++;
                  // add bit for hetero-substituted unsaturated hetero atom if degree>1 is well-defined => catch e.g. nitroso vs. nitro
                  if (ap->color != 6  &&  nmulti > 0  &&
                      (!as_query                 ||
                       ap->sub_desc == SUB_AS_IS ||
                       (ap->sub_desc != NONE     &&
                        ap->sub_desc != SUB_MORE &&
                        ap->sub_desc == degree[i]+SUB_ONE-1)))
                  {
                     ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, nmulti+173*degree[i]));
                     result += 1;
                  }
               }
            }

         if (degree[i] <= 2) continue;
         for (i1=0; i1<nbp[i].n_ligands; i1++)
         {
            if (nbp[i].atoms[i1]+1 == exclude_atom) continue;
            for (i2=i1+1; i2<nbp[i].n_ligands; i2++)
	    {
	       if (nbp[i].atoms[i2]+1 == exclude_atom) continue;
               for (i3=i2+1; i3<nbp[i].n_ligands; i3++)
               {
		  if (nbp[i].atoms[i3]+1 == exclude_atom) continue;
                  seed = AUGMENTED_ATOM_SEED;
                  seed = NEXT_SEED(seed, ap->color);
                  if (mp->atom_array[nbp[i].atoms[i1]].color == 0) continue;
                  if (mp->atom_array[nbp[i].atoms[i2]].color == 0) continue;
                  if (mp->atom_array[nbp[i].atoms[i3]].color == 0) continue;
                  if (mp->bond_array[nbp[i].bonds[i1]].color == 0) continue;
                  if (mp->bond_array[nbp[i].bonds[i2]].color == 0) continue;
                  if (mp->bond_array[nbp[i].bonds[i3]].color == 0) continue;
		  /* count hetero neighbours */
		  nqtmp = 0;
                  if (mp->atom_array[nbp[i].atoms[i1]].color != 6) nqtmp++;
                  if (mp->atom_array[nbp[i].atoms[i2]].color != 6) nqtmp++;
                  if (mp->atom_array[nbp[i].atoms[i3]].color != 6) nqtmp++;

		  /* make sure to add some bits for really odd ones */
		  nmulti = 0;
                  if (mp->bond_array[nbp[i].bonds[i1]].color == 2) nmulti++;
                  if (mp->bond_array[nbp[i].bonds[i1]].color == 3) nmulti++;
                  if (mp->bond_array[nbp[i].bonds[i2]].color == 2) nmulti++;
                  if (mp->bond_array[nbp[i].bonds[i2]].color == 3) nmulti++;
                  if (mp->bond_array[nbp[i].bonds[i3]].color == 2) nmulti++;
                  if (mp->bond_array[nbp[i].bonds[i3]].color == 3) nmulti++;

                  sum = 0;
                  sum += mp->atom_array[nbp[i].atoms[i1]].color *
                         mp->bond_array[nbp[i].bonds[i1]].color;
                  sum += mp->atom_array[nbp[i].atoms[i2]].color *
                         mp->bond_array[nbp[i].bonds[i2]].color;
                  sum += mp->atom_array[nbp[i].atoms[i3]].color *
                         mp->bond_array[nbp[i].bonds[i3]].color;
                  prod = 1;
                  prod *= mp->atom_array[nbp[i].atoms[i1]].color; prod &= 0xFFF;
                  prod *= mp->atom_array[nbp[i].atoms[i2]].color; prod &= 0xFFF;
                  prod *= mp->atom_array[nbp[i].atoms[i3]].color; prod &= 0xFFF;
                  seed = NEXT_SEED(seed, sum);
                  seed = NEXT_SEED(seed, prod);
                  if ((nqtmp > 2  ||  nmulti >= 2)  &&  ap->color == 6)
                  {
                     ADD_BIT(fp_counts, ncounts, seed);
                     ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 73));
		     result += 2;
                  }
                  seed = NEXT_SEED(seed, (sum*prod)&0xFFF);
                  ADD_BIT(fp_counts, ncounts, seed);
		  result++;
		  if (nmulti >= 2  ||  ap->color > 6)   // make sure R-NO2 is covered
                  {
                      ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 53)); result++;
                  }
               }
	    }
	 }
      }
   }

// fprintf(stderr, "USE_AUGMENTED_BOND\n");
   if (which_bits & USE_AUGMENTED_BOND)
   {
      strcpy(prefix, "AB:");
      /* Set bits for all bonds with both end-degrees > 2 */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (bp->atoms[0] == exclude_atom) continue;
         if (bp->atoms[1] == exclude_atom) continue;
         ai1 = bp->atoms[0]-1;
         ai2 = bp->atoms[1]-1;
         if (degree[ai1] <= 2) continue;
         if (degree[ai2] <= 2) continue;
         for (i1=0; i1<nbp[ai1].n_ligands; i1++)
            for (i2=i1+1; i2<nbp[ai1].n_ligands; i2++)
            {
               /* don't reuse current bond */
               if (nbp[ai1].bonds[i1] == i) continue;
               if (nbp[ai1].bonds[i2] == i) continue;
               if (nbp[ai1].atoms[i1]+1 == exclude_atom) continue;
               if (nbp[ai1].atoms[i2]+1 == exclude_atom) continue;
               /* don't use masked atoms/bonds */
               if (mp->atom_array[nbp[ai1].atoms[i1]].color==0) continue;
               if (mp->atom_array[nbp[ai1].atoms[i2]].color==0) continue;
               if (mp->bond_array[nbp[ai1].bonds[i1]].color==0) continue;
               if (mp->bond_array[nbp[ai1].bonds[i2]].color==0) continue;
               sumi = 0;
               sumi += mp->atom_array[nbp[ai1].atoms[i1]].color*
                       mp->bond_array[nbp[ai1].bonds[i1]].color;
               sumi += mp->atom_array[nbp[ai1].atoms[i2]].color*
                       mp->bond_array[nbp[ai1].bonds[i2]].color;
               prodi = 1;
               prodi *= mp->atom_array[nbp[ai1].atoms[i1]].color+
                        mp->bond_array[nbp[ai1].bonds[i1]].color;
               prodi *= mp->atom_array[nbp[ai1].atoms[i2]].color+
                        mp->bond_array[nbp[ai1].bonds[i2]].color;
               sumi  &= 0x0FFF;
               prodi &= 0x0FFF;
               for (j1=0; j1<nbp[ai2].n_ligands; j1++)
                  for (j2=j1+1; j2<nbp[ai2].n_ligands; j2++)
                  {
                     /* don't reuse current bond */
                     if (nbp[ai2].bonds[j1] == i) continue;
                     if (nbp[ai2].bonds[j2] == i) continue;
                     if (nbp[ai2].atoms[j1]+1 == exclude_atom) continue;
                     if (nbp[ai2].atoms[j2]+1 == exclude_atom) continue;
                     /* don't use masked atoms/bonds */
                     if (mp->atom_array[nbp[ai2].atoms[j1]].color==0) continue;
                     if (mp->atom_array[nbp[ai2].atoms[j2]].color==0) continue;
                     if (mp->bond_array[nbp[ai2].bonds[j1]].color==0) continue;
                     if (mp->bond_array[nbp[ai2].bonds[j2]].color==0) continue;
                     sumj = 0;
                     sumj += mp->atom_array[nbp[ai2].atoms[j1]].color*
                             mp->bond_array[nbp[ai2].bonds[j1]].color;
                     sumj += mp->atom_array[nbp[ai2].atoms[j2]].color*
                             mp->bond_array[nbp[ai2].bonds[j2]].color;
                     prodj = 1;
                     prodj *= mp->atom_array[nbp[ai2].atoms[j1]].color+
                              mp->bond_array[nbp[ai2].bonds[j1]].color;
                     prodj *= mp->atom_array[nbp[ai2].atoms[j2]].color+
                              mp->bond_array[nbp[ai2].bonds[j2]].color;
                     sumj  &= 0x0FFF;
                     prodj &= 0x0FFF;
                     seed = AUGMENTED_BOND_SEED;
                     seed = NEXT_SEED(seed, bp->color);
                     seed = NEXT_SEED(seed, prodi+prodj);
                     seed = NEXT_SEED(seed, sumi*sumj);
                     ADD_BIT(fp_counts, ncounts, seed);
		     result++;
                  }
            }
      }
   }

// fprintf(stderr, "USE_HCOUNT_PAIR\n");
   if (1*which_bits & USE_HCOUNT_PAIR)
   {
      strcpy(prefix, "HPR:");
      /* generate bits for hydrogen counted described bonds */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (H_count[bp->atoms[0]] == 0  &&  H_count[bp->atoms[1]] == 0)
            continue;
         if (bp->atoms[0] == exclude_atom) continue;
         if (bp->atoms[1] == exclude_atom) continue;
         /* Don't consider CC single bonds */
         if (mp->atom_array[bp->atoms[0]-1].color == 6  &&
             mp->atom_array[bp->atoms[1]-1].color == 6  &&
             bp->color == 1) continue;
         /* Don't consider explicit AH bonds */
         if (mp->atom_array[bp->atoms[0]-1].color == 0) continue;
         if (mp->atom_array[bp->atoms[1]-1].color == 0) continue;
         for (j1=0; j1<=H_count[bp->atoms[0]]; j1++)
            for (j2=0; j2<=H_count[bp->atoms[1]]; j2++)
            {
               if (j1+j2 == 0) continue;        /* at least one hydrogen */
               // unsaturation triggers bit like a hydrogen
               if (!unsaturated[bp->atoms[0]-1]  &&  !unsaturated[bp->atoms[1]-1]  &&
                   j1*j2 == 0  && bp->color <= 1) continue;
               if (j1+j2 > 3) continue;         /* at most 3 hydrogens */
               // seed = HCOUNT_PAIR_SEED + 53*(j1+j2) + 7*(j1+1)*(j2+1);
               seed = NEXT_SEED(HCOUNT_PAIR_SEED, 7*(j1+1)*(j2+1));
               seed = NEXT_SEED(seed, 53*(j1+j2));
               seed = NEXT_SEED(seed, bp->color);
               seed = NEXT_SEED(seed, mp->atom_array[bp->atoms[0]-1].color +
                                      mp->atom_array[bp->atoms[1]-1].color);
               seed = NEXT_SEED(seed, mp->atom_array[bp->atoms[0]-1].color *
                                      mp->atom_array[bp->atoms[1]-1].color);
               ADD_BIT(fp_counts, ncounts, seed);
	       result++;
               if (bp->color > 1)
               {
                  seed = NEXT_SEED(seed, 83);
                  ADD_BIT(fp_counts, ncounts, seed);
		  result++;
                  /* Add more bits for really special pairs */
                  if (j1+j2 == 1  &&
                      (mp->atom_array[bp->atoms[0]-1].color != 6  ||
                       mp->atom_array[bp->atoms[1]-1].color != 6) &&
                      bp->color <= 3)
                  {
                     seed = NEXT_SEED(seed, 91);
                     ADD_BIT(fp_counts, ncounts, seed);
                     seed = NEXT_SEED(seed, 97);
                     ADD_BIT(fp_counts, ncounts, seed);
                     seed = NEXT_SEED(seed, 103);
                     ADD_BIT(fp_counts, ncounts, seed);
		     result += 3;
                  }
               }
            }
      }
   }

// fprintf(stderr, "USE_HCOUNT_PATH\n");
   if (which_bits & USE_HCOUNT_PATH)
   {
      strcpy(prefix, "HPH:");
      /* generate a short path for each atom that has a hydrogen */
      seed = HCOUNT_PATH_SEED;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
         if (i+1 == exclude_atom) continue;
         if (H_count[i+1] == 0) continue;
         /* don't consider non methyl carbon atoms */
         if (ap->color == 6  &&  H_count[i+1] < 3  && degree[i] < 3) continue;

         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
         if (ap->color != 6)
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     2, 5,  /* path length 2 to 4 */ // EVG
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
// IGNORE_TERM_SYMBOL |
                                     PROCESS_CHAINS,
                                     exclude_atom,
                                     prefix);
         else
         {
            if (degree[i] > 2)  // tertiary hydrogen
            {
               if (0)          // Class disabled to save bit density
               result += SetPathBitsRec(mp, nbp,
                                        fp_counts, ncounts,
                                        NEXT_SEED(seed, 101), touched_indices,
                                        1,
                                        2, 2,  /* path length 2 to 3 */
                                        i, 0, -1,
                                        IGNORE_PATH_SYMBOL |
                                        IGNORE_TERM_SYMBOL |
                                        PROCESS_CHAINS,
                                        exclude_atom,
                                        prefix);
               result += SetPathBitsRec(mp, nbp,
                                        fp_counts, ncounts,
                                        NEXT_SEED(seed, 103), touched_indices,
                                        1,
                                        2, 5,  /* path length 2 to 3 */
                                        i, 0, -1,
                                        IGNORE_PATH_SYMBOL |
                                        FORCED_HETERO_END |  /* i-Pr...Q */
                                        PROCESS_CHAINS,
                                        exclude_atom,
                                        prefix);
            }
            if (H_count[i+1] >= 3)      // methyl
            {
               if (0)          // Class disabled to save bit density
               result += SetPathBitsRec(mp, nbp,
                                        fp_counts, ncounts,
                                        NEXT_SEED(seed,1103), touched_indices,
                                        1,
                                        2, 3,  /* path length 2 to 3 */
                                        i, 0, -1,
                                        IGNORE_PATH_SYMBOL |
                                        PROCESS_CHAINS,
                                        exclude_atom,
                                        prefix);
               result += SetPathBitsRec(mp, nbp,
                                        fp_counts, ncounts,
                                        NEXT_SEED(seed,1105), touched_indices,
                                        1,
                                        3, 6,  /* path length 4 to 4 */
                                        i, 0, -1,
                                        IGNORE_PATH_SYMBOL |
                                        FORCED_HETERO_END |  /* Me...Q */
                                        PROCESS_CHAINS,
                                        exclude_atom,
                                        prefix);
            }
         }
         touched_indices[i] = 0; /* down-dating */
         seed = old_seed;
         if (ap->color == 6) continue;

         if (H_count[i+1] > 1)  /* catch the difference between NH and NH2! */
         {
            touched_indices[i] = 1; /* updating */
            old_seed = seed;
            seed = NEXT_SEED(seed, 113);
            seed = NEXT_SEED(seed, ap->color);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     1, 5,  /* path length 1 to 5 */
                                     i, 0, -1,
// FORCED_HETERO_END |
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_CHAINS,
                                     exclude_atom,
                                     prefix);
            touched_indices[i] = 0; /* down-dating */
            seed = old_seed;
         }

         for (j=1; j<H_count[i+1]; j++)
         {
            seed = NEXT_SEED(seed, 61*j);
            seed = NEXT_SEED(seed, ap->color);
            ADD_BIT(fp_counts, ncounts, seed);
	    result++;
         }
         seed = old_seed;
      }
   }

   /* Compute ring paths */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (atom_status[i] <= 0) ap->color = 0;
      if (i+1 == exclude_atom) ap->color = 0;
      if (ap->color == 0) continue;
   }

   /* remove all bonds with only non-ring atoms from consideration */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bond_status[i] <= 0  &&
          atom_status[bp->atoms[0]-1] == 0  &&
          atom_status[bp->atoms[1]-1] == 0) bp->color = 0;
      if (bp->atoms[0] == exclude_atom) bp->color = 0;
      if (bp->atoms[1] == exclude_atom) bp->color = 0;
   }

// fprintf(stderr, "USE_RING_PATH\n");
   if (which_bits & USE_RING_PATH)
   {
      strcpy(prefix, "RP:");
      seed = RING_PATH_SEED;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
         if (0)         // Class disabled to save bit density
         result += SetPathBitsRec(mp, nbp,
                                  fp_counts, ncounts,
                                  seed, touched_indices,
                                  1,
                                  2, 3,
                                  i, 0, -1,
                                  PROCESS_CHAINS,
				  exclude_atom,
				  prefix);

         if (ap->color > 5  &&  ap->color < 10  &&  atom_status[i]>2)      // only start at common light atoms
         {
            seed = NEXT_SEED(seed, 61);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     3, 8,  /* 3 to 8 */
                                     i, 0, -1,
                                     STOP_AT_HEAVY_ATOM |            // don't cross very heavy atoms
                                     PROCESS_RING_CLOSURES,
				     exclude_atom,
				     prefix);
         }
         seed = old_seed;

         touched_indices[i] = 0; /* down-dating */
      }
   }

   /* Set the color property to represent all different atom types */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      ap->color = StringToInt(periodic_table, ap->atom_symbol);
      if (ap->color == 1  ||  i+1 == exclude_atom)
         ap->color = 0;   /* ignore hydrogens */
      else
         ap->color = ANY_COLOR; /* treat all other atoms alike */
   }

   /* Set the color property to represent the different bond type classes */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->bond_type == SINGLE)        bp->color = 1;
      else if (bp->bond_type == DOUBLE)   bp->color = 2;
      else if (bp->bond_type == TRIPLE)   bp->color = 3;
      else if (bp->bond_type == AROMATIC) bp->color = 4;
      else                                bp->color = 0;
      if (bp->atoms[0] == exclude_atom) bp->color = 0;
      if (bp->atoms[1] == exclude_atom) bp->color = 0;
   }

// fprintf(stderr, "USE_BOND_PATH\n");
   // Here, we have all non-trivial atoms mapped to ANY_COLOR while the
   // bond type is retained.
   if (which_bits & USE_BOND_PATH)
   {
      strcpy(prefix, "BP:");
      seed = BOND_PATH_SEED;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
         if (i+1 == exclude_atom) continue;
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         // start at branch node on ring
         if (degree[i] > 2  &&  atom_status[i] > 1)
         {
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     4, 4,         /* was 4 to 4 */
                                     i, 0, -1,
                                     PROCESS_CHAINS,
                                     exclude_atom,
                                     prefix);
            if (ap->rsize_flags&SPECIAL_RING)
	           result += SetPathBitsRec(mp, nbp,
                                        fp_counts, ncounts,
                                        NEXT_SEED(seed,217), touched_indices,
                                        1,
                                        5, 5,
                                        i, 0, -1,
                                        IGNORE_PATH_SYMBOL |
                                        PROCESS_CHAINS,
                                        exclude_atom,
                                        prefix);
         }
         /* Add other bits to catch poorly specified ring closures */
         seed = old_seed;
         seed = NEXT_SEED(seed, 11);
         result += SetPathBitsRec(mp, nbp,
                                  fp_counts, ncounts,
                                  seed, touched_indices,
                                  1,
                                  4, 6,         /* was 5 to 6 */
                                  i, 0, -1,
                                  // DEBUG_PATH |
                                  FORCED_RING_PATH |
                                  IGNORE_PATH_SYMBOL |
                                  PROCESS_RING_CLOSURES,
                                  exclude_atom,
                                  prefix);
         /* Add bits for paths starting with rare bond orders */
         for (j=0; j<nbp[i].n_ligands; j++)
         {
            bp = &mp->bond_array[nbp[i].bonds[j]];
            ai = nbp[i].atoms[j];
            if (ai+1 == exclude_atom) continue;
            if (bp->color == 0) continue;
            if (bp->bond_type != DOUBLE  &&
                bp->bond_type != TRIPLE)
               continue;
            if (atom_status[i] <= 0  &&  bp->bond_type != TRIPLE)
               continue;
            seed = old_seed;
            seed = NEXT_SEED(seed, bp->color*413);
            touched_indices[ai] = 1;    /* updating */
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     2,
                                     4, 5,
                                     ai, 0, i,
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_RING_CLOSURES |
                                     PROCESS_CHAINS,
                                     exclude_atom,
                                     prefix);
            touched_indices[ai] = 0;    /* down-dating */
         }
         seed = old_seed;
         touched_indices[i] = 0; /* down-dating */
      }
   }

   /* Set the color property to represent the different atom type classes */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (i+1 == exclude_atom)
      {
         ap->color = 0;
	 continue;
      }
      if (0 == strcmp(ap->atom_symbol,"L"))     // Check if the atom list
      {                                         // is completely in one class
         ap->color = 0;
         symbol_lists = mp->symbol_lists;
         while (!IsNULL(symbol_lists))
         {
            if (symbol_lists->atom == i+1)      /* found a match */
            {
               ap->color = 0;
               if (symbol_lists->logic == EXCLUSIVE) break;
               /* Step through all tokens in symbol list */
               strcpy(pat_buf,symbol_lists->string);
               for (tokp = strtok(pat_buf,",");
                    !IsNULL(tokp);
                    tokp = strtok((char *)NULL,","))
               {
#define HYDROGEN_FOUND 1
#define CARBON_FOUND   2
#define ONS_FOUND      4
#define OTHER_FOUND    8
                  if      (0 == strcmp("H",tokp)) ap->color |= HYDROGEN_FOUND;
                  else if (0 == strcmp("C",tokp)) ap->color |= CARBON_FOUND;
                  else if (0 == strcmp("O",tokp)) ap->color |= ONS_FOUND;
                  else if (0 == strcmp("N",tokp)) ap->color |= ONS_FOUND;
                  else if (0 == strcmp("S",tokp)) ap->color |= ONS_FOUND;
                  else                            ap->color |= OTHER_FOUND;
               }
               if      (ap->color == CARBON_FOUND) ap->color = 6;
               else if (ap->color == ONS_FOUND)    ap->color = 8;
               else if (ap->color == OTHER_FOUND)  ap->color = 8;
               else                                ap->color = 0;
            }
            symbol_lists = symbol_lists->next;
         }
      }
      else
      {
         if (0 == strcmp("H", ap->atom_symbol))
            ap->color = 0;      /* ignore hydrogens */
         else if (0 == strcmp("D", ap->atom_symbol))
            ap->color = 0;      /* ignore hydrogens */
         else if (0 == strcmp("T", ap->atom_symbol))
            ap->color = 0;      /* ignore hydrogens */
         else if (0 == strcmp("C", ap->atom_symbol))
            ap->color = 6;      /* carbon second row elements are one class */
         else if (0 == strcmp("N", ap->atom_symbol))
            ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
         else if (0 == strcmp("O", ap->atom_symbol))
            ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
         else if (0 == strcmp("S", ap->atom_symbol))
            ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
         else if (0 == strcmp("Q", ap->atom_symbol))
            ap->color = 8;      /* nitrogen, oxigen, and sulfur are one class */
         else if (0 == strcmp("A", ap->atom_symbol))
            ap->color = 0;
         else 
         {
            tmp = StringToInt(periodic_table, ap->atom_symbol);
            if (1 < tmp  &&  tmp < 115)
               ap->color = 8;     /* non-carbon is the only second class */
            else        /* This could be R atoms or other odd things */
               ap->color = 0;
         }
      }
      if (ap->color > 115) ap->color = 0;       /* ignore special atom types */
   }

// fprintf(stderr, "USE_HCOUNT_CLASS_PATH\n");
   // Here, we unify atom types to carbon/hetero distinction
   if (which_bits & USE_HCOUNT_CLASS_PATH)
   {
      strcpy(prefix, "HCP:");
      /* generate a short path for each atom that has a hydrogen */
      seed = HCOUNT_CLASS_PATH_SEED;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
         if (i+1 == exclude_atom) continue;
         if (H_count[i+1] == 0) continue;
         if (ap->color == 6  &&  H_count[i+1] < 2) continue;
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
         touched_indices[i] = 1; /* updating */
         if (ap->color == 6)
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     2, 4,  /* path length 1 to 3 */
                                     i, 0, -1,
                                     FORCED_HETERO_END |
                                     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);
         else
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     2, 5,
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
                                     FORCED_HETERO_END |
                                     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);
         seed = old_seed;
         touched_indices[i] = 0; /* down-dating */
      }
   }

// fprintf(stderr, "USE_ATOM_CLASS_PATH\n");
   /* Do a small first pass with defined bond types */
   if (which_bits & USE_ATOM_CLASS_PATH)
   {
      strcpy(prefix, "ACP:");
      // seed = ATOM_CLASS_PATH_SEED+117;
      seed = NEXT_SEED(ATOM_CLASS_PATH_SEED,117);
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
         if (i+1 == exclude_atom) continue;
         if (ap->color == 6  &&  degree[i] < 3) continue;
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
         if (ap->color == 6)
         {
            if (0)         // Class disabled to save bit density
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     3, 4,  /* path length 3 to 4 */
                                     i, 0, -1,
                                     FORCED_HETERO_END |
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);
         }
         else
         {
            if (0)         // Class disabled to save bit density
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     2, 3,  /* path length 3 to 3 */
                                     i, 0, -1,
                                     // IGNORE_PATH_SYMBOL |
                                     PROCESS_CHAINS,
                                     exclude_atom,
				     prefix);
            if (0)         // Class disabled to save bit density
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     4, 5,  /* path length 4 to 5 */
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
                                     FORCED_HETERO_END |
                                     PROCESS_CHAINS,
                                     exclude_atom,
				     prefix);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     3, 4,  /* path length 4 to 7 */
                                     i, 0, -1,
                                     FORCED_RING_PATH |
                                     PROCESS_RING_CLOSURES |
                                     PROCESS_CHAINS,
                                     exclude_atom,
				     prefix);
         }
         touched_indices[i] = 0;
         seed = old_seed;
      }
   }

   /* Set the color property to only a single class */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (SINGLE <= bp->bond_type  &&  bp->bond_type <= ANY_BOND)
         bp->color = 5;
      else
         bp->color = 0;
      if (bp->atoms[0] == exclude_atom) bp->color = 0;
      if (bp->atoms[1] == exclude_atom) bp->color = 0;
   }

// fprintf(stderr, "USE_ATOM_CLASS_PATH\n");
   // Here, we've unified atom types to carbon/hetero and made all bond types
   // identical
   if (which_bits & USE_ATOM_CLASS_PATH)
   {
      strcpy(prefix, "ACP:");
      seed = ATOM_CLASS_PATH_SEED;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
         if (i+1 == exclude_atom) continue;
         if (ap->color == 6) continue;
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
         // if (ap->color == 6)
         {
            if (0*degree[i] > 2)         // Class disabled to save bit density
               result += SetPathBitsRec(mp, nbp,
                                        fp_counts, ncounts,
                                        seed, touched_indices,
                                        1,
                                        3, 3,  /* path length 3 to 3 */
                                        i, 0, -1,
                                        FORCED_HETERO_END |
                                        PROCESS_CHAINS,
					exclude_atom,
					prefix);
         }
         // else
         {
            if(0)         // Class disabled to save bit density
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     3, 4,  /* path length 3 to 4 */
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);
            seed = NEXT_SEED(seed, 23+ap->color*19);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     3, 9,  /* path length 3 to 9 */
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_RING_CLOSURES,
                                     exclude_atom,
				     prefix);
         }
         seed = old_seed;
         touched_indices[i] = 0; /* down-dating */
      }
      /* Q-Q and Q-C ring bond count */
      qq_count=0;
      qc_count=0;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (bond_status[i] == 0) continue;
         if (bp->color == 0) continue;
         if (bp->atoms[0] == exclude_atom) continue;
         if (bp->atoms[1] == exclude_atom) continue;
         ai1 = mp->atom_array[bp->atoms[0]-1].color;
         ai2 = mp->atom_array[bp->atoms[1]-1].color;
         if (ai1 == 0  ||  ai2 == 0) continue;
         if (ai1 == 6  &&  ai2 == 6) continue;
         if (ai1 != 6  &&  ai2 != 6)
         {
            qq_count++;
            for (j=3; j<9; j++)        /* set bits for not too large ring size */
               if (bp->rsize_flags&(1<<j))
	       {
                  ADD_BIT(fp_counts, ncounts, NEXT_SEED(ATOM_CLASS_PATH_SEED*17, j*8));
		  result++;
	       }
         }
         else
         {
            qc_count++;
            for (j=3; j<9; j++)        /* set bits for not too large ring size */
               if (bp->rsize_flags&(1<<j))
	       {
                  ADD_BIT(fp_counts, ncounts, NEXT_SEED(ATOM_CLASS_PATH_SEED*19, j*8));
		  result++;
	       }
         }
      }
      seed = 2*ATOM_CLASS_PATH_SEED+3;
      for (i=1; i<=qq_count; i=(int)(1+i*1.5))
      {
         seed = NEXT_SEED(seed, i*153);
         ADD_BIT(fp_counts, ncounts, seed);
	 result++;
         seed = NEXT_SEED(seed, 53);
         if (i <= 1)
         {
            ADD_BIT(fp_counts, ncounts, seed);
            result++;
         }
      }
      seed = 3*ATOM_CLASS_PATH_SEED+5;
      for (i=1; i<=MIN(qc_count,2); i++)
      {
         seed = NEXT_SEED(seed, i*157);
         ADD_BIT(fp_counts, ncounts, seed);
	 result++;
      }
      for (i=3; i<=qc_count; i=(int)(i*1.8))
      {
         seed = NEXT_SEED(seed, i*157);
         ADD_BIT(fp_counts, ncounts, seed);
	 result++;
      }
   }

   /* Compute ring patters with at least one cycle */
   /* remove bonds from consideration that don't have at least one ring atom */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->atoms[0] == exclude_atom) continue;
      if (bp->atoms[1] == exclude_atom) continue;
      bp->color = 5;
      /* ignore non-ring bonds */
      if (atom_status[bp->atoms[0]-1] == 0  &&
          atom_status[bp->atoms[1]-1] == 0)
      {
         bp->color = 0;
      }
      else
      {
         // may be redundant since atom types have already been unified
         ap = &mp->atom_array[bp->atoms[0]-1];
         if (ap->color != 0)
         {
            if (ap->color != 6) ap->color = 8;  /* non-carbons in one class */
            if (ap->rsize_flags == 0) ap->color = 0;
         }
         // may be redundant since atom types have already been unified
         ap = &mp->atom_array[bp->atoms[1]-1];
         if (ap->color != 0)
         {
            if (ap->color != 6) ap->color = 8;  /* non-carbons in one class */
            if (ap->rsize_flags == 0) ap->color = 0;
         }
      }
   }

// fprintf(stderr, "USE_RING_PATTERN\n");
   // Here, we have bond type ignored and atom types mapped to C and Q
   if (which_bits & USE_RING_PATTERN)
   {
      strcpy(prefix, "RPT:");
      /* first process ring bond paths with atom classes */
      seed = RING_PATTERN_SEED;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;
         if (i+1 == exclude_atom) continue;
         /* Don't process fragments starting at carbon in 6-ring only */
         if (ap->color == 6  &&  0 == (ap->rsize_flags&SPECIAL_RING)) continue;
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
         result += SetPathBitsRec(mp, nbp,
                                  fp_counts, ncounts,
                                  seed, touched_indices,
                                  1,
                                  3, 3,  /* ring bond path size 3 to 3 */
                                  i, 0, -1,
                                  PROCESS_CHAINS,
                                  exclude_atom,
                                  prefix);
         seed = old_seed;
         touched_indices[i] = 0; /* down-dating */
      }

      /* Now, we only include complete rings but ignore atom-type */
      /* 'A' atoms are now included nodes */
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         /* add 'A' atom to standard class */
         if (0 == strcmp("A", ap->atom_symbol))
         {
            ap->color = 9;
         }
         if (i+1 == exclude_atom) ap->color = 0;
         if (ap->color == 0) continue;
         ap->color = 9;    /* all ring atoms in same class */
      }
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (bp->atoms[0] == exclude_atom) bp->color = 0;
         if (bp->atoms[1] == exclude_atom) bp->color = 0;
         if (bond_status[i] <= 0) bp->color = 0;
      }
// seed = RING_PATTERN_SEED+23;
      seed = NEXT_SEED(RING_PATTERN_SEED, 23);
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color == 0) continue;
         if (i+1 == exclude_atom) continue;
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
         if (0)         // Class disabled to save bit density
         result += SetPathBitsRec(mp, nbp,
                                  fp_counts, ncounts,
                                  seed, touched_indices,
                                  1,
                                  4, 17,  /* ring size 4 to 17 */
                                  i, 0, -1,
                                  IGNORE_PATH_SYMBOL |
                                  PROCESS_RING_CLOSURES,
                                  exclude_atom,
                                  prefix);

         seed = old_seed;
         if (0)         // Class disabled to save bit density
         if (atom_status[i] > 2)        // start at ring fusion
         {
            seed = NEXT_SEED(seed, 61);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     6, 17,  /* ring path size 6 to 17 */
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_RING_CLOSURES,
                                     exclude_atom,
                                     prefix);
         }

         seed = old_seed;
         if (0)         // Class disabled to save bit density
         if (degree[i] > 2  && atom_status[i] >= 2)     // start at ring substituents
         {
            seed = NEXT_SEED(seed, 67);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     6, 17,  /* ring path size 6 to 17 */
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_RING_CLOSURES,
                                     exclude_atom,
                                     prefix);
         }
         seed = old_seed;

         touched_indices[i] = 0; /* down-dating */
      }
   }

// fprintf(stderr, "USE_RING_SIZE_COUNTS\n");
   if (which_bits & USE_RING_SIZE_COUNTS)
   {
      strcpy(prefix, "RSC:");
      for (j=3; j<10; j++)      /* loop through ring_sizes */
      {
         nrbonds = 0;
         for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
            if (bp->atoms[0] != exclude_atom  &&
                bp->atoms[1] != exclude_atom  &&
                (bp->rsize_flags&(1<<j))) nrbonds++;
         seed = RING_SIZE_SEED;
         seed = NEXT_SEED(seed, j*13);
         for (i=1; i<100; i*=2)
         {
            if (nrbonds >= j*i)
            {
               seed = NEXT_SEED(seed, i);
               /* don't set a bit if just one 5- or 6-ring */
               if ((j != 6  &&  j != 5)  || nrbonds > j*i)
	       {
                  ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
		  result++;
	       }
               if (j != 6  &&  j != 5)    /* set one more bit for odd rings */
               {
                  seed = NEXT_SEED(seed, 17);
                  ADD_BIT_COUNT(fp_counts, ncounts, seed, nrbonds);
		  result++;
               }
            }
            else
               break;
         }
      }
      /* Set bits for different ring sizes connected by a bond */
      for (j=0; j<15; j++)      /* loop through ring_sizes */
         for (k=0; k<15; k++)      /* loop through ring_sizes */
            rscounts[j][k] = 0;
      /* count connections */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         if (bp->atoms[0] == exclude_atom) continue;
         if (bp->atoms[1] == exclude_atom) continue;
         for (j=3; j<15; j++)      /* loop through ring_sizes */
            for (k=3; k<15; k++)      /* loop through ring_sizes */
            {
               if (j==k) continue;
               if ((mp->atom_array[bp->atoms[0]-1].rsize_flags&(1<<j))  &&
                   (mp->atom_array[bp->atoms[1]-1].rsize_flags&(1<<k)))
               {
                  rscounts[j][k]++;
                  rscounts[k][j]++;
               }
            }
      }
      /* set bits */
      for (j=3; j<9; j++)      /* loop through not too large ring_sizes */
         for (k=j+1; k<9; k++)      /* loop through not too large ring_sizes */
         {
            if (rscounts[j][k] == 0) continue;
            seed = 2*RING_SIZE_SEED;
            seed = NEXT_SEED(seed, (j+k)*11);
            seed = NEXT_SEED(seed, j*k*13);
            seed = NEXT_SEED(seed, 19);
            ADD_BIT(fp_counts, ncounts, seed);
	    result++;
            /* more special links => more bits */
            if (rscounts[j][k]==1  ||  j==6  ||  k==6) continue;
            seed = 2*RING_SIZE_SEED+23;
            seed = NEXT_SEED(seed, (j+k)*11);
            seed = NEXT_SEED(seed, j*k*13);
            seed = NEXT_SEED(seed, 19);
            ADD_BIT(fp_counts, ncounts, seed);
	    result++;
         }
   }

   /* Set the color property to represent all different atom types */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      ap->color = StringToInt(periodic_table, ap->atom_symbol);
      if (ap->color <= 1) ap->color = 0;        /* ignore hydrogens */
      /* mark special atom types */
      if (ap->color > 115) ap->color = -1;
      if (0 == strcmp("A", ap->atom_symbol)) ap->color = -1;
      if (i+1 == exclude_atom) ap->color = 0;
      if (ap->color > 1  &&
          (!as_query                 ||
           ap->sub_desc == SUB_AS_IS ||
           (ap->sub_desc != NONE     &&
            ap->sub_desc != SUB_MORE &&
            ap->sub_desc == degree[i]+SUB_ONE-1)))
      {
         ap->color += 32*degree[i];
      }
      else
         ap->color = 0;
   }
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (SINGLE <= bp->bond_type  &&  bp->bond_type <= ANY_BOND)
         bp->color = 5;
      else
         bp->color = 0;
      if (bp->atoms[0] == exclude_atom) bp->color = 0;
      if (bp->atoms[1] == exclude_atom) bp->color = 0;
   }

// fprintf(stderr, "USE_DEGREE_PATH\n");
   /*
    * add special bits for paths starting at atoms with >= 4 neighbours or
    * methyl atoms
    *
    * atom color represent degree bond color identical for all bonds to non-H
    */
   if (which_bits & USE_DEGREE_PATH)
   {
      strcpy(prefix, "DP:");
      seed = DEGREE_PATH_SEED;

      // set bits for degree paths starting with special carbon atoms
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (ap->color <= 0) continue;  // only process if degree defined
         if (i+1 == exclude_atom) continue;
         /* don't start on usual atoms */
         if (degree[i] <= 3  &&		// include if high degree
             atom_status[i] <= 2  &&	// include if ring fusion
             degree[i] != 1)		// include if terminal atom
            continue;
         /* only keep terminals if methyl */
         if (degree[i] == 1  &&  0 != strcmp(ap->atom_symbol, "C"))
            continue;
// fprintf(stderr,"degree of %s atom %d is %d, status is %d\n",
// ap->atom_symbol, i+1, degree[i], atom_status[i]);
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, ap->color);
// fprintf(stderr, "%s(%d) has degree %d(%d)\n",
// ap->atom_symbol, i+1, degree[i], ap->sub_desc);
         result += SetPathBitsRec(mp, nbp,
                                  fp_counts, ncounts,
                                  seed, touched_indices,
                                  1,
                                  2, 4,  /* path length 1 to 3 */
                                  i, 0, -1,
                                  // IGNORE_TERM_SYMBOL |
                                  IGNORE_PATH_SYMBOL |
                                  PROCESS_CHAINS,
                                  exclude_atom,
				  prefix);
         /* special CH fusion atoms */
         if (atom_status[i] > 2  &&  H_count[i+1] >= 1)
         {
            seed = NEXT_SEED(seed, 219);
            result += SetPathBitsRec(mp, nbp,
                                     fp_counts, ncounts,
                                     seed, touched_indices,
                                     1,
                                     2, 5,
                                     i, 0, -1,
                                     IGNORE_PATH_SYMBOL |
                                     PROCESS_CHAINS,
				     exclude_atom,
				     prefix);
         }
         seed = old_seed;
         touched_indices[i] = 0; /* down-dating */
      }

      // set bits for degree paths starting with hetero atoms
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (i+1 == exclude_atom) continue;
         tmp = StringToInt(periodic_table, ap->atom_symbol);
         if (1 >= tmp  ||  tmp >= 115) continue;
         if (tmp == 6) continue;
         if (tmp < 10)	continue; // exclude common hetero atoms
         touched_indices[i] = 1; /* updating */
         old_seed = seed;
         seed = NEXT_SEED(seed, tmp);
// fprintf(stderr, "%s(%d) has degree %d(%d)\n",
// ap->atom_symbol, i+1, degree[i], ap->sub_desc);
         result += SetPathBitsRec(mp, nbp,
                                  fp_counts, ncounts,
                                  seed, touched_indices,
                                  1,
                                  2, 2,
                                  i, 0, -1,
                                  PROCESS_CHAINS,
                                  exclude_atom,
				  prefix);
         if (0) // might overly populate complexes
         result += SetPathBitsRec(mp, nbp,
                                  fp_counts, ncounts,
                                  101+seed, touched_indices,
                                  1,
                                  2, 2,
                                  i, 0, -1,
                                  FORCED_RING_PATH |
                                  PROCESS_CHAINS,
                                  exclude_atom,
				  prefix);
         seed = old_seed;
         touched_indices[i] = 0; /* down-dating */
      }

   }

// fprintf(stderr, "USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS\n");
   if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS))
   {
      /* Collect length_matrix */
      /* allocate storage length_matrix */
      length_tmp = TypeAlloc(mp->n_atoms*mp->n_atoms, int);
      /* allocat indices */
      length_matrix = TypeAlloc(mp->n_atoms, int*);
      for (i=0; i<mp->n_atoms; i++)
      {
         /* set relative pointers */
         length_matrix[i] = length_tmp+i*mp->n_atoms;
         touched_indices[i] = 0;
      }
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
          if (i+1 == exclude_atom) continue;
          touched_indices[i] = 1; /* updating */
// if (FALSE) fprintf(stderr, "starting path search at atom %d(%d)\n", i+1, ap->color);
          SetPathLengthFlags(mp, touched_indices,
                             i, 0, i,
                             12,      /* path perception distance <= 12 */
                             length_matrix, nbp,
                             exclude_atom,
			     prefix);
          touched_indices[i] = 0; /* down-dating */
      }
   }

   /*
    * This screen class will catch non-linear fragments composed of rather
    * frequent linear sub-fragments
    */
   if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS))
   {
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         ap->color = StringToInt(periodic_table, ap->atom_symbol);
         if (0 == strcmp("H", ap->atom_symbol))
            ap->color = 0;      /* ignore hydrogens */
         else if (0 == strcmp("D", ap->atom_symbol))
            ap->color = 0;      /* ignore hydrogens */
         else if (0 == strcmp("T", ap->atom_symbol))
            ap->color = 0;      /* ignore hydrogens */
         else if (0 == strcmp("Q", ap->atom_symbol))
            ap->color = HETERO;
         else if (0 == strcmp("A", ap->atom_symbol))
            ap->color = GENERIC;
         else if (0 == strcmp("L", ap->atom_symbol))
            ap->color = GENERIC;
         else if (0 == strcmp("C", ap->atom_symbol))
         {
            ap->color = 6;      /* carbon second row elements are one class */
            if (cdegree[i] >= 3)
               ap->color = CSP3;
         }
         else if (ap->color > 1  &&  ap->color < 115)
               ap->color = HETERO;
         else        /* This could be R atoms or other odd things */
            ap->color = 0;
         if (i+1 == exclude_atom) ap->color = 0;
// fprintf(stderr, "atom %d has color %d\n", i+1, ap->color);
      }
      /*
       * Bond colors are already set OK, i.e. equal for A-H bonds
       */
      /* NOP */

      /* Now we start setting bits */
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
         if (i+1 == exclude_atom) continue;
         /* Spiders have at least three legs (;-) */
         if (degree[i] < 3) continue;
         /* Spider needs to be special atom or carbon */
         if (ap->color != CSP3  && ap->color != 6) continue;
         touched_indices[i] = 1; /* updating */
         for (j=0; j<=MAX_SPIDER; j++)
            hetero[j] = csp3[j] = 0;
         if (which_bits & USE_CLASS_SPIDERS)
         {
            strcpy(prefix, "CS:");
            SpecialNeighboursRec(mp, touched_indices,
                                 1, i, MAX_SPIDER,
                                 csp3, hetero, nbp,
                                 exclude_atom,
				 prefix);
         }
         touched_indices[i] = 0; /* down-dating */

         /* set bits for spiders with one CSP3 atom and two heteros */
         if (which_bits & USE_CLASS_SPIDERS)
         for (j=1; j<=MAX_SPIDER; j++)
         {
            strcpy(prefix, "CS:");
            if (csp3[j] == 0) continue;
            seed = CLASS_SPIDER_SEED;
            if (ap->color == HETERO)
            {
// seed = CLASS_SPIDER_SEED+HETERO*8+CSP3*11;
               seed = NEXT_SEED(seed, HETERO*8);
               seed = NEXT_SEED(seed, CSP3*11);
            }
            else
            {
// seed = CLASS_SPIDER_SEED+6*8+CSP3*11;
               seed = NEXT_SEED(seed, 6*8);
               seed = NEXT_SEED(seed, CSP3*11);
            }
            for (j1=1; j1<=MAX_SPIDER; j1++)
            {
               tmp1 = hetero[j1];
               if (tmp1 <= 0) continue;
               for (j2=j1; j2<=MAX_SPIDER; j2++)
               {
                  tmp2 = hetero[j2];
                  if (j2 == j1) tmp2--; /* consumed in outer loop */
                  if (tmp2 <= 0) continue;
                  old_seed = seed;
                  seed = NEXT_SEED(seed, j);
                  seed = NEXT_SEED(seed, j1+j2);
                  ADD_BIT(fp_counts, ncounts, seed);
		  result++;
// fprintf(stderr, " setting CSP3 bit for %d: %d|%d|%d\n", ap->color, j, j1, j2);
                  seed = old_seed;
               }
            }
         }

         /* don't hetero-spider normal carbons */
         if (ap->color != CSP3) continue;
         /* set bits for spiders with three defined HETERO atoms */
         if (which_bits & USE_CLASS_SPIDERS)
         for (j=1; j<=MAX_SPIDER; j++)
         {
            strcpy(prefix, "CS:");
            if (hetero[j] == 0) continue;
            seed = CLASS_SPIDER_SEED;
            if (ap->color == HETERO)
            {
// seed = CLASS_SPIDER_SEED+HETERO*8+HETERO*11;
               seed = NEXT_SEED(seed, HETERO*8);
               seed = NEXT_SEED(seed, HETERO*11);
            }
            else
            {
// seed = CLASS_SPIDER_SEED+6*8+HETERO*11;
               seed = NEXT_SEED(seed, 6*8);
               seed = NEXT_SEED(seed, HETERO*11);
            }
            for (j1=j; j1<=MAX_SPIDER; j1++)
            {
               tmp1 = hetero[j1];
               if (j1 == j) tmp1--; /* we've consumed this one in outer loop */
               if (tmp1 <= 0) continue;
               for (j2=j1; j2<=MAX_SPIDER; j2++)
               {
                  tmp2 = hetero[j2];
                  if (j2 == j)  tmp2--; /* consumed in outer loop */
                  if (j2 == j1) tmp2--; /* consumed in outer loop */
                  if (tmp2 <= 0) continue;
                  old_seed = seed;
                  seed = NEXT_SEED(seed, j*j1*j2);
                  ADD_BIT(fp_counts, ncounts, seed);
		  result++;
// fprintf(stderr, " setting HETERO bit for %d: %d|%d|%d\n", ap->color, j, j1, j2);
                  // Additional bits for quarternary centers
                  if (degree[i] > 3) {ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, 501)); result++;}
                  seed = old_seed;
               }
            }
         }
      }

      /**
       * Collect bits that represent feature/path_length/feature triples.
       */
      if (which_bits & USE_FEATURE_PAIRS)
      {
         strcpy(prefix, "FP:");
         /* set feature flags in atom colors */
         for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         {
            if (0 == strcmp(ap->atom_symbol, "C"))      flags = C_FLAG;
            else if (0 == strcmp(ap->atom_symbol, "O")) flags = O_FLAG;
            else if (0 == strcmp(ap->atom_symbol, "N")) flags = N_FLAG;
            else if (0 == strcmp(ap->atom_symbol, "S")) flags = S_FLAG;
            else if (0 == strcmp(ap->atom_symbol, "P")) flags = P_FLAG;
            else if (AtomSymbolMatch(ap->atom_symbol,"F,Cl,Br,I,At"))
               flags = X_FLAG;
            else
               flags = 0;
            if (ap->color == HETERO)
	    {
	       flags |= HETERO_FLAG;
// fprintf(stderr, "HETERO_FLAG set for atom %d\n", i+1);
	    }
            if (cdegree[i] >= 3)     flags |= CSP3_FLAG;
// fprintf(stderr, "%s atom %d has degree %d\n", ap->atom_symbol, i+1, degree[i]);
            if (degree[i] >= 4)
            {
               flags |= QUART_FLAG;
            }
            if (atom_status[i] > 0  &&  degree[i] >= 3)
            {
               flags |= RING_SUBST_FLAG;
               if (0 != (ap->rsize_flags&SPECIAL_RING))
               {
                   flags |= RS_SPECIAL_FLAG;
// fprintf(stderr, "RS_SPECIAL_FLAG set for atom %d\n", i+1);
               }
            }
            if (i+1 == exclude_atom) flags = 0;
            ap->color = flags;
         }
         /* collect bits for selected feature pairs */
         if(0)         // Class disabled to save bit density
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           CSP3_FLAG,
                           HETERO_FLAG,         /* to hetero or ring subst */
                           2, 3,               /* with path length 1 to 9 */
                           FALSE,               /* don't use count */
                           TRUE,                /* use atom type flags */
                           length_matrix,
                           1237,                /* seed = 1237 */
			   exclude_atom,
			   prefix);
         if(0)         // Class disabled to save bit density
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           HETERO_FLAG,         /* from ring substitution */
                           HETERO_FLAG,         /* to hetero or ring subst */
                           1, 12,               /* with path length 1 to 10 */
                           TRUE,                /* use count */
                           TRUE,                /* use atom type flags */
                           length_matrix,
                           1237,                /* seed = 1237 */
			   exclude_atom,
			   prefix);
         if (1)
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           RING_SUBST_FLAG,     /* from ring substitution */
                           RING_SUBST_FLAG,     /* to ring substitution */
                           5, 7,               /* with path length 1 to 12 */
                           FALSE,               /* don't use count */
                           TRUE,                /* use atom type flags */
                           length_matrix,
                           2237,                /* seed = 2237 */
			   exclude_atom,
			   prefix);
         if (0)         // Class disabled to save bit density
         result +=
            SetRingSizePairBits(mp, fp_counts, ncounts,
			    RING_SUBST_FLAG,     /* from ring substitution */
                            HETERO_FLAG,         /* to hetero atom */
                            2, 5,               /* with path length 1 to 12 */
                            length_matrix,
                            3237,                /* seed = 2237 */
			    exclude_atom,
			    prefix);
         if (0)         // Class disabled to save bit density
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
			   RS_SPECIAL_FLAG,/* from special ring substitution */
                           HETERO_FLAG,    /* to hetero atom */
                           2, 4,                /* with path length 1 to 5 */
                           TRUE,               /* don't use count */
                           FALSE,                /* use atom type flags */
                           length_matrix,
                           3237,                /* seed = 3237 */
			   exclude_atom,
			   prefix);
         if (1)
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           QUART_FLAG,          /* from quartenary atom */
                           HETERO_FLAG,         /* to hetero or ring subst */
                           1, 8,                /* with path length 1 to 8 */
                           FALSE,               /* don't use count */
                           TRUE,                /* use atom type flags */
                           length_matrix,
                           4237,                /* seed = 4237 */
			   exclude_atom,
			   prefix);
         if(1)
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           QUART_FLAG,          /* from quartenary atom */
                           RING_SUBST_FLAG,
                           1, 6,                /* with path length 1 to 8 */
                           FALSE,               /* don't use count */
                           TRUE,                /* use atom type flags */
                           length_matrix,
                           5237,                /* seed = 4237 */
			   exclude_atom,
			   prefix);
         if(1)
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           X_FLAG,          /* from halogen atom */
                           CSP3_FLAG,
                           1, 1,                /* with path length 1 to 1 */
                           FALSE,               /* don't use count */
                           TRUE,                /* use atom type flags */
                           length_matrix,
                           15237,                /* seed = 15237 */
			   exclude_atom,
			   prefix);
         if (0)         // Class disabled to save bit density
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
			   RS_SPECIAL_FLAG,/* from special ring substitution */
                           RING_SUBST_FLAG,    /* to hetero atom */
                           2, 4,                /* with path length 1 to 5 */
                           TRUE,               /* don't use count */
                           FALSE,                /* use atom type flags */
                           length_matrix,
                           6237,                /* seed = 6237 */
			   exclude_atom,
			   prefix);
         if(0)         // Class disabled to save bit density
         result +=
            SetFeatureBits(mp, fp_counts, ncounts,
                           HETERO_FLAG,         /* from hetero */
                           RING_SUBST_FLAG,     /* to ring substitution */
                           1, 6,                /* with path length 1 to 8 */
                           FALSE,               /* don't use count */
                           TRUE,                /* use atom type flags */
                           length_matrix,
                           7237,                /* seed = 4237 */
			   exclude_atom,
			   prefix);

         /* Set bits for ring-subst/ring-subst/hetero triples */
         if (0) // too many spurious bits
         {
            for (i1=0, ap1=mp->atom_array; i1<mp->n_atoms; i1++, ap1++)
            {
               if (i1+1 == exclude_atom) continue;
               if (0 == (ap1->color&RING_SUBST_FLAG)) continue;
               /* first atom must be in a non-sixmembered ring */
               if (!(ap1->rsize_flags&SPECIAL_RING)) continue;
               for (i2=0, ap2=mp->atom_array; i2<mp->n_atoms; i2++, ap2++)
               {
                  if (i1 == i2) continue;
                  if (i2+1 == exclude_atom) continue;
                  if (0 == (ap2->color&RING_SUBST_FLAG)) continue;
                  for (i3=0, ap3=mp->atom_array; i3<mp->n_atoms; i3++, ap3++)
                  {
                     if (i1 == i3) continue;
                     if (i2 == i3) continue;
                     if (i3+1 == exclude_atom) continue;
                     if (0 == (ap3->color&HETERO_FLAG)) continue;
                     for (j=2; j<=6; j++)
                     {
                        if (0 == (length_matrix[i1][i2]&(1<<j))) continue;
                        for (j1=2; j1<=5; j1++)
                        {
                           if (0 == (length_matrix[i2][i3]&(1<<j1))) continue;
                           for (j2=2; j2<=5; j2++)
                           {
                              if (0 == (length_matrix[i3][i1]&(1<<j2)))
                                 continue;
                              // Make sure we have a real triangle
                              if (j+j1  == j2) continue;
                              if (j+j2  == j1) continue;
                              if (j1+j1 == j)  continue;
                              if (j+j1+j2 > 11) continue;  // spider too large
                              if (j+j1+j2 <  8) continue;  // spider too small
                              seed = CLASS_SPIDER_SEED;
                              seed = NEXT_SEED(seed, j+j1+j2);
                              seed = NEXT_SEED(seed, j*j1*j2);
                              // distinguish ring sizes of second ring-subst
                              // for (k=3; k<15; k++)
                              for (k=3; k<9; k++)
                              {
                                 if (!(ap2->rsize_flags & (1<<k))) continue;
                                 ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, k*213));
				 result++;
                              }
// fprintf(stderr, "processing tripple (%d|%x)-%d-(%d|%x)-%d-(%d|%s)-%d\n",
// i1+1, ap1->rsize_flags, j,
// i2+1, ap2->rsize_flags, j1,
// i3+1, ap3->atom_symbol, j2);
                           }
                        }
                     }
                  }
               }
            }
         }
      }

   }

   /* Collect bits that describe scaffolds. Those bits cannot (yet?) be used for SSS screening */
   /* The method first collects a variant of the extended connectivity but derived only        */
   /* on the ring bond graph. The bits are then set by combining atom features with the        */
   /* corresponding scaffold extended connectivity of that atom. In addition, we add           */
   /* bits to code connectivity of scaffolds by certain numbers of bonds using the             */
   /* previously collected length_matrix. */
   if (!as_query  &&  (which_bits & (USE_NON_SSS_BITS)))
   {
      strcpy(prefix, "NSB:");
      seed = NON_SSS_SEED;
      extcon = TypeAlloc(mp->n_atoms, int);
      extcon2 = TypeAlloc(mp->n_atoms, int);
      /* initialized extended connectivity */
      for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
      {
          if (atom_status[j] <= 0) continue;
          extcon[j] = ap->rsize_flags;
      }
      /* propagate extended connectivity to neighbours for a few cycles */
      for (i=0; i<32; i++)
      {
          for (j=0; j<mp->n_atoms; j++) extcon2[j] = 0;
          for (j=0; j<mp->n_atoms; j++)
          {
              /* skip non-ring atoms */
              if (atom_status[j] <= 0) continue;
              extcon2[j] = atom_status[j]*3+(extcon[j]*0xF);
              sum = 0; prod = 0;
              for (jj=0; jj<nbp[j].n_ligands; jj++)
              {
                  // only propagate through ring bonds
                  if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
                  sum += extcon[nbp[j].atoms[jj]];
                  prod *= 0xFF&(1+extcon[nbp[j].atoms[jj]]);
              }
              extcon2[j] =
                  ((extcon2[j]+(sum*191)+prod)<<8) +
                  ((extcon2[j]&0xFF0000)>>16);
              extcon2[j] &= 0xFFFFFF;
          }
          for (j=0; j<mp->n_atoms; j++) extcon[j] = extcon2[j];
      }

      /* propagate smallest hash to all members of ring system */
      for (;;)
      {
          changed = FALSE;
          for (j=0; j<mp->n_atoms; j++)
          {
              /* skip non-ring atoms */
              if (atom_status[j] <= 0) continue;
              for (jj=0; jj<nbp[j].n_ligands; jj++)
              {
                  // only propagate through ring bonds
                  if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
                  if (extcon[j] < extcon[nbp[j].atoms[jj]])
                  {
                      changed = TRUE;
                      extcon[nbp[j].atoms[jj]] = extcon[j];
                  }
              }
          }
          if (!changed) break;
      }

      /* Now, use extcon to set bits */
      for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
      {
          if (j+1 == exclude_atom) continue;
          /* skip non-ring atoms */
          if (atom_status[j] <= 0) continue;
          if (which_bits & USE_SCAFFOLD_IDS)
          {
               strcpy(prefix, "SI:");
// ADD_BIT(fp_counts, ncounts, extcon[j]*1013 + seed);
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j]*10013));
// ADD_BIT(fp_counts, ncounts, NEXT_SEED(seed, extcon[j]*40013));
	      result += 1;
          }

          if (which_bits & USE_SCAFFOLD_LINKS)
          {
               strcpy(prefix, "SL:");
              /* needs to be a substitution point */
              if (degree[j] <= atom_status[j]) continue;
              for (jj=0; jj<mp->n_atoms; jj++)
              {
                  if (jj+1 == exclude_atom) continue;
                  /* skip non-ring atoms */
                  if (atom_status[jj] <= 0) continue;
                  if (j == jj) continue;
                  /* needs to be a substitution point */
                  if (degree[jj] <= atom_status[jj]) continue;
                  for (k = 0; k<12; k++)
                      if (length_matrix[j][jj] == (1<<k))
                          break;
                  if (k >= 12) continue;    // multiple paths => no chain connection
                  if (k >= 3) continue;    // only short one counts
                  // bit for ring system
                  ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,3*k), extcon[j]*1013 + extcon[jj]*2003));
                  // bit for position
                  ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,5*k), extcon2[j]*2013 + extcon[jj]*1003));
                  ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,5*k), extcon[j]*2013 + extcon2[jj]*1003));
                  ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,7*k), extcon2[j]*3013 + extcon2[jj]*3003));
		  result += 4;
              }
          }
      }
// fprintf(stderr, "USE_SCAFFOLD_COLORS\n");
      if (which_bits & USE_SCAFFOLD_COLORS)
      {
         strcpy(prefix, "SC:");
          for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
          {
              if (j+1 == exclude_atom) continue;
              if (atom_status[j] <= 0) continue;
              tmp1 = 0;
              if (0 == strcmp(ap->atom_symbol, "C")) tmp1 = 101;
              else if (0 == strcmp(ap->atom_symbol, "O")) tmp1 = 301;
              else if (0 == strcmp(ap->atom_symbol, "N")) tmp1 = 401;
              else if (0 == strcmp(ap->atom_symbol, "S")) tmp1 = 601;
              else if (0 == strcmp(ap->atom_symbol, "P")) tmp1 = 701;
              else if (AtomSymbolMatch(ap->atom_symbol,"F,Cl,Br,I,At"))
                  tmp1 = 901;
              extcon[j] = ap->rsize_flags+tmp1;
          }
          for (i=0; i<32; i++)
          {
              for (j=0; j<mp->n_atoms; j++) extcon2[j] = 0;
              for (j=0; j<mp->n_atoms; j++)
              {
                  /* skip non-ring atoms */
                  if (atom_status[j] <= 0) continue;
                  extcon2[j] = atom_status[j]*3+(extcon[j]*0xF);
                  sum = 0; prod = 0;
                  for (jj=0; jj<nbp[j].n_ligands; jj++)
                  {
                      // only propagate through ring bonds
                      if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
                      sum += extcon[nbp[j].atoms[jj]];
                      prod *= 0xFF&(1+extcon[nbp[j].atoms[jj]]);
                  }
                  extcon2[j] = ((extcon2[j]+(sum*191)+prod)<<8) + ((extcon2[j]&0xFF0000)>>16);
                  extcon2[j] &= 0xFFFFFF;
              }
              for (j=0; j<mp->n_atoms; j++) extcon[j] = extcon2[j];
          }
          /* propagate smallest hash to all members of ring system */
          for (;;)
          {
              changed = FALSE;
              for (j=0; j<mp->n_atoms; j++)
              {
                  /* skip non-ring atoms */
                  if (atom_status[j] <= 0) continue;
                  for (jj=0; jj<nbp[j].n_ligands; jj++)
                  {
                      // only propagate through ring bonds
                      if (bond_status[nbp[j].bonds[jj]] <= 0) continue;
                      if (extcon[j] < extcon[nbp[j].atoms[jj]])
                      {
                          changed = TRUE;
                          extcon[nbp[j].atoms[jj]] = extcon[j];
                      }
                  }
              }
              if (!changed) break;
          }
          for (j=0; j<mp->n_atoms; j++)
          {
              if (j+1 == exclude_atom) continue;
              ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,4), extcon[j]*3013));
// ADD_BIT(fp_counts, ncounts, NEXT_SEED(NEXT_SEED(seed,4), extcon[j]*7013));
// ADD_BIT(fp_counts, ncounts, extcon[j]*16013 + 4*seed);
	      result += 1;
          }
      }
       MyFree((char *)extcon);
       MyFree((char *)extcon2);
   }

   if (which_bits & (USE_CLASS_SPIDERS | USE_FEATURE_PAIRS | USE_NON_SSS_BITS))
   {
       /* de-allocate length_matrix */
       MyFree((char *)length_matrix);
       MyFree((char *)length_tmp);
   }

   for (i=0; i<mp->n_bonds; i++)
      mp->bond_array[i].bond_type = old_bond_types[i];
   MyFree((char *)old_bond_types);
   MyFree((char *)touched_indices);
   MyFree((char *)nbp);
   MyFree((char *)H_count);
   MyFree((char *)atom_status);
   MyFree((char *)degree);
   MyFree((char *)cdegree);
   MyFree((char *)nspecial);
   MyFree((char *)unsaturated);
   MyFree((char *)bond_status);
// seed = 0;
// for (i=0; i<nbytes; i++)
//    seed = (0xFF&fingerprint[i]) ^ ((seed>>8) | ((seed&0xFF)<<16));
// fprintf(stderr, "2: %d bits set by FP %06lX, with %d paths\n",
//         CountBits(fingerprint, nbytes), seed, result);
// FORTIFY    fortifyTest(total_bytes_allocated, "SetFingerprintCountsWithFocus");
   return (result);
}

char *fp_prefixes[] = {
   "RING_PATTERN",
   "RING_PATH",
   "ATOM_SYMBOL_PATH",
   "ATOM_CLASS_PATH",
   "ATOM_COUNT",
   "AUGMENTED_ATOM",
   "HCOUNT_PATH",
   "HCOUNT_CLASS_PATH",
   "HCOUNT_PAIR",
   "BOND_PATH",
   "AUGMENTED_BOND",
   "RING_SIZE_COUNTS",
   "DEGREE_PATH",
   "CLASS_SPIDERS",
   "FEATURE_PAIRS",
   (char *)NULL
};

/*
 * Counts the bit multiplicity in fp_counts[0..ncounts-1] that correspond to the paths
 * with up to some defined size. which_bits is a set of flags that defines
 * which algorithms are triggered. as_query must be set to TRUE if the
 * fp_counts are to be computed for a query, which is mostly for hydrogen
 * counts.
 *
 * exclude_atom is the number of the atom that must not be touched by each
 * fragment. exclude_atom == 0 means all atoms.
 *
 * This function forwards to CountFingerprintPatterns().
 */
int SetFingerprintCountsWithFocus(struct reaccs_molecule_t *mp,
                                  int *fp_counts, int ncounts,
                                  int which_bits,
                                  int as_query,
                                  int fpflags,
                                  int exclude_atom)
{
    return CountFingerprintPatterns(mp, fp_counts, ncounts, which_bits, as_query, fpflags, exclude_atom, (struct string_counter_t *)NULL);
}

int CountBits(char *fingerprint, int nbytes)
{
   int result = 0;
   int i, j;

   for (i=0; i<nbytes; i++)
   {
      for (j=0; j<8; j++)
         if (fingerprint[i]&(1<<j)) result++;
   }
   return (result);
}

/**
 * Removes entries from matches when a previous match hits the same
 * set of atoms. This is a short-cut for symmetry equivalent matches,
 * which is not exact!
 */
ssmatch_t *PruneDuplicateAtomSets(ssmatch_t *matches)
{
    ssmatch_t *result = NULL;
    ssmatch_t *tmp;
    int match_found, i,j;
    while (!IsNULL(matches))
    {
        match_found = FALSE;
        for (tmp=result; !IsNULL(tmp); tmp=tmp->next)
        {
            for (i=0; i<matches->n_match; i++)
            {
                for (j=0; j<tmp->n_match; j++)
                    if (tmp->match_atoms[j] == matches->match_atoms[i]) break;
                if (j >= tmp->n_match) break;
            }
            if (i>=matches->n_match)
            {
                match_found = TRUE;
                break;
            }
        }
        if (match_found)
        {
            tmp = matches; matches = matches->next; tmp->next = NULL;
            FreeSSMatch(tmp);
        }
        else
        {
            tmp = matches->next; matches->next = result;
            result = matches; matches = tmp;
        }
    }
    return result; /* Currently a NOP */
}
