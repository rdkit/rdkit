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
/*    File:           perceive.c                                        */
/*                                                                      */
/*    Purpose:        This file implements perception routines to       */
/*                    classify the bonds and atoms of a chemical        */
/*                    structure (struct reaccs_molecule_t).             */
/*                                                                      */
/************************************************************************/

#include <stdio.h>
#include <string.h>

#include "graph.h"
#include "local.h"
#include "reaccs.h"
#include "set.h"
#include "utilities.h"

#include "perceive.h"

typedef unsigned pair_t[2];

void PerceiveRingBonds(struct reaccs_molecule_t *mp)
/* Sets the topography of the bonds of *mp to the appropriate ring
 * class and flags the ring size in the dummy field.
 */
{
   unsigned int i;
   bond_set_node *ring_list, *plist;
   /* static unsigned graph[255][2]; */
   pair_t *graph;

   graph = TypeAlloc(mp->n_bonds+1, pair_t);
   for (i=0; i<mp->n_bonds; i++)
   {
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
      mp->bond_array[i].topography = CHAIN;
      mp->bond_array[i].dummy      = 0;
   }
   ring_list = RingList(graph,mp->n_bonds);
   ring_list = CombineRings(ring_list);

   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (i=0; i<mp->n_bonds; i++)
         if (IsMember(plist->bond_set,i))
         {
            mp->bond_array[i].topography = RING;
            mp->bond_array[i].dummy     |=
               1<<MIN(8,Cardinality(plist->bond_set));
         }

   DisposeBondSetList(ring_list);

   MyFree((char *)graph);
}

static
void MarkRecursive(struct reaccs_molecule_t *mp,
                   int touched_atoms[],
                   int touched_bonds[],
                   int start_index,
                   int path_length,
                   int current_index,
                   int max_size,
                   neighbourhood_t nbp[])
/*
 * Recursively enumerates the rings in *mp and sets the rsize_flags bits
 * when a ring is found.
 */
{
   int i, j, ai, bi;

   for (i=0; i<nbp[current_index].n_ligands; i++)
   {
      ai = nbp[current_index].atoms[i];
      /* rings only count if started at least atom index */
      if (ai < start_index) continue;
      if (ai == start_index)                    /* ring found */
      {
         if (path_length < 3) continue;
         for (j=0; j<mp->n_atoms; j++)
            if (touched_atoms[j])
               mp->atom_array[j].rsize_flags |= (1<<path_length);
         for (j=0; j<mp->n_bonds; j++)
            if (touched_bonds[j])
               mp->bond_array[j].rsize_flags |= (1<<path_length);
         continue;
      }
      if (touched_atoms[ai]) continue;   /* don't walk backwards */
      /* continue recursion */
      if (path_length+1 > max_size)             /* don't go too far */
         continue;
      if (mp->atom_array[ai].rsize_flags == 0)  /* only ring atoms */
         continue;
      bi = nbp[current_index].bonds[i];
      if (mp->bond_array[bi].rsize_flags == 0)  /* only ring bonds */
         continue;
      touched_atoms[ai] = 1;  /* updating */
      touched_bonds[bi] = 1;
      MarkRecursive(mp, touched_atoms, touched_bonds,
                    start_index, path_length+1, ai, max_size, nbp);
      // touched_bonds[nbp[current_index].bonds[i]] = 0;
      touched_atoms[ai] = 0;  /* down-dating */
      touched_bonds[bi] = 0;
   }
}

void SetRingSizeFlags(struct reaccs_molecule_t *mp, int max_size,
                      neighbourhood_t nbp[])
/*
 * Sets the perceived ring-sizes as flag bits in the bond rsize_flags field.
 * It used a recursive enumeration algorithm pruned to only ring bonds/atoms.
 * nbp[] is the neighbourhood array to speed up processing.
 */
{
   int *atom_status, *bond_status;
   int i;
   int *touched_atoms; /* this array is up- and down-dated during recursion */
   int *touched_bonds; /* this array is up- and down-dated during recursion */

   /* label cyclic parts with unused rsize_flags bit */
   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);
   for (i=0; i<mp->n_atoms; i++)
      if (atom_status[i] > 0) mp->atom_array[i].rsize_flags = 1;
      else                    mp->atom_array[i].rsize_flags = 0;
   for (i=0; i<mp->n_bonds; i++)
      if (bond_status[i] > 0) mp->bond_array[i].rsize_flags = 1;
      else                    mp->bond_array[i].rsize_flags = 0;
   if (atom_status) MyFree((char *)atom_status);
   if (bond_status) MyFree((char *)bond_status);

   touched_atoms = TypeAlloc(mp->n_atoms, int);
   touched_bonds = TypeAlloc(mp->n_bonds, int);
   for (i=0; i<mp->n_atoms; i++)
   {
      if (mp->atom_array[i].rsize_flags == 0) continue;
      touched_atoms[i]=1;     /* updating */
      MarkRecursive(mp, touched_atoms, touched_bonds, i, 1, i, max_size, nbp);
      touched_atoms[i]=0;     /* down-dating */
   }
   MyFree((char *)touched_atoms);
   MyFree((char *)touched_bonds);
}

void PerceiveRingSizes(struct reaccs_molecule_t *mp)
/*
 * Sets the perceived ring-sizes as flag bits in the bond rsize_flags field.
 */
{
   unsigned int i, ai1, ai2;
   int changes, new_flag;
   bond_set_node *ring_pairs, *ring_list, *plist;
   pair_t *graph;

   if (IsNULL(mp)) return;

   graph = TypeAlloc(mp->n_bonds+1, pair_t);
   for (i=0; i<mp->n_bonds; i++)
   {
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
      mp->bond_array[i].rsize_flags      = 0;
   }
   ring_list = RingList(graph,mp->n_bonds);
   ring_list = CombineRings(ring_list);

   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (i=0; i<mp->n_bonds; i++)
         if (IsMember(plist->bond_set,i))
         {
            mp->bond_array[i].rsize_flags     |=
               1<<MIN(10,Cardinality(plist->bond_set));
         }

recombine:
   changes = 0;
   ring_pairs = ProperRingPairs(ring_list, mp->n_atoms, graph);
   /* add the additional rings for the analysis */
   while (ring_pairs!=(bond_set_node *)NULL)
   {
      for (plist = ring_list; plist != (bond_set_node *)NULL; plist=plist->next)
         if (0 == CompareSets(ring_pairs->bond_set,plist->bond_set))
            break;
      if (plist != (bond_set_node *)NULL)
      {
         plist = ring_pairs;
         ring_pairs = ring_pairs->next;
         plist->next = (bond_set_node *)NULL;
         DisposeBondSetList(plist);
      }
      else
      {
         plist = ring_pairs->next;
         ring_pairs->next = ring_list;
         ring_list=ring_pairs;
         ring_pairs = plist;
      }
   }

   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (i=0; i<mp->n_bonds; i++)
         if (IsMember(plist->bond_set,i))
         {
            new_flag = (1<<MIN(10,Cardinality(plist->bond_set)));
            if (0 == (mp->bond_array[i].rsize_flags & new_flag))
            {
               mp->bond_array[i].rsize_flags |= new_flag;
               changes++;
            }
         }
   if (changes > 0) goto recombine;

   for (i=0; i<mp->n_atoms; i++)
      mp->atom_array[i].rsize_flags = 0;
   for (i=0; i<mp->n_bonds; i++)
   {
      ai1 = mp->bond_array[i].atoms[0]-1;
      ai2 = mp->bond_array[i].atoms[1]-1;
      mp->atom_array[ai1].rsize_flags |= mp->bond_array[i].rsize_flags;
      mp->atom_array[ai2].rsize_flags |= mp->bond_array[i].rsize_flags;
   }

   DisposeBondSetList(ring_list);

   MyFree((char *)graph);
}

void PerceiveAromaticBonds(struct reaccs_molecule_t *mp)
/*
 * Sets bond types of the bonds of *mp to "AROMATIC" if in
 * six-ring of sp2 atoms.
 */
{
   unsigned i;
   int ii;
   int changed;
   int ndouble, nsingle, naromatic;
   bond_set_node *ring_list, *ring_pairs, *plist, *plisth;
   pair_t *graph;
   int *bond_is_in_ring;
   int *sp_count;
   int is_cumulene;
   int ring_is_aromatic;

   graph = TypeAlloc(mp->n_bonds+1, pair_t);
   sp_count = TypeAlloc(mp->n_atoms+1, int);
   bond_is_in_ring = TypeAlloc(mp->n_bonds, int);

   for (i=0; i<mp->n_bonds; i++)
   {
      bond_is_in_ring[i] = FALSE;
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
   }
   ring_list = RingList(graph,mp->n_bonds);
   ring_list = CombineRings(ring_list);

   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
         bond_is_in_ring[ii] = TRUE;

   /* add the additional rings for the analysis */
   ring_pairs = ProperRingPairs(ring_list, mp->n_atoms, graph);
   while (ring_pairs!=(bond_set_node *)NULL)
   {
      plisth = ring_pairs->next;
      ring_pairs->next = ring_list;
      ring_list=ring_pairs;
      ring_pairs = plisth;
   }

   do
   {
      changed = FALSE;
      for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      {
         for (i=0; i<=mp->n_atoms; i++) sp_count[i]=0;
         ndouble = nsingle = naromatic = 0;
         for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
         {
            i = (unsigned)ii;
               if (mp->bond_array[i].bond_type == SINGLE)        nsingle++;
               // else if (mp->bond_array[i].bond_type == DOUBLE)   ndouble++;
               else if (mp->bond_array[i].bond_type == AROMATIC) naromatic++;
         }

         /* count number of double bonds per ring atom and mark good ones as sp2 */
         is_cumulene = FALSE;
         for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
         {
            i = (unsigned)ii;
            if (mp->bond_array[i].bond_type == DOUBLE)
            {
               ndouble++;
               sp_count[mp->bond_array[i].atoms[0]]++;
               if (sp_count[mp->bond_array[i].atoms[0]] > 1) is_cumulene = TRUE;
               sp_count[mp->bond_array[i].atoms[1]]++;
               if (sp_count[mp->bond_array[i].atoms[1]] > 1) is_cumulene = TRUE;
            }
            else if (mp->bond_array[i].bond_type == TRIPLE)
            {
               is_cumulene = TRUE;
            }
         }
         // mark proven aromatic atoms as sp2
         for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
         {
            i = (unsigned)ii;
            if (mp->bond_array[i].bond_type == AROMATIC  &&  IsMember(plist->bond_set,i))
            {
               if (sp_count[mp->bond_array[i].atoms[0]] == 0) sp_count[mp->bond_array[i].atoms[0]]++;
               if (sp_count[mp->bond_array[i].atoms[1]] == 0) sp_count[mp->bond_array[i].atoms[1]]++;
            }
         }
         // check if all ring atoms are sp2 */
         ring_is_aromatic =
            !is_cumulene  &&                    // cumulenes are not good
            ((plist->cardinality-2)%4) == 0;    // Hueckle aromatic ring sizes
         for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
         {
            i = (unsigned)ii;
            /* all bonds must connect sp2 atoms */
            if (sp_count[mp->bond_array[i].atoms[0]] != 1) ring_is_aromatic = FALSE;
            if (sp_count[mp->bond_array[i].atoms[1]] != 1) ring_is_aromatic = FALSE;
         }

         if (ring_is_aromatic  &&  (ndouble > 0  ||  nsingle > 0))
         {
            for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
            {
               i = (unsigned)ii;
               if (bond_is_in_ring[i]  &&  mp->bond_array[i].bond_type != AROMATIC)
               {
                  changed = TRUE;
                  mp->bond_array[i].bond_type = AROMATIC;
               }
            }
         }
      }
   } while (changed);

   DisposeBondSetList(ring_list);

   MyFree((char *)bond_is_in_ring);
   MyFree((char *)sp_count);
   MyFree((char *)graph);
}

/**
 * Returns a pointer to the symbol list for atom index or NULL if there is none.
 */
struct symbol_list_t *getList(struct reaccs_molecule_t *mp, int index)
{
   struct symbol_list_t *slp;

   for (slp = mp->symbol_lists; !IsNULL(slp); slp=slp->next)
   {
      if (slp->atom == index+1) break;
   }
   return (slp);
}

void PerceiveDYAromaticity(struct reaccs_molecule_t *mp,
			   neighbourhood_t nbp[])
/*
 * Converts bonds in rings that are perceived as aromatic by Daylight
 * programs into AROMATIC bonds. In addition, hydrogen counts of aromatic
 * non-carbon atoms are remembered in the query_H_count fields.
 */
{
   unsigned i, j;
   int ii;
   int changed;
   bond_set_node *ring_list, *ring_pairs, *plist, *plisth;
   bond_set_node *p, *aromatic_candidates;
   bit_set_t *set;
   pair_t *graph;
//    unsigned graph[MAXBONDS][2];
   int *bond_is_in_ring;
//   int bond_is_in_ring[MAXBONDS];
   int *atom_is_in_ring;
//   int atom_is_in_ring[MAXATOMS];
   int *aromaticity_candidate;
//   int aromaticity_candidate[MAXATOMS];
   int in_ring_double;
   int in_ring_aromatic;
   int exo_pull;
   int conjugated;
   int is_in_ring;
   int local_pi, npi;

   struct reaccs_atom_t *ap, *alp;
   struct reaccs_bond_t *blp;
   neighbourhood_t *nbph;

//   if (mp->n_bonds > MAXBONDS) return;
//   if (mp->n_atoms > MAXATOMS) return;

   bond_is_in_ring = TypeAlloc(mp->n_bonds, int);
   graph = TypeAlloc(mp->n_bonds, pair_t);
   /* Process all rings */
   for (i=0; i<mp->n_bonds; i++)
   {
      bond_is_in_ring[i] = FALSE;
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
   }
   ring_list = RingList(graph,mp->n_bonds);
   if (IsNULL(ring_list))
   {
      MyFree((char *)bond_is_in_ring);
      MyFree((char *)graph);
      return;
   }

   atom_is_in_ring = TypeAlloc(mp->n_atoms, int);
   aromaticity_candidate = TypeAlloc(mp->n_atoms, int);
   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (ii=NextMember(plist->bond_set, 0); ii>=0; ii=NextMember(plist->bond_set, (unsigned)ii+1))
         bond_is_in_ring[ii] = TRUE;

   for (i=0; i<mp->n_atoms; i++)
   {
      atom_is_in_ring[i] = FALSE;
      if (0 != strcmp("C", mp->atom_array[i].atom_symbol))
         aromaticity_candidate[i] = TRUE;
      else
         aromaticity_candidate[i] = FALSE;
   }
   for (i=0; i<mp->n_bonds; i++)
      if (mp->bond_array[i].bond_type >  SINGLE  &&
          mp->bond_array[i].bond_type != TRIPLE)
      {
         aromaticity_candidate[mp->bond_array[i].atoms[0]-1] = TRUE;
         aromaticity_candidate[mp->bond_array[i].atoms[1]-1] = TRUE;
      }
   /* Process just aromaticity candidates */
   DisposeBondSetList(ring_list);
   for (i=0; i<mp->n_bonds; i++)
   {
      graph[i][0] = mp->bond_array[i].atoms[0];
      graph[i][1] = mp->bond_array[i].atoms[1];
      if (!bond_is_in_ring[i])  /* only ring bonds can be aromatic */
      {
         graph[i][0] = 0;
         graph[i][1] = 0;
      }
      else
      {
         atom_is_in_ring[graph[i][0]-1] = TRUE;
         atom_is_in_ring[graph[i][1]-1] = TRUE;
         /* make ring perception ignore bonds with non-candidates */
         if (!aromaticity_candidate[graph[i][0]-1]) graph[i][0] = 0;
         if (!aromaticity_candidate[graph[i][1]-1]) graph[i][1] = 0;
      }
   }
   ring_list = RingList(graph,mp->n_bonds);
   if (IsNULL(ring_list))
   {
      MyFree((char *)bond_is_in_ring);
      MyFree((char *)atom_is_in_ring);
      MyFree((char *)aromaticity_candidate);
      MyFree((char *)graph);
      return;
   }
   ring_list = CombineRings(ring_list);

   ring_pairs = ProperRingPairs(ring_list, mp->n_atoms, graph);

   /* add the additional rings for the analysis */
   while (ring_pairs!=(bond_set_node *)NULL)
   {
      plisth = ring_pairs->next;
      ring_pairs->next = ring_list;
      ring_list=ring_pairs;
      ring_pairs = plisth;
   }

   /* add fused pairs that could be also aromatic */
   set = NewSet(MaxMember(ring_list->bond_set));
   aromatic_candidates = ring_list;
   for (plist=ring_list; plist!=(bond_set_node *)NULL; plist=plist->next)
      for (plisth=plist->next;
           plisth!=(bond_set_node *)NULL;
	   plisth=plisth->next)
      {
         // check first if there is any overlap
         if (IntersectionIsEmpty(plist->bond_set,plisth->bond_set)) continue;;
         /* assumes allocated size is large enough */
         set = SetExclusiveUnion(CopySet(set,plist->bond_set),plisth->bond_set);
	 /* only fused rings are considered */
	 if (Cardinality(set) == plist->cardinality +
	                         plisth->cardinality - 2)
	 {
	    p = NewBondSetNode(MaxMember(plist->bond_set));
	    p->next = aromatic_candidates; 
	    aromatic_candidates = p;
	    CopySet(p->bond_set, set);
            p->cardinality = Cardinality(set);
	 }
      }
   DisposeSet(set);

   do
   {
      changed = FALSE;
      for (plist=aromatic_candidates;
           plist!=(bond_set_node *)NULL;
	   plist=plist->next)
      {
	 npi = 0;
	 conjugated = TRUE;
         /* Look for atoms that are members of this ring */
	 for (i=0, ap=mp->atom_array, nbph=nbp;
	      i<mp->n_atoms;
	      i++, ap++, nbph++)
	 {
            if (!atom_is_in_ring[i]) continue;
            /* Counts # of double bonds in this ring attached to ap */
	    in_ring_double = 0;
            /* Counts # of aromatic bonds in this ring attached to ap */
	    in_ring_aromatic = 0;
            /* becomes true if there is an exocyclic DB that pulls electrons */
	    exo_pull = FALSE;
	    is_in_ring = FALSE;
	    local_pi = 0;       /* count pi electrons of this atom */
	    for (j=0; j<nbph->n_ligands; j++)
	    {
	       alp = &mp->atom_array[nbph->atoms[j]];
	       blp = &mp->bond_array[nbph->bonds[j]];
	       if (bond_is_in_ring[nbph->bonds[j]]  &&  IsMember(plist->bond_set, nbph->bonds[j]))
	       {
		  is_in_ring = TRUE;
		  if (blp->bond_type == AROMATIC) in_ring_aromatic++;
		  if (blp->bond_type == DOUBLE)   in_ring_double++;
	       }
	       else if (mp->bond_array[nbph->bonds[j]].bond_type == DOUBLE)
	       {
		  if (ap->atom_symbol[0] == 'C'  &&  ap->atom_symbol[1] == '\0'  &&
                      // AtomSymbolMatch(ap->atom_symbol, "C")  &&
                      !bond_is_in_ring[nbph->bonds[j]]  && /* patches DY bug */
		      AtomSymbolMatch(alp->atom_symbol, "O,S,P,N,L"))
		     exo_pull = TRUE;
	       }
	    }
	    if (!is_in_ring) continue;  /* no ligands in ring found */
	    if ((in_ring_aromatic >= 1  ||  in_ring_double  == 1)  &&
	        (AtomSymbolMatch(ap->atom_symbol, "C,N,A,*") ||
		 0 == strcmp(ap->atom_symbol, "L"))) /* lists may be aromatic */
	    {
               /* Count this <b>atom</b> for one electron */
	       local_pi = 1;
	    }
	    else if ((in_ring_aromatic == 0  &&  in_ring_double  == 0)  &&
		     ap->charge == 0                                    &&
		     (AtomSymbolMatch(ap->atom_symbol, "N,S,O") ||
                     /* term that catches lists with O,N,S */
                      AtomSymbolMatch(ap->atom_symbol, "L")  &&
                      getList(mp,i) != NULL  &&
                      getList(mp,i)->logic  &&  /* no not-lists */
                      strchr(getList(mp, i)->string,'C') == NULL))
	    {
	       local_pi = 2;
	       if (nbph->n_ligands  == 2  &&
	           in_ring_aromatic == 0  &&
	           in_ring_double   == 0  &&
	           0 == strcmp(ap->atom_symbol, "N"))
		  ap->query_H_count = 1+1;      /* offset 1! */
	    }
	    else if ((in_ring_aromatic == 0  &&  in_ring_double  == 0)  &&
		     ap->charge == 0                                    &&
		     exo_pull						&&
		     0 == strcmp(ap->atom_symbol, "C"))
	    {
	       local_pi = 0;
	    }
	    else
            {
	       conjugated = FALSE;
            }

	    /* catch special cases that Daylight don't consider aromatic */
	    if (ap->charge < 0  &&  AtomSymbolMatch(ap->atom_symbol, "C,N"))
	       conjugated = FALSE;

	    npi += local_pi;
	 }
	 if (conjugated == FALSE) continue;	/* system not conjugated */
	 if (npi%4 != 2)          continue;	/* not a Huckel system */

	 /* now we have an aromatic ring => modify bond_types */
         for (i=0; i<mp->n_bonds; i++)
            if (bond_is_in_ring[i]  &&  IsMember(plist->bond_set,i))
	    {
	       if (mp->bond_array[i].bond_type != AROMATIC)
	       {
		  mp->bond_array[i].bond_type = AROMATIC;
		  changed = TRUE;
	       }
	    }
      }
   } while (changed);

   DisposeBondSetList(aromatic_candidates);
   // add ring bonds between aromatic atoms as aromatic
   for (i=0; i<mp->n_atoms; i++)
      aromaticity_candidate[i] = FALSE;
   for (i=0, blp=mp->bond_array; i<mp->n_bonds; i++, blp++)
   {
      if (blp->bond_type == AROMATIC)
      {
	 aromaticity_candidate[blp->atoms[0]-1] = TRUE;
	 aromaticity_candidate[blp->atoms[1]-1] = TRUE;
      }
   }
   for (i=0, blp=mp->bond_array; i<mp->n_bonds; i++, blp++)
      if (aromaticity_candidate[blp->atoms[0]-1]  &&
          aromaticity_candidate[blp->atoms[1]-1]  &&
	  bond_is_in_ring[i]                      &&
	  blp->bond_type == SINGLE)
	 blp->bond_type = AROMATIC;

   MyFree((char *)atom_is_in_ring);
   MyFree((char *)aromaticity_candidate);
   MyFree((char *)bond_is_in_ring);
   MyFree((char *)graph);
}

void PerceiveMarkush(struct reaccs_molecule_t *mp,
		     neighbourhood_t nbp[])
/*
 * Scans mp for "R" atoms. If there are any, it sets the sub_desc entries
 * of all atoms not attached to an R to s*. This is used to simplify
 * substructure query definition by only specifying free sites.
 */
{
   int i, j;
   struct reaccs_atom_t *ap, *alp;
   neighbourhood_t *nbph;

   for (i=0, ap=mp->atom_array, nbph=nbp;
	i<mp->n_atoms;
	i++, ap++, nbph++)
      if (0 == strcmp(ap->atom_symbol, "R")  &&  nbph->n_ligands == 1) break;

   /* Return if there is no R atom */
   if (i == mp->n_atoms) return;

   for (i=0, ap=mp->atom_array, nbph=nbp;
	i<mp->n_atoms;
	i++, ap++, nbph++)
   {
      if (0 == strcmp(ap->atom_symbol, "R")  &&  nbph->n_ligands == 1)
      {
	 mp->bond_array[nbph->bonds[0]].bond_type = NONE;
	 ap->color = (-1);	/* Note: this is a trick that works,     */
         continue;		/* because SMARTS are assumed not to     */
      }				/* have valid coordinates or attributes. */

      for (j=0; j<nbph->n_ligands; j++)
      {
	 alp = &mp->atom_array[nbph->atoms[j]];
	 if (strcmp(alp->atom_symbol, "R") == 0) break;
      }
      if (j == nbph->n_ligands  &&
          ap->sub_desc == 0)	/* no R attachment => set substitution count */
	 ap->sub_desc = (-2);
   }

   /* Now, we remove the R-atoms from the connection table */
   while (1)
   {
      for (i=0, ap=mp->atom_array, nbph=nbp;
           i<mp->n_atoms;
           i++, ap++, nbph++)
         if (0 == strcmp(ap->atom_symbol, "R")  &&  nbph->n_ligands == 1) break;

      /* Return if there is no R atom */
      if (i == mp->n_atoms)
      {
         return;
      }
      else
      {
         RemoveAtomFromMolecule(mp, i+1);
         SetupNeighbourhood(mp,nbp,mp->n_atoms);
      }
   }

}
