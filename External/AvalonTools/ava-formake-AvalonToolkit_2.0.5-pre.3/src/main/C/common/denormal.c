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
/*                                                                    */
/*    File:           denormal.c                                      */
/*                                                                    */
/*    Purpose:        This file implements denormalization of         */
/*                    connection tables. It changes the bonds marked  */
/*                    as AROMATIC to either single or double bonds.   */
/*                                                                    */
//    History:        25-Apr-1994     Start of development.           */
//                                                                    */
/************************************************************************/

#include <string.h>

#include "casutils.h"
#include "graph.h"
#include "reaccsio.h"
#include "utilities.h"

#include "denormal.h"

static void ColorNormalizedComponents(struct reaccs_molecule_t *mp)
/*
 * 'Disects' the molecule into components connected by AROMATIC bonds.
 * Different components are marked by different colors.
 */
{
   int i, nmerge;
   struct reaccs_bond_t *bp;

   for (i=0; i<mp->n_atoms; i++)       /* all atoms get a different color */
      mp->atom_array[i].color = i+1;

   do  /* merge colors of atoms connected by AROMATIC bonds */
   {
      nmerge = 0;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->bond_type == AROMATIC  &&
           mp->atom_array[bp->atoms[0]-1].color !=
         mp->atom_array[bp->atoms[1]-1].color)
       {
          nmerge ++;
      if (mp->atom_array[bp->atoms[0]-1].color <
           mp->atom_array[bp->atoms[1]-1].color)
         mp->atom_array[bp->atoms[1]-1].color =
             mp->atom_array[bp->atoms[0]-1].color;
     else
               mp->atom_array[bp->atoms[0]-1].color =
             mp->atom_array[bp->atoms[1]-1].color;
  }
   } while (nmerge > 0);
}

int CCTDenormalize(struct reaccs_molecule_t *mp)
/*
 * Denormalizes the molecule *mp which was the result of decoding
 * a CAS CCT entry. It contains the valence and hydrogen range
 * information in the dummy1 and dummy2 fields of the atom structures.
 *
 * The function returns FALSE if anything went wrong in denormalization,
 */
{
   int result;
   int total_order;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int i, j;
   struct atom_constraint_t *acp;
   struct bond_constraint_t *bcp;
   neighbourhood_t *nbp;
   int valence, hco_min, hco_max, charge;
   int nsingle, ndouble, ntriple, naromatic, nbonds;

   if (mp->n_bonds <= 0) return (TRUE);
   acp = TypeAlloc(mp->n_atoms, struct atom_constraint_t);
   bcp = TypeAlloc(mp->n_bonds, struct bond_constraint_t);

   nbp   = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);

   // convert 'AROMATIC' bonds between shortcuts to NONE bonds
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
       if (bp->bond_type != ANY_BOND) continue;
       if (mp->atom_array[bp->atoms[0]-1].atext[0] != '\0') bp->bond_type = NONE;
       if (mp->atom_array[bp->atoms[1]-1].atext[0] != '\0') bp->bond_type = NONE;
   }
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (ap->dummy1 == 0)
         valence = CASValence(ap->atom_symbol);
      else
         valence = ap->dummy1;
// if (valence != CASValence(ap->atom_symbol))
// {
// fprintf(stderr,"DENORMAL>%s|%d: hco_min=%d, hco_max=%d, valence=%d\n",
// mp->name, i+1, HCO_MIN(ap->dummy2), HCO_MAX(ap->dummy2), valence);
// }

      if (HCO_MIN(ap->dummy2) >= HCO_MAX(ap->dummy2)) /* no tautomer */
      {
         hco_min = hco_max = HCO_MIN(ap->dummy2);
      }
      else
         if (valence > 4  &&
             (strcmp(ap->atom_symbol,"S") == 0  ||
//              /* strcmp(ap->atom_symbol,"N") == 0  || */
              strcmp(ap->atom_symbol,"I") == 0  ||
              FALSE))
         {
            hco_min = hco_max = HCO_MIN(ap->dummy2);
// fprintf(stderr,"1: %s(%d): hco_min=%d, hco_max=%d\n",ap->atom_symbol,i+1,hco_min,hco_max);
         }
         else
         {
            hco_min = HCO_MIN(ap->dummy2);
            hco_max = HCO_MAX(ap->dummy2);
// fprintf(stderr,"2: %s(%d): hco_min=%d, hco_max=%d\n",ap->atom_symbol,i+1,hco_min,hco_max);
         }
      nsingle = ndouble = ntriple = naromatic = 0;
      for (j=0; j<nbp[i].n_ligands; j++)
         switch (mp->bond_array[nbp[i].bonds[j]].bond_type)
         {
            case ANY_BOND: break;
            case SINGLE: nsingle++; break;
            case DOUBLE: ndouble++; break;
            case TRIPLE: ntriple++; break;
            default:     naromatic++; break;
         }
      charge = ap->charge;
      if (ap->dummy2 & PLUS_TAUTOMER)  charge++;
      if (ap->dummy2 & MINUS_TAUTOMER) charge--;
      if (valence == 0)             /* estimating valence */
      {
	 valence = (nsingle*2 + ndouble*4 + ntriple*6 + naromatic*3)/2;
	 valence += hco_max;
	 valence += ABS(charge);
      }
      else
        if (ap->dummy1 == 0  &&  ap->dummy2 == 0)  /* estimating H-count */
	{
           // hco_min =
           // hco_max = valence -
                     // (nsingle*2 + ndouble*4 + ntriple*6 + naromatic*3)/2 -
                     // /*ABS*/(charge);
           hco_min = hco_max = ImplicitHydrogens(ap->atom_symbol, nsingle, naromatic, ndouble, ntriple, ap->radical == DOUBLET, charge);
           if (nsingle==0  &&  ndouble==0  &&  ntriple==0  && naromatic==2  &&  hco_max==0  &&
               strcmp("N",ap->atom_symbol)==0)   // catch possibly protonated nitrogen case
           {
              hco_max = 1;
           }
// fprintf(stderr,"3: %s(%d): hco_min=%d, hco_max=%d\n",ap->atom_symbol,i+1,hco_min,hco_max);
	}

      nbonds = nbp[i].n_ligands;
      acp[i].open = FALSE;
      acp[i].is_oxigen = 0 == strcmp(ap->atom_symbol,"O");
      acp[i].dbl_min = valence - ABS(charge) - hco_max - nbonds - 2*ntriple;
      acp[i].dbl_max = valence - ABS(charge) - hco_min - nbonds - 2*ntriple;
      if (acp[i].dbl_min < 0)   /* Patch for tautomeric charges */
      {
          acp[i].dbl_max -= acp[i].dbl_min;
          acp[i].dbl_min = 0;
      }
// if (hco_min != hco_max  ||  valence != CASValence(ap->atom_symbol))
// {
// fprintf(stderr,"%s|%d: hco_min=%d, hco_max=%d, valence=%d, charge=%d, dummy2=%d\n",
// mp->name, i+1, hco_min, hco_max, valence, charge, ap->dummy2);
// fprintf(stderr,
// "%d(%s): dbl_min = %d, dbl_max = %d, hco_min = %d, hco_max = %d, valence = %d\n",
// i+1, ap->atom_symbol,
// acp[i].dbl_min, acp[i].dbl_max,
// hco_min, hco_max, valence);
// }
   }

// for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
// fprintf(stderr,"'%s'(%d):dbl_min=%d, dbl_max=%d\n",
// ap->atom_symbol, i+1, acp[i].dbl_min, acp[i].dbl_max);

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      bcp[i].open = FALSE;
      bcp[i].can_be_single = bp->bond_type == SINGLE  ||
                             bp->bond_type == AROMATIC;
      if (acp[bp->atoms[0]-1].dbl_max == 0  ||  acp[bp->atoms[1]-1].dbl_max == 0)
         bcp[i].can_be_double = FALSE;
      else
         bcp[i].can_be_double = bp->bond_type == DOUBLE  || bp->bond_type == AROMATIC;
      bcp[i].scount = nbp[bp->atoms[0]-1].n_ligands +
                      nbp[bp->atoms[1]-1].n_ligands - 2;
// fprintf(stderr, "(%d-%d): can_be_single=%d, can_be_double=%d, scount=%d\n",
// bp->atoms[0], bp->atoms[1],
// bcp[i].can_be_single, bcp[i].can_be_double, bcp[i].scount);
   }

// PrintREACCSMolecule(stderr,mp,"Before Denormalize()");
   result = Denormalize(mp, acp, bcp);
// PrintREACCSMolecule(stderr,mp,"After Denormalize()");

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      ap->dummy2 &= (PLUS_TAUTOMER|MINUS_TAUTOMER);
      /* Patch back tautomeric charges */
      if (ap->dummy2 != NONE)
      {
         total_order = 0;
         for (j=0; j<nbp[i].n_ligands; j++)
            total_order += mp->bond_array[nbp[i].bonds[j]].bond_type;
         if (total_order == 4  &&
             ap->dummy2&PLUS_TAUTOMER &&
             0==strcmp(ap->atom_symbol,"N"))
            ap->charge += 1;
         if (total_order == 1  &&
             ap->dummy2&MINUS_TAUTOMER &&
             0==strcmp(ap->atom_symbol,"O"))
            ap->charge -= 1;
      }
      ap->color = NONE;
   }

   MyFree((char *)nbp);

   if (bcp) MyFree((char *)bcp);
   MyFree((char *)acp);

   return (result);
}

typedef unsigned atom_pair[2];

static bond_set_node *MoleculeRings(struct reaccs_molecule_t *mp)
/*
 * Returns the list of basis rings of the molecule *mp.
 * Ring with bond_type flag RUBBER_BOND are ignored.
 */
{
   int i;
   atom_pair *bonds;
   bond_set_node *ring_list;

   if (mp->n_bonds == 0) return (bond_set_node *)NULL;

   bonds = TypeAlloc(mp->n_bonds, atom_pair);
   for (i=0; i<mp->n_bonds; i++)
   {
      bonds[i][0] = mp->bond_array[i].atoms[0];
      bonds[i][1] = mp->bond_array[i].atoms[1];
   }

   ring_list = RingList(bonds,mp->n_bonds);
   ring_list = CombineRings(ring_list);

   MyFree((char *)bonds);
   return (ring_list);
}

#define VIOLATION      (-1)

static int Refine(struct atom_constraint_t *acp,
                  struct bond_constraint_t *bcp,
                  struct reaccs_molecule_t *mp,
                  neighbourhood_t          *nbp)
/*
 * Refines the candidate bond and atom constraints bcp[] and acp[] within
 * the molecule context *mp/nbp[]. It returns the number of not yet
 * assigned bonds in the candidate set.
 */
{
   int changed;
   int result;
   int i, j;
   int dbl_min, dbl_max;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;

// for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
// if (bcp[j].open & CANDIDATE) fprintf(stderr, "IN: bond(%d-%d): %s/%s\n", bp->atoms[0], bp->atoms[1], (bcp[j].can_be_single?"SINGLE":""), (bcp[j].can_be_double?"DOUBLE":""));
// for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
// if (acp[i].open & CANDIDATE) fprintf(stderr, "IN: %s(%d) acp.dbl_min=%d, acp.dbl_max=%d\n", ap->atom_symbol, i+1, acp[i].dbl_min, acp[i].dbl_max);
   do
   {
      changed = FALSE;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++) /* refine atoms */
         if (acp[i].open & CANDIDATE  && acp[i].open & TOUCHED)
         {
// fprintf(stderr, "refining '%s' atom %d\n", ap->atom_symbol, i+1);
            acp[i].open = CANDIDATE;
            dbl_min = dbl_max = 0;
            for (j=0; j<nbp[i].n_ligands; j++)
               if (bcp[nbp[i].bonds[j]].can_be_single &&
                   bcp[nbp[i].bonds[j]].can_be_double)
                  dbl_max++;
               else if (bcp[nbp[i].bonds[j]].can_be_single &&
                        !bcp[nbp[i].bonds[j]].can_be_double)
                  /* NOP */;
               else if (!bcp[nbp[i].bonds[j]].can_be_single &&
                        bcp[nbp[i].bonds[j]].can_be_double)
               {
                  dbl_min++;
                  dbl_max++;
               }
               else     /* constraint violation */
               {
// fprintf(stderr, "VIOLATION: dbl_min = %d, dbl_max = %d\n", dbl_min, dbl_max);
                  return (VIOLATION);
               }
// fprintf(stderr, "dbl_min = %d, dbl_max = %d\n", dbl_min, dbl_max);

            if (dbl_min == acp[i].dbl_max)      /* all double bonds assigned */
            {                                   /* other bonds are single */
// fprintf(stderr, "REFINE: %s(%d) dbl_min=%d, dbl_max=%d\n", ap->atom_symbol, i+1, acp[i].dbl_min, acp[i].dbl_max);
               for (j=0; j<nbp[i].n_ligands; j++)
                  if (bcp[nbp[i].bonds[j]].can_be_single &&
                      bcp[nbp[i].bonds[j]].can_be_double)
                  {
                     changed = TRUE;
                     bp=mp->bond_array+nbp[i].bonds[j];
// fprintf(stderr,"setting bond %d-%d to SINGLE\n",bp->atoms[0],bp->atoms[1]);
                     bcp[nbp[i].bonds[j]].can_be_double = FALSE;
                     bcp[nbp[i].bonds[j]].open |= TOUCHED;
                     acp[nbp[i].atoms[j]].open |= TOUCHED;
                  }
               dbl_max = dbl_min;
            }
            else if (dbl_max == acp[i].dbl_min) /* all single bonds assigned */
            {                                   /* other bonds are double */
// fprintf(stderr, "REFINE: %s(%d) dbl_min=%d, dbl_max=%d\n", ap->atom_symbol, i+1, acp[i].dbl_min, acp[i].dbl_max);
               for (j=0; j<nbp[i].n_ligands; j++)
                  if (bcp[nbp[i].bonds[j]].can_be_single &&
                      bcp[nbp[i].bonds[j]].can_be_double)
                  {
                     changed = TRUE;
                     bp=mp->bond_array+nbp[i].bonds[j];
// fprintf(stderr,"setting bond %d-%d to DOUBLE\n",bp->atoms[0],bp->atoms[1]);
                     bcp[nbp[i].bonds[j]].can_be_single = FALSE;
                     bcp[nbp[i].bonds[j]].open |= TOUCHED;
                     acp[nbp[i].atoms[j]].open |= TOUCHED;
                  }
               dbl_min = dbl_max;
            }
            else if (dbl_max < acp[i].dbl_min  ||
                     dbl_min > acp[i].dbl_max)
            {
// fprintf(stderr, "VIOLATION: atom(%d): dbl_min = %d, dbl_max = %d, acp[%d].dbl_min = %d, acp[%d].dbl_max = %d\n", i+1, dbl_min, dbl_max, i, acp[i].dbl_min, i, acp[i].dbl_max);
               return (VIOLATION);
            }
         }
                                                      /* refine bonds */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bcp[i].open & CANDIDATE  && bcp[i].open & TOUCHED)
         {
            bcp[i].open = CANDIDATE;
         }
   } while (changed);

   result = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bcp[i].open & CANDIDATE  &&  bcp[i].can_be_single == bcp[i].can_be_double)
      {
//         fprintf(stderr,"bond %d-%d has can_be_single=%d and can_be_double=%d\n", bp->atoms[0], bp->atoms[1], bcp[i].can_be_single, bcp[i].can_be_double);
         result++;
      }

   return (result);
}

static
int CompareDenormalizations(bond_set_node            *sol1,
                            bond_set_node            *sol2,
                            struct atom_constraint_t *acp,
                            struct bond_constraint_t *bcp,
                            struct reaccs_molecule_t *mp,
                            bond_set_node            *ring_list)
/*
 * Compares the denormalizations *sol1 and *sol2 accoding to the
 * some heuristic priority criteria:
 *
 * It returns 1 if sol1 is better than sol2, 0 if both are identical,
 * and (-1) if sol2 is better than sol1.
 *
 * Note: This is currently a dummy implementation
 */
{
   // DUMMY sol1 = sol1; sol2 = sol2;
   // DUMMY acp = acp; bcp = bcp; mp = mp; ring_list = ring_list;
   return (1);
}

static
bond_set_node *AddNewSolution(bond_set_node            *new_sol,
                              bond_set_node            *old_list,
                              struct atom_constraint_t *acp,
                              struct bond_constraint_t *bcp,
                              struct reaccs_molecule_t *mp,
                              bond_set_node            *ring_list)
/*
 * Adds the denormalization solution *new_sol to the list *old_list.
 * The new one is only added when it either has a new composition
 * (== number of double bonds) or is better than an already existing
 * solution with the same composition.
 * The resulting list is returned.
 */
{
   bond_set_node *bsnp;

   for (bsnp=old_list; bsnp; bsnp=bsnp->next)   /* look for same composition */
      if (new_sol->cardinality == bsnp->cardinality) break;

   if (bsnp)
   {
      if (CompareDenormalizations(new_sol, bsnp, acp, bcp, mp, ring_list) > 0)
      {
	 DisposeSet(bsnp->bond_set); bsnp->bond_set = new_sol->bond_set;
         MyFree((char *)new_sol);
         return (old_list);
      }
      else
      {
	 DisposeBondSetList(new_sol);
	 return (old_list);
      }
   }
   else
   {
      new_sol->next = old_list;
      return (new_sol);
   }
}

bond_set_node *RecDenormalize(struct atom_constraint_t *acp,
                              struct bond_constraint_t *bcp,
                              struct reaccs_molecule_t *mp,
                              neighbourhood_t          *nbp,
                              bond_set_node            *ring_list,
                              bond_set_node            *prev_solution)
/*
 * This function recursively enumerates all denormalizations compatible
 * with the constraints *acp and *bcp in the molecule context *mp, nbp[].
 * It uses the list of rings *ring_list to select the best solution.
 * prev_solution points to the best one so far or NULL if there is none.
 */
{
   int nopen;
   int i, ibranch, itauto, ntauto;
   bond_set_node *new_solution;
   struct atom_constraint_t *acp_new;
   struct bond_constraint_t *bcp_new;
   
// fprintf(stderr, "before Refine()\n");
   nopen = Refine(acp, bcp, mp, nbp);
// fprintf(stderr, "Refine() returned %d\n", nopen);
   if (nopen == VIOLATION)              /* constraints are inconsistent */
   {
      return (prev_solution);
   }
   else if (nopen == 0)                   /* denormalization found */
   {
// fprintf(stderr, "denormalization found\n");
      new_solution = NewBondSetNode(mp->n_bonds);
      for (i=0; i<mp->n_bonds; i++)
      if (bcp[i].open & CANDIDATE)
         if (bcp[i].can_be_double)
         {
            PutMember(new_solution->bond_set, i);
            new_solution->cardinality++;
         }
      new_solution = AddNewSolution(new_solution, prev_solution,
                                    acp, bcp, mp, ring_list);
      return (new_solution);
   }
   else                                    /* branch on open bond */
   {
      ntauto = (-1); ibranch = (-1); /* select best branch point */
      for (i=0; i<mp->n_bonds; i++)  /* prefer those with tautomer ends */
         if (bcp[i].open & CANDIDATE  &&
             bcp[i].can_be_single     &&
             bcp[i].can_be_double)
         {
            if (acp[mp->bond_array[i].atoms[0]-1].dbl_min <
                acp[mp->bond_array[i].atoms[0]-1].dbl_max)
               itauto = 1;
            else
               itauto = 0;
            if (acp[mp->bond_array[i].atoms[1]-1].dbl_min <
                acp[mp->bond_array[i].atoms[1]-1].dbl_max)
               itauto++;
            if (itauto > ntauto)
            {
               ntauto = itauto;
               ibranch = i;
            }
         }
      if (ibranch == (-1))
      {
	 ShowMessage("\ainternal error in RecDenormalize()\n",
		     "RecDenormalize");
	 return (prev_solution);
      }
      acp_new = TypeAlloc(mp->n_atoms, struct atom_constraint_t);
      bcp_new = TypeAlloc(mp->n_bonds, struct bond_constraint_t);
      for (i=0; i<mp->n_atoms; i++) acp_new[i] = acp[i];
      for (i=0; i<mp->n_bonds; i++) bcp_new[i] = bcp[i];
      bcp_new[ibranch].can_be_single = FALSE;
      acp_new[mp->bond_array[ibranch].atoms[0]-1].open |= TOUCHED;
      acp_new[mp->bond_array[ibranch].atoms[1]-1].open |= TOUCHED;
// fprintf(stderr, "DOUBLE branch on open bond %d-%d\n",mp->bond_array[ibranch].atoms[0],mp->bond_array[ibranch].atoms[1]);
      new_solution =
         RecDenormalize(acp_new, bcp_new, mp, nbp, ring_list, prev_solution);
      for (i=0; i<mp->n_atoms; i++) acp_new[i] = acp[i];
      for (i=0; i<mp->n_bonds; i++) bcp_new[i] = bcp[i];
      bcp_new[ibranch].can_be_double = FALSE;
      acp_new[mp->bond_array[ibranch].atoms[0]-1].open |= TOUCHED;
      acp_new[mp->bond_array[ibranch].atoms[1]-1].open |= TOUCHED;
// fprintf(stderr, "SINGLE branch on open bond %d-%d\n",mp->bond_array[ibranch].atoms[0],mp->bond_array[ibranch].atoms[1]);
      new_solution =
         RecDenormalize(acp_new, bcp_new, mp, nbp, ring_list, new_solution);
      MyFree((char *)bcp_new); MyFree((char *)acp_new);
      return (new_solution);
   }
}

int Denormalize(struct reaccs_molecule_t *mp,
                struct atom_constraint_t *acp,
                struct bond_constraint_t *bcp)
/*
 * This function changes the AROMATIC bonds of *mp to single and
 * double bonds. This process is the inverse of the CAS normalization
 * of alternating singe/double bonds and tautomer groups.
 *
 * In the case of tautomeric groups, the hydrogen count fields of the
 * atoms are assumed to contain information on the allowed range in its
 * two bytes.
 *
 * Non-tautomeric groups (i.e. alternating rings) ignore the hydrogen count.
 *
 * The algorithm used works by relaxation and backtracking on each
 * individual connected set of AROMATIC bonds in turn. It uses a number
 * of criteria to select from possible alternative 'denormalizations:'
 *
 *  1: prefer C=O double bonds.
 *  2: prefer alternating six-membered rings.
 *  3: prefer tetrasubstituted double bonds.
 *  4: prefer trisubstituted double bonds.
 * ..
 * 99: prefer higher unsaturation in case of multiple compatible compositions
 *
 * The criteria are listed in decreasing order of importance.
 *
 * The function returns an integer with bits set for different causes
 * of problems.
 *
 * The initial constraints are assumed to be set up in *acp and *bcp.
 */
{
   neighbourhood_t *nbp;
   int result;
   int i;
   int col;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   bond_set_node *ring_list, *solution;

   nbp   = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,nbp,mp->n_atoms);

   ColorNormalizedComponents(mp);
   ring_list = MoleculeRings(mp);
   result = 0;
   for (;;)
   {
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->bond_type == AROMATIC)
         {
            col = mp->atom_array[bp->atoms[0]-1].color;
            break;
         }
      if (i == mp->n_bonds) break;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == col)
	    acp[i].open = CANDIDATE | TOUCHED;
	 else
	    acp[i].open = 0;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->bond_type == AROMATIC  &&
             mp->atom_array[bp->atoms[0]-1].color == col)
         {
            bcp[i].open = CANDIDATE;
         }
         else
            bcp[i].open = 0;

      solution = (bond_set_node *)NULL;
      solution = RecDenormalize(acp, bcp, mp, nbp, ring_list, solution);
      if (solution)
      {
         if (solution->next) result |= MULTIPLE_COMPOSITION;

         for (i=0; i<mp->n_bonds; i++)
            if (bcp[i].open & CANDIDATE)
               if (IsMember(solution->bond_set, i))
                  mp->bond_array[i].bond_type = DOUBLE;
               else
                  mp->bond_array[i].bond_type = SINGLE;
         DisposeBondSetList(solution);
      }
      else
      {
         result |= CONSTRAINT_VIOLATION;
         break;
      }
   }
   DisposeBondSetList(ring_list);

   MyFree((char *)nbp);
   return (result);
}

#ifdef MAIN

FILE *log_file = (FILE *)NULL;

main(int argc, char *argv[])
{
   Fortran_FILE *finp;
   FILE *foutp;
   struct reaccs_molecule_t *mp;

   finp = FortranOpen(argv[1],"r");
   foutp = fopen(argv[2],"w");

   if (argc < 3  ||  (finp && foutp))
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL != ReadREACCSMolecule(finp,mp,""))
         return (EXIT_FAILURE);
      fprintf(stderr,"result = %d\n",CCTDenormalize(mp));
      PrintREACCSMolecule(foutp,mp,"");
      FreeMolecule(mp);
      return (EXIT_SUCCESS);
   }
   else
      return (EXIT_FAILURE);
}
#endif

