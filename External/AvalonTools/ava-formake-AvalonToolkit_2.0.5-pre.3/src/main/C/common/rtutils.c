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

/**********************************************************************/
/*                                                                    */
/*      File:           rtutils.c                                     */
/*                                                                    */
/*      Purpose:        This file implements the functions needed for */
/*                      Reaction Type processing.                     */
/*                                                                    */
//      History:        18-May-94       Start building from pieces of */
//                                      update_rttable.c              */
//                                                                    */
/**********************************************************************/

#include "reaccs.h"
#include "utilities.h"
#include "perceive.h"

void SaveBondOrder(struct reaccs_molecule_t *mp)
/*
 * Save the current bond order in the dummy field of the bond
 * structure.
 */
{
   struct reaccs_bond_t *bp;
   int i;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      bp->dummy = bp->bond_type;
}

void RestoreBondOrder(struct reaccs_molecule_t *mp)
/*
 * Restores the current bond order from the dummy field of the bond
 * structure.
 */
{
   struct reaccs_bond_t *bp;
   int i;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      bp->bond_type = bp->dummy;
      bp->dummy = NONE;
   }
}

void PerceiveReactionCenter(struct reaccs_reaction_t *rp)
/*
 * Uses the atom-atom-mappingto perceive which bonds are
 * unchanged, changed, or broken during the reaction.
 * The procedure uses the color fields of the bonds and atoms.
 */
{
   struct reaccs_molecule_t *r, *p;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bpr, *bpp;
   int ir, ip;

   ApplyToAllMolecules(rp, ResetColors);
   ApplyToAllMolecules(rp, SaveBondOrder);
   ApplyToAllMolecules(rp, PerceiveAromaticBonds);

   /* look for conserved or changed bonds */
   for (r=rp->reactants; !IsNULL(r); r=r->next)
      for (p=rp->products; !IsNULL(p); p=p->next)
      {
         for (bpr=r->bond_array, ir=0; ir<r->n_bonds; bpr++, ir++)
            if (r->atom_array[bpr->atoms[0]-1].mapping != NONE  &&
                r->atom_array[bpr->atoms[1]-1].mapping != NONE)
               for (bpp=p->bond_array, ip=0; ip<p->n_bonds; bpp++, ip++)
                  if ((r->atom_array[bpr->atoms[0]-1].mapping ==
                       p->atom_array[bpp->atoms[0]-1].mapping  &&
                       r->atom_array[bpr->atoms[1]-1].mapping ==
                       p->atom_array[bpp->atoms[1]-1].mapping)   ||
                      (r->atom_array[bpr->atoms[0]-1].mapping ==
                       p->atom_array[bpp->atoms[1]-1].mapping  &&
                       r->atom_array[bpr->atoms[1]-1].mapping ==
                       p->atom_array[bpp->atoms[0]-1].mapping))
                     if (bpr->bond_type == bpp->bond_type)
                     {
                        bpr->color |= UNCHANGED;
                        bpp->color |= UNCHANGED;
                     }
                     else
                     {
                        bpr->color |= CHANGE;
                        bpp->color |= CHANGE;
                     }
      }

   /* look for broken bonds */
   for (r=rp->reactants; !IsNULL(r); r=r->next)
      for (bpr=r->bond_array, ir=0; ir<r->n_bonds; bpr++, ir++)
         if (bpr->color == NONE)
            for (p=rp->products; !IsNULL(p); p=p->next)
               for (ap=p->atom_array, ip=0; ip<p->n_atoms; ap++, ip++)
                  if ((r->atom_array[bpr->atoms[0]-1].mapping == ap->mapping  &&
                       ap->mapping != NONE) ||
                      (r->atom_array[bpr->atoms[1]-1].mapping == ap->mapping  &&
                       ap->mapping != NONE))
                     bpr->color = MAKE_BREAK;

   /* look for made bonds */
   for (p=rp->products; !IsNULL(p); p=p->next)
      for (bpp=p->bond_array, ip=0; ip<p->n_bonds; bpp++, ip++)
         if (bpp->color == NONE)
            for (r=rp->reactants; !IsNULL(r); r=r->next)
               for (ap=r->atom_array, ir=0; ir<r->n_atoms; ap++, ir++)
                  if ((p->atom_array[bpp->atoms[0]-1].mapping == ap->mapping  &&
                       ap->mapping != NONE) ||
                      (p->atom_array[bpp->atoms[1]-1].mapping == ap->mapping  &&
                       ap->mapping != NONE))
                     bpp->color = MAKE_BREAK;

   ApplyToAllMolecules(rp, RestoreBondOrder);

   /* Check result */
   for (r=rp->reactants; !IsNULL(r); r=r->next)
      for (bpr=r->bond_array, ir=0; ir<r->n_bonds; bpr++, ir++)
         if (bpr->color != bpr->reaction_mark)
         {
            bpr->reaction_mark = bpr->color;
            bpr->color = NONE;
         }

   for (p=rp->products; !IsNULL(p); p=p->next)
      for (bpp=p->bond_array, ip=0; ip<p->n_bonds; bpp++, ip++)
         if (bpp->color != bpp->reaction_mark)
         {
            bpp->reaction_mark = bpp->color;
            bpp->color = NONE;
         }
}

