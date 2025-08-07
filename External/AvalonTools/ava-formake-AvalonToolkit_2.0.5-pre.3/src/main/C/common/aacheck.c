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
/*  File:     aacheck.c                                                 */
/*                                                                      */
/*  Purpose:  Implements utilities needed for augmented atom checking   */
/*            and transformation.                                       */
/*                                                                      */
//  Revisions: 1.0.0   11-Sep-92      Creation of module as an extract  */
//                                    from the old struchk.c file.      */
//                                                                      */
/************************************************************************/

#include "aacheck.h"

#include "pattern.h"
#include "reaccs.h"
#include "utilities.h"

extern FILE *aa_log;

int CheckAtoms(struct reaccs_molecule_t *mp,
               augmented_atom_t         good_atoms[],
               int                      ngood)
/*
 * Checks if every atom in *mp matches one of the augmented atoms
 * in godd_atoms[0..ngood-1]. It returns TRUE if all atoms gave a match
 * and FALSE otherwise.
 */
{
   neighbourhood_t *neighbour_array;
   unsigned int match[MAXNEIGHBOURS+1];
   unsigned int i, j;
   unsigned int len;
   unsigned int nmatch;
   char buffer[256];
   char aabuffer[1024];
   int *atom_status, *bond_status;

   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);

   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

   nmatch = 0;
   for (i=0; i<mp->n_atoms; i++)
   {
      for (j=0; j<ngood; j++)
      {
         // check for ring state of central atom
         if (good_atoms[j].topography == RING   &&  atom_status[i] == 0) continue;
         if (good_atoms[j].topography == CHAIN  &&  atom_status[i] != 0) continue;
         if (neighbour_array[i].n_ligands == good_atoms[j].n_ligands &&
             AAMatch(mp,i,match,&good_atoms[j],atom_status,neighbour_array))
         {
            nmatch++;
            break;
         }
      }
      if (j == ngood)
      {
         snprintf(buffer, MAXMSG,
                 "%10s    atom %3d    AA : %s",
                 mp->name,
                 i+1,
                 AAToString(aabuffer, mp, i, &neighbour_array[i]));
         AddMsgToList(buffer);
         if (aa_log)
         {
            fprintf(aa_log,"/* %02d */ \"",neighbour_array[i].n_ligands);
            len = AAPrint(aa_log,mp,i,&neighbour_array[i]);
            fprintf(aa_log,"\",%.*s/* %s */\n",35-len,"",mp->name);
         }
      }
   }

   MyFree((char *)bond_status);
   MyFree((char *)atom_status);
   MyFree((char *)neighbour_array);

   return (nmatch == mp->n_atoms);
}

