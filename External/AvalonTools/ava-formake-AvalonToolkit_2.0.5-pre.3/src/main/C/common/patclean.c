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
/*  File:                     patclean.c                                */
/*                                                                      */
/*  Purpose:                  This module implements the functions      */
/*                            for pattern or template oriented          */
/*                            cleaning.                                 */
/*                                                                      */
//  History:  29-Sep-92       Created module from experimental file     */
//                            sstest.c                                  */
//                                                                      */
/************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "forio.h"
#include "geometry.h"
#include "graph.h"
#include "local.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "ssmatch.h"
#include "utilities.h"

#include "stereo.h"
#include "layout.h"
#include "perceive.h"

#include "patclean.h"

typedef double point[2];

typedef unsigned atom_pair[2];

static bond_set_node *MoleculeRings(struct reaccs_molecule_t *mp)
/*
 * Returns the list of basis rings of the molecule *mp.
 * Rings with bond_type flag RUBBER_BOND are ignored.
 */
{
   int i, idummy;
   atom_pair *bonds;
   bond_set_node *ring_list;

   if (mp->n_bonds == 0) return (bond_set_node *)NULL;

   bonds = TypeAlloc(mp->n_bonds, atom_pair);
   idummy = mp->n_atoms;
   for (i=0; i<mp->n_bonds; i++)
      if (!(mp->bond_array[i].bond_type & RUBBER_BOND))
      {
	 bonds[i][0] = mp->bond_array[i].atoms[0];
	 bonds[i][1] = mp->bond_array[i].atoms[1];
      }
      else
      {
	 idummy++; bonds[i][0] = idummy;
	 idummy++; bonds[i][1] = idummy;
      }
   ring_list = RingList(bonds,mp->n_bonds);
   ring_list = CombineRings(ring_list);

   MyFree((char *)bonds);

   return (ring_list);
}

static double OverlapStrain(struct reaccs_molecule_t *kernel,
                            double new_coords[][2], unsigned nnew)
/*
 * Computes a "strain" value which reflects the overlap of the
 * the coordinates in new_coords[0..nnew-1][1..2] with the atoms
 * of *kernel.
 */
{
   int i, j;
   struct reaccs_atom_t *ap;
   double result;

   result = 0.0;
   for (i=0, ap=kernel->atom_array; i<kernel->n_atoms; i++, ap++)
      for (j=0; j<nnew; j++)
       result +=
             1.0/(0.01 + (ap->x - new_coords[j][0])*(ap->x - new_coords[j][0])
                       + (ap->y - new_coords[j][1])*(ap->y - new_coords[j][1]));
   return (result);
}

static int ColorSubstituents(struct reaccs_molecule_t *mp,
                        int border_color)
/*
 * Colors the substituents of the already colored part of *mp
 * starting with border_color and giving each connected substituent
 * a different color. (border_color values should be as small as possible,
 * because they are used in an array allocation!)
 * The function returns a number >= the largest color used in the
 * coloring.
 */
{
   int *real_colors;
   int max_color;
   int i, j;
   int changed, col1, col2;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;

   max_color = border_color+1;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == NONE)
      {
         ap->color = max_color; max_color++;
      }

   real_colors = TypeAlloc(max_color, int);
   for (i=0; i<max_color; i++)
      real_colors[i] = i;

                /* make bonds joining at one kernel atom connected */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color != NONE && ap->color <= border_color) /* kernel atom */
      {
         col1 = NONE;
         for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
            if (bp->color == border_color &&
                (bp->atoms[0] == i+1 || bp->atoms[1] == i+1))
            {
               if (bp->atoms[0] == i+1)
                  if (col1 == NONE)
                     col1 = mp->atom_array[bp->atoms[1]-1].color;
                  else
                  {
                     col2 = mp->atom_array[bp->atoms[1]-1].color;
                     if (real_colors[col1] < real_colors[col2])
                        real_colors[col2] = real_colors[col1];
                     else
                        real_colors[col1] = real_colors[col2];
                  }
                  if (bp->atoms[1] == i+1)
                     if (col1 == NONE)
                        col1 = mp->atom_array[bp->atoms[0]-1].color;
                     else
                     {
                        col2 = mp->atom_array[bp->atoms[0]-1].color;
                        if (real_colors[col1] < real_colors[col2])
                           real_colors[col2] = real_colors[col1];
                        else
                           real_colors[col1] = real_colors[col2];
                     }
            }
      }

   do
   {
      changed = FALSE;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->color == NONE)
         {
            col1 = mp->atom_array[bp->atoms[0]-1].color;
            col2 = mp->atom_array[bp->atoms[1]-1].color;
            if (real_colors[col1] == real_colors[col2])
               continue;
            if (real_colors[col1] < real_colors[col2])
               real_colors[col2] = real_colors[col1];
            else
               real_colors[col1] = real_colors[col2];
            changed = TRUE;
         }
   } while (changed);

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->color == NONE)
         bp->color = real_colors[mp->atom_array[bp->atoms[0]-1].color];

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color >= border_color)
         ap->color = real_colors[ap->color];

   MyFree((char *)real_colors);

   return (max_color);
}

#define MATCH_COLOR   (NONE+1)
#define BORDER_COLOR  (NONE+2)

static
int TransformColoredSubstituent(struct reaccs_molecule_t *mp,
                                struct reaccs_molecule_t *ssp,
                                ssmatch_t                *match,
                                int icol)
/*
 * Transformes the substituent colored by icol in the molecule
 * *mp. It assumes that the border bonds have been colored by
 * BORDER_COLOR, and that the match is colored by MATCH_COLOR.
 * The coordinates of the kernel are assumed to be set to those
 * of *ssp. *match is the substructure match currently investigated.
 */
{
   int i, j;
   int at1, at2, sat;
   struct reaccs_bond_t *bp;
   double tpoints[MAXATOMS][2], fpoints[MAXATOMS][2], point[2];
   int snumber[MAXATOMS];
   double new_points[MAXATOMS][2];
   int nnew;
   int npoints;

   double strain;

   double coordinates[MAXATOMS][2];
   unsigned edges[MAXBONDS][2];

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      exit (EXIT_FAILURE);
   }

   npoints = 0;
   for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
      if (bp->color == BORDER_COLOR)
      {
         if (mp->atom_array[bp->atoms[0]-1].color == icol)
         {
            at1 = bp->atoms[0]-1; at2 = bp->atoms[1]-1;
         }
         else if (mp->atom_array[bp->atoms[1]-1].color == icol)
         {
            at1 = bp->atoms[1]-1; at2 = bp->atoms[0]-1;
         }
         else
            continue;
         for (sat=0; sat<match->n_match; sat++)
            if (match->match_atoms[sat] == at2) break;

         /* at1 == substituent atom */
         /* at2 == kernel atom */
         /* sat == template atom corresp. to at2 */
fprintf(stderr,"%s(%d)-%s(%d)\n",mp->atom_array[at1].atom_symbol, at1+1, mp->atom_array[at2].atom_symbol, at2+1);
         if (sat == match->n_match)
         {
            fprintf(stderr,
                    "TransformColoredSubstituent: internal error\n");
            exit (EXIT_FAILURE);
         }
                      /* map from kernel to template coords. */
         fpoints[npoints][0] = mp->atom_array[at2].x;
         fpoints[npoints][1] = mp->atom_array[at2].y;
         tpoints[npoints][0] = ssp->atom_array[sat].x;
         tpoints[npoints][1] = ssp->atom_array[sat].y;
                 /* needed to check for spiral connections */
         snumber[npoints] = sat;
         npoints++;
      }

   if (npoints == 0)
      return (FALSE);
   else if (npoints == 1)
   {                /* one connection between kernel and substituent */
fprintf(stderr,"npoints = %d => one connection to kernel\n", npoints);
      for (i=0; i<ssp->n_bonds; i++) /* set up graph */
      {
         edges[i][0] = ssp->bond_array[i].atoms[0]-1;
         edges[i][1] = ssp->bond_array[i].atoms[1]-1;
      }
      for (i=0; i<ssp->n_atoms; i++)
      {
         coordinates[i][0] = ssp->atom_array[i].x;
         coordinates[i][1] = ssp->atom_array[i].y;
      }

      NextSubstituentPoint(point,
                           coordinates, ssp->n_atoms,
                           edges, ssp->n_bonds,
                           sat, USE_INWARDS, (int *)NULL, 0);
      tpoints[npoints][0] = point[0];
      tpoints[npoints][1] = point[1];
      fpoints[npoints][0] = mp->atom_array[at1].x;
      fpoints[npoints][1] = mp->atom_array[at1].y;
      snumber[npoints] = (-1);
      npoints++;
   }
   else /* Check if there is only one attachment atom and map */
   {            /* its environment, too. This is the spiral case. */
fprintf(stderr,"npoints = %d => more connections to kernel\n", npoints);
      for (j=1; j<npoints; j++)
      if (snumber[j] != snumber[j-1]) break;
      if (j == npoints)
      {
         sat = snumber[0];
         for (j=0, bp=ssp->bond_array; j<ssp->n_bonds; j++, bp++)
         {
            if (bp->atoms[0]-1 == sat)
               at1 = bp->atoms[1]-1;
            else if (bp->atoms[1]-1 == sat)
               at1 = bp->atoms[0]-1;
            else
               continue;
            fpoints[npoints][0] =
              mp->atom_array[match->match_atoms[at1]].x;
            fpoints[npoints][1] =
              mp->atom_array[match->match_atoms[at1]].y;
            tpoints[npoints][0] = ssp->atom_array[at1].x;
            tpoints[npoints][1] = ssp->atom_array[at1].y;
            snumber[npoints] = at1;
            npoints++;
         }
      }
   }

   for (j=0, nnew=0; j<mp->n_atoms; j++)
      if (mp->atom_array[j].color == icol)
      {
         new_points[nnew][0] = mp->atom_array[j].x;
         new_points[nnew][1] = mp->atom_array[j].y;
         nnew++;
      }

   PointSetMatchTransformation(new_points, nnew,
                               fpoints, tpoints, npoints,
                               TRUE);

   strain = OverlapStrain(ssp, new_points, nnew);
fprintf(stderr, "strain = %g, nnew=%d\n", strain, nnew);

   for (j=0, nnew=0; j<mp->n_atoms; j++)
      if (mp->atom_array[j].color == icol)
      {
         new_points[nnew][0] = mp->atom_array[j].x;
         new_points[nnew][1] = mp->atom_array[j].y;
         nnew++;
      }

   PointSetMatchTransformation(new_points, nnew,
                               fpoints, tpoints, npoints,
                               FALSE);

fprintf(stderr, "new strain = %g\n", OverlapStrain(ssp, new_points, nnew));
   if (strain < OverlapStrain(ssp, new_points, nnew))
   {
      for (j=0, nnew=0; j<mp->n_atoms; j++)
         if (mp->atom_array[j].color == icol)
         {
            new_points[nnew][0] = mp->atom_array[j].x;
            new_points[nnew][1] = mp->atom_array[j].y;
            nnew++;
         }

      PointSetMatchTransformation(new_points, nnew,
                                  fpoints, tpoints, npoints,
                                  TRUE);
   }

   for (j=0, nnew=0; j<mp->n_atoms; j++)
      if (mp->atom_array[j].color == icol)
      {
         mp->atom_array[j].x = new_points[nnew][0];
         mp->atom_array[j].y = new_points[nnew][1];
         nnew++;
      }

   return (TRUE);
}

static
int AdjustToSubstructureTemplate(struct reaccs_molecule_t *mp,
                                 struct reaccs_molecule_t *ssp,
                                 ssmatch_t *match)
/*
 * This function replaces the coordinates of the match of *ssp in *mp
 * by the corresponding substructure coordinates. The non-matching
 * parts of *mp are moved such that they keep their internal
 * relationships but receive a nicer place with respect to the template.
 *
 * Stereobonds that are unspecified in the source molecule and specified
 * in the match will be set to defined during this process.
 */
{
   int i, j, k;
   int icol;
   int at1, at2;
   int ai1, ai2;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int max_color;

   float oldx, oldy;
   float newx, newy;

   int stereo_ok, swap_flags;
   int parity[MAXATOMS], ph;
   int cistrans[MAXBONDS], cth;
   float x[MAXATOMS], y[MAXATOMS], z[MAXATOMS];
   neighbourhood_t neighbour_array[MAXATOMS], *nbp;
   int numbering[MAXATOMS];

   bond_set_node *ring_list, *plist;
   int *ring_size, *is_ring_atom;

   point p1, p2, p, r12, pp;
   double q;

   // sanity checks
   if (IsNULL(mp)  ||  IsNULL(ssp)) return (FALSE);
   if (mp->n_atoms <= 0  ||  ssp->n_atoms <= 0) return (FALSE);

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      exit (EXIT_FAILURE);
   }

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {                                            /* Save old coordinates */
      x[i] = ap->x; y[i] = ap->y; z[i] = ap->z;
   }
                                    /* Save stereochemistry */
   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);
   nbp = neighbour_array;

   is_ring_atom = TypeAlloc(mp->n_atoms, int);
   ring_size    = TypeAlloc(mp->n_bonds, int);
   for (i=0; i<mp->n_bonds; i++) ring_size[i] = 0;
   ring_list = MoleculeRings(mp);
   for (plist=ring_list; plist; plist=plist->next)
   {
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
	 if (IsMember(plist->bond_set,i))
	 {
	    is_ring_atom[bp->atoms[0]-1] = TRUE;
	    is_ring_atom[bp->atoms[1]-1] = TRUE;
	    if (ring_size[i] == 0  ||  ring_size[i] > plist->cardinality)
	       ring_size[i] = plist->cardinality;
	 }
   }
   DisposeBondSetList(ring_list);

   for (i=0; i<mp->n_atoms; i++)
      parity[i] = AtomParity(mp, i+1, &neighbour_array[i]);

   for (i=0; i<mp->n_atoms; i++) numbering[i] = i;
   CisTransPerception(mp, numbering);
   for (i=0; i<mp->n_bonds; i++)
      cistrans[i] = mp->bond_array[i].color;

   oldx = oldy = 0.0;     /* Save center of match */
   for (i=0; i<match->n_match; i++)
   {
      oldx += mp->atom_array[match->match_atoms[i]].x;
      oldy += mp->atom_array[match->match_atoms[i]].y;
   }
   oldx /= match->n_match; oldy /= match->n_match;

   ResetColors(mp);   /* Draw borderline between match and substituents */
   for (i=0; i<match->n_match; i++)
   {
      at1 = match->match_atoms[i];
      mp->atom_array[at1].color = MATCH_COLOR;
      for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
      {
         if (bp->atoms[0]-1 == at1)
            at2 = bp->atoms[1]-1;
         else if (bp->atoms[1]-1 == at1)
            at2 = bp->atoms[0]-1;
         else
            continue;
	 for (k=0; k<match->n_match; k++)
            if (at2 == match->match_atoms[k])
	    {
	       bp->color = MATCH_COLOR;
               break;
	    }
         if (k == match->n_match)
            bp->color = BORDER_COLOR;
      }
   }

   max_color = ColorSubstituents(mp, BORDER_COLOR);

                                             /* process substituents */
   for (icol=BORDER_COLOR+1; icol<=max_color; icol++)
   {
fprintf(stderr,"transforming substituent with color %d\n", icol);
      TransformColoredSubstituent(mp, ssp, match, icol);
   }

   ResetColors(mp);   /* Remove intermediate coloring */

   for (i=0; i<ssp->n_atoms; i++)
   {
      mp->atom_array[match->match_atoms[i]].x = ssp->atom_array[i].x;
      mp->atom_array[match->match_atoms[i]].y = ssp->atom_array[i].y;
   }

   stereo_ok = TRUE;
   for (i=0; i<mp->n_atoms; i++)
   {
      ph = AtomParity(mp, i+1, &neighbour_array[i]);
      if (ph != parity[i])
      {
// fprintf(stderr, "fixing parity of atom %d that changed from %d to %d\n", i+1, parity[i], ph);
         /* Try all possible combinations of swapping */
         for (swap_flags=0; swap_flags < (1<<4); swap_flags++)
         {
            for (j=0; j<neighbour_array[i].n_ligands; j++)
            {
               bp = &mp->bond_array[neighbour_array[i].bonds[j]];
               if (bp->stereo_symbol == UP  ||  bp->stereo_symbol == DOWN)
               {
                  if (swap_flags & (1<<j))
                     bp->stereo_symbol = UP;
                  else
                     bp->stereo_symbol = DOWN;
               }
            }
            ph = AtomParity(mp, i+1, &neighbour_array[i]);
            if (ph == parity[i]) break;
         }
	 ph = AtomParity(mp, i+1, &neighbour_array[i]);
	 if (ph != parity[i])
	 {
            fprintf(stderr, "parity of atom %d changed from %d to %d\n", i+1, parity[i], ph);
            stereo_ok = FALSE;
	 }
      }
   }

   for (i=0; i<mp->n_atoms; i++) numbering[i] = i;
   CisTransPerception(mp, numbering);

   /* Check if double bonds need to flip and do it if possible */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      ai1 = bp->atoms[0]-1; ai2 = bp->atoms[1]-1;
      cth = mp->bond_array[i].color;
      if (ring_size[i] == 0  &&  cistrans[i] != cth)
      {
      			/* color one part of double bond */
	 mp->atom_array[ai1].color = (-1); mp->atom_array[ai2].color = 1;
	 FloodColor(mp, nbp, ai2, 1); mp->atom_array[ai1].color = 0;

      			/* compute reference vector */
	 p1[0] = mp->atom_array[ai1].x; p1[1] = mp->atom_array[ai1].y;
	 p2[0] = mp->atom_array[ai2].x; p2[1] = mp->atom_array[ai2].y;
	 r12[0] = p2[0] - p1[0]; r12[1] = p2[1]-p1[1];

      			/* flip colored part of molecule */
	 for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
	    if (ap->color == 1)
	    {
	       p[0] = ap->x; p[1] = ap->y;
	       q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
		   (r12[0]*r12[0]       + r12[1]*r12[1]);
	       pp[0] = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
	       pp[1] = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
	       ap->x = pp[0]; ap->y = pp[1];
	    }
	 FlipStereoSymbols(mp, 1);
	 for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
	    ap->color = NONE;
      }
   }

   /* Perceive again with flipped bonds. */
   CisTransPerception(mp, numbering);

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      cth = mp->bond_array[i].color;
      if (cistrans[i] != cth)
      {
         fprintf(stderr, "bond (%d-%d) with order %d changed from %d to %d\n", bp->atoms[0], bp->atoms[1], bp->bond_type, cistrans[i], cth);
         stereo_ok = FALSE;
      }
   }

   if (ring_size) MyFree((char *)ring_size);
   if (is_ring_atom) MyFree((char *)is_ring_atom);

   if (!stereo_ok)
   {
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      {
	 ap->x = x[i]; ap->y = y[i]; ap->z = z[i];
      }
      if (log_file)
	fprintf(log_file,
		  "template clean failed for molecule '%s'\n",mp->name);
      else
        fprintf(stderr,
                 "template clean failed for molecule '%s'\n",mp->name);
      return (FALSE);
   }

   newx = newy = 0.0;   /* Get new origin of match */
   for (i=0; i<match->n_match; i++)
   {
      newx += mp->atom_array[match->match_atoms[i]].x;
      newy += mp->atom_array[match->match_atoms[i]].y;
   }
   newx /= match->n_match; newy /= match->n_match;

   for (i=0; i<match->n_match; i++)       /* color matching fragments */
      mp->atom_array[match->match_atoms[i]].color = MATCH_COLOR;
   FloodFillFragments(mp);

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color != NONE)   /* move colored fragments to old place */
      {
         ap->x += oldx-newx; ap->y += oldy-newy;
         ap->color = NONE;
      }
   ResetColors(mp);

   /* Fix cis/trans stereo for atoms in match */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->stereo_symbol != CIS_TRANS_EITHER) continue;
      if (bp->bond_type     != DOUBLE) 	         continue;
      if (ring_size[i]      == 0)		 continue;

      at1 = bp->atoms[0]-1; at2 = bp->atoms[1]-1;
      for (j=0; j<match->n_match; j++)
      {
	 if (match->match_atoms[j] != at1) continue;
	 for (k=0; k<match->n_match; k++)
	 {
	    if (match->match_atoms[k] != at2) continue;
	    /* We know that the two atoms of the bond belong to the match */
	    /* Just for test: Clear away either spec. */
	    /* We need to check if pattern isn't either too. */
	    bp->stereo_symbol = NONE;
	    break;
	 }
      }
   }

   return (TRUE);
}

static
double AlignMoleculeWithSubstructure(struct reaccs_molecule_t *mp,
                                     struct reaccs_molecule_t *ssp,
                                     ssmatch_t *match,
				     int transform)
/*
 * Moves and rotates *mp such that it 'aligns' with *ssp. *ssp is a
 * substructure of *mp and *match defines the atom-atom correspondence
 * for the match. The function returns the goodness of fit, i.e. the
 * RMS deviation between the rotated molecule and the substructure.
 *
 * The molecule is actually transformed only if transform == TRUE!!
 */
{
   static
   double mpoints[MAXATOMS][2];  /* 2D coordinates of molecule atoms */
   static
   double spoints[MAXATOMS][2];  /* 2D coordinates of substructure atoms */
   static
   double rpoints[MAXATOMS][2];  /* 2D coordinates of molecule atoms that
                                    map to substructure atoms */
   int i;
   double diff, fit;

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      exit (EXIT_FAILURE);
   }

   for (i=0; i<mp->n_atoms; i++)
   {
      mpoints[i][0] = mp->atom_array[i].x;
      mpoints[i][1] = mp->atom_array[i].y;
   }

   for (i=0; i<ssp->n_atoms; i++)
   {
      spoints[i][0] = ssp->atom_array[i].x;
      spoints[i][1] = ssp->atom_array[i].y;
      rpoints[i][0] = mp->atom_array[match->match_atoms[i]].x;
      rpoints[i][1] = mp->atom_array[match->match_atoms[i]].y;
   }

   PointSetMatchTransformation(mpoints, mp->n_atoms,
                               rpoints, spoints, ssp->n_atoms,
                               FALSE);

   for (i=0, fit=0; i<ssp->n_atoms; i++)
   {
      diff = mpoints[match->match_atoms[i]][0]-spoints[i][0];
      fit += diff*diff;
      diff = mpoints[match->match_atoms[i]][1]-spoints[i][1];
      fit += diff*diff;
   }

   if (transform)
      for (i=0; i<mp->n_atoms; i++)
      {
	 mp->atom_array[i].x = mpoints[i][0];
	 mp->atom_array[i].y = mpoints[i][1];
      }

   return (sqrt(fit));
}

ssmatch_t *ClosestMatch(struct reaccs_molecule_t *mp,
                        struct reaccs_molecule_t *ssp)
/*
 * Returns the geometrically closest substructure match of
 * *ssp in *mp. *mp is rotated and flipped such that it
 * best suits the match. NULL is returned if no substructure
 * match is found. Superflous matches are freed.
 */
{
   double fit, best_fit;
   int best_flip;
   ssmatch_t *matches, *hmp, *best_match;

   // check parameters
   if (IsNULL(mp)  ||  IsNULL(ssp)) return ((ssmatch_t *)NULL);

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      exit (EXIT_FAILURE);
   }

   /* find substructure matches */
   matches = SSMatch(mp, ssp, MULTIPLE_MATCHES);

   if (IsNULL(matches)) return (matches);

   /* find geometrically closest match of substructure */
   best_fit = 1.0e17; best_match = matches;
   best_flip = FALSE;
   /* look first for unflipped (proper) matches */
   for (hmp = matches; !IsNULL(hmp); hmp=hmp->next)
   {
      fit = AlignMoleculeWithSubstructure(mp, ssp, hmp, FALSE);
      if (fit < best_fit)
      {
	 best_fit = fit;
         best_match = hmp;
      }
   }
   /* now check flipped ones (i.e. mirror images) */
   FlipMolecule(mp, NONE);
   for (hmp = matches; !IsNULL(hmp); hmp=hmp->next)
   {
      fit = AlignMoleculeWithSubstructure(mp, ssp, hmp, FALSE);
      if (fit < best_fit)
      {
         best_fit = fit;
         best_flip = TRUE;      /* best match is a flipped one */
         best_match = hmp;
      }
   }

   if (best_flip)
      ;                           /* leave molecule flipped */
   else
      FlipMolecule(mp, NONE);     /* retract trial flipping */

   while (matches)
   {
      if (matches == best_match)
      {
         matches = matches->next;
         best_match->next = NULL;
      }
      else
      {
         hmp = matches->next;
         FreeSSMatch(matches);
         matches = hmp;
      }
   }
   FreeSSMatchHeap();

   fit = AlignMoleculeWithSubstructure(mp, ssp, best_match, TRUE);

   return (best_match);
}

void AddStereoHydrogen(struct reaccs_molecule_t *mp,
                       int iatom, int stereo_symbol)
/*
 * Adds a hydrogen atom at atom iatom using a SINGLE
 * bond with stereosymbol stereo_symbol.
 */
{
   int i, hnumber;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;

   double point[2];
   double coordinates[MAXATOMS][2];
   unsigned edges[MAXBONDS][2];

   for (i=0, ap = mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      coordinates[i][0] = ap->x;
      coordinates[i][1] = ap->y;
   }

   for (i=0, bp = mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      edges[i][0] = bp->atoms[0]-1;
      edges[i][1] = bp->atoms[1]-1;
   }

   NextSubstituentPoint(point,
                        coordinates, mp->n_atoms,
                        edges, mp->n_bonds,
                        iatom-1, USE_INWARDS, (int *)NULL, 0);

   hnumber = AddAtomToMolecule(mp, point[0], point[1], 0.0, "H");
   AddBondToMolecule(mp, iatom, hnumber, SINGLE, stereo_symbol);
}

int ForceStereoTemplate(struct reaccs_molecule_t *mp,
                        struct reaccs_molecule_t *ssp)
/*
 * Finds the geometrically closest match of *ssp in *mp and
 * forces the stereosymbols of *ssp onto the matching bonds
 * of *mp. It returns 1 if a match was applied, 0 if no match
 * was found, and (-1) if a possible match could not be applied.
 */
{
   char buffer[255];
   ssmatch_t *best_match;
   int problem_found, match_applied;

   int i, j, k, iatom, jatom;
   int symbol;
   struct reaccs_bond_t *bp, bm;
   neighbourhood_t neighbour_array[MAXATOMS], *nbp;

   // sanity checks
   if (IsNULL(mp)  ||  IsNULL(ssp)) return (FALSE);
   if (mp->n_atoms <= 0  ||  ssp->n_atoms <= 0) return (FALSE);

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      exit (EXIT_FAILURE);
   }

   best_match = ClosestMatch(mp, ssp);
   if (IsNULL(best_match)) return (0);

   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);
   PerceiveRingBonds(mp);

   /* Look for stereoatoms of pattern */
   ResetColors(ssp);
   for (i=0, bp=ssp->bond_array; i<ssp->n_bonds; i++, bp++)
      if (bp->stereo_symbol == UP  || bp->stereo_symbol == DOWN)
         ssp->atom_array[bp->atoms[0]-1].color++;

   problem_found = FALSE;
   match_applied = FALSE;
                                 /* for all matching centers */
   for (i=0; i<best_match->n_match; i++)
   {
      iatom = best_match->match_atoms[i];
      nbp = &neighbour_array[iatom];
      /* if the substructure has stereo and the structure hasn't */
      if (ssp->atom_array[i].color != NONE  &&
          AtomParity(mp, iatom+1, nbp) == UNDEFINED_PARITY)
      {
         /* Look for matching bond and partner atom */
         for (j=0, bp=ssp->bond_array; j<ssp->n_bonds; j++, bp++)
            if (bp->atoms[0] == i+1  &&
                (bp->stereo_symbol == UP  ||  bp->stereo_symbol == DOWN))
            {
               symbol = bp->stereo_symbol;
               jatom  = best_match->match_atoms[bp->atoms[1]-1];
               for (k=0; k<mp->n_bonds; k++)
               {        /* mp->bond_array may be reallocated */
                  bm = mp->bond_array[k];
                       /* fix atom ordering if necessary */
                  if (bm.atoms[1] == iatom+1  &&
                      bm.atoms[0] == jatom+1)
                     if (bm.stereo_symbol == NONE)
                     {   /* only unspecified bonds shall be changed */
                        sprintf(buffer,
                                "ordering of bond bond %d-%d was fixed",
                                iatom+1, jatom+1);
                        AddMsgToList(buffer);
                        mp->bond_array[k].atoms[0] = bm.atoms[1];
                        mp->bond_array[k].atoms[1] = bm.atoms[0];
                        bm = mp->bond_array[k];
                     }
                  else
                  {
                     problem_found = TRUE;
                     sprintf(buffer,"bond %d-%d could not be fixed", iatom+1, jatom+1);
                     AddMsgToList(buffer);
                     break;
                  }
                  /* if this is the matching bond */
                  if (bm.atoms[0] == iatom+1  &&
                      bm.atoms[1] == jatom+1)
                  {
                     /* only unspecified bonds shall be changed */
                     if (bm.stereo_symbol != NONE)
                     {
                        problem_found = TRUE;
                        sprintf(buffer,"bond %d-%d could not be fixed", iatom+1, jatom+1);
                        AddMsgToList(buffer);
                     }
                     else if (bm.topography != RING) /* Easy case */
                     {
                        sprintf(buffer,"bond %d-%d was fixed", iatom+1, jatom+1);
                        AddMsgToList(buffer);
                        mp->bond_array[k].stereo_symbol = symbol;
                        match_applied = TRUE;
                     }
                     else       /* ring stereochemistry */
                     {
                        if (symbol == UP)
                           AddStereoHydrogen(mp, iatom+1, DOWN);
                        else
                           AddStereoHydrogen(mp, iatom+1, UP);
                        sprintf(buffer, "hydrogen was added for ring bond %d-%d", iatom+1, jatom+1);
                        AddMsgToList(buffer);
                        match_applied = TRUE;
                     }
                     break;     /* there can not be more matching bonds */
                  }
               }
               break; /* only one stereobond per center shall be used */
            }
      }
   }

   FreeSSMatch(best_match);
   FreeSSMatchHeap();

   /* Clear ring bond information */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      bp->dummy = bp->topography = 0;

   if (problem_found)       return (-1);
   else if (match_applied)  return (1);
   else                     return (0);
}

#define SQR(x) ((x)*(x))

void ScaleByTemplate(struct reaccs_molecule_t *mp,
                     struct reaccs_molecule_t *ssp,
		     double slimit)
/*
 * Scales the coordinates of molecule *mp according to the average
 * bond length of *ssp, if the scalefactor doesn't exceed slimit.
 */
{
   double mbl, sbl;
   int i;
   struct reaccs_bond_t *bp;

   // sanity checks
   if (IsNULL(mp)  ||  IsNULL(ssp)) return;
   if (mp->n_atoms <= 0  ||  ssp->n_atoms <= 0) return;

   for (i=0, mbl=0.0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      mbl +=
         SQR(mp->atom_array[bp->atoms[0]-1].x -
	     mp->atom_array[bp->atoms[1]-1].x) +
         SQR(mp->atom_array[bp->atoms[0]-1].y -
	     mp->atom_array[bp->atoms[1]-1].y);
   }
   if (mp->n_bonds > 0) mbl = sqrt(mbl/mp->n_bonds);

   for (i=0, sbl=0.0, bp=ssp->bond_array; i<ssp->n_bonds; i++, bp++)
   {
      sbl +=
         SQR(ssp->atom_array[bp->atoms[0]-1].x -
	     ssp->atom_array[bp->atoms[1]-1].x) +
         SQR(ssp->atom_array[bp->atoms[0]-1].y -
	     ssp->atom_array[bp->atoms[1]-1].y);
   }
   if (ssp->n_bonds > 0) sbl = sqrt(sbl/ssp->n_bonds);

   if (slimit > 1.0) slimit = 1.0/slimit;

   if (sbl < mbl  &&  sbl < mbl*slimit) return;
   if (sbl > mbl  &&  sbl*slimit > mbl) return;

   for (i=0; i<mp->n_atoms; i++)
   {
      mp->atom_array[i].x *= sbl/mbl;
      mp->atom_array[i].y *= sbl/mbl;
      mp->atom_array[i].z *= sbl/mbl;
   }
}

int TemplateClean(struct reaccs_molecule_t *mp,
                  struct reaccs_molecule_t *ssp)
/*
 * Checks if the substructure *ssp matches *mp, finds the
 * geometrically closest match, and tries to clean the matching
 * portion of *mp to the coordinates of *ssp, while moving the
 * substituents of the match to reasonable positions.
 * The function returns TRUE if the coordinates of *mp have
 * been modified and FALSE if not.
 *
 * Stereobonds that are unspecified in the source molecule and specified
 * in the match will be set to defined during this process.
 */
{
   ssmatch_t *best_match;
   int adjusted;
   int *H_count;
   struct reaccs_bond_t *bp;
   int i;

   // sanity checks
   if (IsNULL(mp)  ||  IsNULL(ssp)) return (FALSE);
   if (mp->n_atoms <= 0  ||  ssp->n_atoms <= 0) return (FALSE);

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      return (FALSE);
   }

   /* Set up hydrogen count fields in structure for matching */
   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);
   /* Add the explicit hydrogens to the implicit counts */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (0 == strcmp("H", mp->atom_array[bp->atoms[0]-1].atom_symbol))
         H_count[bp->atoms[1]]++;
      else if (0 == strcmp("H", mp->atom_array[bp->atoms[1]-1].atom_symbol))
         H_count[bp->atoms[0]]++;
   }
   /* set the 'query_H_count' field to the correct value */
   for (i=0; i<mp->n_atoms; i++)
      if (H_count[i+1] >= 0) 
         mp->atom_array[i].query_H_count = ZERO_COUNT+H_count[i+1];

   best_match = ClosestMatch(mp, ssp);

   /* return memory and reset *mp's query_H_count fields */
   MyFree((char *)H_count);
   for (i=0; i<mp->n_atoms; i++)
      mp->atom_array[i].query_H_count = NONE;
   for (i=0; i<mp->n_bonds; i++)
      mp->bond_array[i].topography = NONE;

   if (IsNULL(best_match)) return (FALSE);

   /* Scale molecule according to template */
   ScaleByTemplate(mp, ssp, 0.1);
   adjusted = AdjustToSubstructureTemplate(mp, ssp, best_match);

   FreeSSMatch(best_match);
   FreeSSMatchHeap();

   return (adjusted);
}

int TemplateRotate(struct reaccs_molecule_t *mp,
                   struct reaccs_molecule_t *ssp)
/*
 * Checks if the substructure *ssp matches *mp, finds the
 * geometrically closest match, and tries to rotate *mp
 * to most closely resemble the coordinates of the pattern.
 */
{
   ssmatch_t *best_match;
   struct reaccs_bond_t *bp;
   int i;
   struct reaccs_molecule_t *chained_ssp;

   // sanity checks
   if (IsNULL(mp)  ||  IsNULL(ssp)) return (FALSE);
   if (mp->n_atoms <= 0  ||  ssp->n_atoms <= 0) return (FALSE);

   if (mp->n_atoms > MAXATOMS)
   {
      fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
      return (FALSE);
   }

   // first try if template chain bonds would match structure chain bonds only and then perform Clean instead of rotate
   chained_ssp = CopyMolecule(ssp);
   PerceiveRingBonds(chained_ssp);
   for (i=0, bp=chained_ssp->bond_array; i<chained_ssp->n_bonds; i++, bp++)
   {
       if (bp->topography == RING)
           bp->dummy = 0;
       else
           bp->topography = CHAIN;
   }
// PrintREACCSMolecule(stderr,chained_ssp,"");
   best_match = SSMatch(mp, chained_ssp, SINGLE_MATCH);
   if (!IsNULL(best_match))
   {
// fprintf(stderr,"Chain bonds in template all match chain bonds in structure => do a clean instead of rotation\n");
       TemplateClean(mp, chained_ssp);
       /* Clear ring bond information */
       for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
          bp->dummy = bp->topography = 0;
       FreeSSMatch(best_match);
       FreeSSMatchHeap();
       FreeMolecule(chained_ssp);
       return TRUE;
   }
   FreeMolecule(chained_ssp);
// fprintf(stderr,"Just doing template rotate\n");

   best_match = ClosestMatch(mp, ssp); /* implicitly rotates mp */
   /* Clear ring bond information */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      bp->dummy = bp->topography = 0;

   if (IsNULL(best_match)) return (FALSE);

   FreeSSMatch(best_match);
   FreeSSMatchHeap();

   return (TRUE);
}

#if MAIN

#ifdef __TURBOC__
#include <process.h>
unsigned _stklen = 0xFEEE;
#endif

main(int argc, char *argv[])
{
   struct reaccs_molecule_t *mp, *ssp;
   struct reaccs_molecule_t *patterns[10];
   int npat;
   struct data_line_t *data_list; /* variables to store SD-file data */
   struct data_line_t *dph;
   Fortran_FILE *fp, *subfp;
   FILE *outfile;
   ssmatch_t *matches, *hmp, *best_match;
   int i;
   double fit, best_fit;
   int best_flip;
   int is_mol_file;

   if (argc < 3  ||  argc > 4)
   {
      fprintf(stderr,
              "usage: sstest <input File> <pattern file> [<output file>]\n");
      return (EXIT_FAILURE);
   }
   fp = FortranOpen(argv[1],"r");
   if (!fp)
   {
      fprintf(stderr,"Could not open file '%s'\n",argv[1]);
      exit (EXIT_FAILURE);
   }

   mp = TypeAlloc(1,struct reaccs_molecule_t);

   ssp = TypeAlloc(1,struct reaccs_molecule_t);
   subfp = FortranOpen(argv[2],"r");
   if (!subfp)
   {
      fprintf(stderr,"Could not open file '%s'\n",argv[2]);
      exit (EXIT_FAILURE);
   }

   if (argc > 3)
      outfile = CheckedFileOpen(argv[3],"w");
   else
      outfile = stdout;

   npat = 0;
   while (FORTRAN_NORMAL == ReadREACCSMolecule(subfp,ssp,""))
   {
      GetBuffer(subfp);                                /* Skip $$$$ line */
      MakeHydrogensImplicit(ssp);
      for (i=0; i<ssp->n_atoms; i++)
         ssp->atom_array[i].query_H_count = NONE;
      patterns[npat] = ssp; npat++;
      ssp = TypeAlloc(1, struct reaccs_molecule_t);
   }
   if (npat == 0)
   {
      fprintf(stderr,"sstest: could not interpret file '%s'\n",argv[2]);
      return (EXIT_FAILURE);
   }
   FortranClose(subfp);

   while (FORTRAN_NORMAL == ReadREACCSMolecule(fp,mp,""))
   {
      if (mp->n_atoms > MAXATOMS)
      {
         fprintf(stderr,"Too many atoms in molecule '%s'\n",mp->name);
         exit (EXIT_FAILURE);
      }

      data_list = ReadMACCSDataLines(fp);
      if (fp->status == FORTRAN_EOF)
         is_mol_file = TRUE;
      else
      {
         GetBuffer(fp);                           /* Skip $$$$ line */
         is_mol_file = FALSE;
      }

      fprintf(stderr, "aligning molecule '%s' (%ld)\n",
              mp->name, mp->registry_number);

      for (i=0; i<npat; i++)
      {
         ssp = patterns[i];
         changed = TemplateClean(mp, ssp);
         if (changed)
         {
            PrintREACCSMolecule(outfile,mp,"");
            for (dph = data_list; !IsNULL(dph); dph = dph->next)
               fprintf(outfile,"%s\n",dph->data);
            if (!is_mol_file) fprintf(outfile,"$$$$\n");
         }
      }

      while (!IsNULL(data_list))
      {
         dph = data_list->next;
         free((char *)data_list);
         data_list = dph;
      }
      FreeMolecule(mp);
      mp = TypeAlloc(1,struct reaccs_molecule_t);
   }
   FortranClose(fp);

   return (EXIT_SUCCESS);
}
#endif
