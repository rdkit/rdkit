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
/*      File:           layout.c                                        */
/*                                                                      */
/*      Purpose:        This file implements layouting of connection	*/
/*			tables, either from scratch or starting from	*/
/*			prelayouted fragments. It can layout molecules	*/
/*			as well as the components of reactions.		*/
//									*/
//	History:	06-Apr-1994	Start of development.		*/
/*									*/
/************************************************************************/

/* used for debugging support code */
#ifdef __TURBOC__
#include <process.h>
#include <alloc.h>
#endif

#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "local.h"
#include "reaccs.h"
#include "forio.h"
#include "reaccsio.h"
#include "utilities.h"
#include "graph.h"
#include "geometry.h"
#include "stereo.h"

#include "layout.h"

#define STDBOND 1.514
#define PI      3.14159265359
#define SQR(x)  ((x)*(x))

#define FRAME_HEIGHT (1.54*4.0/3.0)

typedef double point[2];

static void LogString(char *message)
{
   FILE *fp;
   fp = fopen("layout.log", "a");
   fprintf(fp, "%s\n", message);
   fclose(fp);
}

void RandomCoordinates(struct reaccs_molecule_t *mp)
/*
 * Sets the X and Y coordinates of the atoms of *mp to random
 * values.
 */
{
   int i;

   for (i=0; i<mp->n_atoms; i++)
      if ((mp->atom_array[i].color & KEEP_POSITION)  == 0)
      {
         mp->atom_array[i].x =
            STDBOND*sqrt((double)mp->n_atoms)*rand()/(double)RAND_MAX;
         mp->atom_array[i].y =
            STDBOND*sqrt((double)mp->n_atoms)*rand()/(double)RAND_MAX;
         mp->atom_array[i].z = 0.0;
      }
}

void RecolorMolecule(struct reaccs_molecule_t *mp)
/*
 * Recolors the atoms and bonds of the molecule *mp, each one
 * with a different color.
 */
{
   int i, ncolor;

   ncolor = 1;
   for (i=0; i<mp->n_atoms; i++) mp->atom_array[i].color = ncolor++;
   for (i=0; i<mp->n_bonds; i++) mp->bond_array[i].color = ncolor++;
}

void PutColorIntoValue(struct reaccs_molecule_t *mp)
/*
 * Makes the colors of the atoms of *mp displayable by putting them
 * into the value fields of the atoms.
 */
{
   int i;
   for (i=0; i<mp->n_atoms; i++)
      mp->atom_array[i].value = mp->atom_array[i].color;
}

int FloodColor(struct reaccs_molecule_t *mp,
               neighbourhood_t *nbp,
               int aindex, int color)
/*
 * Recursively colors the fragment of *mp containing the atom with
 * index aindex with color. It uses a flood fill algorithm.
 * It returns the number of atoms recolored during the process.
 */
{
   int result;
   int i;

   mp->atom_array[aindex].color = color; result = 1;
   for (i=0; i<nbp[aindex].n_ligands; i++)
      if (mp->atom_array[nbp[aindex].atoms[i]].color == NO_COLOR  &&
          !(mp->bond_array[nbp[aindex].bonds[i]].bond_type & RUBBER_BOND))
         result += FloodColor(mp, nbp, nbp[aindex].atoms[i], color);

   return (result);
}

int FloodClearColor(struct reaccs_molecule_t *mp,
                    neighbourhood_t *nbp,
                    int aindex, int color)
/*
 * Recursively clears the color value of those atoms colored with color.
 * The algorithms recursively goes through all neighbours of aindex
 * colored with color. It returns the number of cleared colors.
 */
{
   int result;
   int i;

   if (mp->atom_array[aindex].color != color) return (0);

   mp->atom_array[aindex].color = NONE; result = 1;
   for (i=0; i<nbp[aindex].n_ligands; i++)
      if (mp->atom_array[nbp[aindex].atoms[i]].color == color  &&
          !(mp->bond_array[nbp[aindex].bonds[i]].bond_type & RUBBER_BOND))
         result += FloodClearColor(mp, nbp, nbp[aindex].atoms[i], color);

   return (result);
}

void MakeBondTrans(struct reaccs_molecule_t *mp,
                   neighbourhood_t *nbp,
                   struct reaccs_bond_t *bp)
/*
 * Flips the bond *bp if necessary to make it 'trans'.
 */
{
   int i, j;
   int indexa[MAXNEIGHBOURS], sizea[MAXNEIGHBOURS], na, aia;
   int indexb[MAXNEIGHBOURS], sizeb[MAXNEIGHBOURS], nb, aib;
   int *atom_colors;

   point p1, p2, p, r12;
   double strain1, strain2, q;
   struct reaccs_atom_t *ap;

   atom_colors = TypeAlloc(mp->n_atoms, int);   /* save and clear colors */
   for (i=0; i<mp->n_atoms; i++)
   {
      atom_colors[i] = mp->atom_array[i].color;
      mp->atom_array[i].color = NONE;
   }

   /* collect ligands and sizes */
   aia = bp->atoms[0]-1; aib = bp->atoms[1]-1;
   mp->atom_array[aia].color = 1; mp->atom_array[aib].color = 2;
   na = 0;
   for (i=0; i<nbp[aia].n_ligands; i++)
      if (mp->atom_array[nbp[aia].atoms[i]].color == NONE)
      {
         indexa[na] = nbp[aia].atoms[i];
         sizea[na] = FloodColor(mp, nbp, indexa[na], 3);
         FloodClearColor(mp, nbp, indexa[na], 3);
         na++;
      }
   nb = 0;
   for (i=0; i<nbp[aib].n_ligands; i++)
      if (mp->atom_array[nbp[aib].atoms[i]].color == NONE)
      {
         indexb[nb] = nbp[aib].atoms[i];
         sizeb[nb] = FloodColor(mp, nbp, indexb[nb], 3);
         FloodClearColor(mp, nbp, indexb[nb], 3);
         nb++;
      }

#define DIST(x1,y1,x2,y2) (((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)))

   strain1 = 0;
   for (i=0; i<na; i++)
      for (j=0; j<nb; j++)
         strain1 += sizea[i]*sizeb[j]/(1 + DIST(mp->atom_array[indexa[i]].x,
                                                mp->atom_array[indexa[i]].y,
                                                mp->atom_array[indexb[j]].x,
                                                mp->atom_array[indexb[j]].y));

   p1[0] = mp->atom_array[bp->atoms[0]-1].x;    /* compute bond vector */
   p1[1] = mp->atom_array[bp->atoms[0]-1].y;
   p2[0] = mp->atom_array[bp->atoms[1]-1].x;
   p2[1] = mp->atom_array[bp->atoms[1]-1].y;
   r12[0] = p2[0] - p1[0]; r12[1] = p2[1]-p1[1];
   /* transform ligands a */
   for (i=0; i<na; i++)
   {
      ap = &mp->atom_array[indexa[i]];
      p[0] = ap->x; p[1] = ap->y;
      q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
          (r12[0]*r12[0]       + r12[1]*r12[1]);
      ap->x = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
      ap->y = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
   }

   strain2 = 0;
   for (i=0; i<na; i++)
      for (j=0; j<nb; j++)
         strain2 += sizea[i]*sizeb[j]/(1 + DIST(mp->atom_array[indexa[i]].x,
                                                mp->atom_array[indexa[i]].y,
                                                mp->atom_array[indexb[j]].x,
                                                mp->atom_array[indexb[j]].y));
   /* transform back */
   for (i=0; i<na; i++)
   {
      ap = &mp->atom_array[indexa[i]];
      p[0] = ap->x; p[1] = ap->y;
      q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
          (r12[0]*r12[0]       + r12[1]*r12[1]);
      ap->x = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
      ap->y = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
   }

   if (strain1 > strain2)       /* cis -> trans */
   {
      for (i=0; i<na; i++)
         FloodColor(mp, nbp, indexa[i], 3);
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == 3)
         {
            p[0] = ap->x; p[1] = ap->y;
            q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
                (r12[0]*r12[0]       + r12[1]*r12[1]);
            ap->x = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
            ap->y = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
         }
   }

   for (i=0; i<mp->n_atoms; i++)        /* restore colors */
      mp->atom_array[i].color = atom_colors[i];
   MyFree((char *)atom_colors);
}

int InvertFragmentColor(struct reaccs_molecule_t *mp,
                        neighbourhood_t *nbp,
                        int aindex, int color)
/*
 * Recursively multiplies the color fields in the fragment of *mp containing
 * the atom with index aindex by (-1) when it is already colored with color.
 * The search continues only through positvely colored atoms.
 */
{
   int result;
   int i;

   result = 0;
   for (i=0; i<nbp[aindex].n_ligands; i++)
      if (mp->atom_array[nbp[aindex].atoms[i]].color == color  &&
          !(mp->bond_array[nbp[aindex].bonds[i]].bond_type & RUBBER_BOND))
      {
         result++;
         mp->atom_array[nbp[aindex].atoms[i]].color *= (-1);
         result += InvertFragmentColor(mp, nbp, nbp[aindex].atoms[i], color);
      }

   return (result);
}

struct atom_list
   {
      int atom;
      struct atom_list *next;
   };

typedef unsigned atom_pair[2];
static bond_set_node *MoleculeRings(struct reaccs_molecule_t *mp)
/*
 * Returns the list of basis rings of the molecule *mp.
 * Ring with bond_type flag RUBBER_BOND are ignored.
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

struct ring_node_t
   {
      float xfirst, yfirst;/* coordinates of first atom in fragment path */
      float xlast,  ylast; /* coordinates of last atom in fragment path  */
      int aifirst, ailast; /* atom numbers of first and last atoms 	 */
      int color;           /* color of fragment 			 */
      /* scratch fields */
      float d;             /* distance of aifirst and ailast */
      float angle;         /* angle on ring perimeter */
   };

/* save some typing */
#define ACOLOR(i)       (mp->atom_array[i-1].color)

int FillRingNodeTable(struct ring_node_t 	*rntp,
                      bond_set_node 	   	*ringp,
                      struct reaccs_molecule_t	*mp)
/*
 * Fills the ring node table *rntp with the entries derived
 * from the ring's bond set *ringp and the info from the molecule
 * *mp. It returns the number of entries in the table.
 */
{
   int i, j, nrnt;
   int from, to;
   int h;
   atom_pair *bonds;

   if (mp->n_bonds <= 0) return (0);
   bonds = TypeAlloc(mp->n_bonds, atom_pair);
   for (i=0; i<mp->n_bonds; i++)
   {
      bonds[i][0] = mp->bond_array[i].atoms[0];
      bonds[i][1] = mp->bond_array[i].atoms[1];
   }

   for (i=0, j=0; i<mp->n_bonds; i++)   /* squeeze out irrelevant bonds */
      if (IsMember(ringp->bond_set,i))
      {
         bonds[j][0] = bonds[i][0]; bonds[j][1] = bonds[i][1]; j++;
      }
   nrnt = j;

   for (i=1; i<nrnt; i++)       /* sort bonds into sequence */
      for (j=i; j<nrnt; j++)
         if (bonds[j][0] == bonds[i-1][1])
         {
            h = bonds[i][0]; bonds[i][0] = bonds[j][0]; bonds[j][0] = h;
            h = bonds[i][1]; bonds[i][1] = bonds[j][1]; bonds[j][1] = h;
            break;
         }
         else if (bonds[j][1] == bonds[i-1][1])
         {
            h = bonds[j][0]; bonds[j][0] = bonds[j][1]; bonds[j][1] = h;
            h = bonds[i][0]; bonds[i][0] = bonds[j][0]; bonds[j][0] = h;
            h = bonds[i][1]; bonds[i][1] = bonds[j][1]; bonds[j][1] = h;
            break;
         }

   for (i=1; i<nrnt; i++)
      if (ACOLOR(bonds[i-1][0]) != ACOLOR(bonds[i][0])) break;
   if (i == nrnt)       /* ring has single color -> already layouted */
   {
      MyFree((char *)bonds);
      return (0);
   }
   from = i-1;

   for (i=nrnt-1; i>from; i--)  /* find color sequence of segment */
   {                            /* that might wrap */
      if (ACOLOR(bonds[i][0]) != ACOLOR(bonds[from][0])) break;
   }
   to = i;

   rntp[0].aifirst = bonds[to][1];   rntp[0].ailast  = bonds[from][0];
   rntp[0].color   = ACOLOR(bonds[from][0]);
   nrnt = 1;
   for (i=from+1; i<=to; i++)
      if (ACOLOR(bonds[i][0]) != ACOLOR(bonds[i-1][0]))
      {
         for (j=i; j<to; j++)
            if (ACOLOR(bonds[j][0]) != ACOLOR(bonds[j][1])) break;
         rntp[nrnt].aifirst = bonds[i][0]; rntp[nrnt].ailast  = bonds[j][0];
         rntp[nrnt].color   = ACOLOR(bonds[i][0]); nrnt++;
      }

   for (i=0; i<nrnt; i++)
   {
      rntp[i].xfirst = mp->atom_array[rntp[i].aifirst-1].x;
      rntp[i].yfirst = mp->atom_array[rntp[i].aifirst-1].y;
      rntp[i].xlast  = mp->atom_array[rntp[i].ailast-1].x;
      rntp[i].ylast  = mp->atom_array[rntp[i].ailast-1].y;
   }

   MyFree((char *)bonds);
   return (nrnt);
}
#undef ACOLOR

void LayoutRingSegment(struct reaccs_molecule_t *mp,
                       struct ring_node_t       *rntp,
                       int                       nrnt)
/*
 * Places the atoms of *mp that have colors listed in rntp[0..nrnt-1]
 * relative to each other. It uses the perimeter of a circle to for
 * coordinates of the atoms listed as aifirst and ailast in rntp[].
 */
{
   int i, j, n;
   double xc, yc;
   double x1, y1, x2, y2;
   double peri, len;
   double r, angle;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   point *points, p1, p2, ps, p1p, p2p;
   int npts;

   xc = yc = 0.0; n = 0;                /* compute center of circle */
   for (i=1; i<nrnt; i++)
      for (j=0; j<mp->n_atoms; j++)
         if (mp->atom_array[j].color == rntp[i].color)
         {
            xc += mp->atom_array[j].x;
            yc += mp->atom_array[j].y;
            n++;
         }
   xc /= n; yc /= n;

   len = 0.0;                           /* compute radius of circle */
   for (i=0; i<nrnt; i++)
   {
      len += STDBOND;                   /* bond secante */
      rntp[i].d = sqrt(SQR(rntp[i].xfirst-rntp[i].xlast) +
                       SQR(rntp[i].yfirst-rntp[i].ylast));
      len += rntp[i].d;                                  /* fragment secante */
   }
   peri = 0.0;
   for (i=0; i<nrnt; i++)
   {
      peri += STDBOND*(PI*STDBOND/len)/sin(PI*STDBOND/len);
      rntp[i].angle = 2*PI*rntp[i].d/len;
      if (rntp[i].angle > 0.0001)
         peri += rntp[i].d*(rntp[i].angle/2)/sin(rntp[i].angle/2);
   }
   r = peri/(2*PI);

   points = TypeAlloc(mp->n_atoms, point);
   angle = 0;
   for (i=0; i<nrnt; i++)       /* new coordinates */
   {
      x1 = r*sin(angle)+xc; y1 = r*cos(angle)+yc;
      angle += rntp[i].angle;
      x2 = r*sin(angle)+xc; y2 = r*cos(angle)+yc;
      angle += 2*PI*STDBOND/len;

      p1[0] = rntp[i].xfirst; p1[1] = rntp[i].yfirst;
      p2[0] = rntp[i].xlast;  p2[1] = rntp[i].ylast;
      p1p[0] = x1; p1p[1] = y1; p2p[0] = x2; p2p[1] = y2;
      for (j=0, npts=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
         if (ap->color == rntp[i].color)
         {
            points[npts][0] = ap->x; points[npts][1] =ap->y;
            npts++;
         }

      if (SQR(p1[0]-p2[0])+SQR(p1[1]-p2[1]) > 0.01*SQR(STDBOND))
      {                         /* ring fusion */
         ps[0] = ps[1] = 0.0;   /* for center of gravity of fragment */
         for (j=0,npts=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
            if (ap->color == rntp[i].color)
            {
               ps[0] += ap->x; ps[1] += ap->y; npts++;
            }
         ps[0] /= npts; ps[1] /= npts;
         /* use mirror image of fragment if needed */
         if ((ps[0]-p1[0])*(p2[1]-p1[1]) - (ps[1]-p1[1])*(p2[0]-p1[0]) > 0)
         {
            p1[1] *= (-1); p2[1] *= (-1);
            for (j=0; j<npts; j++) points[j][1] *= (-1);
            FlipStereoSymbols(mp, rntp[i].color);
         }
      }
      else                      /* single atom or spiro */
      {
         ps[0] = ps[1] = 0.0;   /* compute center of gravity of atoms */
         npts = 0;              /* attached to link atom */
         for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
            if ((bp->atoms[0] == rntp[i].aifirst  &&
                 mp->atom_array[bp->atoms[1]-1].color == rntp[i].color) ||
                (bp->atoms[1] == rntp[i].aifirst  &&
                 mp->atom_array[bp->atoms[0]-1].color == rntp[i].color) ||
                (bp->atoms[0] == rntp[i].ailast  &&
                 mp->atom_array[bp->atoms[1]-1].color == rntp[i].color) ||
                (bp->atoms[1] == rntp[i].ailast  &&
                 mp->atom_array[bp->atoms[0]-1].color == rntp[i].color))
            {
               ps[0] += mp->atom_array[bp->atoms[0]-1].x;
               ps[1] += mp->atom_array[bp->atoms[0]-1].y;
               npts++;
               ps[0] += mp->atom_array[bp->atoms[1]-1].x;
               ps[1] += mp->atom_array[bp->atoms[1]-1].y;
               npts++;
            }
         if (npts > 0)
         {
            ps[0] /= npts; ps[1] /= npts;
            p2[0] = ps[0]; p2[1] = ps[1];
            p2p[0] = p1p[0] +
                     (p1p[0]-xc)*sqrt(SQR(ps[0]-p1[0])+SQR(ps[1]-p1[1]))/
                                 sqrt(SQR(p1p[0]-xc)+SQR(p1p[1]-yc));
            p2p[1] = p1p[1] +
                     (p1p[1]-yc)*sqrt(SQR(ps[0]-p1[0])+SQR(ps[1]-p1[1]))/
                                 sqrt(SQR(p1p[0]-xc)+SQR(p1p[1]-yc));
         }
      }

      TransformPoints(points, mp->n_atoms, p1, p2, p1p, p2p);
      for (j=0, npts=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
         if (ap->color == rntp[i].color)
         {
            ap->x = points[npts][0]; ap->y = points[npts][1];
            npts++;
         }
   }
   MyFree((char *)points);
// fprintf(stderr,"peri = %g, len = %g\n",peri,len);
}

void MergeConnectedFragmentColors(struct reaccs_molecule_t *mp)
/*
 * Merge the colors of connected fragments to be the smallest.
 */
{
    int i;
    struct reaccs_bond_t *bp;
    int col1, col2, h;

    for (;;)
    {
        for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
        {
            col1 = mp->atom_array[bp->atoms[0]-1].color;
            col2 = mp->atom_array[bp->atoms[1]-1].color;
            if (col1 != col2) break;
        }
        if (i == mp->n_bonds) break;
        if (col1 > col2) { h = col1; col1 = col2; col2 = h; }
        for (i=0; i<mp->n_atoms; i++)
            if (mp->atom_array[i].color = col2) mp->atom_array[i].color = col1;
    }
}

void MergeColors(struct reaccs_molecule_t *mp,
                 struct ring_node_t       *rntp,
                 int                       nrnt)
/*
 * Merges the atom colors listed in rntp[0..nrnt-1] into one.
 */
{
   int i, j;

   for (i=1; i<nrnt; i++)
      for (j=0; j<mp->n_atoms; j++)
         if (mp->atom_array[j].color == rntp[i].color)
            mp->atom_array[j].color = rntp[0].color;
}

void HideMetalComplexProblems(struct reaccs_molecule_t *mp)
/*
 * Tries to fix the problems with layouting metal complexes by
 * making the bonds to metals RUBBER_BONDs when they are in
 * more than one ring.
 */
{
   struct reaccs_bond_t *bp;
   int i;
   int has_metal;
   int *atom_status;
   int *bond_status;

   if (mp->n_bonds <= 0) return;

   for (i=0, has_metal=FALSE; i<mp->n_atoms; i++)
      if (AtomSymbolMatch(mp->atom_array[i].atom_symbol,"trn,H,D,T"))
         has_metal = TRUE;
   if (!has_metal) return;
// fprintf(stderr, "molecule has metal or hydrogen\n");

   atom_status = TypeAlloc(mp->n_atoms, int);
   bond_status = TypeAlloc(mp->n_bonds, int);
   RingState(mp, atom_status, bond_status);

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bond_status[i] > 1)
      {
         if (AtomSymbolMatch(mp->atom_array[bp->atoms[0]-1].atom_symbol,
                             "trn")    ||
             AtomSymbolMatch(mp->atom_array[bp->atoms[1]-1].atom_symbol,
                             "trn"))
            bp->bond_type |= RUBBER_BOND;
// fprintf(stderr, "ring bond %d-%d\n", bp->atoms[0], bp->atoms[1]);
      }
      else
      { /* Bonds to hydrogen shall be layed out last. */
// fprintf(stderr, "chain bond %d-%d\n", bp->atoms[0], bp->atoms[1]);
         if (AtomSymbolMatch(mp->atom_array[bp->atoms[0]-1].atom_symbol,
                             "H")    ||
             AtomSymbolMatch(mp->atom_array[bp->atoms[0]-1].atom_symbol,
                             "D")    ||
             AtomSymbolMatch(mp->atom_array[bp->atoms[0]-1].atom_symbol,
                             "T")    ||
             AtomSymbolMatch(mp->atom_array[bp->atoms[1]-1].atom_symbol,
                             "H")    ||
             AtomSymbolMatch(mp->atom_array[bp->atoms[1]-1].atom_symbol,
                             "D")    ||
             AtomSymbolMatch(mp->atom_array[bp->atoms[1]-1].atom_symbol,
                             "T"))
         {
            bp->bond_type |= RUBBER_BOND;
// fprintf(stderr, "rubber bond (%d)-(%d)\n", bp->atoms[0], bp->atoms[1]);
         }
      }

   if (atom_status) MyFree((char *)atom_status);
   if (bond_status) MyFree((char *)bond_status);
}

bond_set_node *LayoutRingSort(bond_set_node *ring_list)
/*
 * Orders the rings in *ring_list such that six-membered ones
 * come first follwed by the others by increasing size.
 *
 * If there is a tie, then prefer those with higher total degree.
 */
{
   unsigned       cardinality;
   bit_set_t     *bond_set;
   bond_set_node *rlh;
   int nring;

   nring = 0;
   for (rlh = ring_list; rlh; rlh=rlh->next)
      nring++;
   for (; nring>0; nring--)
      for (rlh = ring_list; rlh  &&  rlh->next; rlh=rlh->next)
         if ((rlh->cardinality != 6  &&
              rlh->next->cardinality == 6)  ||
             (rlh->cardinality != 6  &&
              rlh->next->cardinality != 6  &&
              rlh->cardinality > rlh->next->cardinality))
         {
            cardinality = rlh->next->cardinality;
            rlh->next->cardinality = rlh->cardinality;
            rlh->cardinality = cardinality;
            bond_set = rlh->next->bond_set;
            rlh->next->bond_set = rlh->bond_set;
            rlh->bond_set = bond_set;
         }

   return (ring_list);
}

void LayoutRings(struct reaccs_molecule_t *mp)
/*
 * Set coordinates for the atoms within ring systems. It
 * Recolors atoms and bonds which have been placed relative
 * to each other with the same color.
 */
{
   bond_set_node *ring_list, *plist;
   struct ring_node_t *rnt;
   int nrnt;

   rnt = TypeAlloc(mp->n_bonds+1, struct ring_node_t);
   HideMetalComplexProblems(mp);

   ring_list = MoleculeRings(mp);
   ring_list = LayoutRingSort(ring_list);
   for (plist=ring_list; plist; /* increment in loop */)
   {
      nrnt = FillRingNodeTable(rnt, plist, mp);
      if (nrnt == 0) plist = plist->next;
      else
      {
         LayoutRingSegment(mp, rnt, nrnt);
         MergeColors(mp, rnt, nrnt);
      }
   }
   DisposeBondSetList(ring_list);
   MyFree((char *)rnt);
}

void GetColoredGraph(struct reaccs_molecule_t *mp,
                     atom_pair *edges, int *nedges,
                     point *coords,    int *nnodes,
                     int *numbers, int color)
/*
 * Fetches the graph and the coordinates of those nodes from *mp
 * that are colored with color. The array numbers transforms
 * the indices into mp->atom_array to indices into coords[].
 */
{
   int i;
   struct reaccs_bond_t *bp;

   (*nedges) = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (!(bp->bond_type & RUBBER_BOND)                 &&
          mp->atom_array[bp->atoms[0]-1].color == color  &&
          mp->atom_array[bp->atoms[1]-1].color == color)
      {
         edges[*nedges][0] = bp->atoms[0];
         edges[*nedges][1] = bp->atoms[1];
         (*nedges)++;
      }
   (*nnodes) = 0;
   for (i=0; i<mp->n_atoms; i++)
      if (mp->atom_array[i].color == color)
      {
         coords[*nnodes][0] = mp->atom_array[i].x;
         coords[*nnodes][1] = mp->atom_array[i].y;
         numbers[i] = *nnodes; (*nnodes)++;
      }
   for (i=0; i<*nedges; i++)
   {
      edges[i][0] = numbers[edges[i][0]-1];
      edges[i][1] = numbers[edges[i][1]-1];
   }
}

double ColorStrain(struct reaccs_molecule_t *mp,
                   int                       col1,
                   int                       col2)
/*
 * Computes the strain measure between the two parts of the molecule
 * *mp having the colors col1 and col2.
 *
 * col1 == col2, is also legal
 */
{
   struct reaccs_atom_t *ap1, *ap2;
   int i1, i2;
   double result;

   result = 0;
   for (i1=0, ap1=mp->atom_array; i1<mp->n_atoms; i1++, ap1++)
      if (ap1->color == col1)
         for (i2=0, ap2=mp->atom_array; i2<mp->n_atoms; i2++, ap2++)
            if (ap2->color == col2  &&  i1 != i2)
               result += 1/(0.1 + (ap1->x - ap2->x)* (ap1->x - ap2->x)
                                + (ap1->y - ap2->y)* (ap1->y - ap2->y));

   return (result);
}

int ImproveBondByFlip(struct reaccs_molecule_t *mp,
                      struct reaccs_bond_t     *bp,
                      double threshold)
/*
 * Tries to improve atom collisions by flipping one side of bond bp.
 * The function returnes TRUE if i sucessfully improved the coordinates.
 * Changes are made only if strain improves by at least threshold.
 */
{
   point p1, p2, p, pp, r12;
   double strain1, strain2, q;
   int i;
   struct reaccs_atom_t *ap;

   if (bp->bond_type & DONT_FLIP_BOND)
      return (0);

   strain1 = ColorStrain(mp, mp->atom_array[bp->atoms[0]-1].color,
                         mp->atom_array[bp->atoms[1]-1].color);

   p1[0] = mp->atom_array[bp->atoms[0]-1].x;
   p1[1] = mp->atom_array[bp->atoms[0]-1].y;
   p2[0] = mp->atom_array[bp->atoms[1]-1].x;
   p2[1] = mp->atom_array[bp->atoms[1]-1].y;
   r12[0] = p2[0] - p1[0]; r12[1] = p2[1]-p1[1];

   /* Now flip bond */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == mp->atom_array[bp->atoms[1]-1].color)
      {
         p[0] = ap->x; p[1] = ap->y;
         q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
             (r12[0]*r12[0]       + r12[1]*r12[1]);
         pp[0] = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
         pp[1] = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
         ap->x = pp[0]; ap->y = pp[1];
      }

   strain2 = ColorStrain(mp, mp->atom_array[bp->atoms[0]-1].color,
                         mp->atom_array[bp->atoms[1]-1].color);

   /* transform back if improvement isn't sufficient */
   if (strain1 <= strain2*(1+threshold)+threshold)
   {
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == mp->atom_array[bp->atoms[1]-1].color)
         {
            p[0] = ap->x; p[1] = ap->y;
            q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
                (r12[0]*r12[0]       + r12[1]*r12[1]);
            pp[0] = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
            pp[1] = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
            ap->x = pp[0]; ap->y = pp[1];
         }

      return (FALSE);
   }

   return (TRUE);
}

void ImproveBondByStretch(struct reaccs_molecule_t *mp,
                          struct reaccs_bond_t     *bp)
/*
 * Tries to improve atom collisions by stretching the bond by 1/3rd.
 */
{
   point r12;
   int i, i1, i2;
   int collision;
   struct reaccs_atom_t *ap, *ap1, *ap2;
   int col1, col2;

   col1 = mp->atom_array[bp->atoms[0]-1].color;
   col2 = mp->atom_array[bp->atoms[1]-1].color;

   collision = FALSE;
   for (i1=0, ap1=mp->atom_array; i1<mp->n_atoms; i1++, ap1++)
      if (ap1->color == col1)
         for (i2=0, ap2=mp->atom_array; i2<mp->n_atoms; i2++, ap2++)
            if (ap2->color == col2  &&
                (ap1->x - ap2->x) * (ap1->x - ap2->x) +
                (ap1->y - ap2->y) * (ap1->y - ap2->y) <
                0.05*STDBOND*STDBOND)
            {
               collision = TRUE; break;
            }

   if (!collision) return;

   r12[0] = mp->atom_array[bp->atoms[1]-1].x - mp->atom_array[bp->atoms[0]-1].x;
   r12[1] = mp->atom_array[bp->atoms[1]-1].y - mp->atom_array[bp->atoms[0]-1].y;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == mp->atom_array[bp->atoms[1]-1].color)
      {
         ap->x += 0.3*r12[0]; ap->y += 0.3*r12[1];
      }
}

int ImproveMoleculeByBondFlip(struct reaccs_molecule_t *mp,
                              neighbourhood_t *nbp,
                              int ring_size[],
                              int is_ring_atom[],
                              int frag_color,
                              int only_single_bonds)
/*
 * Tries to improve the coordinates of molecule *mp by flipping
 * around the bonds. It cycles until no further improvement is possible.
 * It returns the number of improvement steps taken.
 * If frag_color != 0, only the corresponding fragment is optimized.
 */
{
   int i, j, changed;
   struct reaccs_bond_t *bp;
   int ai0, ai1;
   int oldcolor, ninvert;
   int ret;
   int nloop;

   ret = 0;
   /* improve bond flipping bonds attached to a ring */
   for (nloop=0; nloop<2; nloop++)
   {
      changed = FALSE;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      {
         ai0 = bp->atoms[0]-1; ai1 = bp->atoms[1]-1;
         if (mp->atom_array[ai0].color & KEEP_POSITION) continue;
         if (mp->atom_array[ai1].color & KEEP_POSITION) continue;
         if (only_single_bonds  &&  bp->bond_type != SINGLE) continue;
         if (!(bp->bond_type & RUBBER_BOND)             &&
             ring_size[i] == 0                          &&
             nbp[ai0].n_ligands > 1  			&&
             nbp[ai1].n_ligands > 1  			&&
             (frag_color == 0  ||
              (mp->atom_array[ai0].color == frag_color  &&
               mp->atom_array[ai1].color == frag_color)))
         {
            oldcolor = mp->atom_array[ai0].color;
            mp->atom_array[ai0].color = 0;
            mp->atom_array[ai1].color *= (-1);
            ninvert = InvertFragmentColor(mp, nbp, ai1, oldcolor);
            if (ninvert > 0)
            {
               mp->atom_array[ai0].color = oldcolor;
               if (is_ring_atom[ai1])	/* temp. swap atoms of bond */
               {
                  bp->atoms[0] = ai1+1; bp->atoms[1] = ai0+1;
               }
               if (ImproveBondByFlip(mp, bp, 0.1))
               { 	/* no cyclyng if only one color is improved */
                  if (frag_color == 0) changed = TRUE;
                  ret++;
               }
               for (j=0; j<mp->n_atoms; j++)
                  if (mp->atom_array[j].color == oldcolor*(-1))
                     mp->atom_array[j].color = oldcolor;
               if (is_ring_atom[ai1])	/* restore bond */
               {
                  bp->atoms[0] = ai0+1; bp->atoms[1] = ai1+1;
               }
            }
            else
            {
               mp->atom_array[ai0].color = oldcolor;
               mp->atom_array[ai1].color *= (-1);
            }
         }
      }
      if (mp->n_bonds > 60)	/* only one try for large molecules */
         changed = FALSE;
      if (!changed) break;
   }

   if (mp->n_atoms > 60) return (ret);/* too expensive for large structures */

   /* improve by flipping chain bonds at atoms with > 2 neighbours */
   for (nloop=0; nloop<2; nloop++)
   {
      changed = FALSE;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (!(bp->bond_type & RUBBER_BOND)       &&
             !is_ring_atom[bp->atoms[0]-1] 	  &&
             !is_ring_atom[bp->atoms[1]-1]        &&
             nbp[bp->atoms[0]-1].n_ligands > 1    &&
             nbp[bp->atoms[1]-1].n_ligands > 1    &&
             (nbp[bp->atoms[1]-1].n_ligands > 2 || nbp[bp->atoms[0]-1].n_ligands > 2)  &&
             (frag_color == 0  ||
              (mp->atom_array[bp->atoms[0]-1].color == frag_color  &&
               mp->atom_array[bp->atoms[0]-1].color == frag_color)))
         {
            if (only_single_bonds  &&  bp->bond_type != SINGLE) continue;
            oldcolor = mp->atom_array[bp->atoms[0]-1].color;
            mp->atom_array[bp->atoms[0]-1].color = 0;
            mp->atom_array[bp->atoms[1]-1].color *= (-1);
            ninvert = InvertFragmentColor(mp, nbp, bp->atoms[1]-1, oldcolor);
            if (ninvert > 0)
            {
               mp->atom_array[bp->atoms[0]-1].color = oldcolor;
               if (ImproveBondByFlip(mp, bp, 1.0))
               { 	/* no cyclyng if only one color is improved */
                  if (frag_color == 0) changed = TRUE;
                  ret++;
               }
               for (j=0; j<mp->n_atoms; j++)
                  if (mp->atom_array[j].color == oldcolor*(-1))
                     mp->atom_array[j].color = oldcolor;
            }
            else
            {
               mp->atom_array[bp->atoms[0]-1].color = oldcolor;
               mp->atom_array[bp->atoms[1]-1].color *= (-1);
            }
         }
      if (!changed) break;
   }

   if (mp->n_atoms > 30) return (ret);/* too expensive for large structures */

   /* improve by flipping chain bonds */
   for (nloop=0; nloop<2; nloop++)
   {
      changed = FALSE;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (!(bp->bond_type & RUBBER_BOND)     &&
             !is_ring_atom[bp->atoms[0]-1] 	&&
             !is_ring_atom[bp->atoms[1]-1]      &&
             nbp[bp->atoms[0]-1].n_ligands > 1  &&
             nbp[bp->atoms[1]-1].n_ligands > 1  &&
             (frag_color == 0  ||
              (mp->atom_array[bp->atoms[0]-1].color == frag_color  &&
               mp->atom_array[bp->atoms[0]-1].color == frag_color)))
         {
            if (only_single_bonds  &&  bp->bond_type != SINGLE) continue;
            oldcolor = mp->atom_array[bp->atoms[0]-1].color;
            mp->atom_array[bp->atoms[0]-1].color = 0;
            mp->atom_array[bp->atoms[1]-1].color *= (-1);
            ninvert = InvertFragmentColor(mp, nbp, bp->atoms[1]-1, oldcolor);
            if (ninvert > 0)
            {
               mp->atom_array[bp->atoms[0]-1].color = oldcolor;
               if (ImproveBondByFlip(mp, bp, 1.0))
               { 	/* no cyclyng if only one color is improved */
                  if (frag_color == 0) changed = TRUE;
                  ret++;
               }
               for (j=0; j<mp->n_atoms; j++)
                  if (mp->atom_array[j].color == oldcolor*(-1))
                     mp->atom_array[j].color = oldcolor;
            }
            else
            {
               mp->atom_array[bp->atoms[0]-1].color = oldcolor;
               mp->atom_array[bp->atoms[1]-1].color *= (-1);
            }
         }
      if (!changed) break;
   }

   return (ret);
}

void ImproveSPAtoms(struct reaccs_molecule_t *mp,
                    neighbourhood_t *nbp,
                    int atom)
/*
 * Makes sp carbon and nitrogen atoms linear.
 * If atom is != 0, only this atom is considered.
 */
{
   struct reaccs_atom_t *ap, *aph;
   struct reaccs_bond_t *bp;
   int i, j, changed;
   atom_pair *edges;
   int nedges;
   point *coords, p1, p2, p1p, p2p;
   int nnodes;
   int *numbers;
   int *colors;

   colors = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++)
      colors[i] = mp->atom_array[i].color;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (nbp[i].n_ligands == 2             &&
          (0 == strcmp(ap->atom_symbol,"C") ||
           0 == strcmp(ap->atom_symbol,"N"))  &&
          (mp->bond_array[nbp[i].bonds[0]].bond_type == SINGLE  &&
           mp->bond_array[nbp[i].bonds[1]].bond_type == TRIPLE    ||
           mp->bond_array[nbp[i].bonds[1]].bond_type == SINGLE  &&
           mp->bond_array[nbp[i].bonds[0]].bond_type == TRIPLE    ||
           mp->bond_array[nbp[i].bonds[0]].bond_type == DOUBLE  &&
           mp->bond_array[nbp[i].bonds[1]].bond_type == DOUBLE)  &&
          (atom == 0  ||  atom == i+1))
      {
         for (j=0, aph=mp->atom_array; j<mp->n_atoms; j++, aph++)
            aph->color = NONE;

         ap->color = 1;
         mp->atom_array[nbp[i].atoms[0]].color = 2;
         mp->atom_array[nbp[i].atoms[1]].color = 3;

         do
         {
            changed = FALSE;
            for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
               if (mp->atom_array[bp->atoms[0]-1].color > 0  &&
                   mp->atom_array[bp->atoms[1]-1].color == 0)
               {
                  changed = TRUE;
                  mp->atom_array[bp->atoms[1]-1].color =
                     mp->atom_array[bp->atoms[0]-1].color;
               }
               else if (mp->atom_array[bp->atoms[1]-1].color > 0  &&
                        mp->atom_array[bp->atoms[0]-1].color == 0)
               {
                  changed = TRUE;
                  mp->atom_array[bp->atoms[0]-1].color =
                     mp->atom_array[bp->atoms[1]-1].color;
               }
         } while (changed);
         if (mp->atom_array[nbp[i].atoms[0]].color !=	/* no ring */
             mp->atom_array[nbp[i].atoms[1]].color)
         {
            numbers = TypeAlloc(mp->n_atoms, int);
            edges = TypeAlloc(mp->n_bonds, atom_pair);
            coords = TypeAlloc(mp->n_atoms, point);

            p1[0] = ap->x;
            p1[1] = ap->y;
            p2[0] = mp->atom_array[nbp[i].atoms[0]].x;
            p2[1] = mp->atom_array[nbp[i].atoms[0]].y;
            p1p[0] = ap->x;
            p1p[1] = ap->y;
            p2p[0] = 2*p1p[0]-mp->atom_array[nbp[i].atoms[1]].x;
            p2p[1] = 2*p1p[1]-mp->atom_array[nbp[i].atoms[1]].y;

            GetColoredGraph(mp,
            		    edges,  &nedges, coords, &nnodes,
                            numbers, 2);
            TransformPoints(coords, nnodes, p1, p2, p1p, p2p);

            for (j=0; j<mp->n_atoms; j++)
               if (mp->atom_array[j].color == 2)
               {
                  mp->atom_array[j].x = coords[numbers[j]][0];
                  mp->atom_array[j].y = coords[numbers[j]][1];
               }

            MyFree((char *)coords);
            MyFree((char *)edges);
            MyFree((char *)numbers);
         }
      }

   for (i=0; i<mp->n_atoms; i++)
      mp->atom_array[i].color = colors[i] ;
   MyFree((char *)colors);
}

int CountColor(struct reaccs_molecule_t *mp, int color)
/*
 * Counts the number of atoms in *mp which arec colored in color.
 */
{
   int result;
   int i;
   struct reaccs_atom_t *ap;

   result = 0;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == color) result++;

   return (result);
}

void LayoutLargeLactams(struct reaccs_molecule_t *mp,
                        int is_ring_atom[], // TRUE if index atom is in some ring
                        int ring_size[],    // smallest size of ring this bond is in or 0 if none
                        neighbourhood_t nbp[])
/*
 * Makes O=C-N of large lactam amid carbonyls trigonal. This shall improve
 * readability of some large-ring structures.
 * It is only applied to lactam carbonyls where the N is not in a small ring.
 */
{
   int i, j, jj, nrings;
   struct reaccs_atom_t *apc, *apo, *apn, *apx;
   neighbourhood_t *nbph;

   for (i=0, apc=mp->atom_array, nbph=nbp; i<mp->n_atoms; i++, apc++, nbph++)
   {
      // check central atom
      if (nbph->n_ligands != 3) continue;               /* need 3 ligands */
      if (!is_ring_atom[i]) continue;                   /* only really large sp2 atoms */
      if (0 != strcmp("C", apc->atom_symbol)) continue; /* central atom needs to be C */
// fprintf(stderr,"found compatible 'C' atom %d\n", i+1);

      // search for compatible '=O' ligand
      for (j=0; j<nbph->n_ligands; j++)
          if (mp->bond_array[nbph->bonds[j]].bond_type == DOUBLE            &&
              0 == strcmp("O", mp->atom_array[nbph->atoms[j]].atom_symbol)  &&
              nbp[nbph->atoms[j]].n_ligands == 1)
              break;
      if (j == nbph->n_ligands) continue;   /* no carbonyl O */
      apo = mp->atom_array + nbph->atoms[j];
// fprintf(stderr,"found compatible 'O' ligand %d\n", nbph->atoms[j]+1);

      // search for compatible '-N' ligand
      for (j=0; j<nbph->n_ligands; j++)
          if (mp->bond_array[nbph->bonds[j]].bond_type == SINGLE            &&
              0 == strcmp("N", mp->atom_array[nbph->atoms[j]].atom_symbol)  &&
              ring_size[nbph->bonds[j]] > 12)
              break;
      if (j == nbph->n_ligands) continue;   /* no amid N */
      // check if this atom has only 2 ring ligands
      nrings = 0;
      for (jj=0; jj<nbp[nbph->atoms[j]].n_ligands; jj++)
          if (is_ring_atom[nbp[nbph->atoms[j]].atoms[jj]]) nrings++;
      if (nrings > 2) continue;
      apn = mp->atom_array + nbph->atoms[j];
// fprintf(stderr,"found compatible 'N' ligand %d\n", nbph->atoms[j]+1);

      // search for compatible '-C' ligand
      for (j=0; j<nbph->n_ligands; j++)
          if (mp->bond_array[nbph->bonds[j]].bond_type == SINGLE            &&
              0 == strcmp("C", mp->atom_array[nbph->atoms[j]].atom_symbol))
              break;
      if (j == nbph->n_ligands) continue;   /* no additional C ligand */
      // check if this atom has only 2 ring ligands
      nrings = 0;
      for (jj=0; jj<nbp[nbph->atoms[j]].n_ligands; jj++)
          if (is_ring_atom[nbp[nbph->atoms[j]].atoms[jj]]) nrings++;
      if (nrings > 2) continue;
      apx = mp->atom_array + nbph->atoms[j];
// fprintf(stderr,"found compatible 'C' ligand %d\n", nbph->atoms[j]+1);

      // now we have located the compatible amid fragment => set its coordinates
      apc->x = 0.0; apc->y = 0.0;
      apo->x =  1.514; apo->y =  0.0; apo->color = apc->color;
      apn->x = -1.514*cos(3.1415926/3); apn->y =  1.514*sin(3.1415926/3); apn->color = apc->color;
      apx->x = -1.514*cos(3.1415926/3); apx->y = -1.514*sin(3.1415926/3); apx->color = apc->color;
   }
}

void FlipRingCarbonyls(struct reaccs_molecule_t *mp,
                       int is_ring_atom[], // TRUE if index atom is in some ring
                       int ring_size[],    // smallest size of ring this bond is in or 0 if none
                       neighbourhood_t nbp[])
/*
 * Makes layed-out O=C-N of large lactam amid carbonyls point into ring. This shall further improve
 * readability of some large-ring structures.
 * It is only applied to lactam carbonyls where the N is not in a small ring.
 */
{
   int i, j, jj, nrings;
   struct reaccs_atom_t *apc, *apo, *apn, *apx;
   double xref, yref;
   neighbourhood_t *nbph;

   for (i=0, apc=mp->atom_array, nbph=nbp; i<mp->n_atoms; i++, apc++, nbph++)
   {
      // check central atom
      if (nbph->n_ligands != 3) continue;               /* need 3 ligands */
      if (!is_ring_atom[i]) continue;                   /* only really large sp2 atoms */
      if (0 != strcmp("C", apc->atom_symbol)) continue; /* central atom needs to be C */
// fprintf(stderr,"found compatible 'C' atom %d\n", i+1);

      // search for compatible '=O' ligand
      for (j=0; j<nbph->n_ligands; j++)
          if (mp->bond_array[nbph->bonds[j]].bond_type == DOUBLE            &&
              0 == strcmp("O", mp->atom_array[nbph->atoms[j]].atom_symbol)  &&
              nbp[nbph->atoms[j]].n_ligands == 1)
              break;
      if (j == nbph->n_ligands) continue;   /* no carbonyl O */
      apo = mp->atom_array + nbph->atoms[j];
// fprintf(stderr,"found compatible 'O' ligand %d\n", nbph->atoms[j]+1);

      // search for compatible '-N' ligand
      for (j=0; j<nbph->n_ligands; j++)
          if (mp->bond_array[nbph->bonds[j]].bond_type == SINGLE            &&
              0 == strcmp("N", mp->atom_array[nbph->atoms[j]].atom_symbol)  &&
              ring_size[nbph->bonds[j]] > 12)
              break;
      if (j == nbph->n_ligands) continue;   /* no amid N */
      // check if this atom has only 2 ring ligands
      nrings = 0;
      for (jj=0; jj<nbp[nbph->atoms[j]].n_ligands; jj++)
          if (is_ring_atom[nbp[nbph->atoms[j]].atoms[jj]]) nrings++;
      if (nrings > 2) continue;
      apn = mp->atom_array + nbph->atoms[j];
// fprintf(stderr,"found compatible 'N' ligand %d\n", nbph->atoms[j]+1);

      // search for compatible '-C' ligand
      for (j=0; j<nbph->n_ligands; j++)
          if (mp->bond_array[nbph->bonds[j]].bond_type == SINGLE            &&
              0 == strcmp("C", mp->atom_array[nbph->atoms[j]].atom_symbol))
              break;
      if (j == nbph->n_ligands) continue;   /* no additional C ligand */
      // check if this atom has only 2 ring ligands
      nrings = 0;
      for (jj=0; jj<nbp[nbph->atoms[j]].n_ligands; jj++)
          if (is_ring_atom[nbp[nbph->atoms[j]].atoms[jj]]) nrings++;
      if (nrings > 2) continue;
      apx = mp->atom_array + nbph->atoms[j];
// fprintf(stderr,"found compatible 'C' ligand %d\n", nbph->atoms[j]+1);
      if (apc->color != apo->color) continue;
      if (apc->color != apn->color) continue;
      if (apc->color != apx->color) continue;

      // now we have located the compatible amid fragment => set its coordinates
      xref = (apx->x+apn->x)/2.0; yref = (apx->y+apn->y)/2.0;
      apc->x = xref-(apc->x-xref); apc->y = yref-(apc->y-yref);
      apo->x = xref-(apo->x-xref); apo->y = yref-(apo->y-yref);
   }
}

void SproutRingSubstituents(struct reaccs_molecule_t *mp,
                            int is_ring_atom[],
                            int ring_size[],
                            neighbourhood_t nbp[])
/*
 * Attaches the first atom of non-ring substituents to already layouted
 * ring fragments in a symmetrical way. This is necessary because
 * step-by-step addition does not produce sufficient results.
 * Only those attachment points are considered which have exclusively
 * non-ring attachments.
 */
{
   atom_pair *edges;
   int nedges;
   point *coords, p1, p2, d, dp, c;
   int nnodes;
   struct reaccs_atom_t *ap;
   int seed;
   int i, j, k, iatom;
   int *numbers;

   double sina, cosa;
   int subst_indices[MAXNEIGHBOURS], h;
   int n_subst, n_rbonds, min_rsize;
   int problem_substituent;

   if (mp->n_bonds < 1) return;	/* trivial molecule */

   /* loop over all ring atoms */
   for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
   {
      if (!is_ring_atom[j]) continue;
      n_subst = 0; n_rbonds = 0;
      problem_substituent = NONE;
      min_rsize = 100;
      for (i=0; i<nbp[j].n_ligands; i++)
         if (!(mp->bond_array[nbp[j].bonds[i]].bond_type & RUBBER_BOND))
            if (ring_size[nbp[j].bonds[i]] != 0)
            {
               if (ring_size[nbp[j].bonds[i]] < min_rsize)
                  min_rsize = ring_size[nbp[j].bonds[i]];
               n_rbonds++;
            }
            else
            {
               iatom = nbp[j].atoms[i];
               if (is_ring_atom[iatom])
               {
                  problem_substituent = MAX(problem_substituent,NONE+1);
               }
               else
               {
                  if (CountColor(mp, mp->atom_array[iatom].color) != 1)
                     problem_substituent = MAX(problem_substituent,NONE+2);
                  subst_indices[n_subst] = iatom;
                  n_subst++;
               }
            }
      if (n_subst == 0)                 /* nothing to be layed out */
      {
//          fprintf(stderr, "atom %d: n_subst = %d\n", j+1, n_subst);
         continue;
      }
      if (n_rbonds != 2)                /* ring fusion atom */
      {
//         fprintf(stderr, "atom %d: n_rbonds = %d\n", j+1, n_rbonds);
         continue;
      }
      if (problem_substituent)          /* substituent is part of other */
      {                                 /* ring or already layed out */
//         if (problem_substituent == NONE+1)
//            fprintf(stderr, "atom %d: ring substituent.\n", j+1);
//         else
//            fprintf(stderr, "atom %d: already layed out.\n", j+1);
         continue;
      }
//      fprintf(stderr, "atom %d: layed out in progress.\n", j+1);

      numbers = TypeAlloc(mp->n_atoms, int);
      edges = TypeAlloc(mp->n_bonds, atom_pair);
      coords = TypeAlloc(mp->n_atoms, point);

      GetColoredGraph(mp,
                      edges,  &nedges, coords, &nnodes,
                      numbers, ap->color);
      seed = numbers[j];
      p1[0] = coords[seed][0];
      p1[1] = coords[seed][1];
// if (j+1 == -1)
// {
// DebugGeometry(TRUE);
// fprintf(stderr,"sprouting at '%s' atom %d\n", ap->atom_symbol, j+1);
// }
      if (min_rsize <= 8)
         NextSubstituentPoint(p2, coords, nnodes, edges, nedges, seed, 0, numbers, mp->n_atoms);
      else
         NextSubstituentPoint(p2, coords, nnodes, edges, nedges, seed, USE_INWARDS, numbers, mp->n_atoms);
// if (j+1 == -1) DebugGeometry(FALSE);
      d[0] = p2[0]-p1[0];
      d[1] = p2[1]-p1[1];
      if (n_subst > 1)
      {
         for (i=1; i<n_subst; i++)	       /* sort substituents by size */
            for (k=i-1; k>=0; k--)
               if (nbp[subst_indices[k]].n_ligands <
                   nbp[subst_indices[k+1]].n_ligands)
               {
                  h = subst_indices[k];
                  subst_indices[k] = subst_indices[k+1];
                  subst_indices[k+1] = h;
               }
               else
                  break;
         cosa = cos(((0.5)*PI/n_subst-PI/2)*80.0/90.0);
         sina = sin(((0.5)*PI/n_subst-PI/2)*80.0/90.0);
         dp[0] = cosa*d[0] - sina*d[1] + p1[0];
         dp[1] = sina*d[0] + cosa*d[1] + p1[1];
         c[0] = c[1] = 0.0;
         for (i=0; i<nnodes; i++)	/* compute center of ring system */
         {
            c[0] += coords[i][0]; c[1] += coords[i][1];
         }
         c[0] /= nnodes; c[1] /= nnodes;
         /* if first substituent is more distant from center than others */
         if ((dp[0]-c[0])*(dp[0]-c[0]) + (dp[1]-c[1])*(dp[1]-c[1]) <
             (p2[0]-c[0])*(p2[0]-c[0]) + (p2[1]-c[1])*(p2[1]-c[1]))
         {	/* swap extreme substituents */
            h = subst_indices[0];
            subst_indices[0] = subst_indices[n_subst-1];
            subst_indices[n_subst-1] = h;
         }
      }
      for (i=0; i<n_subst; i++)
      {			/* scale for 80 deg. with two substituents */
         cosa = cos(((i+0.5)*PI/n_subst-PI/2)*80.0/90.0);
         sina = sin(((i+0.5)*PI/n_subst-PI/2)*80.0/90.0);
         dp[0] = cosa*d[0] - sina*d[1] + p1[0];
         dp[1] = sina*d[0] + cosa*d[1] + p1[1];
         mp->atom_array[subst_indices[i]].x = dp[0];
         mp->atom_array[subst_indices[i]].y = dp[1];
         mp->atom_array[subst_indices[i]].color = ap->color;
      }

      MyFree((char *)numbers);
      MyFree((char *)edges);
      MyFree((char *)coords);
   }
}

void LinkRingFragments(struct reaccs_molecule_t *mp,
                       int is_ring_atom[])      // TRUE if atom is in ring, IO=0
/*
 * Links directly connected ring fragments.
 */
{
   atom_pair *edges;
   int nedges;
   point *coords, p1, p2, p1p, p2p;
   int nnodes;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int col1, col2;
   int n1, n2, h;
   int seed;
   int i, j;
   int *numbers;

   if (mp->n_bonds < 1) return;

   for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
   {
      if (mp->atom_array[bp->atoms[0]-1].color ==
          mp->atom_array[bp->atoms[1]-1].color  ||
              bp->bond_type & RUBBER_BOND           ||
              !is_ring_atom[bp->atoms[0]-1]         ||
              !is_ring_atom[bp->atoms[1]-1])
             continue;

          /* bp is now pointing to a bond linking to ring fragments */
      numbers = TypeAlloc(mp->n_atoms, int);
      edges = TypeAlloc(mp->n_bonds, atom_pair);
      coords = TypeAlloc(mp->n_atoms, point);

      col1 = mp->atom_array[bp->atoms[0]-1].color;
      col2 = mp->atom_array[bp->atoms[1]-1].color;

      for (i=0, n1=0; i<mp->n_atoms; i++)    /* find reference fragment */
             if (mp->atom_array[i].color == col1) n1++;
      for (i=0, n2=0; i<mp->n_atoms; i++)
             if (mp->atom_array[i].color == col2) n2++;
      if (n1 > n2)
      {
             h = bp->atoms[0]; bp->atoms[0] = bp->atoms[1]; bp->atoms[1] = h;
             h = col1; col1 = col2; col2 = h;
             h = n1; n1 = n2; n2 = h;
      }

      /* col2 colors larger/reference fragment */
      GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, col2);
      seed = numbers[bp->atoms[1]-1];
      p2p[0] = coords[seed][0];
      p2p[1] = coords[seed][1];
      NextSubstituentPoint(p1p, coords, nnodes, edges, nedges, seed, NO_RAY_CROSSING|USE_INWARDS, numbers, mp->n_atoms);

      GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, col1);
      seed = numbers[bp->atoms[0]-1];
      p1[0] = coords[seed][0];
      p1[1] = coords[seed][1];
      NextSubstituentPoint(p2, coords, nnodes, edges, nedges, seed, 0, numbers, mp->n_atoms);

      TransformPoints(coords, nnodes, p1, p2, p1p, p2p);

      for (i=0; i<mp->n_atoms; i++)
         if (mp->atom_array[i].color == col1)
         {
            mp->atom_array[i].x = coords[numbers[i]][0];
            mp->atom_array[i].y = coords[numbers[i]][1];
         }

      MyFree((char *)numbers);
      MyFree((char *)edges);
      MyFree((char *)coords);

      ImproveBondByFlip(mp, bp, 0.01);
      ImproveBondByStretch(mp, bp);

      				/* merge colors */
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == col2) ap->color = col1;
   }
}

void LinkRemainingFragments(struct reaccs_molecule_t *mp)
/*
 * Links fragments not yet layed-out.
 *
 * This can happen because of previous layout rules too strict to cover all cases as well as for placing
 * actually implicit hydrogens.
 */
{
   atom_pair *edges;
   int nedges;
   point *coords, p1, p2, p1p, p2p;
   int nnodes;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int col1, col2;
   int n1, n2, h;
   int seed;
   int i, j;
   int *numbers;

   if (mp->n_bonds < 1) return;

   /* make sure hydrogen rubber bond tags are cleared */
   for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
   {
      if (mp->atom_array[bp->atoms[0]-1].color != mp->atom_array[bp->atoms[1]-1].color  &&
          bp->bond_type & RUBBER_BOND                                                   &&
          (AtomSymbolMatch(mp->atom_array[bp->atoms[0]-1].atom_symbol, "H,D,T")  ||
           AtomSymbolMatch(mp->atom_array[bp->atoms[1]-1].atom_symbol, "H,D,T")))
          bp->bond_type &= ~RUBBER_BOND;
   }

   for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
   {
      if (mp->atom_array[bp->atoms[0]-1].color ==
          mp->atom_array[bp->atoms[1]-1].color  ||
          bp->bond_type & RUBBER_BOND)
         continue;
      /* bp is now pointing to a bond linking to fragments */

      numbers = TypeAlloc(mp->n_atoms, int);
      edges = TypeAlloc(mp->n_bonds, atom_pair);
      coords = TypeAlloc(mp->n_atoms, point);

      col1 = mp->atom_array[bp->atoms[0]-1].color;
      col2 = mp->atom_array[bp->atoms[1]-1].color;

      for (i=0, n1=0; i<mp->n_atoms; i++)    /* find reference fragment */
         if (mp->atom_array[i].color == col1)
            n1++;
      for (i=0, n2=0; i<mp->n_atoms; i++)
         if (mp->atom_array[i].color == col2)
            n2++;
      if (n1 > n2)
      {
         h = bp->atoms[0]; bp->atoms[0] = bp->atoms[1]; bp->atoms[1] = h;
         h = col1; col1 = col2; col2 = h;
         h = n1; n1 = n2; n2 = h;
      }
      // bp->atoms[1]/col2/n2 is reference fragment

      GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, col2);
      seed = numbers[bp->atoms[1]-1];
      p2p[0] = coords[seed][0];
      p2p[1] = coords[seed][1];
      if (0 == strcmp(mp->atom_array[bp->atoms[1]-1].atom_symbol, "H"))
      {
// fprintf(stderr,"processing hydrogen atom %d\n", bp->atoms[1]);
         NextSubstituentPoint(p1p, coords, nnodes, edges, nedges, seed, USE_INWARDS|H_GEOMETRY, numbers, mp->n_atoms);
      }
      else
      {
// fprintf(stderr,"processing non-hydrogen atom %d\n", bp->atoms[1]);
         NextSubstituentPoint(p1p, coords, nnodes, edges, nedges, seed, USE_INWARDS, numbers, mp->n_atoms);
      }

      GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, col1);
      seed = numbers[bp->atoms[0]-1];
      p1[0] = coords[seed][0];
      p1[1] = coords[seed][1];
      if (0 == strcmp(mp->atom_array[bp->atoms[0]-1].atom_symbol, "H"))
      {
// fprintf(stderr,"processing hydrogen atom %d\n", bp->atoms[0]);
         NextSubstituentPoint(p2, coords, nnodes, edges, nedges, seed, USE_INWARDS|H_GEOMETRY, numbers, mp->n_atoms);
      }
      else
      {
// fprintf(stderr,"processing non-hydrogen atom %d\n", bp->atoms[0]);
         NextSubstituentPoint(p2, coords, nnodes, edges, nedges, seed, USE_INWARDS, numbers, mp->n_atoms);
      }

      TransformPoints(coords, nnodes, p1, p2, p1p, p2p);

      for (i=0; i<mp->n_atoms; i++)
         if (mp->atom_array[i].color == col1)
         {
            mp->atom_array[i].x = coords[numbers[i]][0];
            mp->atom_array[i].y = coords[numbers[i]][1];
         }

      MyFree((char *)numbers);
      MyFree((char *)edges);
      MyFree((char *)coords);

      				/* merge colors */
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == col1) ap->color = col2;
   }
}

int BranchQuality(struct reaccs_molecule_t *mp, int ai,
                  neighbourhood_t *nbp,
                  int used_atoms[])
/*
 * Computes a measure of the branchedness of atom with index ai
 * or atom number (ai+1). The more smaller substituents the better.
 * The quality function strongly prefers atoms connected to already
 * layouted atoms as indicated by used_atoms[].
 */
{
   int nbranch[MAXNEIGHBOURS];
   int j, k;
   int nsingle, nmore, ntotal;
   int can_be_candidate;
   neighbourhood_t *nbph;
   int placed_neighbour;
   int already_layouted;

   already_layouted = 0;
   can_be_candidate = FALSE; placed_neighbour = FALSE;
   for (j=0; j<nbp[ai].n_ligands; j++)
   {
      nbranch[j] = 0;
      if (mp->bond_array[nbp[ai].bonds[j]].bond_type & RUBBER_BOND)
         continue;
      if (mp->atom_array[ai].color != mp->atom_array[nbp[ai].atoms[j]].color)
         can_be_candidate = TRUE;
      else
         already_layouted++;
      if (used_atoms[nbp[ai].atoms[j]]) placed_neighbour = TRUE;
      nbph = &nbp[nbp[ai].atoms[j]];
      for (k=0; k<nbph->n_ligands; k++)
         if (!(mp->bond_array[nbph->bonds[k]].bond_type & RUBBER_BOND))
         {
            nbranch[j]++;
         }
   }

   if (!can_be_candidate  ||
       already_layouted > 1) return (-1);

   nsingle = nmore = ntotal = 0;
   for (j=0; j<nbp[ai].n_ligands; j++)
   {
      if (nbranch[j] == 1) nsingle++;
      else if (nbranch[j] > 0) nmore++;
      ntotal += nbranch[j];
   }
   if (placed_neighbour)
      return (1000+nsingle*100+nmore*10+ntotal);
   else
      return (nsingle*100+nmore*10+ntotal);
}

int ColorSize(struct reaccs_molecule_t *mp, int col)
/*
 * Counts the number of atoms in mp which have color col.
 */
{
   int i, result;

   result = 0;
   for (i=0; i<mp->n_atoms; i++)
      if (mp->atom_array[i].color == col) result++;

   return (result);
}

void LayoutChainAtoms(struct reaccs_molecule_t *mp, int is_ring_atom[], neighbourhood_t *nbp)
/*
 * Layouts the environment of those atoms which have only chain
 * substituents.
 */
{
   int i, j;
   int col[MAXNEIGHBOURS], oldcolor;
   int size[MAXNEIGHBOURS];
   int index[MAXNEIGHBOURS];
   int bindex[MAXNEIGHBOURS];
   int nneigh;
   int quality, ibest;

   atom_pair *edges;
   int nedges;
   point *coords;
   point p1, p2, p1p, p2p;
   int nnodes;
   struct reaccs_atom_t *ap;
   int h;
   int seed;
   int *numbers;
   int *used_atoms;

   if (mp->n_bonds == 0) return;

   numbers = TypeAlloc(mp->n_atoms, int);
   used_atoms = TypeAlloc(mp->n_atoms, int);
   edges = TypeAlloc(mp->n_bonds, atom_pair);
   coords = TypeAlloc(mp->n_atoms, point);

   for (i=0; i<mp->n_atoms; i++)
      used_atoms[i] = FALSE;

   for (;;)	/* get next best chain atom until all have been used */
   {
      quality = (-1); ibest = (-1);
      for (i=0; i<mp->n_atoms; i++)
         if (!is_ring_atom[i])
         {
            if (quality < BranchQuality(mp, i, nbp, used_atoms))
            {
               ibest = i; quality = BranchQuality(mp, i, nbp, used_atoms);
            }
         }
      if (quality < 0) break;
      used_atoms[ibest] = TRUE;
      /*
       * Here, ibest is the index of the next chain atom around
       * which the ligands are to be arranged.
       */

      /*
       * Now, we find the neighbours of this atom that are not connected
       * by so-called rubber bonds (dative bonds to be layed out last).
       */
      nneigh = 0;
      oldcolor = mp->atom_array[ibest].color;
      mp->atom_array[ibest].color = (-1);
      for (i=0; i<nbp[ibest].n_ligands; i++)
         if (!(mp->bond_array[nbp[ibest].bonds[i]].bond_type & RUBBER_BOND))
         {
            index[nneigh] = nbp[ibest].atoms[i];
            bindex[nneigh] = nbp[ibest].bonds[i];
            col[nneigh]  = mp->atom_array[index[nneigh]].color;
            size[nneigh] = 100*ColorSize(mp, col[nneigh]) +
                           nbp[index[nneigh]].n_ligands;
            nneigh++;
         }

      for (i=1; i<nneigh; i++)	/* already layed out branch is reference */
         if (col[i] == oldcolor)
         {
            h = index[i]; index[i] = index[0]; index[0] = h;
            h = bindex[i]; bindex[i] = bindex[0]; bindex[0] = h;
            h = col[i]; col[i] = col[0]; col[0] = h;
            h = size[i]; size[i] = size[0]; size[0] = h;
         }

      for (;;)			/* optimize for space usage */
      {
         for (i=1; i<nneigh-1; i++)
            if (ABS(size[i]-size[i-1]) + ABS(size[i+1]-size[(i+2)%nneigh]) <
                ABS(size[i]-size[(i+2)%nneigh]) + ABS(size[i+1]-size[i-1]))
            {
               h = index[i]; index[i] = index[i+1]; index[i+1] = h;
               h = bindex[i]; bindex[i] = bindex[i+1]; bindex[i+1] = h;
               h = col[i]; col[i] = col[i+1]; col[i+1] = h;
               h = size[i]; size[i] = size[i+1]; size[i+1] = h;
               break;
            }
         if (i >= nneigh-1) break;
      }

      if (nneigh == 1)		/* Special case for terminal atoms. */
      {
         /*
          * We can make the other fragment the reference and place
          * the terminal atom at the next substitution point from there.
          */
         GetColoredGraph(mp,
                         edges,  &nedges,
                         coords, &nnodes, numbers, col[0]);
         seed = numbers[index[0]];
         NextSubstituentPoint(p1, coords, nnodes, edges, nedges, seed, USE_INWARDS, numbers, mp->n_atoms);
         mp->atom_array[ibest].x = p1[0];
         mp->atom_array[ibest].y = p1[1];
         mp->atom_array[ibest].color = mp->atom_array[index[0]].color;
      }
      else
      {
         /* layout reference branch */
         GetColoredGraph(mp,
                         edges,  &nedges,
                         coords, &nnodes, numbers, col[0]);
         p1p[0] = cos((2*PI*0)/nneigh)*STDBOND;
         p1p[1] = sin((2*PI*0)/nneigh)*STDBOND;
         p2p[0] = 0.0; p2p[1] = 0.0;
         p1[0] = mp->atom_array[index[0]].x;
         p1[1] = mp->atom_array[index[0]].y;
         if (col[0] == oldcolor)	/* real reference branch */
         {
            p2[0] = mp->atom_array[ibest].x;
            p2[1] = mp->atom_array[ibest].y;
         }
         else			/* no reference branch */
         {
            seed = numbers[index[0]];
            NextSubstituentPoint(p2, coords, nnodes, edges, nedges, seed, USE_INWARDS, numbers, mp->n_atoms);
         }
         TransformPoints(coords, nnodes, p1, p2, p1p, p2p);
         for (j=0; j<mp->n_atoms; j++)
            if (mp->atom_array[j].color == col[0])
            {
               mp->atom_array[j].x = coords[numbers[j]][0];
               mp->atom_array[j].y = coords[numbers[j]][1];
            }
         mp->atom_array[ibest].x = 0;
         mp->atom_array[ibest].y = 0;

         /* layout other substituents */
         for (i=1; i<nneigh; i++)
         {
            GetColoredGraph(mp,
                            edges,  &nedges,
                            coords, &nnodes, numbers, col[i]);
            p1p[0] = cos((2*PI*i)/MAX(nneigh,3))*STDBOND;
            p1p[1] = sin((2*PI*i)/MAX(nneigh,3))*STDBOND;
            p2p[0] = 0.0; p2p[1] = 0.0;
            p1[0] = mp->atom_array[index[i]].x;
            p1[1] = mp->atom_array[index[i]].y;
            seed = numbers[index[i]];
            NextSubstituentPoint(p2, coords, nnodes, edges, nedges, seed, USE_INWARDS, numbers, mp->n_atoms);
            TransformPoints(coords, nnodes, p1, p2, p1p, p2p);
            for (j=0; j<mp->n_atoms; j++)
               if (mp->atom_array[j].color == col[i])
               {
                  mp->atom_array[j].x = coords[numbers[j]][0];
                  mp->atom_array[j].y = coords[numbers[j]][1];
               }

         }
                                   /* merge colors */
         for (i=0; i<nneigh; i++)
            for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
               if (ap->color == col[i]) ap->color = oldcolor;
         mp->atom_array[ibest].color = oldcolor;
      }

      if (nneigh == 2  &&
          (mp->bond_array[bindex[0]].bond_type +
           mp->bond_array[bindex[1]].bond_type > SINGLE+SINGLE))
         ImproveSPAtoms(mp, nbp, ibest+1);
      for (i=0; i<nneigh; i++)
         if (mp->bond_array[bindex[i]].bond_type == TRIPLE)
            ImproveSPAtoms(mp, nbp, index[i]+1);

      for (i=0; i<nneigh; i++)
         if (size[i] > 200)	/* fragment was already layed out */
            if (!is_ring_atom[index[i]])
               MakeBondTrans(mp, nbp, &mp->bond_array[bindex[i]]);
   }

   MyFree((char *)used_atoms);
   MyFree((char *)numbers);
   MyFree((char *)edges); MyFree((char *)coords);
}

void LayoutRubberFragments(struct reaccs_molecule_t *mp)
/*
 * Places the fragments connected by rubber bonds in a nice way.
 */
{
   int i, j;

   atom_pair *edges;
   int nedges;
   point *coords;
   int nnodes;
   point p1, p1p;
   point p2, p2p;
   struct reaccs_bond_t *bp, *bph;
   int *numbers;
   double stretch_factor;

   int h;
   int col1, col2;
   int n1, n2;

   if (mp->n_bonds == 0) return;

   for (;;)	/* for all rubber bond fragment connections */
   {
            /* look for rubber connected fragments */
            /* with different color and a root atom i.e. an atom with only RUBBER_BOND connections */
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->bond_type & RUBBER_BOND  &&
             // don't layout shortcut fragments
             (0 != strcmp("R", mp->atom_array[bp->atoms[0]-1].atom_symbol)  ||  mp->atom_array[bp->atoms[0]-1].atext[0] == '\0')  &&
             (0 != strcmp("R", mp->atom_array[bp->atoms[1]-1].atom_symbol)  ||  mp->atom_array[bp->atoms[1]-1].atext[0] == '\0')  &&
             mp->atom_array[bp->atoms[0]-1].color > 0  &&
             mp->atom_array[bp->atoms[1]-1].color > 0  &&
             mp->atom_array[bp->atoms[0]-1].color  !=
             mp->atom_array[bp->atoms[1]-1].color)
         {
            col1 = mp->atom_array[bp->atoms[0]-1].color;
            col2 = mp->atom_array[bp->atoms[1]-1].color;
            n1 = n2 = 0;
            for (j=0, bph=mp->bond_array; j<mp->n_bonds; j++, bph++)
            {
               if (bph->atoms[0] == bp->atoms[0]  &&
                   mp->atom_array[bph->atoms[1]-1].color == col1)
                  n1++;
               if (bph->atoms[1] == bp->atoms[0]  &&
                        mp->atom_array[bph->atoms[0]-1].color == col1)
                  n1++;
               if (bph->atoms[0] == bp->atoms[1]  &&
                        mp->atom_array[bph->atoms[1]-1].color == col2)
                  n2++;
               if (bph->atoms[1] == bp->atoms[1]  &&
                        mp->atom_array[bph->atoms[0]-1].color == col2)
                  n2++;
            }
            if (n2 == 0  &&  n1 > 0)
            {
               break;
            }
            else if (n1 ==0  &&  n2 > 0)
            {
               h = n1; n1 = n2; n2 = h;
               h = col1; col1 = col2; col2 = h;
               h = bp->atoms[0];
               bp->atoms[0] = bp->atoms[1];
               bp->atoms[1] = h;
               break;
            }
         }
      if (i == mp->n_bonds) break;	/* no such fragments -> cleanup */
      // fragment 2 is the root atom
// fprintf(stderr, "n1=%d, n2=%d\n", n1, n2);

      numbers = TypeAlloc(mp->n_atoms, int);
      edges = TypeAlloc(mp->n_bonds, atom_pair);
      coords = TypeAlloc(mp->n_atoms, point);

      GetColoredGraph(mp,
                      edges,  &nedges,
                      coords, &nnodes, numbers, col1);
// fprintf(stderr, "fragment 1 has %d nodes\n", nnodes);
      NextSubstituentPoint(p2p,
                           coords, nnodes,
                           edges, nedges,
                           numbers[bp->atoms[0]-1], USE_INWARDS, numbers, mp->n_atoms);
      p1p[0] = coords[numbers[bp->atoms[0]-1]][0];
      p1p[1] = coords[numbers[bp->atoms[0]-1]][1];

      bp->bond_type &= ~RUBBER_BOND;

      p2[0] = p2[1] = 0.0;
      n2 = 0;
      for (i=0; i<mp->n_atoms; i++)
         if (mp->atom_array[i].color == col2)
         {
            p2[0] += mp->atom_array[i].x;
            p2[1] += mp->atom_array[i].y;
            n2++;
         }
      p2[0] /= n2; p2[1] /= n2;

      GetColoredGraph(mp,
                      edges,  &nedges,
                      coords, &nnodes, numbers, col2);
// fprintf(stderr, "fragment 2 has %d nodes\n", nnodes);
      coords[nnodes][0] = p2[0]; coords[nnodes][1] = p2[1];
      nnodes++;
      NextSubstituentPoint(p1, coords, nnodes, edges, nedges, nnodes-1, USE_INWARDS, numbers, mp->n_atoms);
      stretch_factor = sqrt((double)nnodes/1.5);
      p1[0] = p2[0]   + stretch_factor*(p1[0]-p2[0]);
      p1[1] = p2[1]   + stretch_factor*(p1[1]-p2[1]);
      p2p[0] = p1p[0] + stretch_factor*(p2p[0]-p1p[0]);
      p2p[1] = p1p[1] + stretch_factor*(p2p[1]-p1p[1]);

      TransformPoints(coords, nnodes, p1, p2, p1p, p2p);
      for (j=0; j<mp->n_atoms; j++)
         if (mp->atom_array[j].color == col2)
         {
            mp->atom_array[j].x = coords[numbers[j]][0];
            mp->atom_array[j].y = coords[numbers[j]][1];
            /* Rubber fragments don't count for layout of other atoms */
            mp->atom_array[j].color = -col1;
         }

      MyFree((char *)coords);
      MyFree((char *)edges);
      MyFree((char *)numbers);
   }
   // make sure all freshly layed out ligands of the root atom get their final color
   for (j=0; j<mp->n_atoms; j++)
      if (mp->atom_array[j].color < 0)
         mp->atom_array[j].color *= -1;
}

struct frag_desc_t
   {
      double xll, yll, xur, yur;	/* bounding rectangle */
      double xoffset, yoffset;		/* correction */
      int color;
      int natoms;
      struct frag_desc_t *next;
   };

int BestBondColor(struct reaccs_molecule_t *mp,
                  int *nbestp)
/*
 * Find the color of the best (highest color) bond in mp.
 * Returns this color and sets (*nbestp) to the number of these
 * bonds.
 */
{
   struct reaccs_bond_t *bp;
   int i, best_color;

   if (mp->n_bonds == 0)
   {
      (*nbestp) = 0;
      return (0);
   }

   best_color = mp->bond_array[0].color; (*nbestp) = 1;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->color > best_color)
      {
         best_color = bp->color;
         (*nbestp) = 1;
      }
      else if (bp->color == best_color)
         (*nbestp)++;

   return (best_color);
}

double Diameter(struct reaccs_molecule_t *mp,
                struct reaccs_bond_t *bp,
                point center,
                double *proj)
/*
 * Computes and returns the diameter of the molecule *mp perpendicular
 * to the bond *bp. It also sets *proj to the projection of
 * the direction from the center of *bp to center[] onto the (1,1) direction
 * of the local coordinate system around bp->atom[0] with the bond
 * assumed to point to the y-direction. This projection can be
 * used to rank the bonds being tested.
 * bp is made to point a parallel as possible in the direction of the center.
 */
{
   struct reaccs_atom_t *ap;
   int color, i;
   point p1, p2, e12, e11;
   double x, minx, maxx;
   double l;

   p1[0] = mp->atom_array[bp->atoms[0]-1].x;
   p1[1] = mp->atom_array[bp->atoms[0]-1].y;
   p2[0] = mp->atom_array[bp->atoms[1]-1].x;
   p2[1] = mp->atom_array[bp->atoms[1]-1].y;

   if ((p2[0]-p1[0])*(center[0]-(p1[0]+p2[0])/2) +
       (p2[1]-p1[1])*(center[1]-(p1[1]+p2[1])/2) < 0)
   {
      i = bp->atoms[0];
      bp->atoms[0] = bp->atoms[1];
      bp->atoms[1] = i;
      p1[0] = mp->atom_array[bp->atoms[0]-1].x;
      p1[1] = mp->atom_array[bp->atoms[0]-1].y;
      p2[0] = mp->atom_array[bp->atoms[1]-1].x;
      p2[1] = mp->atom_array[bp->atoms[1]-1].y;
   }

   e12[1] = p2[0] - p1[0]; e12[0] = -(p2[1] - p1[1]);
   e11[0] = (p2[0] - p1[0]) + (p2[1] - p1[1]);
   e11[1] = (p2[1] - p1[1]) - (p2[0] - p1[0]);
   l = sqrt(e12[0]*e12[0] + e12[1]*e12[1]);
   if (l < 0.00001)	/* bond too short */
   {
      (*proj) = 0.0;
      return (0.0);
   }
   e12[0] /= l; e12[1] /= l;
   e11[0] /= l*sqrt(2.0); e11[1] /= l*sqrt(2.0);
   /* e12[] is a unit vector perpendicular to *bp rotated to the right */
   /* e11[] is rotated 45 deg. back towards the bond */

   color = mp->atom_array[bp->atoms[0]-1].color;
   minx = 1.0e7; maxx = -1.0e7;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == color)
      {
         x = (ap->x - p1[0])*e12[0]  +  (ap->y - p1[1])*e12[1];
         if (x < minx) minx = x;
         if (x > maxx) maxx = x;
      }

   (*proj) = (center[0]-(p1[0]+p2[0])/2)*e11[0] +
             (center[1]-(p1[1]+p2[1])/2)*e11[1];
   return (maxx - minx);
}

void OrientFragment(struct reaccs_molecule_t *mp,
                    int color,
                    int ring_size[],
                    int ring_count[])
/*
 * Orient the fragment of *mp with color according to a set of
 * rules. The major goal of this function is to make the orientation
 * as independent as possible of the small details of the structure
 * by refering to the ring systems, to make it as "landscape" as
 * possible, and to make some "preferred bond" vertical.
 * The arrays ring_size[] and ring_count[] are indexed by bond and
 * used to select this preferred bond.
 */
{
   struct reaccs_bond_t *bp, *bpbest;
   struct reaccs_atom_t *ap;
   int i, best_color, nbest;
   double n;
   double bestscore, diameter;
   point center, p;
   double projection;
   double cosa, sina;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (mp->atom_array[bp->atoms[0]-1].color != color  ||
          mp->atom_array[bp->atoms[1]-1].color != color)
         bp->color = 0; 		/* not in fragment */
      else
         bp->color = 1 +
                     ring_size[i] +	/* prefer larger rings */
                     10*ring_count[i];	/* prefer ring (fusion) bonds */

   best_color = BestBondColor(mp, &nbest);

   if (best_color == 0) return;		/* no bond in fragment */

   if (best_color == 1)	/* only chain bonds in fragment */
   {
      n = 0;  /* compute center of fragment */
      center[0] = center[1] = 0.0;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == color)
         {
            n++;
            center[0] += ap->x; center[1] += ap->y;
         }
      center[0] /= n; center[1] /= n;

   }
   else		/* rings in fragment */
   {
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (ring_size[i] > 0)
         {
            if (mp->atom_array[bp->atoms[0]-1].color == color)
               mp->atom_array[bp->atoms[0]-1].color *= (-1);
            if (mp->atom_array[bp->atoms[1]-1].color == color)
               mp->atom_array[bp->atoms[1]-1].color *= (-1);
         }

      n = 0;  /* compute center of gravity of fragment atoms*/
      center[0] = center[1] = 0.0;
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == -color)	/* prefer ring atoms */
         {
            if (0 == strcmp(ap->atom_symbol,"C"))
            {
               n+=200;
               center[0] += 200*ap->x; center[1] += 200*ap->y;
            }
            else	/* give ring hetero atoms more weight */
            {
               n+=1000;
               center[0] += 1000*ap->x; center[1] += 1000*ap->y;
            }
         }
         else if (ap->color == color)	/* use numbering to break symmetry */
         {
            n+=(1.0+0.01*i);
            center[0] += (1.0+0.01*i)*ap->x;
            center[1] += (1.0+0.01*i)*ap->y;
         }
      center[0] /= n; center[1] /= n;
   }


   /* candidate reference bonds are now colored with best_color and   */
   /* the center is the point to be put into the upper right quadrant */

   bestscore = -1;
   for (i=0, bpbest=bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->color == best_color)
      {	/* use bonds only from this fragment */
         diameter = Diameter(mp, bp, center, &projection);
         if (diameter*diameter+0.1*projection*projection > bestscore)
         {
            bestscore = diameter*diameter+0.1*projection*projection;
            bpbest = bp;
         }
      }
   diameter = Diameter(mp, bpbest, center, &projection);

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color < 0) ap->color *= (-1);

   /* Now, bpbest is the reference bond for reorientation. projection */
   /* can be used to fix the orientation of the reference bond. */

   /* First we move bpbest->atoms[0] to the origin */
   p[0] = mp->atom_array[bpbest->atoms[0]-1].x;
   p[1] = mp->atom_array[bpbest->atoms[0]-1].y;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == color)
      {
         ap->x -= p[0]; ap->y -= p[1];
      }

   /* Now, we compute the rotation matrix */
   p[0] = mp->atom_array[bpbest->atoms[1]-1].x;
   p[1] = mp->atom_array[bpbest->atoms[1]-1].y;
   cosa =  p[1]/sqrt(p[0]*p[0]+p[1]*p[1]);
   sina = -p[0]/sqrt(p[0]*p[0]+p[1]*p[1]);
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == color)
      {
         p[0] = ap->x; p[1] = ap->y;
         ap->x =  cosa*p[0] + sina*p[1];
         ap->y = -sina*p[0] + cosa*p[1];
      }

   /* if projection was negative -> make mirror image of fragment */
   if (projection < 0)
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == color)
            ap->x =  -ap->x;

   return;
}

void LayoutFragments(struct reaccs_molecule_t *mp,
                     int ring_size[], int ring_count[])
/*
 * Places the fragments of the molecule with respect to
 * each other. A part of the molecule with the same color
 * is considered a fragment even if it consists of more than
 * one disconnected component.
 */
{
   struct frag_desc_t *flist, *fp, *fph, ftmp;
   int i;
   double xoffset;
   struct reaccs_atom_t *ap;

   			/* compute fragment descriptions */
   flist = (struct frag_desc_t *)NULL;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      for (fp=flist; fp; fp=fp->next)
         if (fp->color == ap->color) break;
      if (fp)	/* fragment already found */
      {
         fp->natoms++;
      }
      else	/* new fragment */
      {
         fp = TypeAlloc(1, struct frag_desc_t);
         fp->natoms = 1; fp->color = ap->color;
         fp->xll = +1.0e7; fp->xur = -1.0e7;
         fp->yll = +1.0e7; fp->yur = -1.0e7;
         fp->next = flist; flist = fp;
      }
   }

                                      /* bubble sort by size of fragment */
   for (fph = flist; fph && fph->next; fph=fph->next)
      for (fp=flist; fp && fp->next; fp=fp->next)
         if (fp->natoms < fp->next->natoms)
         {
            ftmp.natoms = fp->next->natoms; fp->next->natoms = fp->natoms;
            fp->natoms = ftmp.natoms;
            ftmp.color = fp->next->color; fp->next->color = fp->color;
            fp->color = ftmp.color;
         }

   				/* orient fragments */
   for (fp=flist; fp; fp=fp->next)
      if (fp->natoms > 1  &&  (fp->color & KEEP_POSITION) == 0)
         OrientFragment(mp, fp->color, ring_size, ring_count);

                                /* get bounding boxes */
   for (fp=flist; fp; fp=fp->next)
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == fp->color)
         {
            if (ap->x < fp->xll) fp->xll = ap->x;
            if (ap->y < fp->yll) fp->yll = ap->y;
            if (ap->x > fp->xur) fp->xur = ap->x;
            if (ap->y > fp->yur) fp->yur = ap->y;
         }

              /* frame fragment bounding boxes by 0.5 standard bond length */
   for (fp=flist; fp; fp=fp->next)
   {
      fp->xll -= STDBOND/2; fp->yll -= STDBOND/2;
      fp->xur += STDBOND/2; fp->yur += STDBOND/2;
      fp->xoffset = -fp->xll; fp->yoffset = -fp->yll;
   }

   xoffset = 0;
   for (fp=flist; fp; fp=fp->next)
   {
      fp->xoffset += xoffset;
      xoffset += fp->xur - fp->xll;
   }

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      for (fp=flist; fp; fp=fp->next)
         if (fp->color == ap->color) break;

      if (fp)
      {
         ap->x += fp->xoffset; ap->y += fp->yoffset;
      }
   }

   while (flist)
   {
      fp = flist->next; MyFree((char *)flist); flist = fp;
   }
}

#define CHAIN_ONLY		0
#define ALLOW_RING_BONDS	1
#define USE_ALL			2

void LayoutAtomStereo(struct reaccs_molecule_t *mp,
                      neighbourhood_t *nbp,
                      int is_ring_atom[],
                      int ring_size[])
/*
 * Uses the stereoparity information of the atoms of *mp
 * to assign hashed ad wedged bonds.
 */
{
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int ligands[4], parity;
   int i, j, k, h;
   int i0, i1, i2, i3;
   int d[4];
   struct npoint_t tetra[4];
   double x, y;	/* coordinates of root atom */
   double maxvol;
   char buffer[80];
   int irun;
   int nclose;

   ResetColors(mp);

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (ap->stereo_parity == NONE) continue;
      if (nbp[i].n_ligands == 3)	/* not to be expected for Daylight */
      {
         sprintf(buffer, "Unexpected stereochemistry with three ligands at atom %d", i+1);
         ShowMessage(buffer, "LayoutAtomStereo");
         continue;
      }
      else if (nbp[i].n_ligands == 4)
      {
         x = mp->atom_array[i].x;
         y = mp->atom_array[i].y;
         for (j=0; j<4; j++)	/* save coordinate and index information */
         {
            tetra[j].number = nbp[i].bonds[j];
            ligands[j]      = nbp[i].atoms[j];
            tetra[j].x      = mp->atom_array[nbp[i].atoms[j]].x;
            tetra[j].y      = mp->atom_array[nbp[i].atoms[j]].y;
            tetra[j].z      = 0.0;
         }

         parity = 1;
         for (j=1; j<4; j++)	/* sort bonds in order */
            for (k=j-1; k>=0; k--)
               if (ligands[k] > ligands[k+1])
               {
                  parity *= (-1);
                  h = ligands[k]; ligands[k] = ligands[k+1]; ligands[k+1] = h;
               }
               else
                  break;

         for (irun = CHAIN_ONLY; irun <= USE_ALL; irun++)
         {
            maxvol = 0;
            d[0] = d[1] = d[2] = d[3] = 0;
            for (i0=(-1); i0<=1; i0++)  /* i0, i1, i2, i3 == -1, 0, or 1 */
            {
               if (mp->bond_array[tetra[0].number].bond_type != SINGLE &&
                   i0 != 0)
                  continue;
               if (mp->bond_array[tetra[0].number].stereo_symbol != NONE  &&
                   i0 != 0)
                  continue;  /* don't reassign stereobonds */

               bp = mp->bond_array + tetra[0].number;
               if (bp->atoms[0] == i+1) h = bp->atoms[1];
               else                     h = bp->atoms[0];
               if (irun <= CHAIN_ONLY  &&  is_ring_atom[h-1]  &&  i0 != 0)
                  continue;		/* ignore ring ligands for this run */

               if (irun < ALLOW_RING_BONDS         &&
                   ring_size[tetra[0].number] != 0 &&
                   i0 != 0)
                  continue;		/* ignore ring bonds for this run */

               if (irun < USE_ALL                               &&
                   (mp->atom_array[h-1].stereo_parity != NONE ||
                    mp->atom_array[h-1].color != NONE)          &&
                   i0 != 0)
                  continue;		/* don't use confusing stereobonds */

               tetra[0].z = i0*0.1/nbp[nbp[i].atoms[0]].n_ligands;
               for (i1=(-1); i1<=1; i1++)
               {
                  if (mp->bond_array[tetra[1].number].bond_type != SINGLE &&
                      i1 != 0)
                     continue;
                  if (mp->bond_array[tetra[1].number].stereo_symbol != NONE  &&
                      i1 != 0) continue;
                  bp = mp->bond_array + tetra[1].number;
                  if (bp->atoms[0] == i+1) h = bp->atoms[1];
                  else                     h = bp->atoms[0];
                  if (irun <= CHAIN_ONLY && is_ring_atom[h-1] && i1 != 0)
                     continue;
                  if (irun < ALLOW_RING_BONDS          &&
                      ring_size[tetra[1].number] != 0  &&
                      i1 != 0)
                     continue;

                  if (irun < USE_ALL                               &&
                      (mp->atom_array[h-1].stereo_parity != NONE ||
                       mp->atom_array[h-1].color != NONE)          &&
                      i1 != 0)
                     continue;
                  tetra[1].z = i1*0.1/nbp[nbp[i].atoms[1]].n_ligands;
                  for (i2=(-1); i2<=1; i2++)
                  {
                     if (mp->bond_array[tetra[2].number].bond_type != SINGLE &&
                         i2 != 0)
                        continue;
                     if (mp->bond_array[tetra[2].number].stereo_symbol
                            != NONE  &&
                         i2 != 0)
                        continue;
                     bp = mp->bond_array + tetra[2].number;
                     if (bp->atoms[0] == i+1) h = bp->atoms[1];
                     else                     h = bp->atoms[0];
                     if (irun <= CHAIN_ONLY && is_ring_atom[h-1] && i2 != 0)
                        continue;
                     if (irun < ALLOW_RING_BONDS          &&
                         ring_size[tetra[2].number] != 0  &&
                         i2 != 0) continue;

                     if (irun < USE_ALL                               &&
                         (mp->atom_array[h-1].stereo_parity != NONE ||
                          mp->atom_array[h-1].color != NONE)          &&
                         i2 != 0)
                        continue;
                     tetra[2].z = i2*0.1/nbp[nbp[i].atoms[2]].n_ligands;
                     for (i3=(-1); i3<=1; i3++)
                     {
                        if (mp->bond_array[tetra[3].number].bond_type
                              != SINGLE  &&
                            i3 != 0)
                           continue;
                        if (mp->bond_array[tetra[3].number].stereo_symbol
                               != NONE  &&
                            i3 != 0) continue;
// fprintf(stderr, "irun=%d: trying assignment %d|%d|%d|%d\n", irun, i0,i1,i2,i3);
                        bp = mp->bond_array + tetra[3].number;
                        if (bp->atoms[0] == i+1) h = bp->atoms[1];
                        else                     h = bp->atoms[0];
                        if (irun <= CHAIN_ONLY && is_ring_atom[h-1] && i3 != 0)
                           continue;
                        if (irun < ALLOW_RING_BONDS          &&
                            ring_size[tetra[3].number] != 0  &&
                            i3 != 0)
                           continue;

                        if (irun < USE_ALL                               &&
                            (mp->atom_array[h-1].stereo_parity != NONE ||
                             mp->atom_array[h-1].color != NONE)          &&
                            i3 != 0)
                           continue;
                        tetra[3].z = i3*0.1/nbp[nbp[i].atoms[3]].n_ligands;
                        if (i0*i0 + i1*i1 + i2*i2 + i3*i3 < 1)
                           continue;	/* at least one stereo bond */
                        if (i0*i0 + i1*i1 + i2*i2 + i3*i3 > 2)
                           continue;	/* more than 2 stereo bonds */
                        if (i0*i0 + i1*i1 + i2*i2 + i3*i3 == 1)
                        {	/* Check for umbrella atoms */
                           if (i0*i0 == 1) k = 0;	/* which one is it */
                           if (i1*i1 == 1) k = 1;
                           if (i2*i2 == 1) k = 2;
                           if (i3*i3 == 1) k = 3;
                           nclose = 0;
                           for (j = 0; j < 4; j++)
                              if (j != k)
                              {
                                 if ((tetra[k].x-x)*(tetra[j].x-x) +
                                     (tetra[k].y-y)*(tetra[j].y-y)  >
                                     ((tetra[k].x-x)*(tetra[k].x-x) +
                                      (tetra[k].y-y)*(tetra[k].y-y))*0.1)
                                    nclose++; /* points towards stereo bond */
                              }
// fprintf(stderr, "checkking umbrella for atom %d, nclose = %d\n", i+1, nclose);
                           if (nclose<1  ||  (nclose<2 && irun<USE_ALL))
                              continue;
                        }
// fprintf(stderr, "passed umbrella check\n");
                        /* Check if two opposing stereobonds */
                        if (i0*i0 + i1*i1 + i2*i2 + i3*i3 == 2)
                        {
                           if (i0 != 0  &&  i1 != 0) {j=0; k=1;}
                           if (i0 != 0  &&  i2 != 0) {j=0; k=2;}
                           if (i0 != 0  &&  i3 != 0) {j=0; k=3;}
                           if (i1 != 0  &&  i2 != 0) {j=1; k=2;}
                           if (i1 != 0  &&  i3 != 0) {j=1; k=3;}
                           if (i2 != 0  &&  i3 != 0) {j=2; k=3;}
                           /*
                           if (tetra[j].z*tetra[k].z < 0  &&
                               irun >= ALLOW_RING_BONDS)
                               */
                           {
                              if ((tetra[j].x-x)*(tetra[k].x-x) +
                                  (tetra[j].y-y)*(tetra[k].y-y)  <
                                   -0.3*((tetra[j].x-x)*(tetra[j].x-x) +
                                         (tetra[j].y-y)*(tetra[j].y-y)))
                                 continue;
                              if ((tetra[j].x-x)*(tetra[k].x-x) +
                                  (tetra[j].y-y)*(tetra[k].y-y)  <
                                   -0.3*((tetra[k].x-x)*(tetra[k].x-x) +
                                         (tetra[k].y-y)*(tetra[k].y-y)))
                                 continue;
                           }
                        }
// fprintf(stderr, "passed opposition check\n");
                        if (irun <= ALLOW_RING_BONDS  &&
                            i0*i0 + i1*i1 + i2*i2 + i3*i3 > 1) continue;
// fprintf(stderr, "passed one-stereo bond check\n");
                        if (maxvol < Volume(tetra)/(i0*i0+i1*i1+i2*i2+i3*i3))
                        {
                           maxvol = Volume(tetra)/(i0*i0+i1*i1+i2*i2+i3*i3);
                           d[0] = i0; d[1] = i1; d[2] = i2; d[3] = i3;
                        }
                     }
                  }
               }
            }

            /* found biggest volume for current ligand coloring */
            if (maxvol > 0.00001)
            {
               for (j=0; j<4; j++)
               {
                  if (d[j] == 0) continue;
                  if (ap->stereo_parity == 1) d[j] *= (-1);
                  if (parity < 0) d[j] *= (-1);
                  bp = mp->bond_array+tetra[j].number;
                  if (bp->atoms[0] != i+1) /* make bond point towards ligand */
                  {
                     h = bp->atoms[0];
                     bp->atoms[0] = bp->atoms[1];
                     bp->atoms[1] = h;
                  }
                  if (d[j] > 0) bp->stereo_symbol = UP;
                  else	        bp->stereo_symbol = DOWN;
               }

               ap->stereo_parity = NONE; /* This atom has been assigned */
               ap->color = NONE+1;
               break;
            }
         }
         if (maxvol < 0.00001)
            ShowMessageI("Could not find stereo assignment for atom %d",
                        "LayoutAtomStereo", i+1);
      }
      else
      {
         sprintf(buffer,
                 "Unexpected stereochemistry with %d ligands at atom %d",
                 nbp[i].n_ligands, i+1);
         ShowMessage(buffer, "LayoutAtomStereo");
      }
   }

   ResetColors(mp);
}

/*
 * Clears the DONT_FLIP_BOND flags from bond_type fields.
 */
void ClearFlipFlags(struct reaccs_molecule_t *mp)
{
   struct reaccs_bond_t *bp;
   int i;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      bp->bond_type &= ~DONT_FLIP_BOND;
}

void ClearDBStereoInSmallRings(struct reaccs_molecule_t *mp,
                               int ring_size[])
{
   struct reaccs_bond_t *bp;
   int i;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if ((ALL_BOND_TYPES&bp->bond_type) != DOUBLE) continue;

      if (ring_size[i] < 8  &&	/* clear stereo symbols in small rings */
          ring_size[i] > 0  &&
          (ALL_BOND_TYPES&bp->bond_type) == DOUBLE)
      {
         bp->stereo_symbol = NONE;
         continue;
      }
   }
}

void LayoutBondStereo(struct reaccs_molecule_t *mp,
                      neighbourhood_t *nbp,
                      int ring_size[])
/*
 * Uses the stereoparity information of the bonds of *mp
 * to flip double bonds if needed.
 */
{
   struct reaccs_bond_t *bp;
   int i, j;
   struct reaccs_atom_t *ap;
   int l1, ai1, l2, ai2;
   int *atom_colors;

   point p1, p2, p, r12, pp;
   point pdb1, pdb2, pcenter, pdb1_new, pdb2_new;
   double lastq;
   double m00, m01, m10, m11;   // rotation matrix
   double q;

   atom_colors = TypeAlloc(mp->n_atoms, int);	/* save and clear colors */
   for (i=0; i<mp->n_atoms; i++)
   {
      atom_colors[i] = mp->atom_array[i].color;
      mp->atom_array[i].color = NONE;
   }

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      // bond_type needs to be masked because the field is also used to secure coordinate-based placement
      if ((ALL_BOND_TYPES&bp->bond_type) != DOUBLE) continue;

      if (ring_size[i] < 8  &&	/* clear stereo symbols in small rings */
          ring_size[i] > 0  &&
          (ALL_BOND_TYPES&bp->bond_type) == DOUBLE)
      {
         bp->stereo_symbol = NONE;
         continue;
      }
      ai1 = bp->atoms[0]-1; ai2 = bp->atoms[1]-1;

      /* already defined => next bond */
      if (bp->stereo_symbol == NONE  ||
          bp->stereo_symbol == CIS_TRANS_SWAPPED  ||
          bp->stereo_symbol == CIS_TRANS_EITHER) continue;

      if (nbp[ai1].n_ligands < 2  ||
          nbp[ai1].n_ligands > 3  ||
          nbp[ai2].n_ligands < 2  ||
          nbp[ai2].n_ligands > 3)
      {
         ShowMessageI("clearing illegal stereo bond description %d",
                      "LayoutBondStereo", bp->stereo_symbol);
         bp->stereo_symbol = NONE;
         continue;
      }

      /* find reference ligands, i.e. the ones with smallest index */
      l1 = mp->n_atoms;
      for (j=0; j<nbp[ai1].n_ligands; j++)
         if (nbp[ai1].atoms[j] != ai2  &&
             nbp[ai1].atoms[j] < l1) l1 = nbp[ai1].atoms[j];
      l2 = mp->n_atoms;
      for (j=0; j<nbp[ai2].n_ligands; j++)
         if (nbp[ai2].atoms[j] != ai1  &&
             nbp[ai2].atoms[j] < l2) l2 = nbp[ai2].atoms[j];
// fprintf(stderr, "l1 = %d, l2 = %d\n", l1, l2);

      p1[0] = mp->atom_array[l1].x - mp->atom_array[ai1].x;
      p1[1] = mp->atom_array[l1].y - mp->atom_array[ai1].y;
      p2[0] = mp->atom_array[l2].x - mp->atom_array[ai2].x;
      p2[1] = mp->atom_array[l2].y - mp->atom_array[ai2].y;
      r12[0]= mp->atom_array[ai2].x - mp->atom_array[ai1].x;
      r12[1]= mp->atom_array[ai2].y - mp->atom_array[ai1].y;
      q = (p1[0]*r12[1]-p1[1]*r12[0]) * (p2[0]*r12[1]-p2[1]*r12[0]);
      if (((bp->stereo_symbol & TRANS_MASK)  &&  q < 0)  ||
          ((bp->stereo_symbol & CIS_MASK)    &&  q > 0))
      {				/* double bond is OK */
         bp->stereo_symbol = NONE;
         continue;
      }
//      else if (bp->stereo_symbol & (TRANS_MASK|CIS_MASK))
//      {
//fprintf(stderr, "DB %d-%d has stereo_symbol %d and q = %g\n", bp->atoms[0], bp->atoms[1], bp->stereo_symbol, q);
//      }
// fprintf(stderr, "DB %d-%d needs flipping\n", bp->atoms[0], bp->atoms[1]);

      if (ring_size[i] != 0  &&  ring_size[i] < 10)
      {
         bp->stereo_symbol = CIS_TRANS_SWAPPED;
         continue;
      }
      else if (ring_size[i] != 0)
      {
// fprintf(stderr, "1: CIS_TRANS_SWAPPED bond %d-%d: q = %g, ring_size = %d\n", bp->atoms[0], bp->atoms[1], q, ring_size[i]);
          // TODO: try with tilted bonds in case of large rings
          // save original coordinates
          pdb1[0] =  mp->atom_array[ai1].x; pdb1[1] =  mp->atom_array[ai1].y;
          pdb2[0] =  mp->atom_array[ai2].x; pdb2[1] =  mp->atom_array[ai2].y;
          pcenter[0] = 0.5*(pdb1[0]+pdb2[0]); pcenter[1] = 0.5*(pdb1[1]+pdb2[1]);
          // set up rotation
          m00 = cos(PI/6.0);  m01 = sin(PI/6.0);
          m10 = -sin(PI/6.0); m11 = cos(PI/6.0);
          // rotate bond by 30 degrees and check
          mp->atom_array[ai1].x = pcenter[0]+m00*(pdb1[0]-pcenter[0])+m01*(pdb1[1]-pcenter[1]);
          mp->atom_array[ai1].y = pcenter[1]+m10*(pdb1[0]-pcenter[0])+m11*(pdb1[1]-pcenter[1]);
          mp->atom_array[ai2].x = pcenter[0]+m00*(pdb2[0]-pcenter[0])+m01*(pdb2[1]-pcenter[1]);
          mp->atom_array[ai2].y = pcenter[1]+m10*(pdb2[0]-pcenter[0])+m11*(pdb2[1]-pcenter[1]);
          p1[0] = mp->atom_array[l1].x - mp->atom_array[ai1].x;
          p1[1] = mp->atom_array[l1].y - mp->atom_array[ai1].y;
          p2[0] = mp->atom_array[l2].x - mp->atom_array[ai2].x;
          p2[1] = mp->atom_array[l2].y - mp->atom_array[ai2].y;
          r12[0]= mp->atom_array[ai2].x - mp->atom_array[ai1].x;
          r12[1]= mp->atom_array[ai2].y - mp->atom_array[ai1].y;
          lastq = 0;
          q = (p1[0]*r12[1]-p1[1]*r12[0]) * (p2[0]*r12[1]-p2[1]*r12[0]);
          if (((bp->stereo_symbol & TRANS_MASK)  &&  q < 0)  ||
              ((bp->stereo_symbol & CIS_MASK)    &&  q > 0))
          {				/* double bond is OK */
// fprintf(stderr, "2: CIS_TRANS_SWAPPED bond %d-%d: q = %g, ring_size = %d\n", bp->atoms[0], bp->atoms[1], q, ring_size[i]);
             bp->stereo_symbol = NONE;
             pdb1_new[0] =  mp->atom_array[ai1].x; pdb1_new[1] =  mp->atom_array[ai1].y;
             pdb2_new[0] =  mp->atom_array[ai2].x; pdb2_new[1] =  mp->atom_array[ai2].y;
             lastq = fabs(q);
             // continue;
          }
          // set up rotation
          m00 = cos(PI/6.0); m01 = -sin(PI/6.0);
          m10 = sin(PI/6.0); m11 =  cos(PI/6.0);
          // rotate bond by 30 degrees and check
          mp->atom_array[ai1].x = pcenter[0]+m00*(pdb1[0]-pcenter[0])+m01*(pdb1[1]-pcenter[1]);
          mp->atom_array[ai1].y = pcenter[1]+m10*(pdb1[0]-pcenter[0])+m11*(pdb1[1]-pcenter[1]);
          mp->atom_array[ai2].x = pcenter[0]+m00*(pdb2[0]-pcenter[0])+m01*(pdb2[1]-pcenter[1]);
          mp->atom_array[ai2].y = pcenter[1]+m10*(pdb2[0]-pcenter[0])+m11*(pdb2[1]-pcenter[1]);
          p1[0] = mp->atom_array[l1].x - mp->atom_array[ai1].x;
          p1[1] = mp->atom_array[l1].y - mp->atom_array[ai1].y;
          p2[0] = mp->atom_array[l2].x - mp->atom_array[ai2].x;
          p2[1] = mp->atom_array[l2].y - mp->atom_array[ai2].y;
          r12[0]= mp->atom_array[ai2].x - mp->atom_array[ai1].x;
          r12[1]= mp->atom_array[ai2].y - mp->atom_array[ai1].y;
          q = (p1[0]*r12[1]-p1[1]*r12[0]) * (p2[0]*r12[1]-p2[1]*r12[0]);
          if (((bp->stereo_symbol & TRANS_MASK)  &&  q < 0)  ||
              ((bp->stereo_symbol & CIS_MASK)    &&  q > 0))
          {				/* double bond is OK */
// fprintf(stderr, "3: CIS_TRANS_SWAPPED bond %d-%d: q = %g, ring_size = %d\n", bp->atoms[0], bp->atoms[1], q, ring_size[i]);
             if (bp->stereo_symbol != NONE  &&  lastq < fabs(q))
             {
                 pdb1_new[0] =  mp->atom_array[ai1].x; pdb1_new[1] =  mp->atom_array[ai1].y;
                 pdb2_new[0] =  mp->atom_array[ai2].x; pdb2_new[1] =  mp->atom_array[ai2].y;
             }
             bp->stereo_symbol = NONE;
             // continue;
          }
          if (lastq > 0) // did find a fix
          {
             mp->atom_array[ai1].x = pdb1_new[0]; mp->atom_array[ai1].y = pdb1_new[1];
             mp->atom_array[ai2].x = pdb2_new[0]; mp->atom_array[ai2].y = pdb2_new[1];
             continue;
          }
          // restore original coordinates and mark a swapped
          mp->atom_array[ai1].x = pdb1[0]; mp->atom_array[ai1].y = pdb1[1];
          mp->atom_array[ai2].x = pdb2[0]; mp->atom_array[ai2].y = pdb2[1];
          bp->stereo_symbol = CIS_TRANS_SWAPPED;
          continue;
      }

      /* at this point, we have a chain bond */
      /* which has to be flipped to fix it */

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
      ResetColors(mp);
      bp->stereo_symbol = NONE;
   }

   for (i=0; i<mp->n_atoms; i++)	/* restore colors */
   {
      mp->atom_array[i].color = atom_colors[i];
   }
   MyFree((char *)atom_colors);
}

void MakeLandscape(struct reaccs_molecule_t *mp)
/*
 * Rotates the molecule by 90 degrees, if width is smaller than height.
 */
{
   float h, xmin, xmax, ymin, ymax;
   int i;
// BR int prelayed_out;

   if (mp->n_atoms == 0) return;

// BR prelayed_out = FALSE;
   xmin = xmax = mp->atom_array[0].x; ymin = ymax = mp->atom_array[0].y;
   for (i=1; i<mp->n_atoms; i++)
   {
// BR if (mp->atom_array[i].color & KEEP_POSITION) prelayed_out = TRUE;
      if (xmin > mp->atom_array[i].x) xmin = mp->atom_array[i].x;
      if (xmax < mp->atom_array[i].x) xmax = mp->atom_array[i].x;
      if (ymin > mp->atom_array[i].y) ymin = mp->atom_array[i].y;
      if (ymax < mp->atom_array[i].y) ymax = mp->atom_array[i].y;
   }

   if (xmax - xmin  <  ymax - ymin)	/* portrait oreintation? */
      for (i=0; i<mp->n_atoms; i++)	/* then rotate		 */
      {
         h = mp->atom_array[i].x;
         mp->atom_array[i].x = mp->atom_array[i].y;
         mp->atom_array[i].y = -h;
      }
}

#define MAXCOLLISIONS 4    /* maximum number of collisions to be resolved */
#define SYMEPS        0.02 /* collision tolerance if atom has symbol */
#define NOSYMEPS      0.004 /* collision tolerance if atom has no symbol */
// #define NOSYMEPS      0.010 / * C1CC2CCC1CC2 distorted */
// #define NOSYMEPS      0.009 / * C1CC2CCC1CC2 not-distorted */

double BondLengthSQR(struct reaccs_molecule_t *mp,
                     struct reaccs_bond_t     *bp)
/*
 * Computes the square of the bond length of bond *bp within *mp.
 */
{
   double dx, dy;

   dx = mp->atom_array[bp->atoms[0]-1].x - mp->atom_array[bp->atoms[1]-1].x;
   dy = mp->atom_array[bp->atoms[0]-1].y - mp->atom_array[bp->atoms[1]-1].y;
   return (dx*dx + dy*dy);
}

double AtomBondClash(struct reaccs_molecule_t *mp,
                     struct reaccs_atom_t     *ap,
                     struct reaccs_bond_t     *bp,
                     double                    eps,
                     double		      *dist)
/*
 * Returns TRUE if the atom *ap and the bond *bp of molecule *mp clash within
 * eps of the length of *bp and FALSE otherwise.
 */
{
   double bx, by;
   double ax, ay;
   double y2, a2, b2, ab, amb2;

   (*dist) = STDBOND;
   bx = mp->atom_array[bp->atoms[1]-1].x - mp->atom_array[bp->atoms[0]-1].x;
   by = mp->atom_array[bp->atoms[1]-1].y - mp->atom_array[bp->atoms[0]-1].y;
   ax = ap->x - mp->atom_array[bp->atoms[0]-1].x;
   ay = ap->y - mp->atom_array[bp->atoms[0]-1].y;

   a2 = ax*ax + ay*ay; b2 = bx*bx + by*by; ab = ax*bx + ay*by;
   amb2 = SQR(ax-bx) + SQR(ay-by);

   if (b2 <= 0) return (FALSE);		/* Give up on zero length bonds */

   y2 = a2 - ab*ab/b2;

   if (y2 >= eps*b2)	/* no candidate for a clash */
      return (FALSE);
   else
      if (0 <= ab  &&  ab <= b2)	/* clash with bond line */
      {
         (*dist) = y2/eps;
         return (TRUE);
      }
      else
         if (a2 <= eps*b2)		/* clash with first atom */
         {
            (*dist) = a2/eps;
            return (TRUE);
         }
         else if (amb2 <= eps*b2)	/* clash with second atom */
         {
            (*dist) = amb2/eps;
            return (TRUE);
         }

   return (FALSE);
}

int CountCollisions(struct reaccs_molecule_t *mp,
                    struct reaccs_atom_t **apcol,
                    struct reaccs_bond_t **bpcol)
/*
 * Counts the number of collisions in molecule mp and returns this count.
 * *apcol and *bpcol are set to the atom and bond pointers of the
 * last collision detected.
 */
{
   int i, j, ncol;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   double eps;
   double dist, disth;

   ncol = 0;
   dist = 100*STDBOND;
   *apcol = (struct reaccs_atom_t *)NULL;
   *bpcol = (struct reaccs_bond_t *)NULL;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if ((ap->charge                   != 0) ||
          (ap->radical                  != 0) ||
          (ap->mass_difference          != 0) ||
          (strcmp(ap->atom_symbol, "C") != 0))
         eps = SYMEPS;
      else
         eps = NOSYMEPS;

      for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
      {
         if (bp->bond_type & RUBBER_BOND) continue;

         if (bp->atoms[0] == i+1  ||  bp->atoms[1] == i+1) continue;

         /* skip if in different fragments */
         if (ap->color != mp->atom_array[bp->atoms[0]-1].color) continue;
         if (ap->color != mp->atom_array[bp->atoms[1]-1].color) continue;

         if (AtomBondClash(mp, ap, bp, eps, &disth))
         {
            /* only count real clashes */
            if (disth < 100*STDBOND) ncol++;
            if (disth < dist)
            {
               dist = disth;
               *apcol = ap; *bpcol = bp;
            }
         }
      }
   }

   return (ncol);
}

int BondLinksAtoms(struct reaccs_molecule_t *mp,
                   neighbourhood_t *nbp,
                   struct reaccs_bond_t *bp,
                   struct reaccs_atom_t *ap1, int col1,
                   struct reaccs_atom_t *ap2, int col2)
/*
 * Returns TRUE if the bond *bp links the two atoms *ap1 and *ap2
 * via some path. It returns false if this is not the case. It also
 * returns FALSE if *ap1 and *ap2 are also linked via some alternate
 * path not containing *bp.
 *
 * As a side-effect, this routine colors the pieces linked to atoms
 * *ap1 and *ap2 with colors col1 and col2, respectively. This can be
 * used to stretch *bp to remove collisions.
 */
{
   int result;
   neighbourhood_t *nbph;
   int i, j, changed;

   for (i=0; i<mp->n_atoms; i++)
   {
      mp->atom_array[i].color = 0;
   }

   result = FALSE;

   /* flood fill fragment rooted at *ap1 */
   ap1->color = col1;
   do
   {
      changed = FALSE;
      for (i=0, nbph=nbp; i<mp->n_atoms; i++, nbph++)
      {
         if (mp->atom_array[i].color != ap1->color) continue;

         for (j=0; j<nbph->n_ligands; j++)
            if (mp->atom_array[nbph->atoms[j]].color == 0  &&
                (bp-mp->bond_array) != nbph->bonds[j]      &&
                !(mp->bond_array[nbph->bonds[j]].bond_type & RUBBER_BOND))
            {
               mp->atom_array[nbph->atoms[j]].color = ap1->color;
               changed = TRUE;
            }
      }
   } while (changed);
   if (ap2->color == col1) /* *ap1 and *ap2 linked by some other path */
      goto clean_up;
   if (mp->atom_array[bp->atoms[0]-1].color == /* bond within first fragment */
       mp->atom_array[bp->atoms[1]-1].color)   /* or second fragment or not  */
      goto clean_up;			       /* at all connected to atoms  */

   /* flood fill fragment rooted at *ap2 */
   ap2->color = col2;
   do
   {
      changed = FALSE;
      for (i=0, nbph=nbp; i<mp->n_atoms; i++, nbph++)
      {
         if (mp->atom_array[i].color != ap2->color) continue;

         for (j=0; j<nbph->n_ligands; j++)
            if (mp->atom_array[nbph->atoms[j]].color == 0  &&
                (bp-mp->bond_array) != nbph->bonds[j]      &&
                !(mp->bond_array[nbph->bonds[j]].bond_type & RUBBER_BOND))
            {
               mp->atom_array[nbph->atoms[j]].color = ap2->color;
               changed = TRUE;
            }
      }
   } while (changed);
   if (mp->atom_array[bp->atoms[0]-1].color == 0  ||
       mp->atom_array[bp->atoms[1]-1].color == 0) /* bond points off path */
      goto clean_up;

   /* bond really links the atoms */
   result = TRUE;

clean_up:
   return (result);
}

void StretchBond(struct reaccs_molecule_t *mp,
                 struct reaccs_bond_t     *bp,
                 int                       frag_col1,
                 int                       frag_col2,
                 double			   stretch_factor)
/*
 * Stretches the fragments of *mp colored with colors frag_col1 and
 * frag_col2, resp., along the direction defined by *bp. stretch_factor
 * is the extent to which bond *bp was stretched if it connected the two
 * fragments, which is the usual application of this procedure.
 */
{
   point r12;
   struct reaccs_atom_t *ap;
   int i, h;

   if (mp->atom_array[bp->atoms[0]-1].color == frag_col2  &&
       mp->atom_array[bp->atoms[1]-1].color == frag_col1)
   {
      h = frag_col1; frag_col1 = frag_col2; frag_col2 = h;
   }
   r12[0] = mp->atom_array[bp->atoms[1]-1].x -
            mp->atom_array[bp->atoms[0]-1].x;
   r12[1] = mp->atom_array[bp->atoms[1]-1].y -
            mp->atom_array[bp->atoms[0]-1].y;

   /* Give up if bond too short */
   if (SQR(r12[0]) < 0.0000001  &&  SQR(r12[1]) < 0.0000001) return;

   stretch_factor = sqrt(stretch_factor);

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == frag_col1)
      {
         ap->x -= (stretch_factor-1)*r12[0];
         ap->y -= (stretch_factor-1)*r12[1];
      }
      else if (ap->color == frag_col2)
      {
         ap->x += (stretch_factor-1)*r12[0];
         ap->y += (stretch_factor-1)*r12[1];
      }
}

#define STRETCH	1.4

void ImproveCollisions(struct reaccs_molecule_t *mp,
                       neighbourhood_t nbp[],
                       int is_ring_atom[],
                       int ring_count[])
/*
 * Finds collisions between atoms and/or bonds of *mp. Tries to fix
 * them by stretching one of the connecting bonds. It prefers chain
 * bonds * connecting rings over chain bonds sprouted from rings over
 * pure chain bonds. The total number of connections to the atoms
 * of the candidate bonds is also used to order the candidates by
 * decreasing degree. This should be a reasonable order to result
 * in as little distortion as possible.
 *
 * The parameter nbp[] describes the neigbourhood of atoms to speed
 * up access. is_ring_atom[] identifies ring atoms and ring_count[]
 * counts the number of cycles shared by the given bond.
 */
{
   int icol, i, j, k, ncol, ncolh;
   struct reaccs_atom_t *apcol, *apcolh;
   struct reaccs_bond_t *bp, *bpcol, *bpcolh;
   int *atom_colors;

   atom_colors = (int *)NULL;

   for (icol=0; icol <MAXCOLLISIONS; icol++)	/* don't try too hard */
   {
      /* find and count collisions */
      ncol = CountCollisions(mp, &apcol, &bpcol);
      if (ncol == 0) break;	/* no collision => we are done */
// fprintf(stderr, "ncol = %d\n", ncol);

      /* save colors */
      if (!atom_colors)
      {
         atom_colors = TypeAlloc(mp->n_atoms, int);
         for (i=0; i<mp->n_atoms; i++)
            atom_colors[i] = mp->atom_array[i].color;
      }

#define RINGLINKS	1
#define RINGSPROUTS	2
#define CONGESTEDBONDS	3
#define SIMPLECHAIN	4

      for (k=RINGLINKS; k<=SIMPLECHAIN; k++)	/* go through bond classes */
      {					        /* as stretch candidates   */
         for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
         {
            if (ring_count[j] > 0) continue;

            switch (k)	/* select bond class to be used */
            {
               case RINGLINKS:
                  if (!is_ring_atom[bp->atoms[0]-1] ||
                      !is_ring_atom[bp->atoms[1]-1])
                     continue;
                  break;
               case RINGSPROUTS:
                  if (is_ring_atom[bp->atoms[0]-1] ==
                      is_ring_atom[bp->atoms[1]-1])
                     continue;
                  break;
               case CONGESTEDBONDS:
                  if (is_ring_atom[bp->atoms[0]-1] ||
                      is_ring_atom[bp->atoms[1]-1])
                     continue;
                  if (nbp[bp->atoms[0]-1].n_ligands < 3  ||
                      nbp[bp->atoms[1]-1].n_ligands < 3)
                     continue;
                  break;
               case SIMPLECHAIN:
                  if (is_ring_atom[bp->atoms[0]-1] ||
                      is_ring_atom[bp->atoms[1]-1])
                     continue;
                  if (nbp[bp->atoms[0]-1].n_ligands >= 3  &&
                      nbp[bp->atoms[1]-1].n_ligands >= 3)
                     continue;
                  break;
               default:
                  continue;
            }

            if (BondLinksAtoms(mp, nbp, bp,
                               apcol, 1,
                               &mp->atom_array[bpcol->atoms[0]-1], 2))
            {
               StretchBond(mp, bp, 1, 2, STRETCH);
               ncolh = CountCollisions(mp, &apcolh, &bpcolh);
               if (ncol <= ncolh)
               {	/* does not improve => undo change */
                  StretchBond(mp, bp, 1, 2, 1.0/STRETCH);
               }
               else if (ncolh > 0)	/* improvement found => next run */
               {
// fprintf(stderr, "stretching bond %d-%d and continue\n", bp->atoms[0], bp->atoms[1]);
                  break;
               }
               else		/* no further collisions => done */
               {
// fprintf(stderr, "stretching bond %d-%d\n", bp->atoms[0], bp->atoms[1]);
                  goto done_with_it;
               }
            }
         }
         if (j < mp->n_bonds)/* there was a change => don't try other classes*/
            break;
      }

      if (j == mp->n_bonds)   	/* Here we have a hard case. */
      { 	/* Try just moving the atom randomly +/- 1 tenth of a bond */
         apcol->x += STDBOND*(rand()%100-50)/500.0;
         apcol->y += STDBOND*(rand()%100-50)/500.0;
      }
   }

done_with_it:

   /* restore colors */
   if (atom_colors)
   {
      for (i=0; i<mp->n_atoms; i++)
         mp->atom_array[i].color = atom_colors[i];
      MyFree((char *)atom_colors);
   }
}

/*
 * Scales the coordinates of *mp such that the average length of already
 * layed-out bonds (those with the same atom color at both ends)
 * equals STDBOND. This is a NOP if there is no such bond.
 */
void ScaleByFixedFragments(struct reaccs_molecule_t *mp)
{
   struct reaccs_bond_t *bp;
   struct reaccs_atom_t *ap1, *ap2;
   int i, n;
   double scale;

   scale = 0.0; n = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      ap1 = mp->atom_array+bp->atoms[0]-1;
      ap2 = mp->atom_array+bp->atoms[1]-1;
      if (ap1->color == ap2->color)
      {
         scale += SQR(ap1->x - ap2->x) + SQR(ap1->y - ap2->y);
         n++;
      }
   }
   if (scale == 0.0) return;

   scale /= n; scale = STDBOND/sqrt(scale);
   for (i=0, ap1=mp->atom_array; i<mp->n_atoms; i++, ap1++)
   {
      ap1->x *= scale; ap1->y *= scale;
   }
}

void FixLinearBridgeHeads(struct reaccs_molecule_t *mp, int *ring_count, neighbourhood_t *nbp)
/*
 * Scans the molecule *mp for bridgehead atoms that have linearly aligned pairs of bonds and are needed to represent
 * stereochemistry. The bridgehead is deflected along the remaining bond.
 */
{
   int i, j, jj;
   int nring, nh, ih, jh, ibest, nbulk;
   int i1, i2, i3;
   int found;
   struct reaccs_atom_t *ap, *aph;
   double p[2], p1[2], p2[2], p3[2];
   double sprod;
   double lenh, lensum;
   double l1, l2;
   double cosa;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (nbp[i].n_ligands != 4) continue;
      nring = 0;
      nh = 0;
      nbulk = 100;
      for (j=0; j<nbp[i].n_ligands; j++)
      {
         if (strcmp("H",mp->atom_array[nbp[i].atoms[j]].atom_symbol)==0  &&
             ring_count[nbp[i].bonds[j]] == 0  &&    // unlikely case of H atom in ring
             nbp[nbp[i].atoms[j]].n_ligands == 1)    // unlikely case of hydrogen with more than one attachment
         {
            nh++;
            ih = nbp[i].atoms[j];
            jh = j;
         }
         else
         {
            if (ring_count[nbp[i].bonds[j]] > 0)
               nring++;
         }
      }
      p[0] = ap->x; p[1] = ap->y;
      if (nring != 3  ||  nh != 1) continue;
      // look for linear arrangements
      found = FALSE;
      for (j=0; j<nbp[i].n_ligands; j++)
      {
         if (j == jh) continue;
         p1[0] = mp->atom_array[nbp[i].atoms[j]].x;
         p1[1] = mp->atom_array[nbp[i].atoms[j]].y;
         l1 = sqrt((p1[0]-p[0])*(p1[0]-p[0])+(p1[1]-p[1])*(p1[1]-p[1]));
         for (jj=j+1; jj<nbp[i].n_ligands; jj++)
         {
            if (jj == jh) continue;
            p2[0] = mp->atom_array[nbp[i].atoms[jj]].x;
            p2[1] = mp->atom_array[nbp[i].atoms[jj]].y;
            l2 = sqrt((p2[0]-p[0])*(p2[0]-p[0])+(p2[1]-p[1])*(p2[1]-p[1]));
            cosa = (p1[0]-p[0])*(p2[0]-p[0])+(p1[1]-p[1])*(p2[1]-p[1]);
            cosa /= l1*l2;
            if (cosa < -0.99)
            {
               found = TRUE;
               i1 = nbp[i].atoms[j];
               i2 = nbp[i].atoms[jj];
               if (j !=0  &&  jj != 0  &&  jh != 0) i3 = nbp[i].atoms[0];
               if (j !=1  &&  jj != 1  &&  jh != 1) i3 = nbp[i].atoms[1];
               if (j !=2  &&  jj != 2  &&  jh != 2) i3 = nbp[i].atoms[2];
               if (j !=3  &&  jj != 3  &&  jh != 3) i3 = nbp[i].atoms[3];
               break;
            }
         }
      }
      if (!found) continue;
      p3[0] = mp->atom_array[i3].x; p3[1] = mp->atom_array[i3].y;
      ap->x += 0.25*(p3[0]-p[0]); ap->y += 0.25*(p3[1]-p[1]);
   }
}

void FixInwardHydrogens(struct reaccs_molecule_t *mp, int *ring_count, neighbourhood_t *nbp)
/*
 * Scans the molecule *mp for hydrogen atoms pointing into a ring and move them to an outward pointing place. This
 * Helps rendering sterochemistry and consequently canoncalization.
 */
{
   int i, j;
   int nring, nh, ih, ichain, ibest, nbulk;
   struct reaccs_atom_t *ap, *aph;
   double p[2], ph[2], pchain[2], psum[2];
   double sprod;
   double lenh, lensum;

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (nbp[i].n_ligands != 4) continue;
      nring = 0;
      nh = 0;
      nbulk = 100;
      for (j=0; j<nbp[i].n_ligands; j++)
      {
         if (strcmp("H",mp->atom_array[nbp[i].atoms[j]].atom_symbol)==0  &&
             ring_count[nbp[i].bonds[j]] == 0  &&    // unlikely case of H atom in ring
             nbp[nbp[i].atoms[j]].n_ligands == 1)    // unlikely case of hydrogen with more than one attachment
         {
            nh++;
            ih = nbp[i].atoms[j];
            aph = mp->atom_array+ih;
         }
         else
         {
            if (ring_count[nbp[i].bonds[j]] > 0)
            {
               nring++;
               if (nbulk > nbp[nbp[i].atoms[j]].n_ligands)
               {
                  nbulk = nbp[nbp[i].atoms[j]].n_ligands;
                  ibest = nbp[i].atoms[j];
               }
            }
            else
               ichain = nbp[i].atoms[j];
         }
      }
      if (nring != 2  ||  nh != 1) continue;
      p[0] = ap->x; p[1] = ap->y;
      ph[0] = aph->x; ph[1] = aph->y;
      pchain[0] = mp->atom_array[ichain].x; pchain[1] = mp->atom_array[ichain].y;
      sprod = (ph[0]-p[0])*(pchain[0]-p[0]) + (ph[1]-p[1])*(pchain[1]-p[1]);
      if (sprod >= 0) continue;
      psum[0] = mp->atom_array[ibest].x-p[0] + pchain[0]-p[0];
      psum[1] = mp->atom_array[ibest].y-p[1] + pchain[1]-p[1];
      lenh = (ph[0]-p[0])*(ph[0]-p[0]) + (ph[1]-p[1])*(ph[1]-p[1]); lenh = sqrt(lenh);
      lensum = sqrt(psum[0]*psum[0]+psum[1]*psum[1]);
      if (lensum < 0.25*lenh) continue;  // too wide an angle
      aph->x = p[0] + psum[0]*lenh/lensum; aph->y = p[1] + psum[1]*lenh/lensum;
// fprintf(stderr, "FixInwardHydrogens: considering atom %d linking to hydrogen %d and chain atom %d, best neighbour is %d\n", i+1, ih+1, ichain+1, ibest+1);
   }
}

#define CENTER_DISTANCE       1
#define TOP_BOTTOM_SEPARATION 2

double TopBottomSeparation(struct reaccs_molecule_t *mp, int fragment_color, int top_color, int bottom_color)
/*
 * Computes the total separation of the top and bottom attachments measured perpendicular to the diameter of the fragment.
 */
{
   int dia1, dia2;
   double d, dmax, total;
   double p[2], p1[2], edia[2], eper[2];
   int nedges, nnodes, seed;
   int *numbers;
   atom_pair *edges;
   point *coords;
   int i, j;
   struct reaccs_atom_t *ap1, *ap2, *aph;
   struct reaccs_bond_t *bp;

   dmax = 0.0; dia1 = -1; dia2 = -1;
   // first, search for the most distant pair of atoms
   for (i=0, ap1=mp->atom_array; i<mp->n_atoms; i++, ap1++)
   {
      if (ap1->color != fragment_color) continue;
      // ignore hydrogen atom
      if (ap1->atom_symbol[0] == 'H'  &&  ap1->atom_symbol[1] == '\0') continue;
      for (j=i+1, ap2=mp->atom_array+(i+1); j<mp->n_atoms; j++, ap2++)
      {
         if (ap2->color != fragment_color) continue;
         if (ap2->atom_symbol[0] == 'H'  &&  ap2->atom_symbol[1] == '\0') continue;
         d = SQR(ap1->x-ap2->x)+SQR(ap1->y-ap2->y);
         if (d > dmax)
         {
            dia1 = i; dia2 = j;
            dmax = d;
         }
      }
   }
   if (dia1 < 0  || dia2 < 0) return 0.0;
// fprintf(stderr, "diameter extends from %s(%d) to %s(%d)\n", mp->atom_array[dia1].atom_symbol, dia1+1, mp->atom_array[dia2].atom_symbol, dia2+1);
   p[0] = mp->atom_array[dia1].x; p[1] = mp->atom_array[dia1].y;
   // compute the unit vector edia along the diameter
   edia[0] = mp->atom_array[dia2].x-mp->atom_array[dia1].x;
   edia[1] = mp->atom_array[dia2].y-mp->atom_array[dia1].y;
   d = sqrt(SQR(edia[0]) + SQR(edia[1]));
   edia[0] /= d; edia[1] /= d;
   // eper is the perpendicular unit vector
   eper[0] = edia[1]; eper[1] = -edia[0];
   // now, compute the total difference of top and bottom attachment measured along eper.
   total = 0.0;
   numbers = TypeAlloc(mp->n_atoms, int);
   edges = TypeAlloc(mp->n_bonds, atom_pair);
   coords = TypeAlloc(mp->n_atoms, point);
   GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, fragment_color);
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      ap1 = mp->atom_array+(bp->atoms[0]-1); ap2 = mp->atom_array+(bp->atoms[1]-1);
      if (ap1->color != fragment_color  &&  ap2->color != fragment_color) continue;
      if (ap1->color == fragment_color  &&  ap2->color == fragment_color) continue;
      if (ap2->color == fragment_color) // swap to make sure ap1 in in the fragment
      {
         aph = ap1; ap1 = ap2; ap2 = aph;
      }
      if (ap2->color != top_color  &&  ap2->color != bottom_color) continue;
      seed = numbers[ap1-mp->atom_array];
      NextSubstituentPoint(p1, coords, nnodes, edges, nedges, seed, FALSE, numbers, mp->n_atoms);
      if (ap2->color == top_color)
      {
         total += eper[0]*(p1[0]-p[0]) + eper[1]*(p1[1]-p[1]);
      }
      else if (ap2->color == bottom_color)
      {
         total -= eper[0]*(p1[0]-p[0]) + eper[1]*(p1[1]-p[1]);
      }
   }
   MyFree((char *)numbers);
   MyFree((char *)edges);
   MyFree((char *)coords);
   return fabs(total) * dmax;
}

double AttachmentCenterDistance(struct reaccs_molecule_t *mp, int color, int attachment_color)
/*
 * Computes the distance of the atoms bearing the attachment (colored with attachment_color) from the center of the fragment.
 */
{
   double xcenter, ycenter;
   double att_dist;
   double result;
   int i, j, n;
   int ai1, ai2;
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp, *bph;

   att_dist = 0.0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (mp->atom_array[bp->atoms[0]-1].color == color  &&
          mp->atom_array[bp->atoms[1]-1].color != color)
         ai1 = bp->atoms[0]-1;
      else if (mp->atom_array[bp->atoms[1]-1].color == color  &&
          mp->atom_array[bp->atoms[0]-1].color != color)
         ai1 = bp->atoms[1]-1;
      else
         continue;
      // found an attachment at atom ai1;
      for (j=i+1, bph=bp+1; j<mp->n_bonds; j++, bph++)
      {
         if (mp->atom_array[bph->atoms[0]-1].color == color  &&
             mp->atom_array[bph->atoms[1]-1].color != color)
            ai2 = bph->atoms[0]-1;
         else if (mp->atom_array[bph->atoms[1]-1].color == color  &&
             mp->atom_array[bph->atoms[0]-1].color != color)
            ai2 = bph->atoms[1]-1;
         else
            continue;
         // found another attachment at atom ai2;
         // prefer distances of 4 bond length for attachments
         att_dist += SQR(MIN(0.0,4*1.54-sqrt(SQR(mp->atom_array[ai1].x-mp->atom_array[ai2].x) +  SQR(mp->atom_array[ai1].y-mp->atom_array[ai2].y))));
      }
   }

   result = 0.0;
   n = 0; xcenter = ycenter = 0.0;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == color)
      {
         xcenter += ap->x;
         ycenter += ap->y;
         n++;
      }
   if (n == 0) return 0.0;
   xcenter /= n; ycenter /= n;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (mp->atom_array[bp->atoms[0]-1].color == color  &&
          mp->atom_array[bp->atoms[1]-1].color == attachment_color)
      {
         result += SQR(mp->atom_array[bp->atoms[0]-1].x-xcenter);
         result += SQR(mp->atom_array[bp->atoms[0]-1].y-ycenter);
      }
      else if (mp->atom_array[bp->atoms[1]-1].color == color  &&
               mp->atom_array[bp->atoms[0]-1].color == attachment_color)
      {
         result += SQR(mp->atom_array[bp->atoms[1]-1].x-xcenter);
         result += SQR(mp->atom_array[bp->atoms[1]-1].y-ycenter);
      }
      return result-10.0*att_dist;
}

double ImproveFragmentByBondFlip(struct reaccs_molecule_t *mp,
                                 neighbourhood_t *nbp,
                                 struct reaccs_bond_t     *bp,
                                 int improvement_type,
                                 int col1, int col2,
                                 int forced_undo)
/*
 * Tries to improve the layout of the mp fragment with the given color by flipping one side of bond bp.
 * The function returnes TRUE if it sucessfully improved the coordinates.
 * improvement_type indicates which kind of score is to be improved. col1 and col2 are possible parameters for the score.
 */
{
   point p1, p2, p, pp, r12;
   double strain1, strain2, q, color_strain;
   int i;
   int color;
   struct reaccs_atom_t *ap;

   if (bp->bond_type & DONT_FLIP_BOND)
      return (0.0);
   if (mp->atom_array[bp->atoms[0]-1].color != mp->atom_array[bp->atoms[1]-1].color)    // bond crosses fragments
      return (0.0);

   color = mp->atom_array[bp->atoms[0]-1].color;

   color_strain =  -ColorStrain(mp, color, color);
   if (improvement_type == CENTER_DISTANCE)
   {
      strain1 = AttachmentCenterDistance(mp, color, col1) + color_strain;
   }
   else if (improvement_type == TOP_BOTTOM_SEPARATION)
      strain1 = TopBottomSeparation(mp, color, col1, col2) + color_strain;
   else
      return (0.0);

   FloodInvertColor(mp, nbp, bp->atoms[0]-1, color, bp);
   if (mp->atom_array[bp->atoms[1]-1].color == -color)          // loop back => uncolor and return
   {
if (0) fprintf(stderr, "found loop back at bond %s(%d)-%s(%d)\n",
      mp->atom_array[bp->atoms[0]-1].atom_symbol, bp->atoms[0],
      mp->atom_array[bp->atoms[1]-1].atom_symbol, bp->atoms[1]);
      ChangeAtomColors(mp, -color, color);
      return (0.0);
   }
   // now the part of the fragment that contains bp->atoms[0] is -color while the other part has color.
   p1[0] = mp->atom_array[bp->atoms[0]-1].x;
   p1[1] = mp->atom_array[bp->atoms[0]-1].y;
   p2[0] = mp->atom_array[bp->atoms[1]-1].x;
   p2[1] = mp->atom_array[bp->atoms[1]-1].y;
   r12[0] = p2[0] - p1[0]; r12[1] = p2[1]-p1[1];

   /* Now flip bond */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == -color)
      {
         p[0] = ap->x; p[1] = ap->y;
         q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
             (r12[0]*r12[0]       + r12[1]*r12[1]);
         pp[0] = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
         pp[1] = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
         ap->x = pp[0]; ap->y = pp[1];
      }
   ChangeAtomColors(mp, -color, color);

   color_strain =  -ColorStrain(mp, color, color);
   if (improvement_type == CENTER_DISTANCE)
   {
      strain2 = AttachmentCenterDistance(mp, color, col1) + color_strain;
   }
   else if (improvement_type == TOP_BOTTOM_SEPARATION)
   {
      strain2 = TopBottomSeparation(mp, color, col1, col2) + color_strain;
   }

   /* transform back if improvement isn't sufficient */
   if (strain1 > strain2  ||  forced_undo)
   {
      // undoing the flip
      FloodInvertColor(mp, nbp, bp->atoms[0]-1, color, bp);
      for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
         if (ap->color == -color)
         {
            p[0] = ap->x; p[1] = ap->y;
            q = ((p[0]-p1[0])*r12[0] + (p[1]-p1[1])*r12[1])/
                (r12[0]*r12[0]       + r12[1]*r12[1]);
            pp[0] = p1[0] - (p[0]-p1[0]) + 2*r12[0]*q;
            pp[1] = p1[1] - (p[1]-p1[1]) + 2*r12[1]*q;
            ap->x = pp[0]; ap->y = pp[1];
         }
      ChangeAtomColors(mp, -color, color);
   }
   
// if (strain1+0.001 < strain2)
// fprintf(stderr, "%s: Inverting color %d starting at atom %d of bond %s(%d)-%s(%d) improves strain from %g to %g\n",
//       (forced_undo ? "UNDONE" : "FINAL"),
//       color, bp->atoms[0],
//       mp->atom_array[bp->atoms[0]-1].atom_symbol, bp->atoms[0],
//       mp->atom_array[bp->atoms[1]-1].atom_symbol, bp->atoms[1],
//       strain1, strain2);
   return (strain2-strain1);
}

int ImproveFragmentByAllBondFlips(struct reaccs_molecule_t *mp,
                                  neighbourhood_t *nbp,
                                  int color,
                                  int improvement_type,
                                  int col1, int col2)
{
   int i;
   struct reaccs_bond_t *bp, *best_bond;
   int result;
   double improvement, best_improvement;

   result = FALSE;
   
// fprintf(stderr,"Bond Flipping tried with color %d\n", color);
   do
   {
      best_improvement = 0.0;
      best_bond = (struct reaccs_bond_t *)NULL;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->bond_type & DONT_FLIP_BOND)
            continue;
         else if (mp->atom_array[bp->atoms[0]-1].color != color)
            continue;
         else if (mp->atom_array[bp->atoms[1]-1].color != color)
            continue;
         else
         {
            improvement = ImproveFragmentByBondFlip(mp, nbp, bp, improvement_type, col1, col2, TRUE);
            if (improvement > best_improvement+0.001)
            {
               best_improvement = improvement;
               best_bond = bp;
            }
         }
      if (best_bond)
      {
         improvement = ImproveFragmentByBondFlip(mp, nbp, best_bond, improvement_type, col1, col2, FALSE);
         result = TRUE;
      }
   } while (best_bond);
   return result;
}

int FloodInvertColor(struct reaccs_molecule_t *mp,
                     neighbourhood_t          *nbp,
                     int start_index,
                     int fragmentColor,
                     struct reaccs_bond_t *flipBond)
/*
 * Recursively changes the atom colors of the fragment (indicated by fragmentColor) of *mp fragmentColor to -fragmentColor
 * starting at the current start_index but not crossing flipBond using a flood-filling algorithm.
 * The function is used to split the fragment into two parts, one of which can be flipped along the axis defined by flipBond.
 * At the start of the recursion, the caller needs to use an atom of flipBond as the start_index.
 */
{
   int i, nflip, ai1, ai2;
   neighbourhood_t *nbph;

   nbph=nbp+start_index;
   if (mp->atom_array[start_index].color != fragmentColor)      // not in fragment or already flipped
      return 0;
   nflip = 1;
   mp->atom_array[start_index].color *= -1;
   for (i=0; i<nbph->n_ligands; i++)
   {
      if (flipBond-mp->bond_array == nbph->bonds[i])            // don't progress along flipBond
         continue;
      if (mp->atom_array[nbph->atoms[i]].color != fragmentColor) // not in fragment or already flipped, practically redundant
         continue;
      nflip += FloodInvertColor(mp, nbp, nbph->atoms[i], fragmentColor, flipBond);
   }
   return nflip;
}

int ChangeAtomColors(struct reaccs_molecule_t *mp, int fromColor, int toColor)
/*
 * Changes the atom colors fromColor to toColor.
 *
 * The function returned the number of changes made.
 */
{
   int i, result;
   struct reaccs_atom_t *ap;

   result = 0;
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (ap->color == fromColor)
      {
         ap->color = toColor;
         result++;
      }
   return result;
}

static int IsShortcutLikeAtom(struct reaccs_molecule_t *mp, int i, neighbourhood_t *nbp)
/**
 * Labels all R atoms with an atext and all single atom symbol regular atoms attached to just a shortcut atom as shortcut-like.
 */
{
    struct reaccs_atom_t *ap;
    char *cp;
    int j;

    ap = mp->atom_array+i;
    if (ap->atom_symbol[0] == 'R'  &&  ap->atom_symbol[1] == '\0'  &&  0 < strlen(ap->atext))
    {
        // don't render repeat atoms as strings
        for (cp=ap->atext; (*cp) != '\0'; cp++)
            if ('0' <= cp[0]  &&  cp[0] <= '9') return (FALSE);
        // check if this short cut has only real atom neighbours
        for (j=0; j<nbp[i].n_ligands; j++)
        {
            ap = mp->atom_array+nbp[i].atoms[j];
            if (ap->atom_symbol[0] == 'R'  &&  ap->atom_symbol[1] == '\0'  &&  0 < strlen(ap->atext)) return TRUE;
        }
        return (FALSE);  // no other shortcut neighbour
    }
    // check if terminal single atom connected to shortcut
    if (strlen(ap->atom_symbol) > 1  ||  nbp[i].n_ligands > 1) return (FALSE);
    ap = mp->atom_array+nbp[i].atoms[0];
    if (ap->atom_symbol[0] == 'R'  &&  ap->atom_symbol[1] == '\0'  &&  0 < strlen(ap->atext))
    {
        // don't consider repeat atoms as shortcut neighbour
        for (cp=ap->atext; (*cp) != '\0'; cp++)
            if ('0' <= cp[0]  &&  cp[0] <= '9') return (FALSE);
        // check if this short cut neighbour has at least one other shortcut-like neighbours
        for (j=0; j<nbp[nbp[i].atoms[0]].n_ligands; j++)
        {
            ap = mp->atom_array+nbp[nbp[i].atoms[0]].atoms[j];
            if (ap->atom_symbol[0] == 'R'  &&  ap->atom_symbol[1] == '\0'  &&  0 < strlen(ap->atext)) return TRUE;
        }
    }
    return (FALSE);
}

#define SHORTCUT_NONE   0
#define SHORTCUT_START  1
#define SHORTCUT_END    2
#define SHORTCUT_INNER  4
#define SHORTCUT_BRANCH 8
#define SHORTCUT_RING   16
#define SHORTCUT_SINGLE 32

int PerceiveShortcutClasses(struct reaccs_molecule_t *result, int *is_ring_atom,
                            int *shortcut_class, neighbourhood_t *nbp, int max_line_size)
/*
 * Returns TRUE if a shortcut chain was found and FALSE otherwise.
 */
{
    int i, j, jj, next_start;
    struct reaccs_atom_t *ap;
    neighbourhood_t *nbph, *nbphh;
    int nshortcuts;
    int *used_atoms;
    int chain_found;
    int line_size;
// int candidate_found;

    used_atoms = TypeAlloc(result->n_atoms, int);
    // first find special shortcuts
    for (i=0, ap=result->atom_array, nbph=nbp; i<result->n_atoms; i++, ap++, nbph++)
    {
        shortcut_class[i] = SHORTCUT_NONE;
        if (!IsShortcutLikeAtom(result, i, nbp))
        {
            used_atoms[i] = TRUE;
            continue;  // start atom needs to be a shortcut
        }
        nshortcuts = 0;
        for (j=0; j<nbph->n_ligands; j++)
            if (IsShortcutLikeAtom(result, nbph->atoms[j], nbp)) nshortcuts++;
        if (nshortcuts > 2)
        {
            shortcut_class[i] |= SHORTCUT_BRANCH;
            fprintf(stderr, "detected shortcut branching at %s(%d)\n", ap->atext, i+1);
        }
        else if (nshortcuts == 0)
        {
            shortcut_class[i] |= SHORTCUT_SINGLE;
            fprintf(stderr, "detected shortcut single at %s(%d)\n", ap->atext, i+1);
        }
        if (is_ring_atom[i] &&  (shortcut_class[i]&SHORTCUT_BRANCH) != 0)   // don't count ring closures by normal atoms
            shortcut_class[i] |= SHORTCUT_RING;
        if (shortcut_class[i] != SHORTCUT_NONE) used_atoms[i] = TRUE;
    }
// candidate_found = FALSE;
    chain_found = FALSE;
    // find real chains (and rings) of shortcuts
    for (i=0, ap=result->atom_array, nbph=nbp; i<result->n_atoms; i++, ap++, nbph++)
    {
        if (used_atoms[i]) continue;
        nshortcuts = 0;
        for (j=0; j<nbph->n_ligands; j++)
            if (result->bond_array[nbph->bonds[j]].bond_type != NONE  &&
                (!used_atoms[nbph->atoms[j]])  &&  IsShortcutLikeAtom(result, nbph->atoms[j], nbp)) nshortcuts++;
        if (nshortcuts == 2)
        {
// candidate_found = TRUE;
            chain_found = TRUE;
            shortcut_class[i] |= SHORTCUT_INNER;
// fprintf(stderr, "4a:%s[%d]=%d ", result->atom_array[i].atext, i, shortcut_class[i]); fprintf(stderr, "\n");
        }
        else if (nshortcuts == 1)
        {
// candidate_found = TRUE;
            for (j=0; j<nbph->n_ligands; j++)
                if (result->bond_array[nbph->bonds[j]].bond_type != NONE  &&
                    (!used_atoms[nbph->atoms[j]])  &&
                    IsShortcutLikeAtom(result, nbph->atoms[j], nbp))
                {
                    if (nbph->atoms[j] > i) // connected shortcut has higher index => start of sequence
                    {
                        shortcut_class[i] |= SHORTCUT_START;
// fprintf(stderr, "4b:%s[%d]=%d ", result->atom_array[i].atext, i, shortcut_class[i]); fprintf(stderr, "\n");
                    }
                    else                    // end of sequence
                    {
                        shortcut_class[i] |= SHORTCUT_END;
// fprintf(stderr, "4c:%s[%d]=%d ", result->atom_array[i].atext, i, shortcut_class[i]); fprintf(stderr, "\n");
                    }
                    break;
                }
        }
        else
           used_atoms[i] = TRUE; // just to be sure
    }
// if (candidate_found) {for (i=0; i<result->n_atoms; i++) fprintf(stderr, "2:%s[%d]=%d ", result->atom_array[i].atext, i, shortcut_class[i]); fprintf(stderr, "\n");}
    // split long sequences
    for (;;)
    {
        next_start = -1;
        // search for next best start of sequence
        for (i=0, ap=result->atom_array; i<result->n_atoms; i++, ap++)
        {
            if (used_atoms[i]) continue;
            if (shortcut_class[i] & SHORTCUT_START)
            {
                next_start = i;
                break;
            }
            else if (next_start < 0  &&  0 != (shortcut_class[i]&(SHORTCUT_INNER|SHORTCUT_START|SHORTCUT_END)))
                next_start = i;
        }
        if (next_start < 0) break;  // done
        shortcut_class[next_start] = SHORTCUT_START;
// fprintf(stderr, "(0)starting shortcut sequence at %s(%d)\n", result->atom_array[next_start].atext, next_start);
        // trace the long sequence and break it if needed
        nshortcuts = 1;
        line_size = strlen(result->atom_array[next_start].atext);
        for (;;)
        {
            used_atoms[next_start] = TRUE;
            nbph = nbp+next_start;
            // find regular shortcut in sequence
            for (j=0; j<nbph->n_ligands; j++)
            {
                if (used_atoms[nbph->atoms[j]]) continue;
                if (0 != (shortcut_class[nbph->atoms[j]] & (SHORTCUT_INNER|SHORTCUT_START|SHORTCUT_END)))
                {
                    break;
                }
            }
            if (j == nbph->n_ligands)   // the current atom it the last in the sequence
            {
                if (shortcut_class[next_start] == SHORTCUT_START)
                   shortcut_class[next_start] = SHORTCUT_START; // | SHORTCUT_END;
                else
                   shortcut_class[next_start] = SHORTCUT_END;
                break;
            }
            // check if sequence needs to be split here since it's too long
            line_size += strlen(result->atom_array[next_start].atext);
            nshortcuts++;
            if (max_line_size > 1  &&  line_size > max_line_size)
            {
                // check if there is no branch point here or at the next shortcut and avoid splitting in this case
                if (nbph->n_ligands == 2  &&  nbp[nbph->atoms[j]].n_ligands == 2  &&  0 == (shortcut_class[nbph->atoms[j]] & SHORTCUT_END))
                {
                   line_size = strlen(result->atom_array[next_start].atext);
                   // check if next atom would be singleton and don't end sequence if it is
                   nbphh = nbp+nbph->atoms[j];
                   for (jj=0; jj<nbphh->n_ligands; jj++)
                      if ((!used_atoms[nbphh->atoms[jj]]) &&  (shortcut_class[nbphh->atoms[jj]] & (SHORTCUT_INNER|SHORTCUT_START|SHORTCUT_END)))
                         break;
                   if (jj == nbphh->n_ligands)
                   {
// fprintf(stderr, "(1)not braking shortcut sequence at %s(%d)\n", result->atom_array[next_start].atext, next_start);
                   }
                   else
                   {
// fprintf(stderr, "(1)braking shortcut sequence at %s(%d)\n", result->atom_array[next_start].atext, next_start);
                      nshortcuts = 1;
                      shortcut_class[next_start] = SHORTCUT_END;
                      shortcut_class[nbph->atoms[j]] = SHORTCUT_START;
                   }
                }
            }
            next_start = nbph->atoms[j];
// fprintf(stderr, "(2)continuing shortcut sequence at %s(%d) class=%d\n", result->atom_array[next_start].atext, next_start, shortcut_class[next_start]);
        }
    }
// if (candidate_found) {for (i=0; i<result->n_atoms; i++) fprintf(stderr, "3:%s[%d]=%d/%d ", result->atom_array[i].atext, i, shortcut_class[i], used_atoms[i]); fprintf(stderr, "\n");}
    MyFree((char *)used_atoms);
    return chain_found;
}

void LayoutShortcutChains(struct reaccs_molecule_t *result, int *is_ring_atom, int *shortcut_class, neighbourhood_t *nbp)
/*
 * Compute layout of chains of (acyclic) shortcut 'atoms' by arranging them linearly.
 *
 * Returns TRUE if a shortcut chain was found and FALSE otherwise.
 */
{
    struct reaccs_atom_t *ap;
    struct reaccs_bond_t *bp;
    neighbourhood_t *nbph;
    int *used_atoms;
    int i, j, col;
    float x, y;

    used_atoms = TypeAlloc(result->n_atoms, int);

    // layout perceived shortcut chains
    for (;;)
    {
        for (i=0, ap=result->atom_array, nbph=nbp; i<result->n_atoms; i++, ap++, nbph++)
        {
            if (0 == (shortcut_class[i] & (SHORTCUT_START|SHORTCUT_INNER|SHORTCUT_END))) used_atoms[i] = TRUE;
            if (used_atoms[i]) continue;
            if (shortcut_class[i] & SHORTCUT_START)
            {
// fprintf(stderr, "starting shortcut sequence at %s(%d)\n", ap->atext, i+1);
// fprintf(stderr, "%d:(%g,%g) ", i+1, ap->x, ap->y);
                col = ap->color;
                x = ap->x; y = ap->y;
                if (ap->atext[0] == 0)
                    x += 0.5*1.54*(3+strlen(ap->atom_symbol))/4.0;
                else
                    x += 0.5*1.54*(3+strlen(ap->atext))/4.0;
                break;
            }
        }
        if (i == result->n_atoms) break;
        used_atoms[i] = TRUE;
        for (;;)
        {
            nbph = nbp+i;
            // find shortcut in sequence
            for (j=0; j<nbph->n_ligands; j++)
                if (!used_atoms[nbph->atoms[j]]  &&  0 != (shortcut_class[nbph->atoms[j]] & (SHORTCUT_INNER|SHORTCUT_START|SHORTCUT_END)))
                {
                    i = nbph->atoms[j];
                    break;
                }
            if (j == nbph->n_ligands) break;
            used_atoms[i] = TRUE;
            ap = result->atom_array+i;
            ap->color = col;
            if (ap->atext[0] == 0)
                x += 0.5*1.54*(3+strlen(ap->atom_symbol))/4.0;
            else
                x += 0.5*1.54*(3+strlen(ap->atext))/4.0;
            ap->x = x; ap->y = y;
            if (ap->atext[0] == 0)
                x += 0.5*1.54*(3+strlen(ap->atom_symbol))/4.0;
            else
                x += 0.5*1.54*(3+strlen(ap->atext))/4.0;
// fprintf(stderr, "%d:(%g,%g) ", i+1, ap->x, ap->y);
            if (shortcut_class[i] & SHORTCUT_END)
            {
// fprintf(stderr, "\n");
                break;        // catch end-to-end sequences
            }
        }
    }

    // label bonds from shortcut sequences to other fragments as RUBBER_BONDs
    for (i=0, bp=result->bond_array; i<result->n_bonds; i++, bp++)
    {
        if ((shortcut_class[bp->atoms[0]-1] & (SHORTCUT_INNER|SHORTCUT_START|SHORTCUT_END))  &&
            shortcut_class[bp->atoms[1]-1] == SHORTCUT_NONE)
            bp->bond_type |= RUBBER_BOND;
        if ((shortcut_class[bp->atoms[1]-1] & (SHORTCUT_INNER|SHORTCUT_START|SHORTCUT_END))  &&
            shortcut_class[bp->atoms[0]-1] == SHORTCUT_NONE)
            bp->bond_type |= RUBBER_BOND;
        if ((shortcut_class[bp->atoms[1]-1] & (SHORTCUT_START|SHORTCUT_END))  &&
            (shortcut_class[bp->atoms[0]-1] & (SHORTCUT_START|SHORTCUT_END)))
            bp->bond_type |= RUBBER_BOND;
    }

    MyFree((char *)used_atoms);
    return;
}

#define FT_SHORTCUTS    1
#define FT_ATOMS        2

#define MAX_FRAGMENT_NEIGHBOURS 3

struct fragment_t
{
    int color;                                   // color of the fragment to be placed
    int first_index;                             // index of first symbol
    double xcenter, ycenter, radius;             // properties of the fragment's bounding circle
    double xoffset, yoffset;                     // desired position of the fragment center
    int fragment_type;                           // type of fragment to be placed (FT_SHORTCUTS or FT_ATOMS)
    int nsymbols;                                // number of atoms or shortcuts in this fragment
    int nneighbours;                             // number of neighbouring fragments
    int neighbour_cols[MAX_FRAGMENT_NEIGHBOURS]; // colors of up to MAX_FRAGMENT_NEIGHBOURS connected fragments
    int nsneighbours;                            // number of shortcut neighbours
    int used;                                    // label used for sequence threading
    struct fragment_t *next;                     // pointer to next fragment in linked list
};

struct fragment_t **CollectFragmentInfo(struct reaccs_molecule_t *mp, int *shortcut_class, neighbourhood_t *nbp, int *nfragments)
   /*
    * Collects geometric and graph information on the fragments to be layed out and returns an array of shortcut fragments with
    * non-shortcut fragments linked by the *next field.
    * Assumes each fragment to have its own distinct color.
    */
{
    int i, j, k, col1, col2;
    double d;
    struct reaccs_atom_t *ap;
    struct reaccs_bond_t *bp;
    struct fragment_t *fragments, *fp, *fph, *next_fragment, **farray;

    if (IsNULL(mp)  ||  mp->n_atoms == 0) return (struct fragment_t **)NULL;
    // first collect the barebones fragment discriptors
    fragments = (struct fragment_t *)NULL;
    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    {
        // search fragment with matching color
        for (fp = fragments; !IsNULL(fp); fp=fp->next) if (fp->color == ap->color) break;
        if (IsNULL(fp))
        {   // no matching fragment descriptor => make one and initialize
            fp = TypeAlloc(1, struct fragment_t);
            fp->color = ap->color;
            fp->first_index = i;
            fp->xcenter = ap->x; fp->ycenter = ap->y;
            fp->nsymbols = 1;
            fp->nsneighbours = 0;
            if (IsShortcutLikeAtom(mp, i, nbp)  &&  shortcut_class[i]  !=  SHORTCUT_SINGLE)
                fp->fragment_type = FT_SHORTCUTS;
            else
                fp->fragment_type = FT_ATOMS;
            fp->radius = FRAME_HEIGHT*1.0/3.0;
            fp->used = FALSE;
            fp->next = fragments;
            fragments = fp;
        }
        else
        {   // update matching fragment descriptor with atom data
            if (IsShortcutLikeAtom(mp, i, nbp)  &&  shortcut_class[i]  !=  SHORTCUT_SINGLE)
                fp->fragment_type = FT_SHORTCUTS;
            fp->xcenter += ap->x; fp->ycenter += ap->y;
            fp->nsymbols++;
        }
    }
    // finish size calculations
    (*nfragments) = 0;
    for (fp = fragments; !IsNULL(fp); fp=fp->next)
    {
        (*nfragments)++;
        fp->xcenter /= fp->nsymbols;
        fp->ycenter /= fp->nsymbols;
    }
    for (fp = fragments; !IsNULL(fp); fp=fp->next)
    {
        fp->radius = FRAME_HEIGHT*1.0/3.0;
        if (fp->fragment_type == FT_SHORTCUTS) continue;
        for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
        {
            if (fp->color != ap->color) continue;
            d = sqrt(SQR(fp->xcenter-ap->x) + SQR(fp->ycenter-ap->y));
            if (d > fp->radius) fp->radius = d;
        }
    }
// fprintf(stderr, "collect fragment connectivity\n");
    // now, collect fragment connectivity
    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    {
        col1 = mp->atom_array[bp->atoms[0]-1].color;
        col2 = mp->atom_array[bp->atoms[1]-1].color;
        if (col1 == col2) continue;
// fprintf(stderr, "processing bond (%d-%d) connecting colors (%d-%d)\n", bp->atoms[0], bp->atoms[1], col1, col2);
        for (fp=fragments; !IsNULL(fp); fp=fp->next)
        {
            if (fp->color == col1)     // first color matched
            {
                // add neighbour color if needed
                for (j=0; j<fp->nneighbours && j<MAX_FRAGMENT_NEIGHBOURS; j++)
                {
                    if (fp->neighbour_cols[j] == col2) break;     // already connected
                }
                if (j<fp->nneighbours) break;  // already connected
                // look if this neighbour is a shortcut sequence
                for (fph=fragments; !IsNULL(fph); fph=fph->next)
                    if (fph->color == col2)
                    {
                        if (fph->fragment_type == FT_SHORTCUTS) fp->nsneighbours++;
                        break;
                    }
                if (j==MAX_FRAGMENT_NEIGHBOURS) break;  // not enough slots
                fp->neighbour_cols[fp->nneighbours++] = col2;
            }
        }
        for (fp=fragments; !IsNULL(fp); fp=fp->next)
        {
            if (fp->color == col2)     // second color matched
            {
                // add neighbour color if needed
                for (j=0; j<fp->nneighbours && j<MAX_FRAGMENT_NEIGHBOURS; j++)
                {
                    if (fp->neighbour_cols[j] == col1) break;     // already connected
                }
                if (j<fp->nneighbours) break;  // already connected
                // look if this neighbour is a shortcut sequence
                for (fph=fragments; !IsNULL(fph); fph=fph->next)
                    if (fph->color == col1)
                    {
                        if (fph->fragment_type == FT_SHORTCUTS) fp->nsneighbours++;
                        break;
                    }
                if (j==MAX_FRAGMENT_NEIGHBOURS) break;  // not enough slots
                fp->neighbour_cols[fp->nneighbours++] = col1;
            }
        }
    }
    // sort neighbour_cols
    for (fp=fragments; !IsNULL(fp); fp=fp->next)
    {
        for (i=1; i<fp->nneighbours; i++)
            for (j=i-1; j>=0; j--)
            {
                if (fp->neighbour_cols[j] > fp->neighbour_cols[j+1])
                {
                    k = fp->neighbour_cols[j]; fp->neighbour_cols[j] = fp->neighbour_cols[j+1]; fp->neighbour_cols[j+1] = k;
                }
                else
                    break;
            }
    }
// for (next_fragment = fragments; !IsNULL(next_fragment); next_fragment=next_fragment->next)
// {
// fprintf(stderr, "before: col=%d, idx=%d,  center=(%g,%g), radius=%g, type=%s, nsymbols=%d, ncols=(%d,%d,%d), offset=(%g,%g)\n",
// next_fragment->color, next_fragment->first_index, next_fragment->xcenter, next_fragment->ycenter, next_fragment->radius,
// (next_fragment->fragment_type==FT_ATOMS?"FT_ATOMS":"FT_SHORTCUTS"), next_fragment->nsymbols, next_fragment->neighbour_cols[0], next_fragment->neighbour_cols[1], next_fragment->neighbour_cols[2], next_fragment->xoffset, next_fragment->yoffset);
// }
// fprintf(stderr, "thread fragment sequences\n");
    // thread fragment sequences
    farray = TypeAlloc((*nfragments), struct fragment_t *);
    for (i=0; i<(*nfragments); i++) farray[i] = (struct fragment_t *)NULL;
    // first, strip FT_SHORTCUT fragments from list and put them into first slots of farray
    fp = (struct fragment_t *)NULL;
    i=0;
    while (!IsNULL(fragments))
    {
        if (fragments->fragment_type == FT_SHORTCUTS)
        {
            next_fragment = fragments->next;
            farray[i++] = fragments;
            fragments->next = (struct fragment_t *)NULL;
            fragments = next_fragment;
        }
        else
        {
            next_fragment = fragments->next;
            fragments->next = fp;
            fp = fragments;
            fragments = next_fragment;
        }
    }
    fragments = fp;
    // now, all FT_SHORTCUT fragments are in farray[0..i-1].
    // sort them by first_index
    (*nfragments) = i;
// fprintf(stderr, "there are %d FT_SHORTCUT segments\n", (*nfragments));
    for (i=1; i<(*nfragments); i++)
        for (j=i-1; j>=0; j--)
            if (farray[j+1]->first_index < farray[j]->first_index)
            {
                next_fragment = farray[j]; farray[j] = farray[j+1]; farray[j+1] = next_fragment;
            }
            else
                break;
    // connecting the FT_ATOMS fragments to the first FT_SHORTCUT sequence if any or append to array if there is none
    while (!IsNULL(fragments))
    {
        for (i=0; i<(*nfragments); i++)
        {
            fp = farray[i];
            for (j=0; j<fragments->nneighbours; j++)
                if (fragments->neighbour_cols[j] == fp->color)
                    break;
            if (j < fragments->nneighbours) break;
        }
        if (i < (*nfragments))  // matching fragment found => link to this fragment
        {
// fprintf(stderr,"sorting FT_ATOMS fragment with color %d after FT_SHORTCUT fragment colored %d\n", fragments->color, farray[i]->color);
            fp = fragments->next; fragments->next = farray[i]->next; farray[i]->next = fragments; fragments = fp;
        }
        else                    // no matching neighbour found => append to farray
        {
            fp = fragments->next;
            fragments->next = (struct fragment_t *)NULL;
// fprintf(stderr,"appending FT_ATOMS fragment with color %d at position %d\n", fragments->color, i);
            farray[i++] = fragments;
            (*nfragments) = i;
            fragments = fp;
        }
    }
    // for debugging => sort shortcuts to front
    // for (i=1; i<(*nfragments); i++)
        // for (j=i-1; j>=0; j--)
            // if (farray[j]->fragment_type == FT_ATOMS  &&  farray[j+1]->fragment_type == FT_SHORTCUTS)
            // {
                // fph = farray[j]; farray[j] = farray[j+1]; farray[j+1] = fph;
            // }
            // else
                // break;
// fprintf(stderr, "returning %d fragments\n", (*nfragments));
// debug output
// for (i=0; i<(*nfragments); i++)
// {
//     for (j=0, fp=farray[i]; !IsNULL(fp); fp=fp->next, j+=3)
// fprintf(stderr, "%*scol=%d, idx=%d,  center=(%g,%g), radius=%g, type=%s, nsymbols=%d, ncols=(%d,%d,%d), offset=(%g,%g)\n",
//         j, "", fp->color, fp->first_index, fp->xcenter, fp->ycenter, fp->radius,
//         (fp->fragment_type==FT_ATOMS?"FT_ATOMS":"FT_SHORTCUTS"), fp->nsymbols, fp->neighbour_cols[0], fp->neighbour_cols[1], fp->neighbour_cols[2], fp->xoffset, fp->yoffset);
// }
    return farray;
}

double ColorScore(struct reaccs_molecule_t *mp, int moving_color, int fixed_color)
/*
 * Computes a layout score of the atoms colored moving_color wrt. those colored fixed_color.
 * The score considers standard bond lengths of connections and colisions.
 */
{
    int i;
    struct reaccs_bond_t *bp;
    struct reaccs_atom_t *ap1, *ap2, *aph;
    double result, d;

    result = 0.0;
    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
    {
        ap1 = mp->atom_array+(bp->atoms[0]-1);
        ap2 = mp->atom_array+(bp->atoms[1]-1);
        if (ap1->color == ap2->color) continue;
        if (ap1->color == moving_color &&  ap2->color == fixed_color)
        {
            aph = ap1; ap1 = ap2; ap2 = aph;
        }
        else if (ap2->color == moving_color &&  ap1->color == fixed_color)
        {
            // NOP
        }
        else
            continue;
        d = sqrt(SQR(ap1->x-ap2->x) + SQR(ap1->y-ap2->y));
        result += SQR(d-1.54);
    }
    return result;
}

#define MAX_ATTACHMENTS 10

struct attachment_t
{
    int iattachment;    // index of fragment atom where the attachment sprouts
    int iligand;        // index of off-fragment atom to which the attachment links
    int cligand;        // color of ligand
    int nneighbours;    // number of within fragment ligands
    int has_symbol;     // TRUE if this attachment would show with an atom symbol
    double x, y;        // coordinates of this atom
    double xtmp, ytmp;  // test coordinates of this atom 
    double xlig, ylig;  // position of ligand in original molecule
    double xdir, ydir;  // ref. dir. (points to neighbour if nneighbours == 1 and to pref. attachment direction if not)
};

int CollectAttachments(struct reaccs_molecule_t *mp, int color, struct attachment_t *attachments, neighbourhood_t *nbp)
/**
 * Collect the information about the off-fragment attachments of the color fragment in *mp.
 * It puts the information into attachments[] (up to MAX_ATTACHMENTS) and returns the number of attachments found.
 */
{
    struct reaccs_atom_t *ap, *ap1, *ap2;
    struct reaccs_bond_t *bp;
    neighbourhood_t *nbph;
    struct attachment_t *att;
    int na;
    int i, j, seed;
    atom_pair *edges;
    int nedges;
    point *coords;
    int nnodes;
    int *numbers;
    point p;

    na = 0;
    // step through off-fragment bonds
    for (j=0, bp=mp->bond_array; j<mp->n_bonds; j++, bp++)
    {
        ap1 = mp->atom_array + (bp->atoms[0]-1); ap2 = mp->atom_array + (bp->atoms[1]-1);
        if (ap1->color != color  &&  ap2->color != color) continue;
        if (ap1->color == ap2->color) continue;
        if (ap2->color == color) { ap = ap1; ap1 = ap2; ap2 = ap; }
        // ap1 now points to atom of fragment to be placed
        attachments[na].iattachment = (ap1-mp->atom_array);
        attachments[na].iligand = (ap2-mp->atom_array);
        attachments[na].cligand = ap2->color;
        attachments[na].has_symbol = 0 != strcmp("C", ap1->atom_symbol);
        attachments[na].x = ap1->x; attachments[na].y = ap1->y;
        attachments[na].xlig = ap2->x; attachments[na].ylig = ap2->y;
        // init other fields
        attachments[na].xdir = 0.0; attachments[na].ydir = 0.0;
        attachments[na].nneighbours = 0;
        na++;
        if (na >= MAX_ATTACHMENTS) break;
    }
    // fill additional fields
    for (i=0, att=attachments; i<na; i++, att++)
    {
        nbph = nbp+att->iattachment;
        for (j=0; j<nbph->n_ligands; j++)
            if (mp->atom_array[nbph->atoms[j]].color == color)
                att->nneighbours++;
            else
            {
                // just remember in case we have onle one internal neighbour
                att->xdir = mp->atom_array[nbph->atoms[j]].x;
                att->ydir = mp->atom_array[nbph->atoms[j]].y;
            }
        if (att->nneighbours > 0)     // look for direction of next best attachment point
        {
           numbers = TypeAlloc(mp->n_atoms, int);
           edges = TypeAlloc(mp->n_bonds, atom_pair);
           coords = TypeAlloc(mp->n_atoms, point);
           GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, color);
           seed = numbers[att->iattachment];
           NextSubstituentPoint(p, coords, nnodes, edges, nedges, seed, 0, numbers, mp->n_atoms);
           att->xdir = p[0]-att->x; att->ydir = p[1]-att->y;
           MyFree((char *)numbers);
           MyFree((char *)edges);
           MyFree((char *)coords);
        }
    }
    return na;
}

struct tfm_t
{
    double cos_alpha, sin_alpha;  // cosine and sine of the rotation angle alpha
    int do_flip;                  // TRUE if transformation is to flip the coordinate system along the Y axis
};

void TfmPoint(double *xnew, double *ynew,
              double x, double y,
              double xfromcenter, double yfromcenter,
              double xtocenter,   double ytocenter,
              struct tfm_t *tfmp)
/**
 * Transform the point x/y using *tfmp moving xfromcenter/yfromcenter to xtocenter/ytocenter, and put the result into (*newx)/(*newy).
 */
{
    x -= xfromcenter; y -= yfromcenter;
    if (tfmp->do_flip) x = -x;
    (*xnew) = x*tfmp->cos_alpha - y*tfmp->sin_alpha + xtocenter;
    (*ynew) = x*tfmp->sin_alpha + y*tfmp->cos_alpha + ytocenter;
}

double ScoreTransformation(struct tfm_t *tfmp, struct attachment_t *attachments, int na,
                           double xfrom, double yfrom, double xto, double yto,
                           struct reaccs_molecule_t *mp, int color, int bottom_color, int top_color, double target_height)
/**
 * Compute the score of the transformation *tfmp wrt. the attachments[0..na-1] within molecule *mp.
 */
{
    double score;
    double xnew, ynew;
    double xdir, ydir;
    double ymin, ymax;
    double yamin, yamax, xamin, xamax;
    int natoms;
    int i, j;
    struct reaccs_atom_t *ap, *aph;
    struct attachment_t *att;

    score = 0.0;
    yamin = 1.0e10; yamax = -yamin;
    xamin = 1.0e10; xamax = -xamin;
    natoms = 0;
    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    {
        if (ap->color != color) continue;
        natoms++;
        TfmPoint(&xnew, &ynew, ap->x, ap->y, xfrom, yfrom, xto, yto, tfmp);
        if (ynew < yamin) yamin = ynew;
        if (ynew > yamax) yamax = ynew;
        if (xnew < xamin) xamin = xnew;
        if (xnew > xamax) xamax = xnew;
        // collect collision score
        for (j=0, aph=mp->atom_array; j<mp->n_atoms; j++, aph++)
           if (aph->color != color  &&
               (aph->color == top_color  ||  aph->color == bottom_color))
           {
               score += 0.1/(0.1 + (ap->x - aph->x)* (ap->x - ap->x)
                                 + (ap->y - aph->y)* (ap->y - aph->y));
               score += 0.2/(0.1 + (ap->y - aph->y)* (ap->y - aph->y));
           }
    }
    ymin = 1.0e10; ymax = -ymin;
    for (i=0, att=attachments; i<na; i++, att++)
    {
        TfmPoint(&xnew, &ynew, mp->atom_array[att->iligand].x, mp->atom_array[att->iligand].y, xfrom, yfrom, xto, yto, tfmp);
        if (ynew < ymin) ymin = ynew;
        if (ynew > ymax) ymax = ynew;
    }
    for (i=0, att=attachments; i<na; i++, att++)
    {
        TfmPoint(&xnew, &ynew, att->x, att->y, xfrom, yfrom, xto, yto, tfmp);
        if (att->nneighbours <= 1  &&  att->has_symbol  &&  na > 1)
        {
// fprintf(stderr, "x-distance\n");
            score += 1.0*SQR(xnew-att->xlig);
            score += 1.0*SQR(ynew-att->ylig);
        }
        else
        {
            TfmPoint(&xdir, &ydir, att->xdir, att->ydir, 0.0, 0.0, 0.0, 0.0, tfmp);
// fprintf(stderr, "transformed point %g/%g mapping to %g/%g\n", xnew+xdir, ynew+ydir, att->xlig, att->ylig);
            score += SQR(xnew+xdir-att->xlig);
            score += SQR(ynew+ydir-att->ylig);
        }
    }
    if (natoms <= 2)
       return 0.01*score + 1.0*natoms*((yamax-yamin) - 0.00001*(xamax-xamin));
    else if (na <= 1)
       return 1.00*score + 0.0*na*(ymax-ymin) + 0.0*natoms*((yamax-yamin) - 0.1*(xamax-xamin));
    else
       return 1.00*score + 0.0100*SQR(ymax-ymin) + 0.0100*natoms*((yamax-yamin) - (xamax-xamin) + 1.0*SQR(target_height-(yamax-yamin)));
}

#define DEGREE_STEP 15

double FindBestTransformation(struct tfm_t *tfmp, struct attachment_t *attachments, int na,
                              double xfrom, double yfrom, double xto, double yto,
                              struct reaccs_molecule_t *mp, int color, int bottom_color, int top_color, double target_height,
                              int debug)
/**
 * Find the transformation that best rotates (and moves) the attachments[0..na-1] wrt. their connected other atoms.
 */
{
    struct tfm_t tfm;
    double score, new_score;
    int i;

    tfm.cos_alpha = 1.0; tfm.sin_alpha = 0.0; tfm.do_flip = FALSE;
    (*tfmp) = tfm;
    score = ScoreTransformation(&tfm, attachments, na, xfrom, yfrom, xto, yto, mp, color, bottom_color, top_color, target_height);
    for (i=1; i<=360/DEGREE_STEP; i++)
    {
        tfm.cos_alpha =  cos(i*DEGREE_STEP*3.14159265358979/180.0);
        tfm.sin_alpha = -sin(i*DEGREE_STEP*3.14159265358979/180.0);
        new_score = ScoreTransformation(&tfm, attachments, na, xfrom, yfrom, xto, yto, mp, color, bottom_color, top_color, target_height);
        if (new_score < score)
        {
if (debug) fprintf(stderr, "a1:(%d) score=%g > new_score=%g\n", i*DEGREE_STEP, score, new_score);
            (*tfmp) = tfm;
            score = new_score;
        }
        else
        {
if (debug) fprintf(stderr, "a2:(%d) score=%g <= new_score=%g\n", i*DEGREE_STEP, score, new_score);
        }
    }
    tfm.do_flip = TRUE;
    for (i=0; i<=360/DEGREE_STEP; i++)
    {
        tfm.cos_alpha =  cos(i*DEGREE_STEP*3.14159265358979/180.0);
        tfm.sin_alpha = -sin(i*DEGREE_STEP*3.14159265358979/180.0);
        new_score = ScoreTransformation(&tfm, attachments, na, xfrom, yfrom, xto, yto, mp, color, bottom_color, top_color, target_height);
        if (new_score < score)
        {
if (debug) fprintf(stderr, "b1:(%d) score=%g > new_score=%g\n", i*DEGREE_STEP, score, new_score);
            (*tfmp) = tfm;
            score = new_score;
        }
        else
        {
if (debug) fprintf(stderr, "b2:(%d) score=%g <= new_score=%g\n", i*DEGREE_STEP, score, new_score);
        }
    }
    return score;
}

double AlignFragment(struct reaccs_molecule_t *mp, neighbourhood_t *nbp,
                     int color,
                     double top_y,     // y coordinate of sequence string above this fragment (or some estimate if there is none)
                     double bottom_y,  // y-coordinate of sequence string below this fragment
                     int top_color, int bottom_color, double min_spacing,
                     int nsneighbours)
/**
 * Align the atoms colored color to be intercalated between the shortcut sequences colored top_color and bottom_color.
 * Returns the required spacing to fit 4/3 bond-lengths plus the Y dimension of the fragment between the two rows of shortcuts.
 * The subsequent shortcut fragments are supposed to be moved by the caller if needed.
 */
{
    struct attachment_t attachments[MAX_ATTACHMENTS], *att;
    int na;
    int j, n;
    struct reaccs_atom_t *ap, *ap1, *ap2;
    struct reaccs_bond_t *bp;
    double xcenter, ycenter;
    double x, y;
    double xnew, ynew;
    double ymin, ymax;
    struct tfm_t tfm;

    na = CollectAttachments(mp, color, attachments, nbp);

// fprintf(stderr, "sliding fragment with %d attachments and %d sequence neighbours between y=%g and y=%g\n", na, nsneighbours, top_y, bottom_y);

    // find center of gravity of attachment atoms to be used as rotation center and fragment reference point
    xcenter = ycenter = xnew = ynew = 0.0;
    n = 0;
    for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
    {
        if (ap->color != color) continue;
        // ycenter += ap->y;
        n++;
    }
    if (nsneighbours == 1)      // align attachments to single attaching sequence
    {
       for (j=0, att=attachments; j<na; j++, att++)
       {
           xcenter += att->x;
           ycenter += att->y;
           xnew += att->xlig;
           ynew += att->ylig;
// fprintf(stderr, "1: d%d(%g/%g) => (%g/%g), lig=(%g/%g)\n", att->nneighbours, att->x, att->y, att->xdir, att->ydir, att->xlig, att->ylig);
       }
       xcenter /= na; ycenter /= na;   // center of gravity of fragment attachment points
       xnew /= na;                     // center of gravity of sequence atoms to which attachment points are to be connected.
       ynew /= na;                     // y coordinate may be between sequences!
       if (na > 1  &&  n > 2)
          ynew = bottom_y+sqrt(na+0.1*(n-2))*1.54; // align fragment rotation point sqrt(na+0.1*(n-2)) std bond above bottom sequence for bigger multiattached fragments
       else
          ynew = bottom_y+1.54;           // align fragment rotation point 1 std bond above bottom sequence
    }
    else        // splice fragment between top_y/top_color and bottom_y/bottom_color
    {
       for (j=0, att=attachments; j<na; j++, att++)
       {
           xcenter += att->x;
           ycenter += att->y;
           xnew += att->xlig;
           ynew += att->ylig;
// fprintf(stderr, "2: d%d(%g/%g) => (%g/%g), lig=(%g/%g)\n", att->nneighbours, att->x, att->y, att->xdir, att->ydir, att->xlig, att->ylig);
       }
       xcenter /= na; ycenter /= na;   // center of gravity of fragment attachment points
       xnew /= na;                     // center of gravity of sequence atoms to which attachment points are to be connected.
       ynew /= na;                     // y coordinate may be between sequences!
       ynew = (top_y+bottom_y)/2.0;    // put fragment rotation point at the middle of the two sequences
    }
// fprintf(stderr, "na=%d, xcenter=%g, ycenter=%g, xnew=%g, ynew=%g\n", na, xcenter, ycenter, xnew, ynew);
    FindBestTransformation(&tfm, attachments, na, xcenter, ycenter, xnew, ynew, mp, color, bottom_color, top_color, top_y-bottom_y,
          FALSE  ||  bottom_color == 0);
    // now we apply the transformation and collect y extent of resulting fragment
    ymin = ymax = ynew;
    for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
    {
       if (ap->color != color) continue;
       TfmPoint(&x, &y, ap->x, ap->y, xcenter, ycenter, xnew, ynew, &tfm);
       ap->x = x; ap->y = y;
       if (ymin > ap->y) ymin = ap->y;
       if (ymax < ap->y) ymax = ap->y;
    }
// fprintf(stderr, "ymin=%g, ymax=%g\n", ymin, ymax);
    // also consider attachment points
if (TRUE)
    for (j=0, att=attachments; j<na; j++, att++)
    {
       x = att->x+att->xdir; y = att->y+att->ydir;
       TfmPoint(&x, &y, x, y, xcenter, ycenter, xnew, ynew, &tfm);
       if (ymin > y) ymin = y;
       if (ymax < y) ymax = y;
// fprintf(stderr, "ymin=%g, ymax=%g, x=%g, y=%g, xnew=%g, ynew=%g\n", ymin, ymax, att->x+att->xdir, att->y+att->ydir, x, y);
    }
    TfmPoint(&x, &y, xcenter, ycenter, xcenter, ycenter, xnew, ynew, &tfm);
// fprintf(stderr, "after tfm: xcenter=%g, ycenter=%g, xnew=%g, ynew=%g\n", x, y, xnew, ynew);
    // find the y-extent of the newly rotated fragment
// fprintf(stderr, "old-min_spacing=%g, ", min_spacing);
    min_spacing = min_spacing > (FRAME_HEIGHT+ymax-ymin) ? min_spacing : (FRAME_HEIGHT+ymax-ymin);
// fprintf(stderr, "new-min_spacing=%g\n", min_spacing);
    for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
       if (ap->color == color)
       {
          // ap->x += xnew;
          // ap->y -= (top_y-bottom_y-min_spacing)/2.0;
          // if (bottom_color==0)
              // ap->y -= 0.5*min_spacing;
          // else
              // ap->y -= 0.5*min_spacing;
/*
if (0 == strcmp(ap->atom_symbol, "R"))
    fprintf(stderr, "%d:{%s}(%g/%g) ", j+1, ap->atext, ap->x, ap->y);
else
    fprintf(stderr, "%d:%s(%g/%g) ", j+1, ap->atom_symbol, ap->x, ap->y);
*/
       }
// fprintf(stderr, "\n");
    return min_spacing;
}

double FragmentRenderedHeight(struct reaccs_molecule_t *mp, int color)
/**
 * Compute the total height of the already rendered fragment with color.
 */
{
    struct reaccs_atom_t *ap;
    int i, nmatch;
    double ymin, ymax;

    ymin = 1e10; ymax = -1e10; nmatch = 0;
    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
    {
       if (ap->color != color) continue;
       nmatch++;
       if (ymin > ap->y) ymin = ap->y;
       if (ymax < ap->y) ymax = ap->y;
    }
    if (nmatch == 0) return 1.5*FRAME_HEIGHT;
    return (ymax-ymin)+1.5*FRAME_HEIGHT;
}

double FragmentHeight(struct reaccs_molecule_t *mp, int color)
/**
 * Compute the width of the color fragment perpendicular to the largest dimension.
 */
{
    struct reaccs_atom_t *ap1, *ap2;
    double d, dmax;
    double lambda, ymin, ymax;
    int i1, i2;
    point r1, r2, ex, ey;
    
    // search for largest distance
    dmax = -1.0;
    for (i1=0, ap1=mp->atom_array; i1<mp->n_atoms; i1++, ap1++)
    {
        if (ap1->color != color) continue;
        for (i2=i1+1, ap2=mp->atom_array+(i1+1); i2<mp->n_atoms; i2++, ap2++)
        {
            if (ap2->color != color) continue;
            d = SQR(ap1->x-ap2->x)+SQR(ap1->y-ap2->y);
            if (d > dmax)
            {
                dmax = d;
                ex[0] = ap2->x - ap1->x; ex[1] = ap2->y - ap1->y;
                // remember one atom of the diameter
                r1[0] = ap1->x; r1[1] = ap1->y;
            }
        }
    }
    if (dmax <= 0) return 1.5*FRAME_HEIGHT;
    dmax = sqrt(dmax);
    // compute unit vector along largest dimension
    ex[0] /= dmax; ex[1] /= dmax;
    // compute unit vector perpendicular to largest dimension
    ey[0] = ex[1]; ey[1] = -ex[0];
    ymin = 1.0e10; ymax = -1.0e10;
    for (i1=0, ap1=mp->atom_array; i1<mp->n_atoms; i1++, ap1++)
    {
        if (ap1->color != color) continue;
        lambda = (ey[0]*(r1[0]-ap1->x) + ey[1]*(r1[1]-ap1->y));
        if (lambda < ymin) ymin = lambda;
        if (lambda > ymax) ymax = lambda;
    }
    return (ymax-ymin)+1.5*FRAME_HEIGHT;
}

void AlignShortcutsAndFragments(struct reaccs_molecule_t *mp, int *is_ring_atom, int *shortcut_class, neighbourhood_t *nbp)
/*
 * Orients strings of shortcut 'atoms' horizontally, and places them and the connected non-shortcut fragments in a vertical array.
 */
{
    double xoffset, yoffset, xnew, ynew, sin_alpha, cos_alpha, best_score;
    struct reaccs_atom_t *ap, *aph;
    struct reaccs_bond_t *bp;
    int i, j, bestj, nlinks, col, first_color;
    neighbourhood_t *nbph;
    struct fragment_t *fragments, *fp, *fph, **farray;
    int ifrag, jfrag, nfrag;

    atom_pair *edges;
    int nedges;
    point *coords;
    int nnodes;
    int *numbers;
    point p1, p1p;
    point p2, p2p;
    double deltah, new_deltah;

    if (IsNULL(mp)) return;

    // collect information about the fragments to be layed out.
    farray = CollectFragmentInfo(mp, shortcut_class, nbp, &nfrag);

    edges = TypeAlloc(mp->n_bonds, atom_pair);
    coords = TypeAlloc(mp->n_atoms, point);
    numbers = TypeAlloc(mp->n_atoms, int);

    first_color = -1;

    for (ifrag=0; ifrag<nfrag; ifrag++)
        for (fp=farray[ifrag]; !IsNULL(fp); fp=fp->next)
            fp->used = FALSE;

    // first loop through FT_SHORTCUTS fragments and put them on separate lines
    yoffset = 0.0;
    for (ifrag=0; ifrag<nfrag; ifrag++)
    {
        fp = farray[ifrag];
        if (IsNULL(fp)) continue;
        if (ifrag != 0  &&  !IsNULL(farray[ifrag-1])) yoffset += 1.0*FRAME_HEIGHT; // fp->radius+farray[ifrag-1]->radius;
        // skip non-sequence fragments
        if (fp->fragment_type != FT_SHORTCUTS) continue;
// fprintf(stderr, "SHORTCUT_LOOP: ifrag = %d, color = %d, yoffset = %g, ycenter = %g\n", ifrag, fp->color, yoffset, fp->ycenter);
// for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
// {
// if (ap->color != fp->color) continue;
//     if (0 == strcmp(ap->atom_symbol, "R"))
//         fprintf(stderr, "%d:{%s} ", i+1, ap->atext);
//     else
//         fprintf(stderr, "%d:%s ", i+1, ap->atom_symbol);
// }
// fprintf(stderr, "\n");
        // remember first color
        if (first_color == -1) first_color = fp->color;
        fp->used = TRUE;
        // search for start of shortcut chain
        col = -1;
        for (i=0, ap=mp->atom_array, nbph=nbp; i<mp->n_atoms; i++, ap++, nbph++)
        {
            if (ap->color != fp->color) continue;
            if (shortcut_class[i] & SHORTCUT_START)
            {
                col = ap->color;
                // search for shortcut ligand of this atom in same fragment
                for (j=0; j<nbph->n_ligands; j++)
                    if (col == mp->atom_array[nbph->atoms[j]].color  &&
                        0 != (shortcut_class[nbph->atoms[j]] & (SHORTCUT_INNER | SHORTCUT_END)))
                    {
                        break;
                    }
                if (j == nbph->n_ligands)
                {
                    continue;     // no sequence to align
                }
                aph = mp->atom_array + nbph->atoms[j];
                GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, col);
                p1[0] = ap->x; p1[1] = ap->y;
                p2[0] = aph->x; p2[1] = aph->y;
                p1p[0] = 0.0;  p1p[1] = (-1)*yoffset;
                p2p[0] = 1.54; p2p[1] = (-1)*yoffset;
                fp->ycenter = (-1)*yoffset;
// fprintf(stderr, "before: %s(%d) at %g/%g - %s(%d) at %g/%g\n",
// ap->atom_symbol, i+1, ap->x, ap->y,
// aph->atom_symbol, (aph-mp->atom_array)+1, aph->x, aph->y);
                TransformPoints(coords, nnodes, p1, p2, p1p, p2p);
                for (j=0; j<mp->n_atoms; j++)
                   if (mp->atom_array[j].color == col)
                   {
                      mp->atom_array[j].x = coords[numbers[j]][0];
                      mp->atom_array[j].y = coords[numbers[j]][1];
// if (0 == strcmp(mp->atom_array[j].atom_symbol, "R"))
// fprintf(stderr, "%d:{%s}(%g/%g) ", j+1, mp->atom_array[j].atext, mp->atom_array[j].x, mp->atom_array[j].y);
// else
// fprintf(stderr, "%d:%s(%g/%g) ", j+1, mp->atom_array[j].atom_symbol, mp->atom_array[j].x, mp->atom_array[j].y);
                   }
// fprintf(stderr, "\n");
// fprintf(stderr, "after: %s(%d) at %g/%g - %s(%d) at %g/%g\n",
// ap->atom_symbol, i+1, ap->x, ap->y,
// aph->atom_symbol, (aph-mp->atom_array)+1, aph->x, aph->y);
                break;
            }
        }
        if (i == mp->n_atoms)
        {
fprintf(stderr,"## Did not find atoms for fragment color %d\n", fp->color);
            break;
        }
    }

// if (TRUE)
// for (ifrag=0; ifrag<nfrag; ifrag++)
// {
// if (ifrag == 0) fprintf(stderr, "\n");
//     fp = farray[ifrag];
//     if (IsNULL(fp)) continue;
//     if (fp->fragment_type != FT_SHORTCUTS) continue;
// fprintf(stderr, "SHORTCUT_LOOP(1): ifrag = %d, color = %d, yoffset = %g, ycenter = %g\n", ifrag, fp->color, yoffset, fp->ycenter);
// for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
// {
// if (ap->color != fp->color) continue;
//     if (0 == strcmp(ap->atom_symbol, "R"))
//         fprintf(stderr, "%d:{%s} ", i+1, ap->atext);
//     else
//         fprintf(stderr, "%d:%s ", i+1, ap->atom_symbol);
// }
// fprintf(stderr, "\n");
// }
// fprintf(stderr, "\n");

    // move adjacent shortcut sequences further appart if needed
// if (FALSE)
    for (ifrag=0; ifrag<nfrag-1; ifrag++)  // make space to fit interspaced fragments and also improve per fragment layout
    {
        // assuming the neighbouring FT_ATOMS fragments are threaded to the FT_SHORTCUTS fragment.
        if (farray[ifrag]->fragment_type != FT_SHORTCUTS) continue;
        deltah     = farray[ifrag]->ycenter-farray[ifrag+1]->ycenter;
// fprintf(stderr, "y[%d] = %g, y[%d] = %g\n", ifrag, farray[ifrag]->ycenter, ifrag+1, farray[ifrag+1]->ycenter);
        // find biggest fragment
        for (fragments = farray[ifrag]->next; !IsNULL(fragments); fragments = fragments->next)
        {
            if (fragments->nsneighbours == 1)   // only improve single attachment fragments
            {
               ImproveFragmentByAllBondFlips(mp, nbp, fragments->color, CENTER_DISTANCE, farray[ifrag]->color, -1);
               continue;     // to be placed later before this shortcut string
            }
            ImproveFragmentByAllBondFlips(mp, nbp, fragments->color, TOP_BOTTOM_SEPARATION, farray[ifrag]->color, farray[ifrag+1]->color);
            // AlignDoubleAttachmentFragment(mp, fragment->color, farray[ifrag]->color, farray[ifrag+1]->color);
            new_deltah = FragmentHeight(mp, fragments->color);
// fprintf(stderr, "0(%d/%d): deltah=%g, new_deltah=%g, nsymbols=%d, first_atom=%d, nsneighbours=%d\n",
//       ifrag, farray[ifrag]->nsymbols, deltah, new_deltah, fragments->nsymbols, fragments->first_index+1, fragments->nsneighbours);
            if (new_deltah > deltah)
            {
                deltah = new_deltah;
            }
        }
        if (deltah > farray[ifrag]->ycenter-farray[ifrag+1]->ycenter)
        {
            // make room between the two shortcut strings by moving the subsequent ones further down
            deltah -= farray[ifrag]->ycenter-farray[ifrag+1]->ycenter;
// fprintf(stderr, "0: shifting from %d by %g\n", ifrag+1, deltah);
            for (jfrag = ifrag+1; jfrag < nfrag; jfrag++)
            {
                farray[jfrag]->ycenter -= deltah;
// fprintf(stderr, "0: shifting fragment %d/%d to %g\n", jfrag, farray[jfrag]->nsymbols, farray[jfrag]->ycenter);
                for (i=0, ap=mp->atom_array, nbph=nbp; i<mp->n_atoms; i++, ap++, nbph++)
                {
                    if (ap->color != farray[jfrag]->color) continue;
                    ap->y -= deltah;
                }
            }
        }
    }
// fprintf(stderr, "\n");

// if (TRUE)
// for (ifrag=0; ifrag<nfrag; ifrag++)
// {
// if (ifrag == 0) fprintf(stderr, "\n");
//     fp = farray[ifrag];
//     if (IsNULL(fp)) continue;
//     if (fp->fragment_type != FT_SHORTCUTS) continue;
// fprintf(stderr, "SHORTCUT_LOOP(2): ifrag = %d, color = %d, yoffset = %g, ycenter = %g\n", ifrag, fp->color, yoffset, fp->ycenter);
// for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
// {
// if (ap->color != fp->color) continue;
//     if (0 == strcmp(ap->atom_symbol, "R"))
//         fprintf(stderr, "%d:{%s} ", i+1, ap->atext);
//     else
//         fprintf(stderr, "%d:%s ", i+1, ap->atom_symbol);
// }
// fprintf(stderr, "\n");
// }
// fprintf(stderr, "\n");

    // make space before sequence
// if (FALSE)
    for (ifrag=1; ifrag<nfrag; ifrag++) // starts with 1 because first row is not moved
    {
        // assuming the neighbouring FT_ATOMS fragments are threaded to the FT_SHORTCUTS fragment.
        if (farray[ifrag]->fragment_type != FT_SHORTCUTS) continue;
        deltah     = farray[ifrag-1]->ycenter-farray[ifrag]->ycenter;
        // find biggest fragment
        for (fragments = farray[ifrag]->next; !IsNULL(fragments); fragments = fragments->next)
        {
            if (fragments->nsneighbours > 1) continue;     // deal with fragments before the shortcut string if any
            new_deltah = FragmentHeight(mp, fragments->color);
// fprintf(stderr, "1(%d/%d): deltah=%g, new_deltah=%g, nsneighbours=%d, nsymbols=%d\n", ifrag, farray[ifrag]->nsymbols, deltah, new_deltah, fragments->nsneighbours, fragments->nsymbols);
            if (new_deltah > deltah)
            {
                deltah = new_deltah;
            }
        }
        if (deltah > farray[ifrag-1]->ycenter-farray[ifrag]->ycenter)
        {
            // make room between the two shortcut strings by moving the subsequent ones further down
            deltah -= farray[ifrag-1]->ycenter-farray[ifrag]->ycenter;
// fprintf(stderr, "1: shifting from %d by %g\n", ifrag, deltah);
            for (jfrag = ifrag; jfrag < nfrag; jfrag++)
            {
                farray[jfrag]->ycenter -= deltah;
// fprintf(stderr, "1: shifting fragment %d/%d to %g\n", jfrag, farray[jfrag]->nsymbols, farray[jfrag]->ycenter);
                for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
                {
                    if (ap->color != farray[jfrag]->color) continue;
                    ap->y -= deltah;
                }
            }
        }
    }
// fprintf(stderr, "\n");

// if (TRUE)
// for (ifrag=0; ifrag<nfrag; ifrag++)
// {
// if (ifrag == 0) fprintf(stderr, "\n");
//     fp = farray[ifrag];
//     if (IsNULL(fp)) continue;
//     if (fp->fragment_type != FT_SHORTCUTS) continue;
// fprintf(stderr, "SHORTCUT_LOOP(3): ifrag = %d, color = %d, yoffset = %g, ycenter = %g\n", ifrag, fp->color, yoffset, fp->ycenter);
// for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
// {
// if (ap->color != fp->color) continue;
//     if (0 == strcmp(ap->atom_symbol, "R"))
//         fprintf(stderr, "%d:{%s} ", i+1, ap->atext);
//     else
//         fprintf(stderr, "%d:%s ", i+1, ap->atom_symbol);
// }
// fprintf(stderr, "\n");
// }
// fprintf(stderr, "\n");

// fprintf(stderr, "\nnow, we place the FT_ATOMS fragment that connect two neighbouring FT_SHORTCUT fragments\n");
    // now, we place the FT_ATOMS fragment that connect two neighbouring FT_SHORTCUT fragments
if (TRUE)
    for (ifrag=0; ifrag<nfrag; ifrag++)
    {
        // assuming the neighbouring FT_ATOMS fragments are threaded to the FT_SHORTCUTS fragment.
        if (farray[ifrag]->fragment_type != FT_SHORTCUTS) continue;
// fprintf(stderr, "y[%d] = %g\n", ifrag, farray[ifrag]->ycenter);
        fragments = farray[ifrag]->next;
        fp = (struct fragment_t *)NULL;
        for (fragments = farray[ifrag]->next; !IsNULL(fragments); fragments=fragments->next)
        {
// TODO Check: if (fragments->nsymbols > 2) continue;
            // grab the graph of the new fragment
            if (fragments->nsneighbours == 1)  // to be placed next round
            {
               // NOP
            }
            else if (fragments->nsneighbours >= 2  &&  ifrag+1 < nfrag  &&
                     fragments->neighbour_cols[0] == farray[ifrag]->color)
            {                                  // slide between this and next FT_SHORTCUTS fragment
// fprintf(stderr, "2: placing fragment with nsneighbours=%d, nsymbols=%d, nneighbours=%d before string %d\n", fragments->nsneighbours, fragments->nsymbols, fragments->nneighbours, ifrag);
                deltah     = farray[ifrag]->ycenter-farray[ifrag+1]->ycenter;
                GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, fragments->color);
                new_deltah = AlignFragment(mp, nbp, fragments->color,
                                           farray[ifrag]->ycenter, farray[ifrag+1]->ycenter,
                                           fragments->neighbour_cols[0], fragments->neighbour_cols[1],
                                           deltah, fragments->nsneighbours);
                new_deltah = FragmentRenderedHeight(mp, fragments->color);
                if (new_deltah > deltah)    // make room between the two shortcut strings by moving the subsequent ones further down
                {
                    new_deltah -= deltah;
// fprintf(stderr, "2b: shifting aligned fragment starting with atom %d by %g\n", fragments->first_index+1, new_deltah/2.0);
                    // move the current fragment down by 50% of slot change
                    for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
                    {
                        if (ap->color != fragments->color) continue;
                        ap->y -= new_deltah/2.0;
                    }
// fprintf(stderr, "2b: shifting from %d by %g\n", ifrag+1, new_deltah);
                    // move subsequent fragments down by full slot change
                    for (jfrag = ifrag+1; jfrag < nfrag; jfrag++)
                    {
                        farray[jfrag]->ycenter -= new_deltah;
// fprintf(stderr, "2b: shifting fragment %d/%d to %g\n", jfrag, farray[jfrag]->nsymbols, farray[jfrag]->ycenter);
// fprintf(stderr, "  y[%d] = %g\n", jfrag, farray[jfrag]->ycenter);
                        for (i=0, ap=mp->atom_array, nbph=nbp; i<mp->n_atoms; i++, ap++, nbph++)
                        {
                            if (ap->color != farray[jfrag]->color) continue;
                            ap->y -= new_deltah;
                        }
                    }
                }
                else
                   /* NOP since AlignFragment has already placed the fragment*/;
                fragments->used = TRUE;
            }
            else
            {
// fprintf(stderr, "IGNORING fragment with %d symbols\n", fragments->nsymbols);
            }
        }
    }

    // now, we place the FT_ATOMS fragment that are positioned above a shortcut string
// fprintf(stderr, "\nnow, we place the FT_ATOMS fragment that are positioned above a shortcut string\n");
    for (ifrag=0; ifrag<nfrag; ifrag++)
    {
        // assuming the neighbouring FT_ATOMS fragments are threaded to the FT_SHORTCUTS fragment.
        if (farray[ifrag]->fragment_type != FT_SHORTCUTS) continue;
// fprintf(stderr, "y[%d] = %g\n", ifrag, farray[ifrag]->ycenter);
        fragments = farray[ifrag]->next;
        fp = (struct fragment_t *)NULL;
        for (fragments = farray[ifrag]->next; !IsNULL(fragments); fragments=fragments->next)
        {
// TODO Check: if (fragments->nsymbols > 2  &&  fragments->nneighbours > 1) continue;
            // grab the graph of the new fragment
            if (fragments->nsneighbours == 1)  // slide this fragment before current fragment
            {
// fprintf(stderr, "3: placing fragment with nsneighbours=%d, nsymbols=%d, nneighbours=%d before string %d\n", fragments->nsneighbours, fragments->nsymbols, fragments->nneighbours, ifrag);
                if (ifrag > 0) deltah = farray[ifrag-1]->ycenter - farray[ifrag]->ycenter;
                else           deltah = FRAME_HEIGHT;
                GetColoredGraph(mp, edges,  &nedges, coords, &nnodes, numbers, fragments->color);
                new_deltah = AlignFragment(mp, nbp, fragments->color,
                                           (ifrag > 0 ?  farray[ifrag-1]->ycenter : FragmentHeight(mp, fragments->color)-FRAME_HEIGHT/2.0),
                                           farray[ifrag]->ycenter,
                                           0, fragments->neighbour_cols[0],
                                           deltah,
                                           fragments->nsneighbours);
                new_deltah = FragmentRenderedHeight(mp, fragments->color);
// fprintf(stderr, "3: deltah=%g, new_deltah=%g\n", deltah, new_deltah);
                if (ifrag > 0  &&  new_deltah > deltah)    // make room between the two shortcut strings by moving the subsequent ones further down
                {
                    // move the current and subsequent FT_SHORTCUT fragments and the connected FT_ATOMS fragments down by the necessary distance
                    new_deltah -= deltah;
                    for (fph=farray[ifrag]->next; !IsNULL(fph); fph=fph->next)
                    {
// fprintf(stderr, "moving FT_ATOMS fragment colored %d by %g\n", fph->color, new_deltah);
                        for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
                        {
                            if (ap->color != fph->color) continue;
                            // if (ap->color == fragments->color) continue;     // already placed
                            ap->y -= new_deltah;
                        }
                    }
// fprintf(stderr, "3a: shifting from %d by %g\n", ifrag, new_deltah);
                    for (jfrag = ifrag; jfrag < nfrag; jfrag++)
                    {
                        farray[jfrag]->ycenter -= new_deltah;
// fprintf(stderr, "moving FT_SHORTCUT fragment colored %d by %g\n", farray[jfrag]->color, new_deltah);
// fprintf(stderr, "  y[%d] = %g\n", jfrag, farray[jfrag]->ycenter);
                        for (i=0, ap=mp->atom_array, nbph=nbp; i<mp->n_atoms; i++, ap++, nbph++)
                        {
                            if (ap->color != farray[jfrag]->color) continue;
                            ap->y -= new_deltah;
                        }
                    }
                }
                // now, place the fragment between the slots
                // for (i=0, ap=mp->atom_array, nbph=nbp; i<mp->n_atoms; i++, ap++, nbph++)
                // {
                    // if (ap->color != fragments->color) continue;
                    // ap->y += 0.5*(farray[ifrag]->ycenter+farray[ifrag+1]->ycenter);
                // }
                fragments->used = TRUE;
            }
            else if (fragments->nsneighbours >= 2  &&  ifrag+1 < nfrag  &&
                     fragments->neighbour_cols[0] == farray[ifrag]->color)
            {                                  // slide between this and next FT_SHORTCUTS fragment
               // NOP
            }
            else
            {
// fprintf(stderr, "IGNORING fragment with %d symbols\n", fragments->nsymbols);
            }
        }
    }

// if (TRUE)
// for (ifrag=0; ifrag<nfrag; ifrag++)
// {
// if (ifrag == 0) fprintf(stderr, "\n");
//     fp = farray[ifrag];
//     if (IsNULL(fp)) continue;
//     if (fp->fragment_type != FT_SHORTCUTS) continue;
// fprintf(stderr, "SHORTCUT_LOOP(4): ifrag = %d, color = %d, yoffset = %g, ycenter = %g\n", ifrag, fp->color, yoffset, fp->ycenter);
// for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
// {
// if (ap->color != fp->color) continue;
//     if (0 == strcmp(ap->atom_symbol, "R"))
//         fprintf(stderr, "%d:{%s} ", i+1, ap->atext);
//     else
//         fprintf(stderr, "%d:%s ", i+1, ap->atom_symbol);
// }
// fprintf(stderr, "\n");
// }

    // now placing FT_ATOMS fragments aligned to connected FT_SHORTCUTS fragments
    yoffset = 0.0;

    // make sure all layed out bits have same color
    for (ifrag=0; ifrag<nfrag; ifrag++)
        for (fp=farray[ifrag]; !IsNULL(fp); fp=fp->next)
            if (fp->used)
                for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
                    if (fp->color == ap->color) ap->color = first_color;

    for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
    {
       // if (!IsShortcutLikeAtom(mp, j, nbp)) ap->mapping = 1+j;
       if (ap->color == first_color)
          ap->color |= KEEP_POSITION;
    }

    MyFree((char *)numbers);
    MyFree((char *)edges); MyFree((char *)coords);

    // deallocate fragment info
    for (ifrag=0; ifrag<nfrag; ifrag++)
    {
        fragments = farray[ifrag];
        while (!IsNULL(fragments))
        {
            fp = fragments->next;
            MyFree((char *)fragments);
            fragments = fp;
        }
    }
    MyFree((char *)farray);
}

static void SortNeighbourhoodByAtomIndex(struct reaccs_molecule_t *mp, neighbourhood_t *nbp)
{
   int n, i, j, h;
   for (n=0; n<mp->n_atoms; n++, nbp++)
   {
      for (i=1; i<nbp->n_ligands; i++)
         for (j=i-1; j>=0; j--)
            if (nbp->atoms[j] > nbp->atoms[j+1])
            {
               h = nbp->atoms[j]; nbp->atoms[j] = nbp->atoms[j+1]; nbp->atoms[j+1] = h;
               h = nbp->bonds[j]; nbp->bonds[j] = nbp->bonds[j+1]; nbp->bonds[j+1] = h;
            }
            else
               break;
   }
}

struct reaccs_molecule_t *LayoutMolecule(struct reaccs_molecule_t *mp)
/*
 * Imposes a new layout onto the molecule *mp. It assumes that the
 * atoms and bonds of prelayouted fragments have been colored with
 * the same color.
 * The function returns a copy of the molecule with new coordinates.
 * The original molecule is not changed.
 */
{
   struct reaccs_molecule_t *result;
   struct reaccs_bond_t *bp;
   int i;
   int *is_ring_atom, *ring_count;
   neighbourhood_t *nbp;

   bond_set_node *ring_list, *plist;
   int *ring_size;

   int *shortcut_class;
   int process_shortcuts;

   result = CopyMolecule(mp);
   srand(1);    /* make sure the 'random' coordinates are reproducible */
   RandomCoordinates(result);

   ScaleByFixedFragments(result);

   nbp   = TypeAlloc(result->n_atoms, neighbourhood_t);
   SetupNeighbourhood(result,nbp,result->n_atoms);
   SortNeighbourhoodByAtomIndex(result,nbp);

   for (i=0; i<result->n_atoms; i++)
      if (nbp[i].n_ligands == 1)/* terminal atoms don't have DB-stereo */
         if ((ALL_BOND_TYPES&result->bond_array[nbp[i].bonds[0]].bond_type) == DOUBLE)
            result->bond_array[nbp[i].bonds[0]].stereo_symbol = NONE;

   is_ring_atom = TypeAlloc(result->n_atoms, int);
   ring_count   = TypeAlloc(result->n_bonds, int);
   for (i=0; i<result->n_atoms; i++) is_ring_atom[i] = FALSE;
   for (i=0; i<result->n_bonds; i++) ring_count[i] = 0;

   ring_size    = TypeAlloc(result->n_bonds, int);
   for (i=0; i<result->n_bonds; i++) ring_size[i] = 0;
   ring_list = MoleculeRings(result);
   for (plist=ring_list; plist; plist=plist->next)
   {
      for (i=0, bp=result->bond_array; i<result->n_bonds; i++, bp++)
         if (IsMember(plist->bond_set,i))
         {
            is_ring_atom[bp->atoms[0]-1] = TRUE;
            is_ring_atom[bp->atoms[1]-1] = TRUE;
            if (ring_size[i] == 0  ||  ring_size[i] > plist->cardinality)
               ring_size[i] = plist->cardinality;
            bp->rsize_flags |= MAX(9, plist->cardinality);
            ring_count[i] += 1;
         }
   }
   DisposeBondSetList(ring_list);

   // perceive shortcut classes if any
   shortcut_class = TypeAlloc(result->n_atoms, int);
   process_shortcuts = PerceiveShortcutClasses(result, is_ring_atom, shortcut_class, nbp, 30);

   if (process_shortcuts) LayoutShortcutChains(result, is_ring_atom, shortcut_class, nbp);

   LayoutLargeLactams(result, is_ring_atom, ring_size, nbp); /* trial implementation */

   LayoutRings(result);

   FlipRingCarbonyls(result, is_ring_atom, ring_size, nbp); /* trial implementation */

// FixTransDBsInRings(result, nbp); /* not yet implemented */
//

   LinkRingFragments(result, is_ring_atom);

   SproutRingSubstituents(result, is_ring_atom, ring_size, nbp);

   LayoutChainAtoms(result, is_ring_atom, nbp);

   FixLinearBridgeHeads(result, ring_count, nbp);

   LinkRemainingFragments(result);

   FixInwardHydrogens(result, ring_count, nbp);

   ImproveMoleculeByBondFlip(result, nbp, ring_size, is_ring_atom, 0, FALSE);

   ClearFlipFlags(result);

   ImproveCollisions(result, nbp, is_ring_atom, ring_count);

   LayoutRubberFragments(result);  // testing

   if (process_shortcuts) AlignShortcutsAndFragments(result, is_ring_atom, shortcut_class, nbp);

   MergeConnectedFragmentColors(result);

   LayoutFragments(result, ring_size, ring_count);

   LayoutBondStereo(result, nbp, ring_size);

   // try flipping single bonds to improve layout
   ImproveMoleculeByBondFlip(result, nbp, ring_size, is_ring_atom, 0, TRUE);
   ClearFlipFlags(result);

   LayoutAtomStereo(result, nbp, is_ring_atom, ring_size);

   ClearDBStereoInSmallRings(result, ring_size);

// MakeLandscape(result); */

   for (i=0, bp=result->bond_array; i<mp->n_bonds; i++, bp++)
      bp->bond_type &= ALL_BOND_TYPES;	/* upper nibble was used for internal flags */

   if (shortcut_class) MyFree((char *)shortcut_class);
   if (ring_size)    MyFree((char *)ring_size);
   if (ring_count)   MyFree((char *)ring_count);
   if (is_ring_atom) MyFree((char *)is_ring_atom);
   if (nbp) 	     MyFree((char *)nbp);

   return (result);
}

#ifdef DEBUG
main(int argc, char *argv[])
{
   Fortran_FILE *finp;
   FILE *foutp;
   struct reaccs_molecule_t *mp, *mph;

   finp = FortranOpen(argv[1],"r");
   foutp = fopen(argv[2],"w");

   if (argc == 3  &&  (finp && foutp))
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL != ReadREACCSMolecule(finp,mp,""))
         return (EXIT_FAILURE);
      RecolorMolecule(mp);
      mph = LayoutMolecule(mp);
//      PutColorIntoValue(mph);
      PrintREACCSMolecule(foutp,mph,"");
      FreeMolecule(mp); FreeMolecule(mph);
      return (EXIT_SUCCESS);
   }
   else
   {
      fprintf(stderr,"usage: layout <input-file> <output-file>\n");
      return (EXIT_FAILURE);
   }
}
#endif
