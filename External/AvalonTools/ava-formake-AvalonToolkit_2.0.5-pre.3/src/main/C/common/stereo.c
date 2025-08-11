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
/*  File:     stereo.c                                                  */
/*                                                                      */
/*  Purpose:  Implements utilities needed for stereochemistry           */
/*            checking.                                                 */
/*                                                                      */
/*  History:  11-Sep-92       Creation of module as an extract from     */
/*                            the old struchk.c file.                   */
/*                                                                      */
/************************************************************************/

#include <math.h>
#include <string.h>

#include "local.h"
#include "mdl.h"
#include "utilities.h"

#include "stereo.h"

struct stereo_bond_t
   {
      double x, y;   /* relative 2D coordinates */
      int symbol;    /* stereo symbol */
      int number;    /* atom number of ligand atom */
      double angle;  /* angle in radiants rel. to first bond */
                     /* in array (counted counter clockwise) */
   };

double Volume(struct npoint_t tetra[4])
/*
 * Computes the signed volume of the tetrahedron defined by
 * the four points in tetra[].
 */
{
   double ax, ay, az, bx, by, bz, cx, cy, cz;

   ax = tetra[1].x - tetra[0].x;
   ay = tetra[1].y - tetra[0].y;
   az = tetra[1].z - tetra[0].z;
   bx = tetra[2].x - tetra[0].x;
   by = tetra[2].y - tetra[0].y;
   bz = tetra[2].z - tetra[0].z;
   cx = tetra[3].x - tetra[0].x;
   cy = tetra[3].y - tetra[0].y;
   cz = tetra[3].z - tetra[0].z;

   return (ax*(by*cz-bz*cy) +
           ay*(bz*cx-bx*cz) +
           az*(bx*cy-by*cx));
}

char *stereo_error;

#define PI     3.14159265359

#define ANGLE_EPSILON      (PI/180*5)      /* 5 degrees */
#define EPS	0.0000001

int Atom3Parity(struct stereo_bond_t ligands[3])
/*
 * Computes the stereo parity defined by three ligands.
 */
{
   int i, j, a, b;
   int reference;
   struct npoint_t tetrahedron[4], h;
   double angle;
   int maxnum;

   maxnum = ligands[0].number;
   for (i=1; i<3; i++)
      if (maxnum < ligands[i].number) maxnum = ligands[i].number;

   reference = (-1);
   for (i=0; i<3; i++)
      if (ligands[i].symbol != NONE)
         if (reference == (-1))
            reference = i;
         else
         {
            stereo_error = "three attachments with more than 2 stereobonds";
            return (ILLEGAL_REPRESENTATION);
         }

   if (reference == (-1)) return (UNDEFINED_PARITY);

   if (reference == 0)      {a = 1; b = 2;}
   else if (reference == 1) {a = 0; b = 2;}
   else                     {a = 0; b = 1;}

   angle = Angle(ligands[a].x, ligands[a].y, ligands[b].x, ligands[b].y);

   if (angle < ANGLE_EPSILON  || ABS(PI-angle) < ANGLE_EPSILON)
   {
      stereo_error = "three attachments: colinearity violation";
      return (ILLEGAL_REPRESENTATION);
   }

   tetrahedron[0].x = 0.0;
   tetrahedron[0].y = 0.0;
   tetrahedron[0].z = 0.0;
   tetrahedron[0].number = maxnum+1;
   for (i=0; i<3; i++)
   {
      tetrahedron[i+1].x = ligands[i].x;
      tetrahedron[i+1].y = ligands[i].y;
      if (ligands[i].symbol == UP)
         tetrahedron[i+1].z = 1.0;
      else if (ligands[i].symbol == DOWN)
         tetrahedron[i+1].z = -1.0;
      else if (ligands[i].symbol == NONE)
         tetrahedron[i+1].z = 0.0;
      else
      {
         stereo_error = "three attachments: illegal bond symbol";
         return (ILLEGAL_REPRESENTATION);
      }
      tetrahedron[i+1].number = ligands[i].number;
   }

   for (i=1; i<4; i++)
      for (j=i; j>0; j--)
    if (tetrahedron[j].number < tetrahedron[j-1].number)
    {
       h = tetrahedron[j];
       tetrahedron[j] = tetrahedron[j-1];
       tetrahedron[j-1] = h;
    }
    else
       break;

   if (Volume(tetrahedron) > 0.0) return (EVEN_PARITY);
   else                           return (ODD_PARITY);
}

int Atom4Parity(struct stereo_bond_t ligands[4])
/*
 * Computes the stereo parity defined by four ligands.
 * Assumes central atom at 0/0/0.
 */
{
   int i, j;
   struct npoint_t tetrahedron[4], h;
   int nup, ndown, nopposite;
   double angle;

   nup = ndown = 0;
   for (i=0; i<4; i++)
   {
      tetrahedron[i].x = ligands[i].x;
      tetrahedron[i].y = ligands[i].y;
      tetrahedron[i].z = 0.0;
      tetrahedron[i].number = ligands[i].number;
      if      (ligands[i].symbol == UP)
      {
         nup++;
         tetrahedron[i].z = 1.0;
      }
      else if (ligands[i].symbol == DOWN)
      {
         ndown++;
         tetrahedron[i].z = (-1.0);
      }
      else if (ligands[i].symbol != NONE)
      {
         stereo_error = "illegal bond symbol";
         return (ILLEGAL_REPRESENTATION);
      }
   }
   if (nup == 0  &&  ndown == 0)
   {
      return (UNDEFINED_PARITY);
   }

   if (nup > 2 || ndown > 2)
   {
      stereo_error = "too many stereobonds";
      return (ILLEGAL_REPRESENTATION);
   }

   if (nup + ndown == 1)    // check for 'umbrellas'
   {
       for (i=0; i<4; i++)
           if (ligands[i].symbol == UP  ||  ligands[i].symbol == DOWN) break;
       nopposite = 0;
       for (j=0; j<4; j++)
           if (i==j) continue;
           else
           {
               if (ligands[i].x*ligands[j].x + ligands[i].y*ligands[j].y < 0) nopposite++;
           }
       if (nopposite > 2)
       {
           stereo_error = "UMBRELLA: all non-stereo bonds opposite to single stereo bond";
           return (ILLEGAL_REPRESENTATION);
       }
   }

   for (i=0; i<2; i++)
      if (ligands[i].symbol == UP  &&  ligands[i+2].symbol == DOWN   ||
          ligands[i].symbol == DOWN  &&  ligands[i+2].symbol == UP)
      {
         stereo_error = "UP/DOWN opposition";
         return (ILLEGAL_REPRESENTATION);
      }

   for (i=0; i<4; i++)
      if (ligands[i].symbol == UP  &&  ligands[(i+1)%4].symbol == UP   ||
          ligands[i].symbol == DOWN  &&  ligands[(i+1)%4].symbol == DOWN)
      {
         stereo_error = "Adjacent like stereobonds";
         return (ILLEGAL_REPRESENTATION);
      }

   for (i=0; i<4; i++)
      if (ligands[i].symbol       == NONE  &&
          ligands[(i+1)%4].symbol == NONE  &&
          ligands[(i+2)%4].symbol == NONE)
      {
         angle = Angle(ligands[i].x - ligands[(i+1)%4].x,
                       ligands[i].y - ligands[(i+1)%4].y,
                       ligands[(i+2)%4].x - ligands[(i+1)%4].x,
                       ligands[(i+2)%4].y - ligands[(i+1)%4].y);
         if (angle < (185*PI/180))
         {
            stereo_error = "colinearity or triangle rule violation";
            return (ILLEGAL_REPRESENTATION);
         }
      }

   for (i=1; i<4; i++)
      for (j=i; j>0; j--)
         if (tetrahedron[j].number < tetrahedron[j-1].number)
         {
            h = tetrahedron[j];
            tetrahedron[j] = tetrahedron[j-1];
            tetrahedron[j-1] = h;
         }
         else
            break;

   if (Volume(tetrahedron) > 0.0) return (EVEN_PARITY);
   else                           return (ODD_PARITY);
}

int AtomParity(struct reaccs_molecule_t *mp,
               int iatom,
               neighbourhood_t *nbp)
/*
 * Computes the stereo parity of atom number iatom in *mp relative
 * to its numbering. The immediate neighbours are defined by *nbp
 * to speed up processing.
 *
 * Pro-memoriam copy. Should be deleted some time.
 */
{
   int i, j, jatom, ndb;
   int multiple, allene, stereo;

   struct stereo_bond_t stereo_ligands[4], h;

   if (nbp->n_ligands < 3  ||  nbp->n_ligands > 4)
   {
      return (ILLEGAL_REPRESENTATION);
   }

   multiple = FALSE; allene = FALSE; stereo = FALSE;
   for (i=0; i<nbp->n_ligands; i++)
   {
      if (mp->bond_array[nbp->bonds[i]].bond_type != SINGLE)
      {
          multiple = TRUE;
          // check if the multiple bond is part of an allene
          jatom = nbp->atoms[i]+1;
          ndb = 0;
          for (j=0; j<mp->n_bonds; j++)
              if (mp->bond_array[j].atoms[0] == jatom  ||  mp->bond_array[j].atoms[1] == jatom)
                  if (mp->bond_array[j].bond_type == DOUBLE) ndb++;
          if (ndb == 2) allene = TRUE;
      }
      stereo_ligands[i].x = mp->atom_array[nbp->atoms[i]].x- mp->atom_array[iatom-1].x;
      stereo_ligands[i].y = mp->atom_array[nbp->atoms[i]].y- mp->atom_array[iatom-1].y;
      stereo_ligands[i].number = nbp->atoms[i]+1;
      if (mp->bond_array[nbp->bonds[i]].atoms[0] == iatom)
      {
         stereo_ligands[i].symbol =
         mp->bond_array[nbp->bonds[i]].stereo_symbol;
         if (stereo_ligands[i].symbol == UP  ||
             stereo_ligands[i].symbol == DOWN)
            stereo = TRUE;
      }
      else
         stereo_ligands[i].symbol = NONE;
   }

   if (multiple &&
       stereo &&
       0 != strcmp(mp->atom_array[iatom-1].atom_symbol,"P")  &&
       0 != strcmp(mp->atom_array[iatom-1].atom_symbol,"S"))
   {
      if (allene)
          return (ALLENE_PARITY);
      else
      {
          stereo_error = "AtomParity: Stereobond at unsaturated atom";
          return (ILLEGAL_REPRESENTATION);
      }
   }
   else if (multiple && 0 != strcmp(mp->atom_array[iatom-1].atom_symbol,"S"))
      return (UNDEFINED_PARITY);

   stereo_ligands[0].angle = 0.0;   /* comp. angle rel. to first ligand */
   for (i=1; i<nbp->n_ligands; i++)
      stereo_ligands[i].angle =
         Angle(stereo_ligands[0].x, stereo_ligands[0].y,
               stereo_ligands[i].x, stereo_ligands[i].y);
   for (i=2; i<nbp->n_ligands; i++)    /* sort ligands */
      for (j=i; j>1; j--)
         if (stereo_ligands[j].angle < stereo_ligands[j-1].angle)
         {
            h = stereo_ligands[j];
            stereo_ligands[j] = stereo_ligands[j-1];
            stereo_ligands[j-1] = h;
         }
         else
            break;

   if (nbp->n_ligands == 3) return (Atom3Parity(stereo_ligands));
   else                     return (Atom4Parity(stereo_ligands));
}

int CheckStereo(struct reaccs_molecule_t *mp)
/*
 * Checks if all potential stereocenters are either completely undefined
 * or attributed with hashes and wedges according to MDL rules.
 */
{
   int parity;
   int center_defined;
   int result = TRUE;
   int i, j;

   neighbourhood_t *neighbour_array, *nbp;

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

   center_defined = FALSE;
   for (i=0, nbp=neighbour_array; i<mp->n_atoms; i++, nbp++)
      if (0 == strcmp(mp->atom_array[i].atom_symbol, "C")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4         ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "S")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "N")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "O")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "P")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "Si") &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)
      {
         parity = AtomParity(mp, i+1, nbp);
	 if (parity == ILLEGAL_REPRESENTATION)
	 {
            snprintf(msg_buffer, MAXMSG, "%10s    atom %3d : %s",
                    mp->name,i+1,stereo_error);
            AddMsgToList(msg_buffer);
	    result = FALSE;
         }
         else if (parity == EVEN_PARITY  || parity == ODD_PARITY  ||  parity == ALLENE_PARITY)
            center_defined = TRUE;
      }
      else
      {
         for (j=0; j<nbp->n_ligands; j++)
            if (mp->bond_array[nbp->bonds[j]].atoms[0] == i+1  &&
                (mp->bond_array[nbp->bonds[j]].stereo_symbol == UP  ||
                 mp->bond_array[nbp->bonds[j]].stereo_symbol == DOWN))
            {
	       snprintf(msg_buffer, MAXMSG,
                       "%10s    atom %3d : %s",
                       mp->name,
                       i+1,
                       "stereobond to non-stereogenic atom");
               AddMsgToList(msg_buffer);
               result = FALSE;
            }
      }

   if (mp->chiral_flag && !center_defined)
   {
      AddMsgToList("chiral flag set but no stereocenter defined");
      result = FALSE;
   }

   MyFree((char *)neighbour_array);
   return (result);
}

int FixDubious3DMolecule(struct reaccs_molecule_t *mp)
/*
 * Checks if the structure has 3D coordinates and/or flat sp3-carbons with stereo-bonds and
 * converts the designation to 2D, clearing away any Z-component of the coordinates.
 * Real 3D structures without stereo designations go through untouched.
 */
{
    int non_zero_z;
    int i, j, i1, i2, i3;
    int nstereo;
    neighbourhood_t *neighbour_array, *nbp;
    struct npoint_t tetra[4];
    int nflat_sp3;
    int stereo_triple;
    double length, vol;
    struct reaccs_bond_t *bp;
    int result = 0;

    // first check if this is a trivial case i.e. designated '2D' and all Z-coordinates == 0.0
    non_zero_z = FALSE;
    for (i=0; i<mp->n_atoms; i++)
        if (mp->atom_array[i].z != 0.0) non_zero_z = TRUE;
    if (!non_zero_z  &&  0 == strcmp(mp->dimensionality, "2D")) return (0);

    // count the number of stereo bonds
    nstereo = 0;
    for (i=0; i<mp->n_bonds; i++)
       if (mp->bond_array[i].stereo_symbol == UP  ||  mp->bond_array[i].stereo_symbol == DOWN)
           nstereo++;
    // compute average bond length to use in Volume significance testing
    length = 0.0;
    for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
        length += (mp->atom_array[bp->atoms[0]-1].x-mp->atom_array[bp->atoms[1]-1].x)*(mp->atom_array[bp->atoms[0]-1].x-mp->atom_array[bp->atoms[1]-1].x)
               +  (mp->atom_array[bp->atoms[0]-1].y-mp->atom_array[bp->atoms[1]-1].y)*(mp->atom_array[bp->atoms[0]-1].y-mp->atom_array[bp->atoms[1]-1].y);
    length /= mp->n_bonds;
    length = sqrt(length);

    // check if there is a flat sp3 carbon
    neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);
    SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);
    nflat_sp3 = 0;
    for (i=0, nbp=neighbour_array; i<mp->n_atoms; i++, nbp++)
    {
        if (nbp->n_ligands < 3) continue;
        if (0 != strcmp(mp->atom_array[i].atom_symbol, "C")  &&
            0 != strcmp(mp->atom_array[i].atom_symbol, "N")  &&
            0 != strcmp(mp->atom_array[i].atom_symbol, "P")  &&
            0 != strcmp(mp->atom_array[i].atom_symbol, "S"))
            continue;
        for (j=0; j<mp->n_bonds; j++)
            if ((mp->bond_array[j].bond_type == UP  ||  mp->bond_array[j].bond_type == DOWN)  &&
                i+1 == mp->bond_array[j].atoms[0])
                break;
        if (j < mp->n_bonds) continue;  // no stereo designation
        tetra[0].x = mp->atom_array[i].x; tetra[0].y = mp->atom_array[i].y; tetra[0].z = mp->atom_array[i].z;
        for (i1=0; i1<nbp->n_ligands; i1++)
            if (mp->bond_array[nbp->bonds[i1]].bond_type != SINGLE  &&
                0 != strcmp(mp->atom_array[i].atom_symbol, "P")     &&
                0 != strcmp(mp->atom_array[i].atom_symbol, "S")) break;
        if (i1 >= nbp->n_ligands) continue;     // multiple bond found => no sp3 carbon
        stereo_triple = 0;
        for (i1=0; i1<nbp->n_ligands; i1++)
        {
            tetra[1].x = mp->atom_array[nbp->atoms[i1]].x; tetra[1].y = mp->atom_array[nbp->atoms[i1]].y; tetra[1].z = mp->atom_array[nbp->atoms[i1]].z;
            if (mp->bond_array[nbp->bonds[i1]].bond_type == UP  ||  mp->bond_array[nbp->bonds[i1]].bond_type == DOWN) stereo_triple |= 1;
            for (i2=i1+1; i2<nbp->n_ligands; i2++)
            {
                tetra[2].x = mp->atom_array[nbp->atoms[i2]].x; tetra[2].y = mp->atom_array[nbp->atoms[i2]].y; tetra[2].z = mp->atom_array[nbp->atoms[i2]].z;
                if (mp->bond_array[nbp->bonds[i2]].bond_type == UP  ||  mp->bond_array[nbp->bonds[i2]].bond_type == DOWN) stereo_triple |= 2;
                for (i3=i2+1; i3<nbp->n_ligands; i3++)
                {
                    tetra[3].x = mp->atom_array[nbp->atoms[i3]].x; tetra[3].y = mp->atom_array[nbp->atoms[i3]].y; tetra[3].z = mp->atom_array[nbp->atoms[i3]].z;
                    if (mp->bond_array[nbp->bonds[i3]].bond_type == UP  ||  mp->bond_array[nbp->bonds[i3]].bond_type == DOWN) stereo_triple |= 4;
                    vol = Volume(tetra); if (vol < 0) vol = -vol;
                    if (!stereo_triple) continue;
                    if (vol < 0.01*length*length*length)
                    {
                        nflat_sp3++;
                        break;
                    }
                    stereo_triple &= ~4;
                }
                stereo_triple &= ~2;
                if (i3 >= nbp->n_ligands) break;
            }
            stereo_triple &= ~1;
            if (i2 >= nbp->n_ligands) break;
        }
    }

    if (non_zero_z  &&  0 == strcmp(mp->dimensionality, "2D"))
    {
        for (i=0; i<mp->n_atoms; i++)
            mp->atom_array[i].z = 0.0;
        result |= ZEROED_Z_COORDINATES;
        AddMsgToList("Cleared z-coordinates in 2D MOL file");
    }
    else if (non_zero_z  &&  nstereo > 0  &&  nflat_sp3 > 0)
    {
        strcpy(mp->dimensionality, "2D");
        result |= CONVERTED_TO_2D;
        for (i=0; i<mp->n_atoms; i++)
            mp->atom_array[i].z = 0.0;
        result |= ZEROED_Z_COORDINATES;
        AddMsgToList("Cleared z-coordinates of tilted 2D MOL file");
    }

    MyFree((char *)neighbour_array);
    return (result);
}

int DubiousStereochemistry(struct reaccs_molecule_t *mp)
/*
 * Checks if there is some ill-defined stereochemistry in the
 * molecule *mp. The function returns a bit set integer which defines
 * the problems encountered.
 */
{
   int i, j;
   int nmulti;
   int result = 0;
   struct reaccs_bond_t *bp;
   int is_allene, ndb, jatom, jj;

   neighbourhood_t *neighbour_array, *nbp;

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);

   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

                                             /* look for EITHER bonds */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->stereo_symbol == EITHER)
      {
	snprintf(msg_buffer, MAXMSG,
	        "%10s    EITHER bond %d-%d",
                mp->name,
                bp->atoms[0], bp->atoms[1]);
        AddMsgToList(msg_buffer);
        result |= EITHER_BOND_FOUND;
      }

                 /* look for stereo bonds to non-stereogenic atoms */
   for (i=0, nbp=neighbour_array; i<mp->n_atoms; i++, nbp++)
   {
      nmulti = 0;
      is_allene = FALSE;
      for (j=0; j<nbp->n_ligands; j++)
      {
          if (mp->bond_array[nbp->bonds[j]].bond_type != SINGLE)
          {
             nmulti++;
             if (mp->bond_array[nbp->bonds[j]].bond_type == DOUBLE)
             {
                 jatom = nbp->atoms[j];
                 ndb = 0;
                 for (jj=0; jj<neighbour_array[jatom].n_ligands; jj++)
                     if (mp->bond_array[neighbour_array[jatom].bonds[jj]].bond_type == DOUBLE) ndb++;
             }
             if (ndb == 2) is_allene = TRUE;
          }
      }
   
      if (!((0 == strcmp(mp->atom_array[i].atom_symbol, "C")  &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4  &&  nmulti == 0) ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "C")  &&
             nbp->n_ligands >= 2  &&  nbp->n_ligands <= 3  &&  nmulti == 1  &&  is_allene) ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "S")  &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)            ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "N")  &&
             nbp->n_ligands > 3  &&  nbp->n_ligands <= 4  &&  nmulti == 0) ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "P")  &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)                ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "Si") &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)  &&  nmulti == 0))
        for (j=0; j<nbp->n_ligands; j++)
           if (mp->bond_array[nbp->bonds[j]].atoms[0] == i+1  &&
               (mp->bond_array[nbp->bonds[j]].stereo_symbol == UP  ||
                mp->bond_array[nbp->bonds[j]].stereo_symbol == DOWN))
           {
// fprintf(stderr, "atom %d: nmulti = %d, is_allene = %d\n", i+1, nmulti, is_allene);
              snprintf(msg_buffer, MAXMSG,
                      "%10s    atom %3d : %s",
                      mp->name,
                      i+1,
                      "stereobond to non-stereogenic atom");
              AddMsgToList(msg_buffer);
              result |= STEREO_BOND_AT_NON_STEREO_ATOM;
           }
   }

   MyFree((char *)neighbour_array);
   return (result);
}

void RemoveDubiousStereochemistry(struct reaccs_molecule_t *mp)
/*
 * Removes ill-defined stereodescriptors.
 */
{
   int i, j;
   int nmulti;
   struct reaccs_bond_t *bp;
   neighbourhood_t *neighbour_array, *nbp;

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);

   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

                       /* remove EITHER marks */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->stereo_symbol == EITHER) bp->stereo_symbol = NONE;

                   /* remove stereo marks to non-stereogenic atoms */
   for (i=0, nbp=neighbour_array; i<mp->n_atoms; i++, nbp++)
   {
      nmulti = 0;
      for (j=0; j<nbp->n_ligands; j++)
         if (mp->bond_array[nbp->bonds[j]].bond_type != SINGLE)
            nmulti++;
   
      if (!((0 == strcmp(mp->atom_array[i].atom_symbol, "C")  &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4  &&  nmulti == 0) ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "S")  &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)            ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "N")  &&
             nbp->n_ligands > 3  &&  nbp->n_ligands <= 4  &&  nmulti == 0) ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "P")  &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)                ||
            (0 == strcmp(mp->atom_array[i].atom_symbol, "Si") &&
             nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)  &&  nmulti == 0))
        for (j=0; j<nbp->n_ligands; j++)
           if (mp->bond_array[nbp->bonds[j]].atoms[0] == i+1  &&
               (mp->bond_array[nbp->bonds[j]].stereo_symbol == UP  ||
                mp->bond_array[nbp->bonds[j]].stereo_symbol == DOWN))
              mp->bond_array[nbp->bonds[j]].stereo_symbol = NONE;
   }
   MyFree((char *)neighbour_array);
}

static double clash_limit = 0.15;

void SetCollisionLimit(int percent)
/*
 * Sets the limit for atom/atom and atom/bond collision in percent
 * of the average bond length.
 */
{
   clash_limit = percent/100.0;
}

int AtomClash(struct reaccs_molecule_t *mp)
/*
 * Checks if any two atoms in *mp come closer than 10% of the
 * average bond length or if an atom is too close the line
 * between two bonded atoms.
 */
{
   int i, j;
   struct reaccs_bond_t *bp;
   struct reaccs_atom_t *ap, *ap1, *ap2;
   double bond_square_median, dist, min_dist;
   double rr, rb, bb;
   double *blengths, h;

                  /* compute median of square of bond lenght (quick/dirty) */
   if (mp->n_bonds <= 0) return (FALSE);
   blengths = TypeAlloc(mp->n_bonds, double);
   blengths[0] = 1.0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      ap1 = &mp->atom_array[bp->atoms[0]-1];
      ap2 = &mp->atom_array[bp->atoms[1]-1];
      blengths[i] = (ap1->x - ap2->x)*(ap1->x - ap2->x) +
                    (ap1->y - ap2->y)*(ap1->y - ap2->y);
   }
   for (i=1; i<mp->n_bonds; i++)
      for (j=i-1; j>=0; j--)
         if (blengths[j] > blengths[j+1])
         {
            h = blengths[j]; blengths[j] = blengths[j+1]; blengths[j+1] = h;
         }
	 else
	    break;
   bond_square_median = blengths[mp->n_bonds/2];
   MyFree((char *)blengths);

                  /* Check if two atoms get too close to each other */
   min_dist = bond_square_median;
   for (i=0, ap1=mp->atom_array; i<mp->n_atoms; i++, ap1++)
      for (j=i+1, ap2 = ap1+1; j<mp->n_atoms; j++, ap2++)
      {
         dist = (ap1->x - ap2->x)*(ap1->x - ap2->x) +
                (ap1->y - ap2->y)*(ap1->y - ap2->y);
         if (dist < clash_limit*clash_limit*bond_square_median)
         {
            snprintf(msg_buffer, MAXMSG,
                    "%10s    atom %3d : only %.2g%% of avg. bond away from atom %d",
                    mp->name, j+1,
                    100*sqrt(dist/(bond_square_median+EPS)), i+1);
	    AddMsgToList(msg_buffer);
            return (TRUE);
         }
         if (dist < min_dist) min_dist = dist;
      }

                   /* check if atom lies on top of some bond */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      for (j=0, ap=mp->atom_array; j<mp->n_atoms; j++, ap++)
      {
         if (bp->atoms[0]-1 == j  ||  bp->atoms[1]-1 == j) continue;

         ap1 = &mp->atom_array[bp->atoms[0]-1];
         ap2 = &mp->atom_array[bp->atoms[1]-1];
         rr = (ap->x - ap1->x)*(ap->x - ap1->x) +
              (ap->y - ap1->y)*(ap->y - ap1->y);
         bb = (ap2->x - ap1->x)*(ap2->x - ap1->x) +
              (ap2->y - ap1->y)*(ap2->y - ap1->y);
         rb = (ap->x - ap1->x)*(ap2->x - ap1->x) +
              (ap->y - ap1->y)*(ap2->y - ap1->y);
         if (0 <= rb  &&   /* cos alpha > 0 */
             rb <= bb  &&  /* projection of r onto b does not exceed b */
             (rr*bb - rb*rb)/(bb+EPS) <       /* distance from bond < limit */
                clash_limit*clash_limit*bond_square_median)
         {
            snprintf(msg_buffer, MAXMSG,
                    "%10s    atom %3d : only %.2g%% %s from bond %d-%d",
                    mp->name,
                    j+1,
                    100*sqrt((rr*bb-rb*rb)/(bb*bond_square_median+EPS)),
                    "of average bond length",
                    bp->atoms[0], bp->atoms[1]);
            AddMsgToList(msg_buffer);
            return (TRUE);
         }
      }

   return (FALSE);     /* no clash */
}

int CisTransPerception(struct reaccs_molecule_t *mp,
                       int                       numbering[])
/*
 * Sets the color field of the defined double bonds in *mp to CIS,
 * TRANS, or NONE depending on the ligands with the lowest numbering[].
 * It returns the number of defined double bonds found.
 */
{
   neighbourhood_t *nba;
   int i, j1, j2, k, nmin, maxnum;
   int equal;
   int result;
   struct reaccs_bond_t *bp;
   int at1, at2;
   double x21, y21;
   double x23, y23;
   double x32, y32;
   double x34, y34;
   double sign;

   nba = TypeAlloc(mp->n_atoms, neighbourhood_t);

   SetupNeighbourhood(mp,nba,mp->n_atoms);

   result = 0;
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      bp->color = NONE;

   maxnum = (-1);
   for (i=0; i<mp->n_atoms; i++)
      if (numbering[i] > maxnum) maxnum = numbering[i];
   maxnum++;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->bond_type == DOUBLE && bp->stereo_symbol != CIS_TRANS_EITHER)
      {
         j1 = bp->atoms[0]-1; j2 = bp->atoms[1]-1;
         if (!AtomSymbolMatch(mp->atom_array[j1].atom_symbol,"C,N"))
            continue;
         if (!AtomSymbolMatch(mp->atom_array[j2].atom_symbol,"C,N"))
            continue;
         if (nba[j1].n_ligands <= 1  || /* no subst. */
             nba[j2].n_ligands <= 1) continue;
         if (nba[j1].n_ligands > 3  ||  /* valence error in *mp */
             nba[j2].n_ligands > 3) continue;

         equal = FALSE;     /* find lowest numbered neighbour of j1 */
         for (k=0, nmin=maxnum; k<nba[j1].n_ligands; k++)
            if (nba[j1].atoms[k] != j2) /* no loop back */
            {
               if (numbering[nba[j1].atoms[k]] < nmin)
               {
                  at1 = nba[j1].atoms[k];
                  nmin = numbering[at1];
               }
               else if (numbering[nba[j1].atoms[k]] == nmin)
                  equal = TRUE;
            }
	 if (equal)     /* identical substituents */
            continue;   /* no stereochemistry */

         equal = FALSE;     /* find lowest numbered neighbour of j1 */
         for (k=0, nmin=maxnum; k<nba[j2].n_ligands; k++)
            if (nba[j2].atoms[k] != j1) /* no loop back */
            {
               if (numbering[nba[j2].atoms[k]] < nmin)
               {
                  at2 = nba[j2].atoms[k];
                  nmin = numbering[at2];
               }
               else if (numbering[nba[j2].atoms[k]] == nmin)
                  equal = TRUE;
            }
	 if (equal)     /* identical substituents */
	    continue;   /* no stereochemistry */

	 /* Now, bp points to a double bond, at1 and at2 are the */
	 /* indices (not numbers) of the atoms with lowest numbering */
	 /* attached to each end of the bond, and the bond is */
	 /* guaranteed to be either CIS or TRANS */

	 x21 = mp->atom_array[at1].x - mp->atom_array[bp->atoms[0]-1].x;
	 y21 = mp->atom_array[at1].y - mp->atom_array[bp->atoms[0]-1].y;
	 x23 = mp->atom_array[bp->atoms[1]-1].x -
	 mp->atom_array[bp->atoms[0]-1].x;
	 y23 = mp->atom_array[bp->atoms[1]-1].y -
	       mp->atom_array[bp->atoms[0]-1].y;
	 x32 = (-x23); y32 = (-y23);
	 x34 = mp->atom_array[at2].x - mp->atom_array[bp->atoms[1]-1].x;
	 y34 = mp->atom_array[at2].y - mp->atom_array[bp->atoms[1]-1].y;
	 sign = (x21*y23-x23*y21)*(x32*y34-x34*y32);
	 if (-0.001 < sign  &&  sign < 0.001) continue;
	 result++;
	 if (sign > 0.0) bp->color = CIS;
	 else            bp->color = TRANS;
      }

   MyFree((char *)nba);
   return (result);
}

int NoParityDefined(struct reaccs_molecule_t *mp)
/*
 * Returns TRUE if the molecule *mp has no defined stereocenter.
 */
{
   neighbourhood_t *neighbour_array, *nbp;
   int parity;
   int i;

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

   for (i=0, nbp=neighbour_array; i<mp->n_atoms; i++, nbp++)
      if (0 == strcmp(mp->atom_array[i].atom_symbol, "C")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4         ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "S")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "N")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "O")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "P")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "Si") &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)
      {
         parity = AtomParity(mp, i+1, nbp);
         if (parity == EVEN_PARITY  ||  parity == ODD_PARITY)
	 {
	    MyFree((char *)neighbour_array);
	    return (TRUE);
	 }
      }

   MyFree((char *)neighbour_array);
   return (FALSE);
}

int AllCentersRefined(struct reaccs_molecule_t *mp,
                      int                       numbering[])
/*
 * Checks if all defined stereocenters in *mp are unambiguously
 * numbered by numbering, i.e. all four ligands get different
 * numbers. It returns TRUE if this is the case and FALSE otherwise.
 */
{
   neighbourhood_t *neighbour_array, *nbp;
   int parity;
   int i, j, k, h;
   int refnum[4], nref;

   neighbour_array = TypeAlloc(mp->n_atoms, neighbourhood_t);
   SetupNeighbourhood(mp,neighbour_array,mp->n_atoms);

   for (i=0, nbp=neighbour_array; i<mp->n_atoms; i++, nbp++)
      if (0 == strcmp(mp->atom_array[i].atom_symbol, "C")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4         ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "S")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "N")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "O")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "P")  &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4            ||
          0 == strcmp(mp->atom_array[i].atom_symbol, "Si") &&
          nbp->n_ligands > 2  &&  nbp->n_ligands <= 4)
      {
         parity = AtomParity(mp, i+1, nbp);
         if (parity == EVEN_PARITY  ||  parity == ODD_PARITY)
         {
                         /* extraxt numbers of stereoligands */
            nref = nbp->n_ligands;
            for (j=0; j<nref; j++)
               refnum[j] = numbering[nbp->atoms[j]];
                      /* sort ligands */
            for (j=nref; j>0; j--)
               for (k=j; k>0; k--)
                  if (refnum[k] < refnum[k-1])
                  {
		     h = refnum[k];
		     refnum[k] = refnum[k-1];
		     refnum[k-1] = h;
		  }
                       /* check if there is a duplicate */
	     for (j=1; j<nref; j++)
		if (refnum[j] == refnum[j-1])
	        {
                   snprintf(msg_buffer, MAXMSG, "%10s    unrefined parity at atom %d",
                           mp->name,i+1);
                   AddMsgToList(msg_buffer);
		   return (FALSE);
                }
         }
      }

   MyFree((char *)neighbour_array);
   return (TRUE);
}

int PermutationParity(int *indices, int size)
/**
 * Computes the parity of the list of distinct numbers indices[0..size-1], i.e. if the number of exchanges to sort the numbers is
 * even (EVEN_PARITY) or odd (ODD_PARITY). The function destroys the input order!
 */
{
    int nswap, i, j, h;

    nswap = 0;
    for (i=1; i<size; i++)
        for (j=i-1; j>=0; j--)
            if (indices[j] > indices[j+1])
            {
                nswap++;
                h = indices[j]; indices[j] = indices[j+1]; indices[j+1] = h;
            }
            else
                break;
    
    return 0 == (nswap%2) ? EVEN_PARITY : ODD_PARITY;
}

int IsStereoMatch(struct reaccs_molecule_t *mp, neighbourhood_t *nbp,
                  struct reaccs_molecule_t *qp, neighbourhood_t *qnbp,
                  ssmatch_t *matches, int invert_stereo)
/**
 * Checks if the match matches into mp obeys the center stereochemistry of qp.
 * Uses the mirror image of *qp if invert_stereo == TRUE.
 */
{
    int i, j, is_match;
    int indices[4], qindices[4], size;
    int qperm, mperm;
    int test_parity;

    if (IsNULL(matches)) return FALSE;
    is_match = TRUE;
    for (i=0; i<qp->n_atoms; i++)
        if (qp->atom_array[i].stereo_parity != NONE)    // only stereo query atom can cause match to fail
        {
            for (j=0; j<qnbp[i].n_ligands; j++)
            {
                indices[j] = matches->match_atoms[qnbp[i].atoms[j]];
                qindices[j] = qnbp[i].atoms[j];
            }
            qperm = PermutationParity(qindices, qnbp[i].n_ligands);
            mperm = PermutationParity(indices, qnbp[i].n_ligands);
            test_parity = qp->atom_array[i].stereo_parity;
            if (invert_stereo) test_parity = INVERT_PARITY(test_parity);
            if (qperm != mperm) test_parity = INVERT_PARITY(test_parity);
            if (AtomParity(mp, matches->match_atoms[i]+1, nbp+matches->match_atoms[i]) !=  test_parity)
            {
                is_match = FALSE;
                break;
            }
// else
// {
// fprintf(stderr, "stereo match found between query atom %d with %d ligands and target atom %d with %d ligands\n", i+1, qnbp[i].n_ligands, matches->match_atoms[i]+1, nbp[matches->match_atoms[i]].n_ligands);
// }
        if (FALSE &&  qp->atom_array[i].stereo_parity == EVEN_PARITY  &&  !invert_stereo)
        {
            /*
            for (j=0; j<qnbp[i].n_ligands; j++)
                indices[j] = matches->match_atoms[qnbp[i].atoms[j]];
            fprintf(stderr, "testing query atom %d with EVEN_PARITY(%d) mapped to atom %d with parity %d and sort-parity %d\n",
                    i+1, qp->atom_array[i].stereo_parity,
                    matches->match_atoms[i]+1, AtomParity(mp, matches->match_atoms[i]+1, nbp+matches->match_atoms[i]),
                                               PermutationParity(indices, qnbp[i].n_ligands));
            */
            if (AtomParity(mp, matches->match_atoms[i]+1, nbp+matches->match_atoms[i]) !=  PermutationParity(indices, qnbp[i].n_ligands))
            {
                is_match = FALSE;
                break;
            }
        }
        else if (FALSE && qp->atom_array[i].stereo_parity == ODD_PARITY   &&   invert_stereo)
        {
            /*
            for (j=0; j<qnbp[i].n_ligands; j++)
                indices[j] = matches->match_atoms[qnbp[i].atoms[j]];
            fprintf(stderr, "testing query atom %d with EVEN_PARITY(%d) mapped to atom %d with parity %d and sort-parity %d\n",
                    i+1, qp->atom_array[i].stereo_parity,
                    matches->match_atoms[i]+1, AtomParity(mp, matches->match_atoms[i]+1, nbp+matches->match_atoms[i]),
                                               PermutationParity(indices, qnbp[i].n_ligands));
            */
            for (j=0; j<qnbp[i].n_ligands; j++)
                indices[j] = matches->match_atoms[qnbp[i].atoms[j]];
            if (AtomParity(mp, matches->match_atoms[i]+1, nbp+matches->match_atoms[i]) !=  PermutationParity(indices, qnbp[i].n_ligands))
            {
                is_match = FALSE;
                break;
            }
        }
        else if (FALSE && qp->atom_array[i].stereo_parity == EVEN_PARITY  &&   invert_stereo)
        {
            for (j=0; j<qnbp[i].n_ligands; j++)
                indices[j] = matches->match_atoms[qnbp[i].atoms[j]];
            if (AtomParity(mp, matches->match_atoms[i]+1, nbp+matches->match_atoms[i]) ==  PermutationParity(indices, qnbp[i].n_ligands))
            {
                is_match = FALSE;
                break;
            }
        }
        else if (FALSE && qp->atom_array[i].stereo_parity == ODD_PARITY   &&  !invert_stereo)
        {
            for (j=0; j<qnbp[i].n_ligands; j++)
                indices[j] = matches->match_atoms[qnbp[i].atoms[j]];
            if (AtomParity(mp, matches->match_atoms[i]+1, nbp+matches->match_atoms[i]) ==  PermutationParity(indices, qnbp[i].n_ligands))
            {
                is_match = FALSE;
                break;
            }
        }
        }

    return is_match;
}
