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

/****************************************************************************/
/*                                                                          */
/*   File:   	didepict.c                                                  */
/*		                                                            */
/*   Purpose:	Generate device independent graphics commands from an MDL   */
/*		data structure. These commands are primaryly intended to be */
/*		intepreted by a Java class.                                 */
/*                                                                          */
//   Histroy:	1996-06-21	Start of development from the Windows-DLL   */
//				version depict.c.                           */
//                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <string.h>

   /* assuming scale mode == TWIPS => 1/1440 inch per logical unit */
   /* Std. bond = 0.35, one point = 1/72 inch, font = 10 points */
static int stdbnd   = (int)(0.35*1440);
static int  fontsize = 10*1440/72;

#define HELVETICA	"Helvetica"
#define TIMES		"TimesRoman"
#define COURIER		"Courier"
#define DIALOG		"Dialog"

static void ChooseFont(char *buffer, double scale, char *fontname)
/*
 * Select the font for depiction of atom symbols.
 */
{
   sprintf(buffer,"f %s %d\n", fontname, (int)(scale*(fontsize+0.5)));
}

#include "local.h"

#include "utilities.h"

#include "forio.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "perceive.h"
#include "layout.h"
#include "smi2mol.h"
#include "didepict.h"
#include "rtutils.h"

#define SNG           1
#define BND_RIGHT     2
#define BND_LEFT      4
#define DBL_RIGHT     (1|2)
#define DBL_LEFT      (1|4)
#define DBL           DBL_RIGHT
#define TRP           (1|2|4)

#define SYM_DBL       8
#define WDG           16
#define HSH           32
#define UNK           64
#define DSH           128

#define CUT          1
#define CUT_SEC      2

static
void PlotBond(char *buffer,
              int x1, int y1, int x2, int y2,
              int d, int type, int flags1, int flags2)
/*
 * Prints the commands to draw a bond from (x1,y1) to (x2,y2) to *fp.
 * It uses the font size d as the separation characteristic, the symbol type,
 * and the shortening flags flags1 and flags2 to determine the appearance.
 */
{
   int i, len;
   int dbx, dby;       /* multiple bond spacing */
   double bf = 0.45;
   int dwx, dwy;       /* wedge width */
   double wf = 0.35;
   double dhx, dhy;    /* hash line distance */
   int nhash;
   double hf = 0.40;
   int dcx, dcy;       /* line end cutting */
   double cf = 0.55;
   int ddx, ddy;       /* multiple bond cutting */
   double df = 0.45;
   int wedge[3][2];

   char buf[100];	/* sprintf buffer */

   buffer[0] = '\0';

   len = (int)sqrt((double)(x1-x2)*(double)(x1-x2) +
                   (double)(y1-y2)*(double)(y1-y2));
   dbx = bf*((double)d*(x2-x1))/len; dby = bf*((double)d*(y2-y1))/len;
   dwx = wf*((double)d*(x2-x1))/len; dwy = wf*((double)d*(y2-y1))/len;
   dcx = cf*((double)d*(x2-x1))/len; dcy = cf*((double)d*(y2-y1))/len;
   ddx = df*((double)d*(x2-x1))/len; ddy = df*((double)d*(y2-y1))/len;

   if (flags1 & CUT) { x1 += dcx; y1 += dcy; }
   if (flags2 & CUT) { x2 -= dcx; y2 -= dcy; }

   if (type & SNG)
   {
      sprintf(buf, "l %d %d %d %d\n", x1, y1, x2, y2);
      strcat(buffer,buf);
   }

   if (type & BND_RIGHT)
   {
      if (flags1 & CUT_SEC)
         sprintf(buf, "l %d %d ", x1-dby+ddx, y1+dbx+ddy);
      else
         sprintf(buf, "l %d %d ", x1-dby,     y1+dbx);
      strcat(buffer,buf);
      if (flags2 & CUT_SEC)
         sprintf(buf, "%d %d\n", x2-dby-ddx, y2+dbx-ddy);
      else
         sprintf(buf, "%d %d\n", x2-dby,     y2+dbx);
      strcat(buffer,buf);
   }

   if (type & BND_LEFT)
   {
      if (flags1 & CUT_SEC)
         sprintf(buf, "l %d %d ", x1+dby+ddx, y1-dbx+ddy);
      else
         sprintf(buf, "l %d %d ", x1+dby,     y1-dbx);
      strcat(buffer,buf);
      if (flags2 & CUT_SEC)
         sprintf(buf, "%d %d\n", x2+dby-ddx, y2-dbx-ddy);
      else
         sprintf(buf, "%d %d\n", x2+dby,     y2-dbx);
      strcat(buffer,buf);
   }

   if (type & SYM_DBL)
   {
      sprintf(buf, "l %d %d ", x1-dby/2, y1+dbx/2);
      strcat(buffer,buf);
      sprintf(buf, "%d %d\n", x2-dby/2, y2+dbx/2);
      strcat(buffer,buf);
      sprintf(buf, "l %d %d ", x1+dby/2, y1-dbx/2);
      strcat(buffer,buf);
      sprintf(buf, "%d %d\n", x2+dby/2, y2-dbx/2);
      strcat(buffer,buf);
   }

   if (type & WDG)
   {
      wedge[0][0] = x1; wedge[0][1] = y1;
      wedge[1][0] = x2-dwy/2; wedge[1][1] = y2+dwx/2;
      wedge[2][0] = x2+dwy/2; wedge[2][1] = y2-dwx/2;
      sprintf(buf, "pf 3 %d %d %d %d %d %d\n",
              wedge[0][0], wedge[0][1],
              wedge[1][0], wedge[1][1],
              wedge[2][0], wedge[2][1]);
      strcat(buffer,buf);
   }

   if (type & HSH)
   {
      nhash = (int)(len/(hf*d)+0.5);
      if (nhash < 2) nhash = 2;
      dhx = wf*((double)d*(x2-x1))/len; dhy = wf*((double)d*(y2-y1))/len;
      for (i=0; i<=nhash; i++)
      {
         sprintf(buf, "l %d %d ", (int)(x1+i*(x2-x1-dhy/2)/nhash),
                                  (int)(y1+i*(y2-y1+dhx/2)/nhash));
	 strcat(buffer,buf);
         sprintf(buf, "%d %d\n", (int)(x1+i*(x2-x1+dhy/2)/nhash),
                                 (int)(y1+i*(y2-y1-dhx/2)/nhash));
	 strcat(buffer,buf);
      }
   }

   if (type & DSH)
   {
      x1 += -dby+ddx; y1 -= -dbx-ddy;
      x2 += -dby-ddx; y2 -= -dbx+ddy;
      sprintf(buf, "l %d %d %d %d", x1, y1, x2, y2);
      strcat(buffer,buf);
   }
}

static void GetWorldWindow(struct reaccs_molecule_t *mp,
                           float *xminp, float *xmaxp,
                           float *yminp, float *ymaxp,
                           float *lenp)
/*
 * Computes the drawing window of the molecule *mp in world coordinates.
 */
{
   struct reaccs_atom_t *ap;
   struct reaccs_bond_t *bp;
   int i;
   double d1, d2;

   *xminp = *yminp = 1e10; *xmaxp = *ymaxp = -1e10;      /* World Window */
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
   {
      if (*xminp > ap->x) *xminp = ap->x;
      if (*xmaxp < ap->x) *xmaxp = ap->x;
      if (*yminp > ap->y) *yminp = ap->y;
      if (*ymaxp < ap->y) *ymaxp = ap->y;
   }
   if (*xminp == *xmaxp) {*xminp -= 0.5; *xmaxp += 0.5;}
   if (*yminp == *ymaxp) {*yminp -= 0.5; *ymaxp += 0.5;}

   *lenp = 0.0;
   for (i=0, bp= mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      d1 = mp->atom_array[bp->atoms[0]-1].x -
           mp->atom_array[bp->atoms[1]-1].x;
      d2 = mp->atom_array[bp->atoms[0]-1].y -
           mp->atom_array[bp->atoms[1]-1].y;
      *lenp += d1*d1 + d2*d2;
   }
   if (*lenp == 0.0) *lenp = 1.5;
   else              *lenp = sqrt(*lenp/mp->n_bonds);
}

static double xoffset, yoffset;
static double scale = 1.0, fontscale = 1.0;

static void SetTransformation(int metalen,
                              double xmin,  double xmax,
                              double ymin,  double ymax,
                              double len)
{
   scale = metalen/len;
   fontscale = 1.0;
   xoffset = (xmin+xmax)/2;
   yoffset = (ymin+ymax)/2;
   if ((xmax-xmin)*scale  > 16000)    /* prevent integer overflow */
      scale = 16000/(xmax-xmin);
   if ((ymax-ymin)*scale  > 16000)        /* prevent integer overflow */
      scale = 16000/(ymax-ymin);
}

static void FixScale(double xext_best, double yext_best,
                     double xext_new,  double yext_new)
/*
 * Sets the global scale varaible such that the picture will now
 * fit into the box defined by xext_new and yext_new.
 */
{
   if (xext_best < xext_new  &&  yext_best < yext_new)
   {
      return;     /* best size would already fix into new box => NOOP */
   }

   if (xext_new/xext_best <= yext_new/yext_best)
   {
      /* xscale defines scale */
      scale *= xext_new/xext_best;
      fontscale *= xext_new/xext_best;
   }
   else
   {
      /* yscale defines scale */
      scale *= yext_new/yext_best;
      fontscale *= yext_new/yext_best;
   }
   if (fontscale < 0.3) fontscale = 0.3;
}

static int WorldToMetaX(double x)
{
   return ((int)((x-xoffset)*scale));
}

static int WorldToMetaY(double y)
{
   return (-(int)((y-yoffset)*scale));
}

static neighbourhood_t nba[MAXATOMS];

static void SetDrawFlags(struct reaccs_molecule_t *mp)
/*
 * This procedure looks for each atom and bond how they should be drawn.
 * The information is stored in the color fields of atoms and bonds.
 */
{
   int i;
   struct reaccs_atom_t *ap, *ap1, *ap2;
   struct reaccs_bond_t *bp;
   int at1, at2;

   MakeRingsClockwise(mp);

   ResetColors(mp);
   SetupNeighbourhood(mp,nba,mp->n_atoms);
   PerceiveRingBonds(mp);

   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (0 != strcmp(ap->atom_symbol, "C")  ||
          ap->charge  != NONE                ||
          ap->radical != NONE                ||
          ap->query_H_count != NONE          ||
	  nba[i].n_ligands == 0		     ||
          ap->mass_difference != NONE)
          ap->color |= CUT;
      else if (nba[i].n_ligands != 1)
          ap->color |= CUT_SEC;

   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      bp->color = NONE;
      switch (bp->bond_type)
      {
         case SINGLE: if (bp->stereo_symbol == UP)
                         bp->color |= WDG;
                      else if (bp->stereo_symbol == DOWN)
                         bp->color |= HSH;
                      else
                         bp->color |= SNG;
                      break;
         case DOUBLE: at1 = bp->atoms[0]; at2 = bp->atoms[1];
                      ap1 = &mp->atom_array[at1-1];
                      ap2 = &mp->atom_array[at2-1];
                      if (bp->topography != RING)
                      {         /* dbl bond should have common direction */
                         if ((ap2->x-ap1->x)*2 + (ap2->y-ap1->y)*3  < 0)
                         {
                            bp->atoms[0] = at2; bp->atoms[1] = at1;
                         }
                      }
                      if (nba[bp->atoms[0]-1].n_ligands < 2   ||
                          nba[bp->atoms[1]-1].n_ligands < 2   ||
                          (nba[bp->atoms[0]-1].n_ligands > 3 &&
			   nba[bp->atoms[1]-1].n_ligands > 3 &&
                           bp->topography != RING)            ||
                          (nba[bp->atoms[0]-1].n_ligands == 2 &&
                           nba[bp->atoms[1]-1].n_ligands == 2 &&
                           0 != strcmp(ap1->atom_symbol,"C")  &&
                           0 != strcmp(ap2->atom_symbol,"C")  &&
                           bp->topography != RING))
                         bp->color |= SYM_DBL;
		      else
                         bp->color |= DBL;
                      break;
         case TRIPLE: bp->color |= TRP;
                      at1 = bp->atoms[0]; at2 = bp->atoms[1];
                      ap1 = &mp->atom_array[at1-1];
                      ap2 = &mp->atom_array[at2-1];
                      ap1->color &= ~CUT_SEC; ap2->color &= ~CUT_SEC;
                      break;
         case AROMATIC: bp->color |= DSH;
                        bp->color |= SNG;
                        break;
         default: bp->color |= SNG;
                  break;
      }
   }
}

static void AttributesToString(char *buffer,
                               int charge,
                               int radical,
                               int mass_difference)
{
   if (mass_difference != 0)
   {
      buffer += sprintf(buffer,"*");
   }

   if (abs(charge) > 1)
   {
      buffer += sprintf(buffer,"%d", abs(charge));
   }

   if (charge > 0)
   {
      buffer[0] = '+'; buffer++;
   }
   else if (charge < 0)
   {
      buffer[0] = '-'; buffer++;
   }

   if (radical == SINGLET)
   {
      buffer[0] = '|'; buffer++;
   }
   else if (radical == DOUBLET)
   {
      buffer[0] = '.'; buffer++;
   }
   else if (radical == TRIPLET)
   {
      buffer[0] = ':'; buffer++;
   }

   buffer[0] = '\0';
}


#define NCOLORS               17
static long colortable[NCOLORS] =
   {
      0x000000,        /* 0 black */
      0xFF0000,        /* 1 red */
      0x0000FF,        /* 2 blue */
      0x007F00,        /* 3 dark green */
      0x7F3F00,        /* 4 brown */
      0x3F007F,        /* 5 dark violet */
      0x003F3F,        /* 6 dark cyan */
      0x7F7F00,        /* 7 dark yellow */
      0x7F007F,        /* 8 dark pink */
      0x7F7F7F,        /* 9 grey */
      0x1CC6E7,        /* 10 light blue */
      0x0063FF,        /* 11 middle blue */
      0x0000B5,        /* 12 dark blue */
      0x000000,        /* 13 black */
      0xFF00FF,        /* 14 pink */
      0xFF0000,        /* 15 red */
      0xFF8C42,        /* 16 orange */
   };

#define MAXWRITEBUFFER 2000

static
void BufferedWrite(FILE *fp, char buffer[], char app[])
{
   if (strlen(buffer) + strlen(app) + 1 > MAXWRITEBUFFER)
   {					/* flush buffer */
      fprintf(fp, "%s%s", buffer, app);
      buffer[0] = '\0';
   }
   else
   {					/* append app */
      strcat(buffer, app);
   }
}

void MoleculeToPlotFile(FILE *fp,
			struct reaccs_molecule_t *mp,
			int xext, int yext,
                        int flags,
                        int colors[],
			char *labels[])
/*
 * Chemical structure depiction workhorse. Takes the MDL formatted
 * chemical structure *mp and writes device independent drawing
 * commands to *fp.
 * 
 * The flags parameter currently tells if color is to be used.
 *
 * The array colors[] is used to color the different parts of the molecule.
 * colors[i] == 0 means draw in black. 1 <= colors[i] <= 9 means use one of
 * the nine most dark colors to get reasonable readablility. If colors ==
 * NULL, the standard atom label coloring with black bonds is used.
 *
 * xext and yext are used for the desired x- and y- dimensions of the resulting
 * image in twips.
 *
 * The array labels[] contains alternate atom labels to be used instead
 * of the ones defined in the atom_symbol fields.
 * Labels which are "\0" strings are not considered. Use " " to get
 * a blank label.
 */
{
   int current_color;
   int bcolor1;
   int bcolor2;

   struct reaccs_bond_t *bp;
   struct reaccs_atom_t *ap;

   float xmin, xmax, ymin, ymax, len;
   int i;
   int color_changed;
   int xExt, yExt;
   int x1, y1, x2, y2;
   int flags1, flags2;
   char chargestring[10];
   char atom_symbol[10], *symbolp;
   int charge;
   int *H_count;

   char buf[MAXWRITEBUFFER], buffer[MAXWRITEBUFFER];

   if (!mp) return;

   GetWorldWindow(mp, &xmin, &xmax, &ymin, &ymax, &len);

   SetTransformation(stdbnd, xmin,xmax, ymin,ymax, len);

   xExt = WorldToMetaX(xmax)-WorldToMetaX(xmin);
   yExt = WorldToMetaY(ymin)-WorldToMetaY(ymax);
   xExt += 1.5*stdbnd; yExt += 1.5*stdbnd;

   if ((xext) > 0  &&  (yext) > 0)
   {
      FixScale((double)xExt,   (double)yExt,
	       (double)(xext), (double)(yext));
      xExt = xext; yExt = yext;
   }

   buffer[0] = '\0';	/* clear buffer */

   sprintf(buf, "e %d %d\n", xExt, yExt);
   BufferedWrite(fp, buffer, buf);
   sprintf(buf, "o %d %d\n", -xExt/2, -yExt/2);
   BufferedWrite(fp, buffer, buf);

   current_color = 0;

   SetDrawFlags(mp);

   /* Draw bonds */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      x1 = WorldToMetaX(mp->atom_array[bp->atoms[0]-1].x);
      y1 = WorldToMetaY(mp->atom_array[bp->atoms[0]-1].y);
      x2 = WorldToMetaX(mp->atom_array[bp->atoms[1]-1].x);
      y2 = WorldToMetaY(mp->atom_array[bp->atoms[1]-1].y);
      flags1 = mp->atom_array[bp->atoms[0]-1].color;
      flags2 = mp->atom_array[bp->atoms[1]-1].color;
      if (colors)    /* use coded colors */
      {
	 if (colors[bp->atoms[0]-1] < NCOLORS && colors[bp->atoms[1]-1] < NCOLORS )
	 {
            int xmid, ymid;
            xmid = (x1+x2)/2;
            ymid = (y1+y2)/2;

            /* use this colors */
            if (current_color != colors[bp->atoms[0]-1])
            {
               current_color = colors[bp->atoms[0]-1];
               sprintf(buf, "c 0x%lX\n", colortable[current_color % NCOLORS]);
               BufferedWrite(fp, buffer, buf);
            }

            PlotBond(buf, x1, y1, xmid, ymid,
	             fontscale*fontsize, bp->color, flags1, 0);
            BufferedWrite(fp, buffer, buf);

            if (current_color != colors[bp->atoms[1]-1])
            {
               current_color = colors[bp->atoms[1]-1];
               sprintf(buf, "c 0x%lX\n", colortable[current_color % NCOLORS]);
               BufferedWrite(fp, buffer, buf);
            }

            PlotBond(buf, xmid, ymid, x2, y2,
                     fontscale*fontsize, bp->color, 0, flags2);
            BufferedWrite(fp, buffer, buf);

	 }
	 else
	 {
	    if (bp->reaction_mark != NONE      				      &&
	        bp->reaction_mark != UNCHANGED 				      &&
		strcmp(mp->atom_array[bp->atoms[0]-1].atom_symbol, "H") != 0  &&
		strcmp(mp->atom_array[bp->atoms[1]-1].atom_symbol, "H") != 0)
	    {
	       /* use red */
	       if (current_color != 1)
	       {
		  current_color = 1;
		  sprintf(buf, "c 0x%lX\n", colortable[current_color % NCOLORS]);
		  BufferedWrite(fp, buffer, buf);
	       }
	    }
	    else
	    {
	       /* use black */
	       if (current_color != 0)
	       {
		  current_color = 0;
		  sprintf(buf, "c 0x%lX\n", colortable[current_color % NCOLORS]);
		  BufferedWrite(fp, buffer, buf);
	       }
	    }
	 }
      }
      else
      {
	 if (bp->reaction_mark != NONE      				   &&
	     bp->reaction_mark != UNCHANGED 				   &&
	     strcmp(mp->atom_array[bp->atoms[0]-1].atom_symbol, "H") != 0  &&
	     strcmp(mp->atom_array[bp->atoms[1]-1].atom_symbol, "H") != 0)
	 {
	    /* use red */
	    if (current_color != 1)
	    {
	       current_color = 1;
	       sprintf(buf, "c 0x%lX\n", colortable[current_color % NCOLORS]);
	       BufferedWrite(fp, buffer, buf);
	    }
	 }
	 else
	 {
	    /* use black */
	    if (current_color != 0)
	    {
	       current_color = 0;
	       sprintf(buf, "c 0x%lX\n", colortable[current_color % NCOLORS]);
	       BufferedWrite(fp, buffer, buf);
	    }
	 }
         PlotBond(buf, x1, y1, x2, y2,
	          fontscale*fontsize, bp->color, flags1, flags2);
         BufferedWrite(fp, buffer, buf);
      }
   }

   /* draw atom symbol strings */
   ChooseFont(buf, fontscale, HELVETICA);
   BufferedWrite(fp, buffer, buf);

   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);

   color_changed = FALSE;
   sprintf(buf, "c 0x%lX\n", colortable[0]);
   BufferedWrite(fp, buffer, buf);
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (0 != strcmp(ap->atom_symbol,"C")  ||
	  (ap->color & CUT)      	    ||
	  ap->charge != NONE                ||
	  ap->radical != NONE               ||
	  ap->mass_difference != NONE       ||
	  mp->n_atoms == 1                  ||
	  (labels && labels[i] && labels[i][0]))
      {

	 if (ap->charge == NONE  &&
             ap->radical == NONE &&
             ap->mass_difference == NONE)
         {
            if (0 == strcmp(ap->atom_symbol, "Cl")  && H_count[i+1] == 1)
               strcpy(atom_symbol, "HCl");
            else if (0 == strcmp(ap->atom_symbol, "Br")  && H_count[i+1] == 1)
               strcpy(atom_symbol, "HBr");
            else if (0 == strcmp(ap->atom_symbol, "I")  && H_count[i+1] == 1)
               strcpy(atom_symbol, "HI");
            else if (0 == strcmp(ap->atom_symbol, "F")  && H_count[i+1] == 1)
               strcpy(atom_symbol, "HF");
            else if (0 == strcmp(ap->atom_symbol, "O")  && H_count[i+1] == 2)
               strcpy(atom_symbol, "H2O");
            else if (0 == strcmp(ap->atom_symbol, "S")  && H_count[i+1] == 2)
               strcpy(atom_symbol, "H2S");
            else if (H_count[i+1] == 0)
	       strcpy(atom_symbol, ap->atom_symbol);
	    else if (H_count[i+1] == 1)
	       sprintf(atom_symbol, "%sH", ap->atom_symbol);
	    else if (H_count[i+1] > 1)
	       sprintf(atom_symbol, "%sH%d", ap->atom_symbol, H_count[i+1]);
         }
         else
         {
            if (H_count[i+1] == 0)
	       strcpy(atom_symbol, ap->atom_symbol);
	    else if (H_count[i+1] == 1)
	       sprintf(atom_symbol, "%sH", ap->atom_symbol);
	    else if (H_count[i+1] > 1)
	       sprintf(atom_symbol, "%sH%d", ap->atom_symbol, H_count[i+1]);
         }

	 charge = ap->charge;
	 if (0 == strcmp(atom_symbol, "R#"))
	 {
	    if (!GetNumProperty(mp->prop_lines, "M  RGP", i+1, &charge)  ||  charge < 0)
	       strcpy(atom_symbol, "R?");
	    else
	       sprintf(atom_symbol, "R%d", charge);
	    charge = 0;
	 }
#define RGB(r, g, b) 	((((long)r) << 16) | (((long)g) << 8) | (((long)b) << 0))
	 if ((flags & USE_COLORS)  &&  (colors == NULL))
	 {
	    color_changed = TRUE;
	    if (strcmp(ap->atom_symbol, "N") == 0)          /* blue */
	    {
	       sprintf(buf, "c 0x%lX\n", RGB(0, 0, 255));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "O") == 0)     /* red */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(255, 0, 0));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "S") == 0)     /* brown */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(127, 63, 0));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "P") == 0)     /* dark pink */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(127, 0, 127));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "I") == 0)     /* violet */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(255, 0, 255));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "F") == 0)     /* green */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(0, 127, 0));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "Cl") == 0)    /* green */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(0, 127, 0));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "Br") == 0)    /* green */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(0, 127, 0));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else if (strcmp(ap->atom_symbol, "C") != 0)     /* blue green */
	    {
	       sprintf(buf, "c 0x%lX\n",  RGB(0, 63, 127));
	       BufferedWrite(fp, buffer, buf);
	    }
	    else
	       color_changed = FALSE;
	 }
	 else if ((flags & USE_COLORS)  &&  (colors != NULL))
	 {
	    if (colortable[colors[i]] != 0) color_changed = TRUE;
	    else                            color_changed = FALSE;
	    sprintf(buf, "c 0x%lX\n",  colortable[colors[i]]);
	    BufferedWrite(fp, buffer, buf);
	 }

	 if (labels  &&  labels[i] && labels[i][0])	/* there is a label */
	    symbolp = labels[i];
	 else
	    symbolp = atom_symbol;

	 if (charge              != NONE  ||	/* tricky case */
	     ap->radical         != NONE  ||
	     ap->mass_difference != NONE)
	 {
	    AttributesToString(chargestring,
			       charge, ap->radical,
			       ap->mass_difference);
	    sprintf(buf, "t 01 %d %d %s %s\n",
		    WorldToMetaX(ap->x), WorldToMetaY(ap->y),
		    symbolp, chargestring);
	 }
	 else					/* ordinary case */
	    sprintf(buf, "t 0 %d %d %s\n", WorldToMetaX(ap->x),
	                                  WorldToMetaY(ap->y),
					  symbolp);
	 BufferedWrite(fp, buffer, buf);
	 if (color_changed  &&  flags & USE_COLORS)
	 {
	    sprintf(buf, "c 0x%lX\n", RGB(0, 0, 0));
	    BufferedWrite(fp, buffer, buf);
	 }
      }

   MyFree((char *)H_count);

   /* draw stereo uncertainty of double bonds */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->stereo_symbol == CIS_TRANS_EITHER)
      {
	 x1 = WorldToMetaX(mp->atom_array[bp->atoms[0]-1].x);
	 y1 = WorldToMetaY(mp->atom_array[bp->atoms[0]-1].y);
	 x2 = WorldToMetaX(mp->atom_array[bp->atoms[1]-1].x);
	 y2 = WorldToMetaY(mp->atom_array[bp->atoms[1]-1].y);
	 sprintf(buf, "t 0 %d %d %s\n", (x1+x2)/2, (y1+y2)/2, "?");
	 BufferedWrite(fp, buffer, buf);
      }

   fprintf(fp, "%s", buffer);
}

int DepictSmilesToPlotFile(FILE *fp, char *smiles, int xext, int yext)
/*
 * Writes a device independent metafile of the depiction of the structure
 * defined by the SMILES string smiles[] to the file *fp.
 *
 * The fuction returns TRUE if a depiction was created and FALSE in case of error.
 *
 * xext and yext are desired deminsions of the picture in twips or 0 if the caller
 * doesn't care.
 */
{
   struct reaccs_molecule_t *mp;

   if (!smiles)
   {
      fprintf(stderr, "DepictSmilesToPlotFile: SMILES string was NULL\n");
      return (FALSE);
   }

   mp = SMIToMOL(smiles, DO_LAYOUT);

   if (!mp)
   {
      fprintf(stderr, "MetaDepictSmiles: SMILES conversion failed for '%s'\n", smiles);
      return (FALSE);
   }

   MoleculeToPlotFile(fp, mp, xext, yext, USE_COLORS, (int *)NULL, (char **)NULL);

   FreeMolecule(mp);
   return (TRUE);
}

int DepictRxnSmilesToPlotFile(FILE *fp, char *smiles)
/*
 * Writes a device independent metafile of the depiction of the reaction
 * defined by the reaction SMILES string smiles[] to the file *fp.
 *
 * The fuction returns TRUE if a depiction was created and FALSE in case of error.
 */
{
   struct reaccs_reaction_t *rxn;
   int xext, yext;
   char *reac_stop, *prod_start;

   if (!smiles)
   {
      fprintf(stderr, "DepictRxnSmilesToPlotFile: SMILES string was NULL\n");
      return (FALSE);
   }

   reac_stop  = strchr(smiles, '>');	/* end of reactant SMILES */
   prod_start = strrchr(smiles, '>');	/* start of product SMILES */

   if (!reac_stop || !prod_start || prod_start <= reac_stop)
   {
      fprintf(stderr, "DepictRxnSmilesToPlotFile: SMILES string not a reaction\n");
      return (FALSE);
   }

   rxn = SMIToRXN(smiles);

   xext = yext = 0;
   MoleculeToPlotFile(fp,
                      rxn->reactants, xext, yext,
		      NONE, (int *)NULL, (char **)NULL);
   xext = yext = 0;
   MoleculeToPlotFile(fp,
                      rxn->products, xext, yext,
		      NONE, (int *)NULL, (char **)NULL);

   FreeReaction(rxn);
   return (TRUE);
}

void SmilesToMWMF(char *smiles, double *mwp, char *mfbuffer, int bufsize)
/*
 * Computes the molecular weight of the molecule defined
 * by smiles and puts it into *mwp. The buffer mfbuffer is filled with
 * a '\0' delimited string of the molecular formula. The caller is
 * responsible for providing enough space for this string.
 * bufsize is the usable size of the buffer including the '\0' character.
 */
{
   struct reaccs_molecule_t *mp;

   (*mwp) = 0;	mfbuffer[0]= '\0';

   mp = SMIToMOL(smiles, 0);	/* Don't do a layout */

   if (mp)
   {
      (*mwp) = MolecularWeight(mp);
      MolecularFormula(mp, mfbuffer, bufsize);
      FreeMolecule(mp);
      return;
   }
   else
      return;
}

int AttributedSmilesToPlotFile(FILE *fp,
			       char *smiles,
			       int xext, int yext,
			       char *coordinatestring,
			       char *colorstring,
			       char *labelstring)
/*
 * Writes a device independent depiction of the structure defined by the
 * SMILES string smiles[] to the file *fp.
 *
 * xext and yext are desired deminsions of the picture twips or 0 if the caller
 * doesn't care.
 *
 * The fuction returns TRUE on success and FALSE in case of failure.
 *
 * The string coordinatestring[] contains a comma separated list of
 * 2D coordinates, i.e. 2 values per atom. They are used to place the
 * the atoms of the SMILE that are not implicit hydrogen atoms.
 * The implicit ones are then positioned with the layout algorithm.
 * If coordinatestring == NULL, all atoms will be layed out.
 * Atoms with missing coordinates are also layed out. Only relative
 * positioning is considered and only connected atoms are layed out
 * as one unit.
 * E.g.: "1.54,1.54,,,3.08,3.08"
 *
 * The string colorstring[] lists the colors to be used for the different
 * atoms of the SMILES. If NULL, then a standard atom symbol dependent
 * scheme is used. Only 16 colors are supported.
 * E.g.: "0,1,1,1,2,0,0,0,3,4,5"
 *
 * The string labelstring[] lists comma separated values for
 * alternate atom labels. If this string is NULL, then standard atom
 * labels will be used. This string can also contain empty values in which
 * case the atom symbol will be used, too.
 * E.g.: "CH3,,X,Y,Z,,"
 */
{
   struct reaccs_molecule_t *mp, *mph;
   struct reaccs_bond_t *bp;
   int i, atno, *atnos;
   char *cp;
   double x, y;
   int *colors, col;
   char **labels, label[20], *cph;
   int novalue;

   if (!smiles)
   {
      fprintf(stderr, "AttributedSmilesToPlotFile: SMILES string was NULL\n");
      return (FALSE);
   }

   mp = SMIToMOL(smiles, 0);
   if (!mp)
   {
      fprintf(stderr, "AttributedSmilesToPlotFile: SMILES conversion failed\n");
      return (FALSE);
   }

   /* Save atom numbers in SMILES order */
   atnos = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++)
      atnos[i] = mp->atom_array[i].color;
   RecolorMolecule(mp);

   if (coordinatestring)
   {
      fprintf(stderr, "AttributedSmilesToPlotFile: Using coordinate string\n");
      cp = coordinatestring;
      atno = 0;
      /* parse coordinate string */
      for (;;)
      {
	 atno++; novalue = FALSE;
	 if (1 != sscanf(cp, "%lf", &x))	/* no value found */
	 {
	    x = 0.0; novalue = TRUE;
	    /* skip field */
	    while ((*cp) && (*cp) != ',') cp++; if (*cp) cp++;
	 }
	 if (1 != sscanf(cp, "%lf", &y))	/* no value found */
	 {
	    y = 0.0; novalue = TRUE;
	    /* skip field */
	    while ((*cp) && (*cp) != ',') cp++; if (*cp) cp++;
	 }

	 if (!novalue)
	    for (i=0; i<mp->n_atoms; i++)
	       if (atnos[i] == atno)
	       {
		   mp->atom_array[i].x = x;
		   mp->atom_array[i].y = y;
		   mp->atom_array[i].color = KEEP_POSITION | atno;
		   break;
	       }

	 if (*cp == '\0') break;	/* end of coordinate string */
	 if (atno > mp->n_atoms) break;	/* safeguard agains too long input */
      }

//    /* fuse connected fragments with coordinates */
//    /* TO BE IMPLEMENTED */
//
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
	 if (mp->atom_array[bp->atoms[0]-1].color & KEEP_POSITION  &&
	     mp->atom_array[bp->atoms[1]-1].color & KEEP_POSITION)
	    bp->bond_type |= DONT_FLIP_BOND;
   }
   mph = LayoutMolecule(mp);
   FreeMolecule(mp);
   mp = mph;

   if (colorstring)
   {
      colors = TypeAlloc(mp->n_atoms, int);
      atno = 0; cp = colorstring;
      while (1 == sscanf(cp, "%d", &col))
      {
	 atno++;
	 /* skip fields */
	 while ((*cp) && (*cp) != ',') cp++; if (*cp) cp++;
	 for (i=0; i<mp->n_atoms; i++)
	    if (atnos[i] == atno)
	    {
		if (col < NCOLORS) colors[i] = col;
		else               colors[i] = 0;
	        break;
	    }
	 if (atno > mp->n_atoms) break;	/* safeguard agains too long input */
      }
   }
   else
      colors = (int *)NULL;

   if (labelstring)
   {
      labels = TypeAlloc(mp->n_atoms, char *);
      atno = 0; cp = labelstring;
      while (*cp)
      {
	 atno++;
	 /* copy label */
	 cph = label;
	 while ((*cp) && (*cp) != ',' && (cph-label) < 19)
	 {
	    (*cph) = (*cp);
	    cp++; cph++;
	 }
	 (*cph) = '\0'; if (*cp) cp++;
	 /* attach label to atom */
	 for (i=0; i<mp->n_atoms; i++)
	    if (atnos[i] == atno  &&  label[0] != '\0')
	    {
		labels[i] = TypeAlloc(strlen(label)+1, char);
		strcpy(labels[i], label);
	        break;
	    }
	 if (atno > mp->n_atoms) break;	/* safeguard agains too long input */
      }
   }
   else
      labels = (char **)NULL;

   MoleculeToPlotFile(fp, mp, xext, yext, USE_COLORS, colors, labels);
   if (labels)
   {
      for (i=0; i<mp->n_atoms; i++)
	 if (labels[i]) MyFree(labels[i]);
      MyFree((char *)labels);
   }
   FreeMolecule(mp);
   if (colors) MyFree((char *)colors);
   MyFree((char *)atnos);
   return (TRUE);
}
