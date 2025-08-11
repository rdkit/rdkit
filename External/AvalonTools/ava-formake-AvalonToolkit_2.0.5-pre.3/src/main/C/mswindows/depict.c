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
/*   File:         depict.c                                             */
/*                                                                      */
/*   Purpose:      Dynamic link library to depict MDL MOL-files.        */
/*                 It provides functions to convert MOL-files into      */
/*                 plot commands and plot commands into meta files.     */
/*                                                                      */
/************************************************************************/

/* This set of macros is necessary to make a source file that exports */
/* the DLL entry points in uppercase on gcc, VC++ and Borland C++     */
#if defined(__BORLANDC__)
#else
#define _export __declspec(dllexport)  __stdcall
#endif

#include <windows.h>
#include <stdlib.h>
#include <stdio.h>

#include "reaccs.h"
#include "depict.h"
#include "sketch.h"

#define DEBUG 1
#define DEBUG_HEAP FALSE

FILE *log_file = (FILE *)NULL; /* global file */

/* Series of mappings between WIN16 and WIN32 funtions */
#ifdef WIN32

#define SetWindowExt(hdc, nXExtent, nYExtent) \
        SetWindowExtEx(hdc, nXExtent, nYExtent, NULL)

#define SetWindowOrg(hdc, X, Y) SetWindowOrgEx(hdc, X, Y, NULL)

#define MoveTo(hdc, X, Y) MoveToEx(hdc, X, Y, NULL)

DWORD GetTextExtent(HDC hdc, LPCSTR lpszString, int cbString)
{
   SIZE size;
   BOOL bool;

   bool = GetTextExtentPoint(hdc, lpszString, cbString, &size);
   return (MAKELONG(size.cx, size.cy));
}

#endif

static void checkHeap()
{
   HGLOBAL heap;
   int result = 0;
   if (!DEBUG_HEAP) return;

   heap = GetProcessHeap();
   result = HeapValidate(heap,0,NULL);
   if (!result)
      fprintf(stderr, "HeapValidate(%d, %d, %d) returns %d\n",
         heap, 0, NULL, HeapValidate(heap,0,NULL));
}

static void WMFLog(char *message, int natoms, int nbonds, int length)
{
   FILE *fp;

   fp = fopen("wmf.log", "a+");
   fprintf(fp, "%15s %3d %3d %5d\n", message, natoms, nbonds, length);
   fclose(fp);
}

#ifdef WIN32

static int stdbnd = (int)(0.35*1440);
static int penwidth = 0.9*1440/72;
static int fontsize = 10*1440/72;

static int min_penwidth = 1.0*1440/72;

#else

static int stdbnd, penwidth, fontsize;
static int stdbnd, min_penwidth;

int _export FAR PASCAL LibMain(HANDLE hInstance, WORD wDataSeg, WORD wHeapSize,
                              LPSTR lpszCmdLine)
{
   hInstance = hInstance;     /* used to satisfy the compiler */
   wDataSeg = wDataSeg;
   wHeapSize = wHeapSize;
   lpszCmdLine = lpszCmdLine;

   /* assuming MM_TWIPS => 1/1440 inch per logical unit */
   /* Std. bond = 0.35, one point = 1/72 inch, font = 10 points */

   stdbnd = (int)(0.35*1440); penwidth = 1440/72; fontsize = 10*1440/72;
   min_penwidth = 1.0*1440/72;

   return (1);
}
#endif

HFONT DepictFont(double scale)
{
   HFONT hFont;

   if (scale < 0.5) scale = 0.5;
   hFont = CreateFont((-1)*           /* use character height not cell height */
                      (scale*         /* scale if needed */
                      fontsize+0.5),  /* correct for rounding */
                      0,              /* use best guess for width */
                      0,              /* parallel to bottom of page */
                      0,              /* parallel to bottom of page */
                      800,            /* extra bold weight */
                      0,              /* not slanted */
                      0,              /* not underlined */
                      0,              /* don't strike out */
                      ANSI_CHARSET,
                      OUT_DEFAULT_PRECIS,
                      CLIP_DEFAULT_PRECIS,
                      ANTIALIASED_QUALITY,
                      VARIABLE_PITCH|FF_SWISS,
                      "Arial");

   return (hFont);
}

#include <math.h>
#include <string.h>

static void SendMetaFileToClipboard(HANDLE hmf, int xExt, int yExt)
{
   GLOBALHANDLE   hGMem;
   LPMETAFILEPICT lpMFP;

   hGMem = GlobalAlloc(GHND, (DWORD) sizeof (METAFILEPICT));
   lpMFP = (LPMETAFILEPICT) GlobalLock(hGMem);

   lpMFP->mm   = MM_ANISOTROPIC;
   lpMFP->xExt = xExt;
   lpMFP->yExt = yExt;
   lpMFP->hMF  = hmf;

   GlobalUnlock(hGMem);

   OpenClipboard(NULL);
   EmptyClipboard();
   SetClipboardData(CF_METAFILEPICT, hGMem);
   CloseClipboard();
}

HGLOBAL FileToMemoryHandle(LPSTR fname, int *fsizep)
/*
 * Reads *fname and returns a handle to a global shared memory block
 * which contains the lines of the file as (line length)+data
 * records. This is the format used by the MDLCT clipboard entry.
 * *fsizep is set to the size of the used part of the allocated memory.
 * The two sizes differ because space has been allocated on the basis
 * of the file size including CR+LF sequences, but only one byte of
 * length information is used in the memory coding of the file.
 */
{
   char buffer[256];      /* line buffer */
   FILE *fp;
   long fsize;             /* size of data block to being allocated */
   int lsize;
   HGLOBAL returnhnd;
   LPSTR molbuffer, cp;
   int data_length;

   fp = fopen(fname,"r");    /* get size of file */
   fseek(fp, 0L, SEEK_END);
   fsize = ftell(fp);
   fseek(fp, 0L, SEEK_SET);

   returnhnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, fsize);
   molbuffer = GlobalLock(returnhnd);

   data_length = 0;
   cp = molbuffer;
   while (fgets(buffer, 255, fp))
   {
      buffer[255] = '\0';
      lsize = strlen(buffer);
      lsize--; buffer[lsize] = '\0';     /* get rid of '\n' */
      data_length += strlen(buffer)+1;
      (*cp) = lsize; cp++;
      strcpy(cp, buffer); cp += lsize;
   }
   GlobalUnlock(returnhnd);
   fclose(fp);
   (*fsizep) = data_length;
   (*fsizep) = (int)fsize;

   return (returnhnd);
}

#include "local.h"
#include "forio.h"
#include "reaccs.h"
#include "reaccsio.h"
#include "utilities.h"
#include "perceive.h"
#include "layout.h"
#include "smi2mol.h"

#include "mymacros.h"

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
#define S_D           256
#define S_A           512
#define D_A           1024

#define CUT          1
#define CUT_SEC      2

void DrawBond(HDC hdcMeta,
              int x1, int y1, int x2, int y2,
              int d, int type, int flags1, int flags2)
/*
 * Draws a bond from (x1,y1) to (x2,y2). It uses the font size d
 * as the separation characteristic, the symbol type, and the
 * shortening flags flags1 and flags2 to determine the appearance.
 */
{
   int i, len;
   int dbx, dby;      /* multiple bond spacing */
   double bf = 0.45;
   int dwx, dwy;      /* wedge width */
   double wf = 0.35;
   double dhx, dhy;   /* hash line distance */
   int nhash;
   double hf = 0.40;
   int dcx, dcy;      /* line end cutting */
   double cf = 0.55;
   int ddx, ddy;      /* multiple bond cutting */
   double df = 0.45;
   POINT wedge[4];

   len = (int)sqrt((double)(x1-x2)*(double)(x1-x2) +
                      (double)(y1-y2)*(double)(y1-y2));
   if (len == 0) return;

   dbx = bf*((double)d*(x2-x1))/len; dby = bf*((double)d*(y2-y1))/len;
   dwx = wf*((double)d*(x2-x1))/len; dwy = wf*((double)d*(y2-y1))/len;
   dcx = cf*((double)d*(x2-x1))/len; dcy = cf*((double)d*(y2-y1))/len;
   ddx = df*((double)d*(x2-x1))/len; ddy = df*((double)d*(y2-y1))/len;

   if (flags1 & CUT) { x1 += dcx; y1 += dcy; }
   if (flags2 & CUT) { x2 -= dcx; y2 -= dcy; }

   if (type & SNG)
   {
      MoveTo(hdcMeta, x1, y1);
      LineTo(hdcMeta, x2, y2);
   }

   if (type & BND_RIGHT)
   {
      if (flags1 & CUT_SEC)
         MoveTo(hdcMeta, x1-dby+ddx, y1+dbx+ddy);
      else
         MoveTo(hdcMeta, x1-dby,     y1+dbx);
      if (flags2 & CUT_SEC)
         LineTo(hdcMeta, x2-dby-ddx, y2+dbx-ddy);
      else
         LineTo(hdcMeta, x2-dby,     y2+dbx);
   }

   if (type & BND_LEFT)
   {
      if (flags1 & CUT_SEC)
         MoveTo(hdcMeta, x1+dby+ddx, y1-dbx+ddy);
      else
         MoveTo(hdcMeta, x1+dby,     y1-dbx);
      if (flags2 & CUT_SEC)
         LineTo(hdcMeta, x2+dby-ddx, y2-dbx-ddy);
      else
         LineTo(hdcMeta, x2+dby,     y2-dbx);
   }

   if (type & SYM_DBL)
   {
      MoveTo(hdcMeta, x1-dby/2, y1+dbx/2);
      LineTo(hdcMeta, x2-dby/2, y2+dbx/2);
      MoveTo(hdcMeta, x1+dby/2, y1-dbx/2);
      LineTo(hdcMeta, x2+dby/2, y2-dbx/2);
   }

   if (type & WDG)
   {
      wedge[0].x = x1; wedge[0].y = y1;
      wedge[1].x = x2-dwy/2; wedge[1].y = y2+dwx/2;
      wedge[2].x = x2+dwy/2; wedge[2].y = y2-dwx/2;
      wedge[3].x = x1; wedge[3].y = y1;
      Polygon (hdcMeta, wedge, 4);
   }

   if (type & HSH)
   {
      nhash = (int)(len/(hf*d)+0.5);
      if (nhash < 2) nhash = 2;
      dhx = wf*((double)d*(x2-x1))/len; dhy = wf*((double)d*(y2-y1))/len;
      for (i=0; i<=nhash; i++)
      {
         MoveTo(hdcMeta, (int)(x1+i*(x2-x1-dhy/2)/nhash),
                         (int)(y1+i*(y2-y1+dhx/2)/nhash));
         LineTo(hdcMeta, (int)(x1+i*(x2-x1+dhy/2)/nhash),
                         (int)(y1+i*(y2-y1-dhx/2)/nhash));
      }
   }

   if (type & DSH)
   {
      if ((type & SNG))
      {
         x1 += -dby+ddx; y1 -= -dbx-ddy;
         x2 += -dby-ddx; y2 -= -dbx+ddy;
      }
      nhash = (int)(len/(hf*d)+0.5);
      if (nhash < 2) nhash = 2;
      for (i=0; i<nhash; i+=2)
      {
         MoveTo(hdcMeta, (int)(x1+i*(x2-x1)/nhash),
                         (int)(y1+i*(y2-y1)/nhash));
         LineTo(hdcMeta, (int)(x1+(i+1)*(x2-x1)/nhash),
                         (int)(y1+(i+1)*(y2-y1)/nhash));
      }
   }
}

static void GetWorldWindow(struct reaccs_molecule_t *mp,
                           double *xminp, double *xmaxp,
                           double *yminp, double *ymaxp,
                           double *lenp)
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
 * Sets the global scale variable such that the picture will now
 * fit into the box defined by xext_new and yext_new.
 */
{
   if (xext_new*0.5 < xext_best  &&  xext_best < xext_new  &&
       yext_new*0.5 < yext_best  &&  yext_best < yext_new)
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

   if (fontscale < 0.4) fontscale = 0.4;
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
          ap->query_H_count != NONE          ||
          ap->radical != NONE                ||
          nba[i].n_ligands == 0              ||
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
         case SINGLE_AROMATIC: bp->color |= DSH;
                               bp->color |= SNG;
                               bp->color |= S_A;
                               break;
         case DOUBLE_AROMATIC: bp->color |= DSH;
                               bp->color |= SNG;
                               bp->color |= D_A;
                               break;
         case SINGLE_DOUBLE: bp->color |= DSH;
                             bp->color |= SNG;
                             bp->color |= S_D;
                             break;
         default: bp->color |= DSH;
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

struct placer_t
   {
      short key1, key2;         /* 0xCDD7 and 0x9AC6 */
      short hmf;                /* == 0 for files */
      short x1; short y1;       /* lower left and upper right corner */
      short x2; short y2;       /* in logical coordinates */
                                /* use SetWindowExt/SetWindowOrg to set */
                                /* window coordinates */
                                /* mapping mode is always MM_ANISOTROPIC */
      short wInch;              /* dots per inch */
      short res1, res2;         /* == 0 */
      short checksum;           /* XOR of all other components */
   };

static void MakePlacer(struct placer_t *pp, int xext, int yext)
/*
 * Puts the placer *pp with recomputed check sum to the file *fp.
 */
{
   pp->x1 = -xext/2; pp->y1 = -yext/2;
   pp->x2 =  xext/2; pp->y2 =  yext/2;
   pp->wInch = 1440;    /* 1440 twips == 1 inch */

   pp->key1 = 0xCDD7; pp->key2 = 0x9AC6;
   pp->hmf = 0;
   pp->res1 = 0; pp->res2 = 0;
   pp->checksum = 0;

   pp->checksum = pp->checksum ^ pp->key1 ^ pp->key2
                               ^ pp->hmf
                               ^ pp->x1   ^ pp->y1
                               ^ pp->x2   ^ pp->y2
                               ^ pp->wInch
                               ^ pp->res1 ^ pp->res2;
}

#define USE_COLORS 0x001
#define BINARY_WMF 0x010
#define NO_PLACER  0x020   /* suppress placer for binary metafiles */

#define NCOLORS               10
static long colortable[NCOLORS] =
   {
      RGB(  0,   0,   0),        /* black */
      RGB(255,   0,   0),        /* red */
      RGB(  0,   0, 255),        /* blue */
      RGB(  0, 128,   0),        /* dark green */
      RGB(128,  64,   0),        /* brown */
      RGB( 64,   0, 128),        /* dark violet */
      RGB(  0,  64,  64),        /* dark cyan */
      RGB(128, 128,   0),        /* dark yellow */
      RGB(128,   0, 128),        /* dark pink */
      RGB(128, 128, 128),        /* grey */
   };

int  MoleculeToWMFBuffer(struct reaccs_molecule_t *mp,
                         char *buffer, int bufsize,
                         int *xextp, int *yextp,
                         int flags,
                         int colors[],
                         char *labels[])
/*
 * Chemical structure depiction workhorse. Takes the MDL formatted
 * chemical structure *mp and converts it into WMF format. This
 * WMF is either posted to the clipboard (buffer == NULL) or copied
 * to buffer.
 * 
 * The flags parameter provides information
 * on how to render the WMF (color or b/w) and which kind of output
 * to produce (binary or hex-code), and whether a placable header
 * shall be prepended.
 *
 * *xextp and *yextp receive the x- and y- dimensions of the resulting
 * image in twips.
 *
 * The array colors[] is used to color the different parts of the molecule.
 * colors[i] == 0 means draw in black. 1 <= colors[i] <= 9 means use one of
 * the nine most dark colors to get reasonable readablility. If colors ==
 * NULL, the standard atom label coloring with black bonds is used.
 *
 * The array labels[] contains alternate atom labels to be used instead
 * of the ones defined in the atom_symbol fields.
 * Labels which are "\0" strings are not considered. Use " " to get
 * a blank label.
 *
 * The function returns <=0 in case of failure and the used size of buffer
 * when succeeded. 0 means general failure, while < 0 is the negative size
 * of the required buffer.
 */
{
   struct placer_t placer;       /* header for ALDUS placable metafile */
   int             result;
   char           *lpwmf;             /* pointer to metafile memory */

   static HANDLE   hmf;
   HDC             hdcMeta, hdc;
   HPEN            hPen;
   HPEN            hPenArray[NCOLORS];
   int             current_color;
   HBRUSH          hBrush;
   HFONT           hFont, hFontMeta, hFontScr;
   TEXTMETRIC      tm;
   DWORD           txtExt;
   int             txtHeight, txtWidth;

   struct reaccs_bond_t     *bp;
   struct reaccs_atom_t     *ap;

   double xmin, xmax, ymin, ymax, len;
   int i, nout;
   int color_changed;
   int xExt, yExt;
   int x1, y1, x2, y2;
   int flags1, flags2;
   char chargestring[10];
   char atom_symbol[10], *symbolp;
   int charge;
   int *H_count;

   if (!mp) return (0);
   if (mp->n_atoms > MAXATOMS)
   {
      ShowMessageI("Too many (%d) atoms in molecule",
                   "MoleculeToWMFBuffer", mp->n_atoms);
      return (0);
   }

   GetWorldWindow(mp, &xmin, &xmax, &ymin, &ymax, &len);

   SetTransformation(stdbnd, xmin,xmax, ymin,ymax, len);

   xExt = WorldToMetaX(xmax)-WorldToMetaX(xmin);
   yExt = WorldToMetaY(ymin)-WorldToMetaY(ymax);

   xExt += 1.5*stdbnd; yExt += 1.5*stdbnd;

   if ((*xextp) > 0  &&  (*yextp) > 0)
   {
      FixScale((double)xExt,     (double)yExt,
               (double)(*xextp), (double)(*yextp));
      xExt = *xextp; yExt = *yextp;
   }
   else
   {
      (*xextp) = xExt; (*yextp) = yExt;              /* give size info */
   }

   hdcMeta = CreateMetaFile(NULL);

   SetMapMode(hdcMeta, MM_ANISOTROPIC); /* Test */
   SetWindowExt(hdcMeta, xExt, yExt);
   SetWindowOrg(hdcMeta, -xExt/2, -yExt/2);

   if (colors)
   {
       for (i=0; i<NCOLORS; i++)
          hPenArray[i] = CreatePen(PS_SOLID,
                                   MAX(min_penwidth, fontscale*penwidth),
                                   colortable[i]);
       SelectObject(hdcMeta, hPenArray[0]);
       current_color = 0;
   }
   else
   {
      hPen   = CreatePen(PS_SOLID,
                         MAX(min_penwidth, fontscale*penwidth),
                         RGB(0, 0, 0));
      SelectObject(hdcMeta, hPen);
   }

   hBrush = CreateSolidBrush(RGB(0,0,0));
   SelectObject(hdcMeta, hBrush);
   SetBkMode(hdcMeta, OPAQUE);

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
         if (colors[bp->atoms[0]-1]  ==  colors[bp->atoms[1]-1]   &&
             colors[bp->atoms[0]-1] < NCOLORS)
         {   /* use this color */
            if (current_color != colors[bp->atoms[0]-1])
            {
               current_color = colors[bp->atoms[0]-1];
               SelectObject(hdcMeta, hPenArray[current_color]);
            }
         }
         else
         {   /* use black */
            if (current_color != 0)
            {
               current_color = 0;
               SelectObject(hdcMeta, hPenArray[current_color]);
            }
         }
      }
      DrawBond(hdcMeta,
               x1, y1, x2, y2,
               fontscale*fontsize,
               bp->color, flags1, flags2);
   }

   /* draw atom symbol strings */
   hFontMeta = DepictFont(fontscale);
   SelectObject(hdcMeta, hFontMeta);
   SetTextAlign(hdcMeta, TA_BASELINE);

   H_count = TypeAlloc(mp->n_atoms+1, int);
   ComputeImplicitH(mp, H_count);

   hFont = DepictFont(fontscale);
   hdc = CreateDC("DISPLAY", NULL, NULL, NULL);
   hFontScr = SelectObject(hdc, hFont);
   GetTextMetrics(hdc, &tm);
   SetTextColor(hdcMeta, RGB(0, 0, 0));
   for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
      if (0 != strcmp(ap->atom_symbol,"C")  ||
          ap->charge != NONE                ||
          ap->radical != NONE               ||
          ap->mass_difference != NONE       ||
          (ap->color & CUT)                 ||
          mp->n_atoms == 1                  ||
          (labels && labels[i] && labels[i][0]))
      {

         if (ap->charge          == NONE  &&
             ap->radical         == NONE &&
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
         if (0 == strncmp(atom_symbol, "R#", 2))
         {
            if (!GetNumProperty(mp->prop_lines, "M  RGP", i+1, &charge)  ||
                charge < 0)
               strcpy(atom_symbol, "R?");
            else
            {
               atom_symbol[1] = '0'+charge;
            }
            charge = 0;
         }
         if ((flags & USE_COLORS)  &&  (colors == NULL))
         {
            color_changed = TRUE;
            if (strcmp(ap->atom_symbol, "N") == 0)          /* blue */
               SetTextColor(hdcMeta, RGB(0, 0, 255));
            else if (strcmp(ap->atom_symbol, "O") == 0)     /* red */
               SetTextColor(hdcMeta, RGB(255, 0, 0));
            else if (strcmp(ap->atom_symbol, "S") == 0)     /* brown */
               SetTextColor(hdcMeta, RGB(127, 64, 0));
            else if (strcmp(ap->atom_symbol, "P") == 0)     /* dark pink */
               SetTextColor(hdcMeta, RGB(127, 0, 127));
            else if (strcmp(ap->atom_symbol, "I") == 0)     /* violet */
               SetTextColor(hdcMeta, RGB(255, 0, 255));
            else if (strcmp(ap->atom_symbol, "F") == 0)     /* green */
               SetTextColor(hdcMeta, RGB(0, 128, 0));
            else if (strcmp(ap->atom_symbol, "Cl") == 0)    /* green */
               SetTextColor(hdcMeta, RGB(0, 128, 0));
            else if (strcmp(ap->atom_symbol, "Br") == 0)    /* green */
               SetTextColor(hdcMeta, RGB(0, 128, 0));
            else if (strcmp(ap->atom_symbol, "C") != 0)     /* blue green */
               SetTextColor(hdcMeta, RGB(0, 64, 128));
           else
               color_changed = FALSE;
         }
         else if ((flags & USE_COLORS)  &&  (colors != NULL))
         {
            if (colortable[colors[i]] != 0) color_changed = TRUE;
            else                            color_changed = FALSE;
            SetTextColor(hdcMeta, colortable[colors[i]]);
         }

         if (labels  &&  labels[i] && labels[i][0])  /* there is a label */
            symbolp = labels[i];
         else
            symbolp = atom_symbol;
         txtExt = GetTextExtent(hdc, symbolp, 1);
         txtWidth  = LOWORD(txtExt);
         txtHeight = HIWORD(txtExt) -
                     tm.tmInternalLeading -
                     tm.tmDescent;

         if (charge              != NONE  ||
             ap->radical         != NONE  ||
             ap->sub_desc        != NONE  ||
             ap->mass_difference != NONE)
         {
            SetTextAlign(hdcMeta, TA_TOP|TA_UPDATECP);
            MoveTo(hdcMeta, WorldToMetaX(ap->x) -txtWidth/2,
                     - tm.tmAscent +
                       WorldToMetaY(ap->y)+txtHeight/2);
            TextOut(hdcMeta, 0, 0, /* coords. ignored with TA_UPDATECP */
                    symbolp,
                    strlen(symbolp));
            AttributesToString(chargestring,
                               charge, ap->radical,
                               ap->mass_difference);
            SetTextAlign(hdcMeta, TA_TOP|TA_UPDATECP);
            MoveTo(hdcMeta,
                   WorldToMetaX(ap->x) -txtWidth/2
                    + LOWORD(GetTextExtent(hdc,
                                           symbolp,
                                           strlen(symbolp))),
                    - tm.tmAscent + WorldToMetaY(ap->y));
            TextOut(hdcMeta, 0, 0, /* coords. ignored with TA_UPDATECP */
                    chargestring, strlen(chargestring));
            if (ap->sub_desc != 0)
            {
               if (ap->sub_desc == (-1)) TextOut(hdcMeta, 0, 0, " (s0)", strlen(" (s0)"));
               else if (ap->sub_desc == (-2)) TextOut(hdcMeta, 0, 0, " (s*)", strlen(" (s*)"));
               else if (ap->sub_desc == 1) TextOut(hdcMeta, 0, 0, " (s1)", strlen(" (s1)"));
               else if (ap->sub_desc == 2) TextOut(hdcMeta, 0, 0, " (s2)", strlen(" (s2)"));
               else if (ap->sub_desc == 3) TextOut(hdcMeta, 0, 0, " (s3)", strlen(" (s3)"));
               else if (ap->sub_desc == 4) TextOut(hdcMeta, 0, 0, " (s4)", strlen(" (s4)"));
               else if (ap->sub_desc == 5) TextOut(hdcMeta, 0, 0, " (s5)", strlen(" (s5)"));
               else TextOut(hdcMeta, 0, 0, " (s?)", strlen(" (s?)"));
            }
            SetTextAlign(hdcMeta, TA_TOP|TA_NOUPDATECP);
         }
         else
         {
            SetTextAlign(hdcMeta, TA_TOP|TA_NOUPDATECP);
            TextOut(hdcMeta, WorldToMetaX(ap->x) -txtWidth/2,
                            - tm.tmAscent + WorldToMetaY(ap->y)+txtHeight/2,
                            symbolp,
                            strlen(symbolp));
         }
         if (color_changed && flags & USE_COLORS)
            SetTextColor(hdcMeta, RGB(0, 0, 0));
      }

   // label query bonds if needed
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
   {
      if (bp->color & (S_D|S_A|D_A))         // Only draw bond label if needed
      {
         txtExt = GetTextExtent(hdc, "s/d", 3);
         txtWidth  = LOWORD(txtExt);
         txtHeight = HIWORD(txtExt) -
                     tm.tmInternalLeading -
                     tm.tmDescent;
         x1 = WorldToMetaX(mp->atom_array[bp->atoms[0]-1].x);
         y1 = WorldToMetaY(mp->atom_array[bp->atoms[0]-1].y);
         x2 = WorldToMetaX(mp->atom_array[bp->atoms[1]-1].x);
         y2 = WorldToMetaY(mp->atom_array[bp->atoms[1]-1].y);
         SetTextAlign(hdcMeta, TA_CENTER|TA_BASELINE|TA_NOUPDATECP);
         SetTextColor(hdcMeta, RGB(127, 127, 255));
         if (bp->color & S_D)
            TextOut(hdcMeta, (x1+x2)/2, (y1+y2)/2+txtHeight/2, "s/d", strlen("s/d"));
         else if (bp->color & S_A)
            TextOut(hdcMeta, (x1+x2)/2, (y1+y2)/2+txtHeight/2, "s/a", strlen("s/a"));
         else if (bp->color & D_A)
            TextOut(hdcMeta, (x1+x2)/2, (y1+y2)/2+txtHeight/2, "d/a", strlen("d/a"));
         else
            TextOut(hdcMeta, (x1+x2)/2, (y1+y2)/2+txtHeight/2, "?/?", strlen("?/?"));
         SetTextColor(hdcMeta, RGB(0, 0, 0));
      }
   }

   MyFree((char *)H_count);

   /* draw stereo uncertainty of double bonds */
   for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
      if (bp->stereo_symbol == CIS_TRANS_EITHER  &&
          !(bp->dummy & 252)) /* Ignore in small rings */
      {
         x1 = WorldToMetaX(mp->atom_array[bp->atoms[0]-1].x);
         y1 = WorldToMetaY(mp->atom_array[bp->atoms[0]-1].y);
         x2 = WorldToMetaX(mp->atom_array[bp->atoms[1]-1].x);
         y2 = WorldToMetaY(mp->atom_array[bp->atoms[1]-1].y);
         SetTextAlign(hdcMeta, TA_TOP|TA_NOUPDATECP);
         TextOut(hdcMeta, (x1+x2)/2 -txtWidth/2,
                           - tm.tmAscent +
                            (y1+y2)/2 +txtHeight/2,
                           "?", strlen("?"));
      }
   SelectObject(hdc, hFontScr);
   DeleteObject(hFont);
   DeleteDC(hdc);

   SelectObject(hdcMeta, GetStockObject(DEVICE_DEFAULT_FONT));
   DeleteObject(hFontMeta);

   SelectObject(hdcMeta, GetStockObject(BLACK_PEN));
   if (colors)
      for (i=0; i<NCOLORS; i++)
         DeleteObject(hPenArray[i]);
   else
      DeleteObject(hPen);

   SelectObject(hdcMeta, GetStockObject(BLACK_BRUSH));
   DeleteObject(hBrush);

   hmf = CloseMetaFile(hdcMeta);

   /* Metafile was created in global memory. */
   /* Now put it into buffer.                */

#ifdef WIN32
   result = (int)GetMetaFileBitsEx(hmf, 0, NULL);
#else
   result = (int)GlobalSize(hmf);
#endif

   if (!buffer)
   {
      SendMetaFileToClipboard(hmf, xExt, yExt);
      return (result);
   }

   if (result + sizeof(struct placer_t) > bufsize)
   {                              /* Not enough space in buffer */
      DeleteMetaFile(hmf);
       return (-(result + sizeof(struct placer_t)));
   }

   if (flags & BINARY_WMF)
   {
#ifdef WIN32
      result = GetMetaFileBitsEx(hmf, 0, NULL);
      lpwmf = TypeAlloc(result+10, char);
      result = GetMetaFileBitsEx(hmf, result+10, lpwmf);
      DeleteMetaFile(hmf);
#else
      result = (int)GlobalSize(hmf);
      lpwmf = GlobalLock(hmf);
#endif
      if (flags & NO_PLACER)
      {
         memcpy(buffer, lpwmf, result);
      }
      else
      {
         MakePlacer(&placer, xExt, yExt);
         memcpy(buffer, &placer, sizeof(struct placer_t));
         memcpy(buffer + sizeof(struct placer_t), lpwmf, result);
         result += sizeof(struct placer_t);
      }
#ifdef WIN32
      MyFree((char *)lpwmf);
#else
      GlobalUnlock(hmf);
      DeleteMetaFile(hmf);
#endif
      return (result);
   }
   else
   {
#ifdef WIN32
      result = GetMetaFileBitsEx(hmf, 0, NULL);
      lpwmf = TypeAlloc(result, char);
      result = GetMetaFileBitsEx(hmf, result, lpwmf);
      DeleteMetaFile(hmf);
#else
      hmf = GetMetaFileBits(hmf);
      result = (int)GlobalSize(hmf);
      if (result*2.1+2 > bufsize)    /* ASCII version to big for buffer */
      {
         DeleteMetaFile(hmf);
         return (-(result*2.1+2));
      }
      lpwmf = GlobalLock(hmf);
#endif

      for (i=0; i<4; i++)    /* patch WMF size into WMF */
      {
         lpwmf[6+i] = (((result/2) >> (8*i))) & 0xFF;
      }
      buffer[0] = '\0';
      nout = 0;
      for (i=0; i<result; i++)
      {
         sprintf(buffer, "%02X", 0xFF & (unsigned)lpwmf[i]); buffer+=2;
         nout += 2;
         if (i % 20  == 19)
         {
            sprintf(buffer, "\r\n"); buffer+=2;
            nout += 2;
         }
      }
      if (i % 20  !=  0)
      {
         sprintf(buffer, "\r\n"); buffer+=2;
         nout += 2;
      }
#ifdef WIN32
      MyFree((char *)lpwmf);
#else
      GlobalUnlock(hmf);
      GlobalFree(hmf);
#endif
      return (nout);
   }
}

int _export FAR PASCAL WMFDepictSmiles(LPSTR smiles, LPSTR fname,
                                       int *xextp, int *yextp)
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] to the file named fname[].
 * The fuction returns TRUE if a depiction was created and FALSE otherwise.
 */
{
   struct reaccs_molecule_t *mp;
   struct reaccs_bond_t     *bp;
   int i;
   int result;
   int nstereo;

   HFILE fpout;
   char *buffer;

checkHeap();

   if (!smiles  ||  !fname)
   {
      MessageBox(NULL, "SMILES string or file name was NULL",
                       "WMFDepictSmiles", MB_OK);
      return (FALSE);
   }


   mp = SMIToMOL(smiles, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);

   if (mp)
   {
      /* Conservetive estimate of required buffer size 950+65*n_atoms */
      nstereo = 0;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->stereo_symbol != NONE) nstereo++;
      buffer = TypeAlloc(950+65*mp->n_atoms+100*nstereo, char);
      result = MoleculeToWMFBuffer(mp,
                                   buffer, 950+65*mp->n_atoms+100*nstereo,
                                   xextp, yextp,
                                   USE_COLORS | BINARY_WMF,
                                   (int *)NULL,
                                   (char **)NULL);
      if (result <= 0)
         result = FALSE;
      else
      {
         fpout = _lcreat(fname, 0);
         _lwrite(fpout, buffer, result);
         _lclose(fpout);
         result = TRUE;
      }
      MyFree(buffer);

      FreeMolecule(mp);
      return (result);
   }
   else
      return (FALSE);
}

_export void FAR PASCAL SmilesToMWMF(LPSTR smiles,
                                     double *mwp,
                                     LPSTR mfbuffer, int nbuf)
/*
 * Computes the molecular weight of the molecule defined
 * by smiles and puts it into *mwp. The buffer mfbuffer is filled with
 * a '\0' delimited string of the molecular formula. The caller is
 * responsible for providing enough space for this string.
 */
{
   struct reaccs_molecule_t *mp;

   (*mwp) = 0; mfbuffer[0]= '\0';
   if (!smiles  ||  !mfbuffer)
   {
      MessageBox(NULL, "SMILES string or formula buffer was NULL", "SmilesToMWMF", MB_OK);
      return;
   }

   mp = SMIToMOL(smiles, 1); /* Don't do a layout */

   if (mp)
   {
      (*mwp) = MolecularWeight(mp);
      MolecularFormula(mp, mfbuffer, nbuf);
      FreeMolecule(mp);
      return;
   }
   else
      return;
}

int _export FAR PASCAL AttributedSmilesToWMFBuffer(LPSTR smiles,
                                                   LPSTR buffer, int bufsize,
                                                   int *xextp, int *yextp,
                                                   LPSTR coordinatestring,
                                                   LPSTR colorstring,
                                                   LPSTR labelstring)
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] to the buffer buffer[0..bufsize-1].
 *
 * The fuction returns the number of bytes written to buffer[] or 0
 * in case of failure.
 *
 * The string coordinatestring[] contains a comma separated list of
 * 2D coordinates, i.e. 2 values per atom. They are used to place the
 * the atoms of the SMILES that are not implicit hydrogen atoms.
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
 *
 * The function prepends the ALDUS placeable header to the buffer.
 *
 * It returns the size of the buffer used if successful, and a number <= 0
 * on failure. 0 means general error, a number < 0 is the size of the necessary
 * buffer to hold all the data. buffer[] is not changed in this case.
 */
{
   struct reaccs_molecule_t *mp, *mph;
   struct reaccs_bond_t *bp;
   int result;
   int i, atno, *atnos;
   char *cp;
   double x, y;
   double xmin, xmax, ymin, ymax;
   int *colors, col;
   char **labels, label[20], *cph;
   int novalue;
   int ncoordinates;

checkHeap();

   if (!smiles  ||  !buffer)
   {
      MessageBox(NULL, "SMILES string or buffer was NULL",
                       "AttributedSmilesToWMFBuffer", MB_OK);
      return (FALSE);
   }

   mp = SMIToMOL(smiles, 0);
   if (!mp) return(FALSE);

   /* Save atom numbers in SMILES order */
   atnos = TypeAlloc(mp->n_atoms, int);
   for (i=0; i<mp->n_atoms; i++)
      atnos[i] = mp->atom_array[i].color;
   RecolorMolecule(mp);

   ncoordinates = 0;
   if (coordinatestring)
   {
      xmin = ymin =  1.0e+7;
      xmax = ymax = -1.0e+7;
      cp = coordinatestring;
      atno = 0;
      /* parse coordinate string */
      for (;;)
      {
         atno++; novalue = FALSE;
         if (1 != sscanf(cp, "%lf", &x)) /* no value found */
         {
            x = 0.0; novalue = TRUE;
            /* skip field */
         }
         while ((*cp) && (*cp) != ',') cp++; if (*cp) cp++;
         if (1 != sscanf(cp, "%lf", &y)) /* no value found */
         {
            y = 0.0; novalue = TRUE;
            /* skip field */
         }
         while ((*cp) && (*cp) != ',') cp++; if (*cp) cp++;

         if (!novalue)
         {
            for (i=0; i<mp->n_atoms; i++)
               if (atnos[i] == atno)
               {
                   mp->atom_array[i].x = x;
                   mp->atom_array[i].y = y;
                   if (xmin > x) xmin = x; if (xmax < x) xmax = x;
                   if (ymin > y) ymin = y; if (ymax < y) ymax = y;
                   mp->atom_array[i].color = KEEP_POSITION;
                   break;
               }
            ncoordinates++;
         }

         if (*cp == '\0') break;        /* end of coordinate string */
         if (atno > mp->n_atoms) break; /* safeguard agains too long input */
      }

      /* fuse connected fragments with coordinates */
      /* TO BE IMPLEMENTED */

      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (mp->atom_array[bp->atoms[0]-1].color & KEEP_POSITION  &&
             mp->atom_array[bp->atoms[1]-1].color & KEEP_POSITION)
            bp->bond_type |= DONT_FLIP_BOND;

      if (xmax <= xmin  ||  ymax <= ymin) /* invalid coordinates */
         RecolorMolecule(mp);
   }
   mph = LayoutMolecule(mp);

   /*
    * Assume that hydrogens are superflous when there are coordinates
    */
   if (ncoordinates > 2)
      MakeSomeHydrogensImplicit(mph, NO_QUERY_H_COUNT | ANCILLARY_STEREO);

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
         if (atno > mp->n_atoms) break; /* safeguard agains too long input */
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
         if (atno > mp->n_atoms) break; /* safeguard agains too long input */
      }
   }
   else
      labels = (char **)NULL;

   result = MoleculeToWMFBuffer(mp, buffer, bufsize, xextp, yextp,
                                USE_COLORS | BINARY_WMF,
                                colors,
                                labels);
   if (labels)
   {
      for (i=0; i<mp->n_atoms; i++)
         if (labels[i]) MyFree(labels[i]);
      MyFree((char *)labels);
   }
   FreeMolecule(mp);
   if (colors) MyFree((char *)colors);
   MyFree((char *)atnos);
   return (result);
}

_export int FAR PASCAL SmilesToWMFBuffer(LPSTR smiles,
                                         LPSTR buffer, int bufsize,
                                         int *xextp, int *yextp)
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] to the buffer buffer[0..bufsize-1].
 *
 * The fuction returns the number of bytes written to buffer[], 0
 * in case of failure, or -n where n is the required size of the buffer
 * if the buffer is too small.
 *
 * The function prepends the ALDUS placeable header to the buffer.
 */
{
   int result;

   result = AttributedSmilesToWMFBuffer(smiles,
                                       buffer, bufsize,
                                       xextp, yextp,
                                       (char *) NULL,
                                       (char *) NULL,
                                       (char *) NULL);
   return (result);
}

int _export FAR PASCAL RTFDepictSmiles(LPSTR smiles, LPSTR fname,
                                       int *xextp, int *yextp)
/*
 * Writes a WMF depiction of the structure defined by the SMILES string
 * smiles[] in hex-code form to the file named fname[].
 * This form can directly be used for WMF fields in RTF files.
 */
{
   struct reaccs_molecule_t *mp;
   struct reaccs_bond_t     *bp;
   int nstereo, i;
   int result;

   HFILE fpout;
   char *buffer;
   int bufsize;

checkHeap();

   if (!smiles  ||  !fname)
   {
      MessageBox(NULL, "SMILES string or file name was NULL",
                       "WMFDepictSmiles", MB_OK);
      return (FALSE);
   }

   mp = SMIToMOL(smiles, DO_LAYOUT | DROP_TRIVIAL_HYDROGENS);

   result = TRUE;

   if (mp)
   {
      /* Conservative estimate of required buffer size 700+90*n_atoms */
      nstereo = 0;
      for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         if (bp->stereo_symbol != NONE) nstereo++;
      bufsize = 700+90*mp->n_atoms+200*nstereo;
      buffer = TypeAlloc(bufsize, char);
      result = MoleculeToWMFBuffer(mp,
                                   buffer, bufsize,
                                   xextp, yextp,
                                   0,
                                   (int *)NULL,
                                   (char **)NULL);
      if (result < -bufsize) /* second try */
      {
         MyFree(buffer);
         bufsize = (-result)+16;
         buffer = TypeAlloc(bufsize, char);
         result = MoleculeToWMFBuffer(mp,
                                      buffer, bufsize,
                                      xextp, yextp,
                                      0,
                                      (int *)NULL,
                                      (char **)NULL);
      }
      if (result <= 0)
      {
         result = FALSE;
      }
      else
      {
         fpout = _lcreat(fname, 0);
         _lwrite(fpout, buffer, result);
         _lclose(fpout);
         result = TRUE;
      }
      MyFree(buffer);

      FreeMolecule(mp);
      return (result);
   }
   else
      return (FALSE);
}

_export int FAR PASCAL WMFDepictCTString(LPSTR ct_string, LPSTR fname,
                                         int *xextp, int *yextp)
/*
 * Writes a WMF depiction of the structure defined by the connection
 * table string ct_string[] to the file named fname[].
 * The fuction returns TRUE if a depiction was created and FALSE otherwise.
 * *xextp and *yextp are set to the dimensions of the resulting picture.
 */
{
   struct reaccs_molecule_t *mp;
   struct reaccs_bond_t     *bp;
   int nstereo, i;
   Fortran_FILE *fp;
   int result;

   char *buffer;
   HFILE fpout;

   if (!ct_string  ||  !fname)
   {
      MessageBox(NULL, "SMILES string or file name was NULL",
                       "WMFDepictCTString", MB_OK);
      return (FALSE);
   }

   result = FALSE;

   fp = FortranStringOpen(ct_string);
   if (fp)
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL == ReadREACCSMolecule(fp,mp,"") && mp->n_atoms > 0)
      {
         /* Conservetive estimate of required buffer size 950+65*n_atoms */
         nstereo = 0;
         for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
            if (bp->stereo_symbol != NONE) nstereo++;
         buffer = TypeAlloc(950+65*mp->n_atoms+100*nstereo, char);
         result = MoleculeToWMFBuffer(mp,
                                      buffer, 950+65*mp->n_atoms+100*nstereo,
                                      xextp, yextp,
                                      USE_COLORS | BINARY_WMF,
                                      (int *)NULL,
                                      (char **)NULL);
         if (result <= 0)
            result = FALSE;
         else
         {
            fpout = _lcreat(fname, 0);
            _lwrite(fpout, buffer, result);
            _lclose(fpout);
            result = TRUE;
         }
         MyFree(buffer);
      }
   }

   FreeMolecule(mp);
   FortranClose(fp);
   return (result);
}

int _export FAR PASCAL CTStringToSmiles(LPSTR ct_string,
                                        LPSTR smiles, int nbuf)
/*
 * Write the isomeric SMILES corresponding to the MOL file in ct_string into
 * the buffer smiles[0..nbuf-1].
 *
 * The function returns the size of the SMILES written, 0 on error, and
 * the negative required size of the buffer if not enough space was provided.
 *
 * Note: The buffer must also provide space for the terminal '\0'.
 */
{
   struct reaccs_molecule_t *mp;
   char *tmp_smiles;
   Fortran_FILE *fp;
   int result;

checkHeap();

   if (!ct_string  ||  !smiles)
   {
      MessageBox(NULL, "CT string or buffer was NULL", "CTStringToSmiles", MB_OK);
      return (FALSE);
   }

   result = FALSE;

   fp = FortranStringOpen(ct_string);
   if (fp)
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL == ReadREACCSMolecule(fp,mp,"") && mp->n_atoms > 0)
      {
         tmp_smiles = MOLToSMI(mp, ISOMERIC_SMILES);
         if (nbuf >= strlen(tmp_smiles)+1)
         {
            strcpy(smiles, tmp_smiles);
            result = strlen(tmp_smiles);
         }
         else
         {
            result = -strlen(tmp_smiles)-1;
         }
         MyFree(tmp_smiles);
      }
   }

   FreeMolecule(mp);
   FortranClose(fp);
   return (result);
}

int _export FAR PASCAL CTStringToSmilesExt(LPSTR ct_string,
                                           LPSTR smiles, int nbuf)
/*
 * Write the isomeric SMILES corresponding to the MOL file in ct_string into
 * the buffer smiles[0..nbuf-1] including the coordinate list after a '|'.
 *
 * The function returns the size of the SMILES written, 0 on error, and
 * the negative required size of the buffer if not enough space was provided.
 *
 * Note: The buffer must also provide space for the terminal '\0'.
 */
{
   struct reaccs_molecule_t *mp;
   char *tmp_smiles;
   Fortran_FILE *fp;
   int result;
   char *coordp;

checkHeap();

   coordp = (char *)NULL;

   if (!ct_string  ||  !smiles)
   {
      MessageBox(NULL, "CT string or buffer was NULL",
                       "CTStringToSmiles", MB_OK);
      return (FALSE);
   }

   result = FALSE;

   fp = FortranStringOpen(ct_string);
   if (fp)
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL == ReadREACCSMolecule(fp,mp,"") && mp->n_atoms > 0)
      {
         // tmp_smiles = MOLToSMI(mp, ISOMERIC_SMILES);
         tmp_smiles = MOLToSMIExt(mp,
                                  ISOMERIC_SMILES|USE_Z,
                                  (int *)NULL,
                                  &coordp);
             if (nbuf >= strlen(tmp_smiles)+1+strlen(coordp)+1)
             {
                strcpy(smiles, tmp_smiles);
                strcat(smiles, "|");
                strcat(smiles, coordp);
                result = strlen(tmp_smiles);
             }
             else
             {
                result = -strlen(tmp_smiles)-1;
         }
         MyFree(tmp_smiles);
      }
   }

   if (coordp)
   {
      MyFree((char *)coordp);
      coordp = NULL;
   }
   FreeMolecule(mp);
   FortranClose(fp);
   return (result);
}

_export int FAR PASCAL DrawMOLFile(LPSTR ct_string)
{
   static HANDLE hmf;
   HDC           hdcMeta, hdc;
   HPEN          hPen;
   HBRUSH        hBrush;
   HFONT         hFont, hFontMeta, hFontScr;
   TEXTMETRIC    tm;
   DWORD         txtExt;
   int           txtHeight, txtWidth;

   Fortran_FILE *fp;
   struct reaccs_molecule_t *mp;
   struct reaccs_bond_t     *bp;
   struct reaccs_atom_t     *ap;

   double xmin, xmax, ymin, ymax, len;
   int i;
   int xExt, yExt;
   int x1, y1, x2, y2;
   int flags1, flags2;
   char chargestring[10];

   fp = FortranStringOpen(ct_string);
   if (fp)
   {
      mp = TypeAlloc(1,struct reaccs_molecule_t);
      if (FORTRAN_NORMAL == ReadREACCSMolecule(fp,mp,"") &&
          mp->n_atoms > 0)
      {
         GetWorldWindow(mp, &xmin, &xmax, &ymin, &ymax, &len);
         SetTransformation(stdbnd, xmin,xmax, ymin,ymax, len);

         hdcMeta = CreateMetaFile(NULL);
         xExt = WorldToMetaX(xmax)-WorldToMetaX(xmin);
         yExt = WorldToMetaY(ymin)-WorldToMetaY(ymax);
         xExt += 1.5*stdbnd; yExt += 1.5*stdbnd;
         SetWindowExt(hdcMeta, xExt, yExt);
         SetWindowOrg(hdcMeta, -xExt/2, -yExt/2);

         hPen   = CreatePen(PS_SOLID,
                            MAX(min_penwidth, fontscale*penwidth),
                            RGB(0, 0, 0));
         SelectObject(hdcMeta, hPen);
         hBrush = CreateSolidBrush(RGB(0,0,0));
         SelectObject(hdcMeta, hBrush);
         SetBkMode(hdcMeta, OPAQUE);

         SetDrawFlags(mp);

         for (i=0, bp=mp->bond_array; i<mp->n_bonds; i++, bp++)
         {
            x1 = WorldToMetaX(mp->atom_array[bp->atoms[0]-1].x);
            y1 = WorldToMetaY(mp->atom_array[bp->atoms[0]-1].y);
            x2 = WorldToMetaX(mp->atom_array[bp->atoms[1]-1].x);
            y2 = WorldToMetaY(mp->atom_array[bp->atoms[1]-1].y);
            flags1 = mp->atom_array[bp->atoms[0]-1].color;
            flags2 = mp->atom_array[bp->atoms[1]-1].color;
            DrawBond(hdcMeta,
                     x1, y1, x2, y2,
                     fontscale*fontsize,
                     bp->color, flags1, flags2);
         }

         hFont = DepictFont(fontscale);
         hFontMeta = DepictFont(fontscale);
         SelectObject(hdcMeta, hFontMeta);
         hdc = CreateDC("DISPLAY", NULL, NULL, NULL);
         hFontScr = SelectObject(hdc, hFont);
         SetTextAlign(hdcMeta, TA_BASELINE);
         GetTextMetrics(hdc, &tm);
         for (i=0, ap=mp->atom_array; i<mp->n_atoms; i++, ap++)
            if (0 != strcmp(ap->atom_symbol,"C")  ||
                ap->charge != NONE                ||
                (ap->color & CUT)                 ||
                ap->radical != NONE               ||
                ap->mass_difference != NONE)
            {
               txtExt = GetTextExtent(hdc, ap->atom_symbol, 1);
               txtWidth  = LOWORD(txtExt);
               txtHeight = HIWORD(txtExt) -
                           tm.tmInternalLeading -
                           tm.tmDescent;
               if (ap->charge          != NONE  ||
                   ap->radical         != NONE ||
                   ap->mass_difference != NONE)
               {
                   SetTextAlign(hdcMeta, TA_BASELINE|TA_UPDATECP);
                   MoveTo(hdcMeta, WorldToMetaX(ap->x) -txtWidth/3.5,
                                   WorldToMetaY(ap->y)+txtHeight/2);
                   TextOut(hdcMeta, WorldToMetaX(ap->x) -txtWidth/3.5,
                                    WorldToMetaY(ap->y)+txtHeight/2,
                                    ap->atom_symbol,
                                    strlen(ap->atom_symbol));
                   AttributesToString(chargestring,
                                      ap->charge, ap->radical,
                                      ap->mass_difference);
                   SetTextAlign(hdcMeta, TA_BOTTOM|TA_UPDATECP);
                   TextOut(hdcMeta, WorldToMetaX(ap->x) -txtWidth/3.5,
                                    WorldToMetaY(ap->y)+txtHeight/2,
                                    chargestring,
                                    strlen(chargestring));
                   SetTextAlign(hdcMeta, TA_BASELINE|TA_NOUPDATECP);
               }
               else
                  TextOut(hdcMeta, WorldToMetaX(ap->x) -txtWidth/3.5,
                                   WorldToMetaY(ap->y)+txtHeight/2,
                                   ap->atom_symbol,
                                   strlen(ap->atom_symbol));
            }
         SelectObject(hdcMeta, GetStockObject(DEVICE_DEFAULT_FONT));
         SelectObject(hdc,     hFontScr);
         DeleteObject(hFontMeta); DeleteObject(hFont);
         DeleteDC(hdc);

         SelectObject(hdcMeta, GetStockObject(BLACK_PEN));
         DeleteObject(hPen);
         SelectObject(hdcMeta, GetStockObject(BLACK_BRUSH));
         DeleteObject(hBrush);
         hmf = CloseMetaFile(hdcMeta);

         SendMetaFileToClipboard(hmf, xExt, yExt);
         LocalCompact(0); GlobalCompact(0);
      }

      if (mp) FreeMolecule(mp);
      if (fp) FortranClose(fp);
      return (TRUE);
   }
   else
   {
      MessageBox(NULL, "Could not open ct_string!", "BOX", MB_OK);
      return (FALSE);
   }
}

static LPSTR data_buffer = NULL;        /* Current data buffer pointer */
static HGLOBAL hbuffer = 0;             /* Current data buffer handle */

LPSTR _export FAR PASCAL GetDataFromClipboard(LPINT sizep, LPSTR format)
/*
 * Return a pointer to the data of the given format or NULL if
 * the format is no available.
 */
{
   HWND hwnd;
   int fmt;
   LPSTR result;
   char buffer[100];
   HGLOBAL hCMem;
   LPSTR lpCMem, cp;

   int data_size, line_size;

checkHeap();

   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   fmt = 0; result = NULL; (*sizep) = 0;
   for (;;)
   {
      fmt = EnumClipboardFormats(fmt);
      if (fmt == 0) break;
      GetClipboardFormatName(fmt, buffer, 99);
      if (strcmp(buffer, format) == 0)    /* format was found */
      {
         if (data_buffer)
         {
            GlobalUnlock(hbuffer);
            GlobalFree(hbuffer);
            hbuffer = 0; data_buffer = NULL;
         }
         hCMem = GetClipboardData(fmt);
         lpCMem = GlobalLock(hCMem);    /* make a pointer to handle's memory */
         data_size = (int)GlobalSize(hCMem);
         hbuffer = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, data_size);
         if (hbuffer == 0)
         {
            hbuffer = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, data_size);
         }
         cp = data_buffer = GlobalLock(hbuffer);
         while (data_size > 0)
         {
            line_size = lpCMem[0];
            lpCMem++; data_size--;
            if (line_size >= 6  &&  0 == strcmp(lpCMem, "M  END"))
            {
               sprintf(cp,"%.*s\n", strlen("M  END"), lpCMem); cp += strlen(cp);
               break;
            }
            sprintf(cp,"%.*s\n",line_size, lpCMem);
            cp += strlen(cp);
            lpCMem += line_size; data_size -= line_size;
         }

         GlobalUnlock(hCMem);
         result = data_buffer; (*sizep) = strlen(data_buffer);
         break;
      }
   }
   CloseClipboard();
   return (result);
}

_export void FAR PASCAL RemoveAllButOne(LPSTR format)
/*
 * Removes all but one data format from the clipboard.
 */
{
   HWND hwnd;
   int fmt;
   HGLOBAL hCMem, hCMem1;
   LPSTR lpCMem, lpCMem1;
   char buffer[100];
   int data_size;
   
   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   fmt = 0;
   for (;;)
   {
      fmt = EnumClipboardFormats(fmt);
      if (fmt == 0) break;
      GetClipboardFormatName(fmt, buffer, 99);
      MessageBox(NULL,
               buffer,
                 "Current Format",
               MB_OK);
      if (strcmp(buffer, format) == 0)  /* format was found */
      {
         MessageBox(NULL, buffer, "Found Format", MB_OK);
         hCMem = GetClipboardData(fmt);
         lpCMem = GlobalLock(hCMem);    /* make a pointer to handle's memory */
         data_size = (int)GlobalSize(hCMem);

         hCMem1 = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, data_size);
         lpCMem1 = GlobalLock(hCMem1);

         while (data_size > 0)
         {
            (*lpCMem1) = (*lpCMem);
            lpCMem1++; lpCMem++; data_size--;
         }
         GlobalUnlock(hCMem); GlobalUnlock(hCMem1);
         EmptyClipboard();
         SetClipboardData(fmt, hCMem1);
break;
      }
   }
   CloseClipboard();
}

void _export FAR PASCAL PostMOLFileToClipboard(LPSTR fname)
/*
 * Reads *fname and posts the contents to the clipboard in
 * MDLCT format.
 */
{
   HWND hwnd;
   int fsize;              /* size of data block to be posted */
   HANDLE molhnd;
   UINT fmt;

   molhnd = FileToMemoryHandle(fname, &fsize);

   fmt = RegisterClipboardFormat("MDLCT");  /* post to the clipboard */
   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   EmptyClipboard();
   SetClipboardData(fmt, molhnd);
   CloseClipboard();

   return;
}

_export void FAR PASCAL PostDataToClipboard(LPSTR data, int size, LPSTR fmtstr)
/*
 * Posts the data contained in data[0..size-1] to the clipboard as format fmtstr.
 * Registers fmt if necessary. The procedure is designed to work with
 * null terminated strings.
 */
{
   HWND hwnd;
   HANDLE datahnd;
   LPSTR databuffer;
   UINT fmt;

   datahnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, size+1);
   databuffer = GlobalLock(datahnd);
   memcpy(databuffer, data, size);
   databuffer[size] = '\0';
   GlobalUnlock(datahnd);

   fmt = RegisterClipboardFormat(fmtstr);
   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   SetClipboardData(fmt, datahnd);
   CloseClipboard();

   return;
}

void _export FAR PASCAL SaveCTFileFromClipboard(LPSTR fname)
/*
 * Fetches the connection table from the clipboard and places it
 * into file *fname. The file is cleared if there is no MDLCT format
 * available.
 */
{
   HWND hwnd;
   char buffer[100];
   int fmt;
   HGLOBAL hCMem;
   LPSTR lpCMem;
   FILE *fp;

   int data_size, line_size;

   fp = fopen(fname,"w");
   if (!fp) return;

   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   fmt = 0;
   for (;;)                  /* loop over all available formats */
   {
      fmt = EnumClipboardFormats(fmt);
      if (fmt == 0) break;
      GetClipboardFormatName(fmt, buffer, 99);
      if (strcmp(buffer, "MDLCT") == 0)  /* MDLCT format was found */
      {
         hCMem = GetClipboardData(fmt);
         lpCMem = GlobalLock(hCMem);    /* make a pointer to handle's memory */
         data_size = (int)GlobalSize(hCMem);

         while (data_size > 1)
         {
            line_size = lpCMem[0];
            lpCMem++; data_size--;
            if (data_size < line_size) break;   /* safeguard */
            strncpy(buffer, lpCMem, line_size); buffer[line_size] = '\0';
            strcat(buffer,"\n");
            fputs(buffer, fp);
            if (data_size > strlen("M  END")+1 &&
                0 == strncmp(lpCMem, "M  END", strlen("M  END"))) break;
            lpCMem += line_size; data_size -= line_size;
         }

         GlobalUnlock(hCMem);
         break;
      }
   }
   CloseClipboard();

   fclose(fp);

   return;
}

int _export FAR PASCAL FetchDataFromClipboard(LPSTR buffer,
                                              int   bufsize,
                                              LPSTR fmtname)
/*
 * Fetches data of the format fmtname from the clipboard and places them
 * into buffer. Returns the actual number of bytes retrieved if successful.
 * If bufsize is not big enough, -1 is returned.
 * If fmtname is not on the clipboard, FetchDataFromClipboard returns 0.
 */
{
   char tmpname[100];
   int result;
   HWND hwnd;
   int fmt;
   HGLOBAL hCMem;
   LPSTR lpCMem;
   long data_size;

   result = 0;
   hwnd = GetActiveWindow();
   OpenClipboard(hwnd);
   fmt = 0;
   for (;;)                               /* loop over all available formats */
   {
      fmt = EnumClipboardFormats(fmt);
      if (fmt == 0) break;
      GetClipboardFormatName(fmt, tmpname, 99);
      if (strcmp(tmpname, fmtname) == 0)  /* requested format was found */
      {
         hCMem = GetClipboardData(fmt);
         lpCMem = GlobalLock(hCMem);
         data_size = GlobalSize(hCMem);
         if (data_size > bufsize)
         {
            result = (-1);
            break;
         }
         memcpy(buffer, lpCMem, (int)data_size);
         result = (int)data_size;

         GlobalUnlock(hCMem);
         break;
      }
   }
   CloseClipboard();

   return (result);
}

void _export FAR PASCAL SMILESStringToMOLFile(LPSTR smiles, LPSTR fname)
/*
 * Converts the SMILES string *smiles to a MOL-file named *fname.
 *
 * The file is cleared in case of error.
 */
{
   FILE *fp;
   struct reaccs_molecule_t FAR *mp;
   int i;

   mp = SMIToMOL(smiles, DO_LAYOUT);
   /* The following is a patch to get the correct sizes into ISIS */
   fp = fopen(fname,"w");
   if (mp)
   {
      for (i=0; i<mp->n_atoms; i++)
      {
         mp->atom_array[i].x *= 0.5;
         mp->atom_array[i].y *= 0.5;
      }

      PrintREACCSMolecule(fp,mp,"");
      FreeMolecule(mp);
   }
   fclose(fp);
   return;
}

void _export FAR PASCAL PostMOLFileEditable(LPSTR fname)
/*
 * Posts (adds) the MOL File fname to the Clipboard in the following formats:
 * - MDLSK as the named data format for ISISDraw, and
 * - OwnerLink ISISServer to inform receivers about the appropriate editor
 * - Native as the data content in a compound document.
 * The presentation picture (and the MDLCT if needed) have to be put onto
 * the Clipboard by the caller.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;

   char *smiles = NULL;
   static LPSTR   smibuffer   = NULL;
   static HGLOBAL smihnd      = NULL;
#ifdef __WIN32__
#define NBUFFER 40000
#else
#define NBUFFER 10000
#endif
   static char buffer[NBUFFER];
   int size;
   int tmp;

   ffp = FortranOpen(fname,"r");  /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
   {
      if (mp) FreeMolecule(mp);
      return;
   }
   FortranClose(ffp);

   size = MoleculeToSKCBuffer(mp, buffer, NBUFFER);
#define ISISLink "ISISServer\0\0NONE\0\0"
#define link_size 18
   if (size > 0)        // Could be rendered => post to Clipboard
   {
      PostDataToClipboard(buffer, size, "MDLSK");
      PostDataToClipboard(buffer, size, "Native");
      PostDataToClipboard(ISISLink, link_size, "OwnerLink");
   }

   if (mp) FreeMolecule(mp);
}

LPSTR _export FAR PASCAL MOLFileToSMILESString(LPINT sizep, LPSTR fname)
/*
 * Converts the MOL file stored in *fname to the corresponding
 * SMILES string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * In case of error, the function returns a NULL pointer.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;

   char *smiles = NULL;
   static LPSTR   smibuffer   = NULL;
   static HGLOBAL smihnd      = NULL;

checkHeap();

   ffp = FortranOpen(fname,"r");  /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
   {
      if (mp) FreeMolecule(mp);
      return (NULL);
   }
   FortranClose(ffp);

   smiles = MOLToSMI(mp, ISOMERIC_SMILES);          /* Convert MOL to SMILES */
   FreeMolecule(mp);

   if (smihnd != NULL)                   /* allocate persitent memory */
   {
      GlobalUnlock(smihnd);
      GlobalFree(smihnd);
      smihnd = NULL;
   }
   if (smiles == NULL)
      smihnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, 1);
   else
      smihnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, strlen(smiles)+1);
   smibuffer = GlobalLock(smihnd);

   if (smiles != NULL)
   {
      strcpy(smibuffer, smiles);            /* copy SMILES */
      (*sizep) = strlen(smiles);
   }
   else
   {
      strcpy(smibuffer, "");
      (*sizep) = 0;
   }

   if (smiles != NULL)
      MyFree((char *) smiles);               /* free temporay storage */

   return (smibuffer);
}

#define UNIQUESMILES   1
#define ADDCOORDINATES 2

LPSTR _export FAR PASCAL MOLFileToSMILESExt(LPINT sizep,
                                            LPSTR fname,
                                            int flags)
/*
 * Converts the MOL file stored in *fname to the corresponding
 * SMILES string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * The parameter flags defines the kind processing required. It
 * can request that stereochemistry be ignored (UNIQUESMILES) and
 * that a comma separated list of coordinates for each atom in the
 * SMILES be appended.
 *
 * In case of error, the function returns a NULL pointer.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;

   char *smiles = NULL;
   static LPSTR   smibuffer   = NULL;
   static HGLOBAL smihnd      = NULL;

   char *coords = NULL;

checkHeap();

   ffp = FortranOpen(fname,"r");  /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,""))
   {
      if (mp) FreeMolecule(mp);
      return (NULL);
   }
   FortranClose(ffp);

   if (flags & UNIQUESMILES)
      smiles = MOLToSMIExt(mp, 0, (int *)NULL, &coords);   /* ignore stereo */
   else
      smiles = MOLToSMIExt(mp, ISOMERIC_SMILES, (int *)NULL, &coords);
   FreeMolecule(mp);

   if (smihnd != NULL) /* Free current persistent memory */
   {
      GlobalUnlock(smihnd);
      GlobalFree(smihnd);
      smihnd = NULL;
   }
   /* allocate persitent memory */
   if (coords && smiles)
   {
      smihnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT,
                           strlen(smiles)+1 +
                           strlen(coords)+1);
      smibuffer = GlobalLock(smihnd);

      strcpy(smibuffer, smiles);            /* copy SMILES */
      (*sizep) = strlen(smiles);

      if (flags & ADDCOORDINATES)
      {
         strcat(smibuffer, " ");
         strcat(smibuffer, coords);            /* append coordinates */
         (*sizep) += 1+strlen(coords);
      }

      MyFree((char *) smiles);               /* free temporay storage */
      MyFree((char *) coords);               /* free temporay storage */
   }
   else
      smibuffer = NULL;

   return (smibuffer);
}

LPSTR _export FAR PASCAL MOLFileToSMARTSString(LPINT sizep, LPSTR fname)
/*
 * Converts the MOL file query stored in *fname to the corresponding
 * SMARTS string and returns a pointer to it. The data becomes
 * private to the DLL. The caller must therefore take care of
 * Copying it to his data area.
 *
 * In case of error, the function returns a NULL pointer.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;
   char *smiles = NULL;

   static LPSTR   smibuffer = NULL;
   static HGLOBAL smihnd = NULL;

checkHeap();

   ffp = FortranOpen(fname,"r");    /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,"")  ||
       mp->n_atoms == 0)
   {
      if (mp) FreeMolecule(mp);
      (*sizep) = 0;
      return (NULL);
   }
   FortranClose(ffp);

   /* Convert MOL to SMARTS */
   smiles = MOLToSMI(mp, ISOMERIC_SMILES | SMARTS_PERCEPTION);
   FreeMolecule(mp);
   if (smiles == NULL  ||  strlen(smiles) == 0)
   {
      (*sizep) = 0;
      return (NULL);
   }

   if (smihnd != NULL)                   /* allocate persitent memory */
   {
      GlobalUnlock(smihnd);
      GlobalFree(smihnd);
      smihnd = NULL;
   }
   smihnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, strlen(smiles)+1);
   smibuffer = GlobalLock(smihnd);

   strcpy(smibuffer, smiles);            /* copy SMILES */
   (*sizep) = strlen(smiles);

   MyFree((char *) smiles);               /* free temporay storage */

   return (smibuffer);
}

void _export FAR PASCAL MOLFileToSMARTSBuffer(LPINT sizep, LPSTR buffer, LPSTR fname)
/*
 * Converts the MOL file query stored in *fname to the corresponding
 * SMARTS string and writes the result to buffer[]. Temporary storage is
 * removed. (*sizep) shall contain the buffer size on input and will
 * become the actual size of the copied SMARTS on output.
 *
 * In case of error, the function sets (*sizep) to -1 or -required_buffer_size
 * if the size of the buffer wasn't sufficient.
 */
{
   Fortran_FILE *ffp;
   struct reaccs_molecule_t *mp;
   char *smiles = NULL;

checkHeap();

   ffp = FortranOpen(fname,"r");    /* read molecule */
   mp = TypeAlloc(1,struct reaccs_molecule_t);
   if (FORTRAN_NORMAL != ReadREACCSMolecule(ffp,mp,"")  ||
       mp->n_atoms == 0)
   {
      if (mp)
      {
         FreeMolecule(mp);
         (*sizep) = -1;
      }
      else
         (*sizep) = -1;
      return;
   }
   FortranClose(ffp);

   /* Convert MOL to SMARTS */
   smiles = MOLToSMI(mp, ISOMERIC_SMILES | SMARTS_PERCEPTION);
   FreeMolecule(mp);
   if (smiles == NULL  ||  strlen(smiles) == 0)
   {
      (*sizep) = -1;
      return;
   }

   if ((*sizep) <= strlen(smiles))
   {
      (*sizep) = -strlen(smiles);
   }
   else
   {
      strcpy(buffer, smiles);            /* copy SMILES */
      (*sizep) = strlen(smiles);
   }

   MyFree((char *) smiles);               /* free temporay storage */

   return;
}

#define CONVERSION_FLOAT      0
#define CONVERSION_DOUBLE      1
#define CONVERSION_INT         2
#define CONVERSION_LONG                3

union number_t
   {
      unsigned char bytestring[8];
      float         single_num;
      double        double_num;
      int           integer_num;
      long          long_num;
   } number_union;

LPSTR _export FAR PASCAL NumberBytesToString(LPINT sizep,
                                             LPSTR bytestring,
                                             int conversion)
/*
 * Converts the number stored as a binary representation in bytestring[]
 * into the corresponding string representation. It uses conversion
 * to select what kind of number is to be converted.
 */
{
   static LPSTR   numbuffer = NULL;
   static HGLOBAL numhnd = NULL;
   char buffer[30];
   int i;

checkHeap();

   switch (conversion)
   {
      case CONVERSION_FLOAT:
         for (i=0; i<4; i++) number_union.bytestring[i] = bytestring[3-i];
         sprintf(buffer,"%g",number_union.single_num);
         break;
     case CONVERSION_DOUBLE:
        for (i=0; i<8; i++) number_union.bytestring[i] = bytestring[7-i];
        sprintf(buffer,"%g",number_union.double_num);
        break;
     case CONVERSION_INT:
        for (i=0; i<2; i++) number_union.bytestring[i] = bytestring[i];
        sprintf(buffer,"%g",number_union.integer_num);
        break;
     case CONVERSION_LONG:
        for (i=0; i<4; i++) number_union.bytestring[i] = bytestring[i];
        sprintf(buffer,"%g",number_union.long_num);
        break;
     default:
        MessageBox(NULL, "Conversion not known", "NumberBytesToString", MB_OK);
        sprintf(buffer,"error (%d)",conversion);
        break;
   }

   if (numhnd != NULL)                        /* allocate persitent memory */
   {
      GlobalUnlock(numhnd);
      GlobalFree(numhnd);
      numhnd = NULL;
   }
   numhnd = GlobalAlloc(GMEM_SHARE | GMEM_ZEROINIT, strlen(buffer)+1);
   numbuffer = GlobalLock(numhnd);

   strcpy(numbuffer, buffer);            /* copy number string */
   (*sizep) = strlen(buffer);

   return (numbuffer);
}
