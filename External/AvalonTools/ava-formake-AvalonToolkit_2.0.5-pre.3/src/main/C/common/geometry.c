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
/*   File:           geometry.c                                         */
/*                                                                      */
/*  Purpose:        Implements the functions to handle geometric        */
/*                  objects like point sets.                            */
/*                                                                      */
/************************************************************************/

#include "geometry.h"

#include <math.h>

#include "local.h"

void TransformPoints(double r[][2], int n,
                     double p1[2],  double p2[2],
                     double p1p[2], double p2p[2])
/*
 * Transforms the 2D points in r[0...n-1][2] such that the points p1[] and
 * p2[] would be transformed to p1p[] and p2p[], resp. in a best fit
 * manner by just rotating and moving them without scaling.
 */
{
   double rp[2];              /* r prime                              */
   double dp[2], dpp[2];      /* delta p and delta p prime            */
   double mp[2], mpp[2];      /* mean p and mean p prime              */
   double a, b, c;            /* parameters of rotation matrix        */
   double ls, lsp;            /* distances between source points and  */
                              /* destination points, resp.            */
   int i, j;

   ls = lsp = 0.0;
   for (i=0; i<2; i++)          /* compute transformation parameters    */
   {
      dp[i] = p1[i]-p2[i];       dpp[i] = p1p[i]-p2p[i];
      mp[i] = (p1[i]+p2[i])/2.0; mpp[i] = (p1p[i]+p2p[i])/2.0;
      ls += dp[i]*dp[i];         lsp += dpp[i]*dpp[i];
   }
   if (ls < 0.00001  ||  lsp < 0.00001)
   {
      a = 1.0; b = 0.0; c = 1.0;
   }
   else
   {
      a = (dp[0]*dpp[0] + dp[1]*dpp[1])/ls;
      b = (dpp[0]*dp[1] - dp[0]*dpp[1])/ls;
      c = sqrt(lsp/ls);
   }

   for (j=0; j<n; j++)              /* transform points                     */
   {
      for (i=0; i<2; i++) /* move center of gravity */
         rp[i] = r[j][i]-mp[i];

      r[j][0] =  a*rp[0] + b*rp[1];  /* rotate point */
      r[j][1] = -b*rp[0] + a*rp[1];
      for (i=0; i<2; i++) rp[i] = r[j][i]/c;

      for (i=0; i<2; i++)      /* restore new center of gravity */
      r[j][i] = rp[i]+mpp[i];
   }
}

void PointSetMatchTransformation(double points[][2], unsigned npoints,
                                 double from[][2],
                                 double to[][2],     unsigned nmatch,
                                 int reflection)
/*
 * Computes the transformation which rotates and translates from[0..nmatch-1]
 * such that the points with the same index in from[] and to[] come as
 * close together as possible.  The resulting transformation is applied to
 * points[0..npoints-1]. If reflection is TRUE, the transformation will be
 * improper.
 */
{
   double dx1, dy1, dx2, dy2;
   double calpha, salpha;
   double sxxp, syyp, sxyp, sxpy;
   double denominator;
   double flipx;
   int i;
   double xh, yh, xhh, yhh;

   if (reflection) flipx = (-1.0);
   else                flipx = 1.0;

   dx1 = dy1 = dx2 = dy2 = 0.0;
   for (i=0; i<nmatch; i++)
   {
      dx1 += flipx*from[i][0]; dy1 += from[i][1];
      dx2 += to[i][0];   dy2 += to[i][1];
   }
   dx1 /= nmatch; dy1 /= nmatch;
   dx2 /= nmatch; dy2 /= nmatch;

   sxxp = syyp = sxyp = sxpy = 0.0;
   for (i=0; i<nmatch; i++)
   {
      sxxp += (flipx*from[i][0]-dx1)*(to[i][0]-dx2);
      syyp += (from[i][1]-dy1)*(to[i][1]-dy2);
      sxyp += (flipx*from[i][0]-dx1)*(to[i][1]-dy2);
      sxpy += (from[i][1]-dy1)*(to[i][0]-dx2);
   }

   denominator = sqrt((sxxp+syyp)*(sxxp+syyp) + (sxpy-sxyp)*(sxpy-sxyp));
   if (denominator < 1.0e-5)
   {
      calpha = 1.0; salpha = 0.0;
   }
   else
   {
      calpha = (sxxp + syyp)/denominator;
      salpha = (sxpy - sxyp)/denominator;
   }
   if (ABS(salpha) < 1.0e-10) salpha = 0.0;
   if (ABS(calpha) < 1.0e-10) calpha = 0.0;

   for (i=0; i<npoints; i++)
   {
      xh = flipx*points[i][0]-dx1; yh = points[i][1]-dy1;
      xhh = calpha*xh + salpha*yh;
      yhh = -salpha*xh + calpha*yh;
      points[i][0] = xhh+dx2; points[i][1] = yhh+dy2;
   }
}

typedef double coords[2];
#define LEN(t) sqrt((t)[0]*(t)[0]+(t)[1]*(t)[1])
#define DIST(x1,y1,x2,y2) \
        sqrt(((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)))
#define DIST2(x1,y1,x2,y2) \
        (((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)))

int NextGraphNeighbour(int node,
                       unsigned edges[][2], unsigned nedges,
                       int last_neighbour)
/*
 * Scans the graph defined by edges[0..nedges-1] to find the next
 * neighbour of node. The node last_neighbour was the last node reported.
 * If there is no such node, (-1) is reported.
 */
{
   int result, i, j, maxnode;

   maxnode = 0;
   for (i=0; i<nedges; i++)
      for (j=0; j<2; j++)
         if (edges[i][j] > maxnode) maxnode = edges[i][j];
   result = maxnode+1;
   for (i=0; i<nedges; i++)
   {
      if (edges[i][0] == node  &&  edges[i][1] < result)
         result = edges[i][1];
      if (edges[i][1] == node  &&  edges[i][0] < result)
         result = edges[i][0];
   }

   if (result == maxnode+1 || result <= last_neighbour)
      return (-1);
   else
      return (result);
}

int HasCusp(double m[2], double neigh[][2], int nneigh)
/*
 * Tests if the neighbour coordinates neigh[0..nneigh-1][0..1] form a cusp
 * or umbrella around m[0..1];
 */
{
   int i, j;
   double ti[2], tj[2];      /* temporary points */
   double tb[2];             /* direction of bisection of widest angle */
   double len, cosa, mincosa;

// fprintf(stderr, "nneigh = %d\n", nneigh);
   /* need at least three bonds for real cusp */
   if (nneigh <= 1) return (FALSE);
   if (nneigh == 2) return (TRUE);

   mincosa = 2;
   /* first, search for pair of neighbours with largest angle between them */
   for (i=0; i<nneigh; i++)
   {
      ti[0] = neigh[i][0]-m[0]; ti[1] = neigh[i][1]-m[1];
      len = LEN(ti); ti[0]/=len; ti[1]/=len;
      for (j=i+1; j<nneigh; j++)
      {
         tj[0] = neigh[j][0]-m[0]; tj[1] = neigh[j][1]-m[1];
         len = LEN(tj); tj[0]/=len; tj[1]/=len;
         cosa = ti[0]*tj[0] + ti[1]*tj[1];
// fprintf(stderr, "cosa = %g\n", cosa);
         if (cosa < -0.999) return (TRUE);
         if (cosa < mincosa)
         {
            mincosa = cosa;
// fprintf(stderr, "mincosa = %g\n", mincosa);
            tb[0] = (ti[0]+tj[0])/2.0;
            tb[1] = (ti[1]+tj[1])/2.0;
            len = LEN(tb); tb[0]/=len; tb[1]/=len;
         }
      }
   }

   /* Check if all bonds point along bisection of widest angle */
   for (i=0; i<nneigh; i++)
   {
      ti[0] = neigh[i][0]-m[0]; ti[1] = neigh[i][1]-m[1];
      len = LEN(ti); ti[0]/=len; ti[1]/=len;
      cosa = (ti[0]*tb[0]+ti[1]*tb[1]);
// fprintf(stderr, "cosa = %g\n", cosa);
      if (cosa < 0.0) return (FALSE);
   }
   return (TRUE);
}

static int debug_flag = FALSE;
void DebugGeometry(int flag)
{
   debug_flag = flag;
}

int GraphToMolecule(int graphIndex, int numbers[], int natoms)
/*
 * Translates from graphIndex to molecule atom number using the index translation table numbers[].
 * This function is used for debugging purposed.
 *
 * Note: numbers uses index origin 0, while atom numbers are IO 1.
 */
{
    int i;
    for (i=0; i<natoms; i++)
        if (graphIndex == numbers[i])
            return (i+1);
    return (-1);
}

int PointsIntoRing(int seed, double p1[2], 
                   double coordinates[][2], unsigned nnodes,
                   unsigned edges[][2], unsigned nedges,
                   int numbers[], int natoms)
/*
 * Checks if the vector coordinates[seed] -> p1 crosses one of the non-adjacent edges.
 * Usually, this means that it points into a ring.
 *
 * numbers[] contains the mapping from molecule to graph numbering.
 */
{
    int i;
    double *p0, *pp0, *pp1;
    double denominator, lambda, mu;

    p0 = coordinates[seed];

    for (i=0; i<nedges; i++)
    {
        if (edges[i][0] == seed  ||  edges[i][1] == seed) continue;
        pp0 = coordinates[edges[i][0]];
        pp1 = coordinates[edges[i][1]];
        denominator = (pp1[0]-pp0[0])*(p0[1]-p1[1])-(pp1[1]-pp0[1])*(p0[0]-p1[0]);
        if (ABS(denominator) < 1.0e-7) continue;
        lambda = (p0[0]-pp0[0])*(p0[1]-p1[1])-(p0[1]-pp0[1])*(p0[0]-p1[0]);
        lambda /= denominator;
        if (0.0 > lambda  ||  lambda > 1.0) continue;   // ray does not cross this bond
        mu = (pp1[0]-pp0[0])*(p0[1]-pp0[1])-(p0[0]-pp0[0])*(pp1[1]-pp0[1]);
        mu /= denominator;
        if (mu <= 0.0) continue;                        // crosses at back end of ray
// fprintf(stderr, "ray starting at atom %d crosses bond %d-%d\n", GraphToMolecule(seed, numbers, natoms), GraphToMolecule(edges[i][0], numbers, natoms), GraphToMolecule(edges[i][1], numbers, natoms));
        return TRUE;
    }

    return FALSE;
}

void NextSubstituentPoint(double point[2],
                          double coordinates[][2], unsigned nnodes,
                          unsigned edges[][2], unsigned nedges,
                          unsigned seed,
	                  int flags,
                          int numbers[], int natoms)
/*
 * Finds the next best point to substitute the graph defined by
 * coords[0..nnodes-1][2] and edges[0..nedges-1][2] at node seed.
 * The coordinates found are stored in points[0..1].
 * The parameter flags tells whether candidated positions may 
 * also be put between existing substituents and if a hydrogen needs to be layed out.
 *
 * numbers[0..natoms-1] contains the mapping from molecule to graph numbering.
 */
{
   int i, j;
   double value, dist, len, vtmp, vtmp_local, cos_alpha;
   double m[2];        /* coordinates of seed */
   double (*neigh)[2]; /* coordinates of neighbours of seed */
   int nneigh;
   int *neigh_degree, *neigh_index;
   double (*cand_points)[2];       /* coordinates of candidates for point */
   int ncand;
   double t[2], th[2];            /* temporary points */
   int has_cusp, points_into_ring;

   static int    nphi=2;  /* directions of sprouting */
   //                           120 deg   -120 deg   45 deg    -45 deg    60 deg    -60 deg
   static double sinphi[2*2] = {0.866025, -0.866025, 0.707107, -0.707107};
   static double cosphi[2*2] = {-0.5,     -0.5,      0.707107,  0.707107};
   int nphi_used = FALSE;

                                       /* get # of neighbours of seed */
   m[0] = coordinates[seed][0];
   m[1] = coordinates[seed][1];
   for (i=0, nneigh=0; i<nedges; i++)
      if (edges[i][0] == seed  ||  edges[i][1] == seed) nneigh++;

   if (nneigh == 0)                        /* trivial boundary cases */
   {
      neigh = NULL;
      neigh_degree = NULL;
      cand_points = TypeAlloc(4*nphi+2, coords);
      ncand=0;
      th[0] = 1.514;
      // shorten hydrogen bonds a bit
      if (flags & H_GEOMETRY) th[0] *= 0.666;
      th[1] = 0.0;
      cand_points[ncand][0] = th[0] + m[0];
      cand_points[ncand][1] = th[1] + m[1];
      ncand++;
      cand_points[ncand][0] = -th[0] + m[0];
      cand_points[ncand][1] = -th[1] + m[1];
      ncand++;

      for (i=0; i<2*nphi; i++)
      {
         cand_points[ncand][0] = m[0] + th[0]*cosphi[i] - th[1]*sinphi[i];
	 cand_points[ncand][1] = m[1] + th[0]*sinphi[i] + th[1]*cosphi[i];
         ncand++;
      }

      for (i=0; i<2*nphi; i++)
      {
         cand_points[ncand][0] = m[0] - th[0]*cosphi[i] + th[1]*sinphi[i];
	 cand_points[ncand][1] = m[1] - th[0]*sinphi[i] - th[1]*cosphi[i];
         ncand++;
      }
      nphi_used = TRUE;
   }
   else
   {
      neigh = TypeAlloc(nneigh, coords);     /* set-up neighbourhood */
      neigh_degree = TypeAlloc(nneigh, int);
      neigh_index = TypeAlloc(nneigh, int);
      for (i=0, nneigh=0; i<nedges; i++)
      {
        if (edges[i][0] == seed)
        {
           neigh[nneigh][0] = coordinates[edges[i][1]][0];
           neigh[nneigh][1] = coordinates[edges[i][1]][1];
	   neigh_index[nneigh] = edges[i][1];
           nneigh++;
        }
        if (edges[i][1] == seed)
        {
           neigh[nneigh][0] = coordinates[edges[i][0]][0];
           neigh[nneigh][1] = coordinates[edges[i][0]][1];
	   neigh_index[nneigh] = edges[i][1];
           nneigh++;
        }
      }
      for (i=0; i<nedges; i++)
	 for (j=0; j<nneigh; j++)
	 {
	    if (edges[i][0] == neigh_index[j]) neigh_degree[j]++;
	    if (edges[i][1] == neigh_index[j]) neigh_degree[j]++;
	 }

                                                  /* set-up candidates */
      cand_points = TypeAlloc(nphi+(nneigh-1)*nneigh, coords);
      ncand=0;

      for (i=0, len=0.0; i<nneigh; i++)   /* compute average edge length */
      {
         len += DIST(m[0],m[1], neigh[i][0], neigh[i][1]);
      }
      len = len/nneigh;
      dist = DIST(m[0],m[1],neigh[0][0],neigh[0][1]);

      /* add angle candidates */
      if (dist > 0.001  &&  nneigh <= 2)
      {
         th[0] = (neigh[0][0]-m[0])*len/dist;
         th[1] = (neigh[0][1]-m[1])*len/dist;
         for (i=0; i<nphi; i++)
         {
            cand_points[ncand][0] = th[0]*cosphi[i] - th[1]*sinphi[i] + m[0];
            cand_points[ncand][1] = th[0]*sinphi[i] + th[1]*cosphi[i] + m[1];
            ncand++;
         }
         nphi_used = TRUE;
      }

      has_cusp = HasCusp(m, neigh, nneigh);
      for (i=0; i<nneigh; i++)   /* add vector sum and diff. candidates */
         for (j=i+1; j<nneigh; j++)
         {
            t[0] = neigh[i][0]-m[0]; t[1] = neigh[i][1]-m[1];
            dist = sqrt(t[0]*t[0]+t[1]*t[1]);
            if (dist > 0.001) {t[0] /= dist; t[1] /= dist;}
            else              continue;
            th[0] = neigh[j][0]-m[0]; th[1] = neigh[j][1]-m[1];
            dist = sqrt(th[0]*th[0]+th[1]*th[1]);
            if (dist > 0.001) {th[0] /= dist; th[1] /= dist;}
            else              continue;
            cos_alpha = t[0]*th[0]+t[1]*th[1];
            t[0] = t[0]+th[0];
            t[1] = t[1]+th[1];
            dist = sqrt(t[0]*t[0]+t[1]*t[1]);
	    if (dist > 0.001)
            {
               if (has_cusp)
               {
                  /* pointing outwards */
                  cand_points[ncand][0] = m[0]-t[0]*len/dist;
                  cand_points[ncand][1] = m[1]-t[1]*len/dist;
// if (debug_flag) fprintf(stderr, "%d: Cusp candidate %g/%g\n", ncand, cand_points[ncand][0], cand_points[ncand][1]);
                  ncand++;
               }
	       /* pointing inwards */
	       if (( (flags & USE_INWARDS) != 0) &&
                   nneigh>2 &&
                   cos_alpha > -0.9)    /* only if not almost 180 degrees */
	       {
		  cand_points[ncand][0] = m[0]+t[0]*len/dist;
		  cand_points[ncand][1] = m[1]+t[1]*len/dist;
// if (debug_flag) fprintf(stderr, "%d: Non-cusp candidate %g/%g, dist=%g\n", ncand, cand_points[ncand][0], cand_points[ncand][1], dist);
		  ncand++;
	       }
	       else
	       {
		  cand_points[ncand][0] = m[0]+t[0]*len/dist;
		  cand_points[ncand][1] = m[1]+t[1]*len/dist;
// if (debug_flag) fprintf(stderr, "%d: Non candidate %g/%g, dist=%g, cos_alpha=%g\n", ncand, cand_points[ncand][0], cand_points[ncand][1], dist, cos_alpha);
	       }
            }
         }

// if (debug_flag) fprintf(stderr, "%d candidates\n", ncand);
   }

   value = 1.0e7;             /* find best candidate */
   for (i=0; i<ncand; i++)
   {
// if (debug_flag) fprintf(stderr, "0: nnodes = %d\n", nnodes);
      points_into_ring = FALSE;
      has_cusp = HasCusp(m, neigh, nneigh);
      if (flags|NO_RAY_CROSSING  && !has_cusp)
         points_into_ring = PointsIntoRing(seed, cand_points[i], coordinates, nnodes, edges, nedges, numbers, natoms);
      // repulsion by neighbour degree weighted force
      vtmp_local = 0.0;
      if (neigh != NULL)
	 for (j=0; j<nneigh; j++)
	    vtmp_local +=
	       neigh_degree[j] *
	       1.0/(0.1+DIST2(cand_points[i][0],cand_points[i][1],
			      neigh[j][0],neigh[j][1]));

      for (j=0, vtmp=0.0; j<nnodes; j++)
         vtmp +=
            1.0/(0.1+DIST2(cand_points[i][0],cand_points[i][1],
                           coordinates[j][0],coordinates[j][1]));

// if (debug_flag) fprintf(stderr, "1: Candidate %g/%g has value %g|%g|%d\n", cand_points[i][0], cand_points[i][1], vtmp, vtmp_local, has_cusp);
      if (i < nphi  &&  nphi_used)
         vtmp *= 1.85;        /* prefer bisection candidates */

      /* make sure that cusps are resolved */
      if (has_cusp) vtmp += 100*vtmp_local;
      else          vtmp += 10*vtmp_local;

      if (points_into_ring) vtmp *= 10.0;
// if (debug_flag) fprintf(stderr, "2: Candidate %g/%g has value %g\n", cand_points[i][0], cand_points[i][1], vtmp);
      if (vtmp < value)
      {
         value = vtmp;
         point[0] = cand_points[i][0];
         point[1] = cand_points[i][1];
      }
   }
   if (neigh != NULL)
   {
      MyFree((char *)neigh);
      MyFree((char *)neigh_degree);
      MyFree((char *)neigh_index);
   }

   MyFree((char *)cand_points);
}

void TransformGraphSubstituent(double scoords[][2], unsigned nsnodes,
                               unsigned sedges[][2], unsigned nsedges,
                               unsigned apo,
                               double mcoords[][2], unsigned nmnodes,
                               unsigned medges[][2], unsigned nmedges,
                               unsigned seed)
/*
 * The subtituent graph definded by the node positions scoords[0..nsnodes][]
 * and the edges sedges[0..nsedges-1][] is rotated and moved such that it
 * when linked to the graph (mcoords/medges) between the nodes sapo and
 * mapo, respectively, does not overlap too much with it.
 */
{
   int i, j, k, h;
   int igood;
   double value, dist, vtmp;
   double m1[2], m2[2];     /* coordinates of seed and its neighbour  */
   double s1[2], s2[2];     /* coordinates of apo and one neighbour   */
   double (*r)[2];          /* coordinatesof subst. to be transformed */
   double t[2], th[2];      /* temporary points */

   static int    nphi=4;  /* directions of sprouting */
   static double sinphi[4] = {0.866025, -0.866025, 0.866025, -0.866025};
   static double cosphi[4] = {-0.5,     -0.5,      0.5,      0.5};

   r = TypeAlloc(nsnodes,coords);
   
   if (nsnodes <= 1  ||  nmnodes <= 1)       /* trivial boundary cases */
   {
      for (i=0; i<nsnodes; i++)
      {
   scoords[i][0] -= scoords[apo][0] - mcoords[seed][0];
    scoords[i][1] -= scoords[apo][1] - mcoords[seed][1];
      }
      return;
   }

                                   /* get molecule reference points */
   m1[0] = mcoords[seed][0];
   m1[1] = mcoords[seed][1];
   h = NextGraphNeighbour(seed,medges,nmedges,-1);
   if (h == (-1))
   {
      m2[0] = m1[0]+1.5; m2[1] = m1[1];
   }
   else
   {
      m2[0] = mcoords[h][0];
      m2[1] = mcoords[h][1];
   }
                                    /* get substituent reference points */
   s1[0] = scoords[apo][0];
   s1[1] = scoords[apo][1];
   h = NextGraphNeighbour(apo,sedges,nsedges,-1);
   if (h == (-1))
   {
      s2[0] = s1[0]+1.5; s2[1] = s1[1];
   }
   else
   {
      s2[0] = scoords[h][0];
      s2[1] = scoords[h][1];
   }
   
   igood = 0; value = 1.0e7;
   dist = DIST2(s1[0],s1[1],s2[0],s2[1]);
   for (i=0; i<nphi; i++)            /* find good orientation */
   {
      th[0] = s2[0]-s1[0]; th[1] = s2[1]-s1[1];
      t[0] = th[0]*cosphi[i] - th[1]*sinphi[i] + s1[0];
      t[1] = th[0]*sinphi[i] + th[1]*cosphi[i] + s1[1];
      
      for (j=0; j<nsnodes; j++)
         if (j != apo  &&
             dist*0.01 > DIST2(t[0],t[1],scoords[j][0],scoords[j][1]))
     break;
      if (j != nsnodes) continue;     /* collision within substituent */
      
      for (j=0; j<nsnodes; j++) /* determine value wrt. molecule */
      {
       r[j][0] = scoords[j][0]; r[j][1] = scoords[j][1];
      }
      TransformPoints(r,nsnodes,s1,t,m1,m2);
      vtmp = 0.0;
      for (j=0; j<nsnodes; j++)
         for (k=0; k<nmnodes; k++)
            vtmp += 1.0/(0.1+DIST2(r[j][0], r[j][1],
                                mcoords[k][0], mcoords[k][0]));
      if (vtmp < value)
      {
         igood = i; value = vtmp;
      }
   }
                                    
   th[0] = s2[0]-s1[0]; th[1] = s2[1]-s1[1];    /* apply good transformation */
   t[0] = th[0]*cosphi[igood] - th[1]*sinphi[igood] + s1[0];
   t[1] = th[0]*sinphi[igood] + th[1]*cosphi[igood] + s1[1];
   for (j=0; j<nsnodes; j++)
   {
      r[j][0] =  scoords[j][0]; r[j][1] = scoords[j][1];
   }
   TransformPoints(r,nsnodes,s1,t,m1,m2);
   for (j=0; j<nsnodes; j++)
   {
      scoords[j][0] = r[j][0];
      scoords[j][1] = r[j][1];
   }

   MyFree((char *)r);
}
