/*
 * This file is part of the RingDecomposerLib, licensed
 * under BSD New license (see LICENSE in the root directory).
 * Copyright (c) 2016
 * University of Hamburg, ZBH - Center for Bioinformatics
 * Niek Andresen, Florian Flachsenberg, Matthias Rarey
 * 
 * Please cite:
 * 
 * Kolodzik, A.; Urbaczek, S.; Rarey, M.
 * Unique Ring Families: A Chemically Meaningful Description
 * of Molecular Ring Topologies.
 * J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
 * 
 * Flachsenberg, F.; Andresen, N.; Rarey, M.
 * RingDecomposerLib: An Open-Source Implementation of
 * Unique Ring Families and Other Cycle Bases.
 * J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
 */

#include "RDLcycleFams.h"

#include <stdlib.h>
#include <limits.h>
#include <stdio.h>

#include "RDLapsp.h"
#include "RDLgraph.h"

/** returns 1 if the intersection of P(r,y) and P(r,z) is equal to {r};
0 otherwise */
static int RDL_pathsShareOnlyStart(unsigned r, unsigned y, unsigned z, RDL_graph *gra, RDL_sPathInfo *spi)
{
  unsigned result = 0, i, pnt, count=0;
  unsigned *vertInRY, *vertInRZ; /*edges in P(r,y) and P(r,z)*/
  vertInRY = malloc(gra->V * sizeof(*vertInRY));
  vertInRZ = malloc(gra->V * sizeof(*vertInRZ));

  for(i=0; i<gra->V; ++i)
  {
    vertInRY[i] = 0;
    vertInRZ[i] = 0;
  }
  vertInRY[y] = 1;
  vertInRZ[z] = 1;
  /*find out which vertices are on the path P(r,y)*/
  pnt = y;
  do
  {
    pnt = spi->pred[r][pnt];
    vertInRY[pnt] = 1;
  }while(pnt != r);
  /*find out which vertices are on the path P(r,z)*/
  pnt = z;
  do
  {
    pnt = spi->pred[r][pnt];
    vertInRZ[pnt] = 1;
  }while(pnt != r);
  /*find out if more than one vertex is shared by the paths*/
  for(i=0; i<gra->V && count<2; ++i)
  {
    if(vertInRY[i] == 1 && vertInRZ[i] == 1)
    {
      ++count;
    }
  }
  /*if r is the only shared vertex of the two paths return 1*/
  if(count == 1 && (vertInRY[r] == 1) && (vertInRZ[r] == 1))
  {
    result = 1;
  }

  free(vertInRY);
  free(vertInRZ);
  return result;
}

/** returns a cycle vector (element of {0,1}^m). odd (x=UINT_MAX) or even cycle.*/
static char *RDL_findPrototype(unsigned r, unsigned y, unsigned z, unsigned x,
    RDL_graph *gra, RDL_sPathInfo *spi)
{
  unsigned i, vert1, vert2;
  char *proto;

  proto = malloc(gra->E * sizeof(*proto));
  for(i=0; i<gra->E; ++i)
  {
    proto[i] = 0;
  }
  /*path from r to y*/
  vert1 = y;
  do
  {
    vert2 = vert1;
    vert1 = spi->pred[r][vert1];
    proto[RDL_edgeId(gra, vert1, vert2)] = 1;
  }while(vert1 != r);
  /*path from r to z*/
  vert1 = z;
  do
  {
    vert2 = vert1;
    vert1 = spi->pred[r][vert1];
    proto[RDL_edgeId(gra, vert1, vert2)] = 1;
  }while(vert1 != r);
  if(x == UINT_MAX)/*odd cycle*/
  {
    proto[RDL_edgeId(gra,y,z)] = 1;
  }
  else /*even cycle*/
  {
    proto[RDL_edgeId(gra,y,x)] = 1;
    proto[RDL_edgeId(gra,z,x)] = 1;
  }
  return proto;
}

/** fills the rc datastructure with the odd cycle r-y-z-r */
static void RDL_addOdd(unsigned r, unsigned y, unsigned z,
    RDL_graph *gra, RDL_sPathInfo *spi, RDL_cfURF *rc)
{
  RDL_cfam *new;
  if (rc->alloced == rc->nofFams) {
    rc->alloced *= 2;
    rc->fams = realloc(rc->fams, rc->alloced * sizeof(*rc->fams));
  }
  /* skip if we're out of memory... */
  if (rc->fams) {
    new = malloc(sizeof(*new));
    new->r = r;
    new->p = y;
    new->q = z;
    new->x = UINT_MAX; /*odd cycle*/
    new->mark = 0;
    new->prototype = RDL_findPrototype(r, y, z, UINT_MAX, gra, spi);
    new->weight = spi->dist[r][y] + spi->dist[r][z] + 1;
    rc->fams[rc->nofFams++] = new;
  }
}

/** fills the rc datastructure with the even cycle r-y-x-z-r */
static void RDL_addEven(unsigned r, unsigned y, unsigned x,
    unsigned z, RDL_graph *gra, RDL_sPathInfo *spi, RDL_cfURF *rc)
{
  RDL_cfam *new;
  if (rc->alloced == rc->nofFams) {
    rc->alloced *= 2;
    rc->fams = realloc(rc->fams, rc->alloced * sizeof(*rc->fams));
  }
  /* skip if we're out of memory... */
  if (rc->fams) {
    new = malloc(sizeof(**rc->fams));
    new->r = r;
    new->p = y;
    new->q = z;
    new->x = x; /*even cycle*/
    new->mark = 0;
    new->prototype = RDL_findPrototype(r, y, z, x, gra, spi);
    new->weight = spi->dist[r][y] + spi->dist[r][z] + 2;
    rc->fams[rc->nofFams++] = new;
  }
}

/** finds a number of cycle families that contain at least all RELEVANT cycle
families - just like in Vismara's pseudocode */
static void RDL_vismara(RDL_cfURF *rc, RDL_graph *gra, RDL_sPathInfo *spi)
{
  unsigned i,j;
  unsigned rv,yv,zv,pv,qv; /*variables as in Vismara's algorithm, extended by a 'v'*/
  unsigned *evenCand; /*'S' in Vismara's algorithm*/
  unsigned nofCandidates = 0;
  evenCand = malloc(gra->V * sizeof(*evenCand));

  for(rv = 0; rv < gra->V; ++rv) {
    for(yv = 0; yv < gra->V; ++yv) {
      /*all yv reachable from rv respecting the ordering*/
      if(spi->reachable[rv][yv] == 1) {
        nofCandidates = 0;
        for(i = 0; i < gra->degree[yv]; ++i) {
          zv = gra->adjList[yv][i][0];
          /*all zv reachable from rv respecting the ordering and adjacent to yv*/
          if(spi->reachable[rv][zv] == 1) {
            if(spi->dist[rv][zv] + 1 == spi->dist[rv][yv]) {
              evenCand[nofCandidates] = zv;
              ++nofCandidates;
            }
            else if(spi->dist[rv][zv] != spi->dist[rv][yv] + 1
                && (gra->degree[zv] < gra->degree[yv] ||
                    (gra->degree[zv] == gra->degree[yv] && zv<yv))
                && RDL_pathsShareOnlyStart(rv, yv, zv, gra, spi) == 1) {
              /*add odd cycle rv-yv rv-zv zv-yv*/
              RDL_addOdd(rv, yv, zv, gra, spi, rc);
            }
            /*to fill dPaths in sPathInfo with the edges to r*/
            if(spi->dist[rv][zv] == 1) {
              RDL_addEdge(spi->dPaths[rv], zv, rv);
            }
          }
        }
        /*any pair in evenCand*/
        for(i = 0; i < nofCandidates; ++i) {
          pv = evenCand[i];
          for(j = i+1; j < nofCandidates; ++j) {
            qv = evenCand[j];
            if((RDL_pathsShareOnlyStart(rv, pv, qv, gra, spi) == 1)) {
              /*add even cycle rv-pv rv-qv pv-yv-qv*/
              RDL_addEven(rv, pv, yv, qv, gra, spi, rc);
            }
          }
        }
        /*to fill dPaths in sPathInfo/fill U_r (see Vismara)*/
        for(i = 0; i < nofCandidates; ++i) {
          pv = evenCand[i];
          RDL_addEdge(spi->dPaths[rv], yv, pv);
        }
      }
    }
  }

  free(evenCand);
}

static int RDL_cycleFamsComp(const void *cf1, const void *cf2)
{
    if((*((RDL_cfam **)cf1))->weight < (*((RDL_cfam **)cf2))->weight)
    {
        return -1;
    }
    else if((*((RDL_cfam **)cf1))->weight > (*((RDL_cfam **)cf2))->weight)
    {
        return 1;
    }
    else return 0;
}

RDL_cfURF *RDL_findCycleFams(RDL_graph *gra, RDL_sPathInfo *spi)
{
  RDL_cfURF *rc = malloc(sizeof(*rc));
  /*number of CFs is at most 2m^2+vn (Vismara Lemma 3)*/

  rc->alloced = 64;
  rc->fams = malloc(rc->alloced * sizeof(*rc->fams));
  rc->nofFams = 0;
  RDL_vismara(rc, gra, spi);

  if (!rc->fams) {
    RDL_outputFunc(RDL_ERROR, "Graph is too large, can't allocate memory!\n");
    free(rc);
    return NULL;
  }

  rc->alloced = rc->nofFams;
  rc->fams = realloc(rc->fams, rc->alloced * sizeof(*rc->fams));

  /* sort by size */
  qsort(rc->fams, rc->nofFams, sizeof(*rc->fams), RDL_cycleFamsComp);
  return rc;
}

void RDL_deleteCycleFams(RDL_cfURF *rc)
{
  unsigned i;
  for(i=0; i<rc->nofFams; ++i)
  {
    free(rc->fams[i]->prototype);
    free(rc->fams[i]);
  }
  free(rc->fams);
  free(rc);
}

