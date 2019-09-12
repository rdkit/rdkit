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

#include "RDLapsp.h"

#include <stdlib.h>
#include <limits.h>
#include <stdio.h>

#include "RDLgraph.h"
#include "RDLutility.h"

static void RDL_initializeSPathInfo(RDL_sPathInfo *info, RDL_graph *gra)
{
  unsigned i;
  info->pred = RDL_alloc2DUIntArray(gra->V, gra->V);
  info->dist = RDL_alloc2DUIntArray(gra->V, gra->V);
  info->reachable = RDL_alloc2DCharArray(gra->V, gra->V);
  info->dPaths = malloc(gra->V * sizeof(*info->dPaths));
  for(i=0; i<gra->V; ++i)
  {
    info->dPaths[i] = RDL_initNewGraph_g(gra->V, 0);
  }
}

static void RDL_findpaths(RDL_sPathInfo *spi, RDL_graph *gra)
{
  const unsigned INFINITY = UINT_MAX;
  unsigned i,j,w,adj,run;
  unsigned q_head, q_nextfree, q_size;
  char *color = malloc(gra->V * sizeof(*color));
  unsigned *queue = malloc(gra->V * sizeof(*queue));

  for(run=1; run<3; ++run) /*2 run throughs to get Vr correct*/
  {
    for(i=0; i<gra->V; ++i)
    {
      /*initialization*/
      q_head=0;
      q_nextfree=0;
      q_size=0;
      for(j=0; j<gra->V; ++j)
      {
        if(run==1)
        {
          spi->dist[i][j] = INFINITY;
          spi->pred[i][j] = INFINITY;
          spi->reachable[i][j] = 0;
        }
        color[j] = 'w'; /*white*/
      }
      spi->dist[i][i] = 0;
      spi->pred[i][i] = i;
      color[i] = 'b'; /*black*/
      queue[q_nextfree]=i; /*enqueue*/
      ++q_nextfree; /*enqueue*/
      ++q_size; /*enqueue*/

      while(q_size > 0)
      {
        w = queue[q_head]; /*deqeue*/
        ++q_head; /*dequeue*/
        --q_size; /*dequeue*/
        for(adj=0; adj<gra->degree[w]; ++adj) /*for each node adj to w*/
        {
          j=gra->adjList[w][adj][0];
          /*if j precedes i in order (only first run)*/
          if(run==2 || ((gra->degree[j] == gra->degree[i]) && j<i) ||
          gra->degree[j] < gra->degree[i])
          if(color[j] == 'w') /*unvisited*/
          {
            if(run==2)
            {/*if in the 2nd run a dist gets shorter*/
              if(spi->dist[i][w]+1 < spi->dist[i][j])
              {/*not element of Vr*/
                spi->reachable[i][j] = 0;
              }
            }
            if(run==1 || (spi->dist[i][w]+1 < spi->dist[i][j]))
            {/*if 2nd run and dist stays the same, pred shouldn't
            change to keep the shortest path along ordering*/
              /*predecessor of j on a shortest path from i to j*/
              spi->pred[i][j] = w;
            }
            if(spi->dist[i][j] > spi->dist[i][w]+1)
            {
              spi->dist[i][j] = spi->dist[i][w] + 1;
            }
            color[j] = 'b';
            queue[q_nextfree] = j; /*enqueue*/
            ++q_nextfree; /*enqueue*/
            ++q_size; /*enqueue*/
            if(run==1)
            {/*reachable should not change to 1 in the 2nd run*/
              spi->reachable[i][j] = 1;
            }
          }
        }
      }
    }
  }

  free(color);
  free(queue);
}

RDL_sPathInfo *RDL_AllPairsShortestPaths(RDL_graph *gra)
{
  RDL_sPathInfo *info = malloc(sizeof(*info));

  RDL_initializeSPathInfo(info, gra);
  RDL_findpaths(info, gra);

  return info;
}

void RDL_deleteAPSP(RDL_sPathInfo *info, unsigned V)
{
  unsigned i;
  RDL_delete2DUIntArray(info->pred, V);
  RDL_delete2DUIntArray(info->dist, V);
  RDL_delete2DCharArray(info->reachable, V);
  for(i=0; i<V; ++i)
  {
    RDL_deleteGraph(info->dPaths[i]);
  }
  free(info->dPaths);
  free(info);
}
