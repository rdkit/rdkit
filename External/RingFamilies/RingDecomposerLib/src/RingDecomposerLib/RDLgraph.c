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

#include "RDLgraph.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "RDLutility.h"

/** initializes a new graph with the value |V| and the array degree, then
  allocates enough space for 'adjList', so that they can be filled.
  E is set to 0, edges are added later.*/
static void RDL_initGraph(RDL_graph *gra, unsigned V, unsigned *degree, char owns_edges)
{
  unsigned i;
  RDL_edge **adjList;

  gra->V = V;
  gra->E = 0;
  gra->degree = degree;

  adjList = malloc(V * sizeof(*adjList));
  for(i=0; i<V; ++i)
  {
    if(degree[i] > 0)
    {
      adjList[i] = malloc(degree[i] * sizeof(**adjList));
    }
  }
  gra->adjList = adjList;

  if (owns_edges) {
    gra->edgesAlloced = 64;
    gra->edges = malloc(gra->edgesAlloced * sizeof(*gra->edges));
  }
  else {
    gra->edgesAlloced = 0;
    gra->edges = NULL;
  }
}

static void RDL_DFSvisit(const RDL_graph *gra, unsigned vertex, char *visited)
{
  unsigned i,j;
  for(i=0; i<gra->degree[vertex]; ++i)
  {
    j = gra->adjList[vertex][i][0];
    if(visited[j] == 0)
    {
      visited[j] = 1;
      RDL_DFSvisit(gra,j,visited);
    }
  }
}

/** returns 1 if the graph is connected, 0 otherwise. uses DFS */
char RDL_checkGraphConnected(const RDL_graph *gra)
{
  unsigned i;
  char *visited;
  char result;
  result = 1;
  visited = malloc(gra->V * sizeof(*visited));
  visited[0] = 1; /*start vertex*/
  for(i=1; i<gra->V; ++i)
  {
    visited[i] = 0; /*unvisited*/
  }
  RDL_DFSvisit(gra, 0, visited);
  for(i=0; i<gra->V; ++i)
  {
    if(visited[i] == 0) /*one was unvisited*/
    {
      result = 0;
      break;
    }
  }

  free(visited);
  return result;
}

/* if nodes are adjacent */
int RDL_isAdj(const RDL_graph *graph, unsigned i, unsigned j)
{
  unsigned idx;
  for(idx=0; idx<graph->degree[i]; ++idx)
  {
    if(graph->adjList[i][idx][0] == j)
    {
      return 1;
    }
  }
  return 0;
}

void RDL_deleteGraph(RDL_graph *gra)
{
  unsigned i;
  assert(gra != NULL);
  if(gra->edges)
  {
    for(i=0; i<gra->E; ++i)
    {
      free(gra->edges[i]);
    }
    free(gra->edges);
  }
  for(i=0; i<gra->V; ++i)
  {
    if(gra->degree[i] > 0)
    {
      free(gra->adjList[i]);
    }
  }
  free(gra->adjList);
  free(gra->degree);
  free(gra);
}

void RDL_printGraph(const RDL_graph *graph)
{
  unsigned i,j;

  printf("|V|=%d, |E|=%d\n",graph->V,graph->E);
  for(i=0; i<graph->V; ++i)
  {
    printf("%d:  ",i);
    for(j=0; j<graph->degree[i]; ++j)
    {
      printf("%d ",graph->adjList[i][j][0]);
    }
    printf("\n");
  }
  if(graph->edges)/*undirected*/
  {
    printf("edges:\n");
    for(i=0; i<graph->E; ++i)
    {
      printf("%d: [%d,%d]\n", i, graph->edges[i][0], graph->edges[i][1]);
    }
  }
}

RDL_graph *RDL_initNewGraph_g(unsigned V, char owns_edges)
{
  RDL_graph *graph = malloc(sizeof(*graph));
  unsigned *degree;
  unsigned i;
  degree = malloc(V * sizeof(*degree));
  for(i=0; i<V; ++i)
  {
    degree[i] = 0;
  }
  RDL_initGraph(graph,V,degree,owns_edges);

  return graph;
}

void RDL_addEdge(RDL_graph *gra, unsigned from, unsigned to)
{
  unsigned i;

  /* loops */
  if (from == to) {
    return;
  }

  for(i=0; i<gra->degree[from]; ++i) {
    if(gra->adjList[from][i][0] == to) {
      /*edge already exists*/
      return;
    }
  }
  ++gra->E;
  ++gra->degree[from];
  if(gra->degree[from] == 1)/*was 0, has never been initialized*/
  {
    gra->adjList[from] = malloc(gra->degree[from] * sizeof(*gra->adjList[from]));
  }
  else
  {
    gra->adjList[from] = realloc(gra->adjList[from], gra->degree[from] *
                                 sizeof(*gra->adjList[from]));
  }
  gra->adjList[from][ gra->degree[from]-1 ][0] = to;
}

static unsigned RDL_addToEdgeArray(RDL_graph *gra, unsigned from, unsigned to)
{
  unsigned temp;
  if(from > to)
  {
    temp = from;
    from = to;
    to = temp;
  }
  if(gra->E == gra->edgesAlloced)
  {
    gra->edgesAlloced *= 2;
    gra->edges = realloc(gra->edges, gra->edgesAlloced * sizeof(*gra->edges));
  }
  gra->edges[gra->E-1] = malloc(2 * sizeof(**gra->edges));
  gra->edges[gra->E-1][0] = from;
  gra->edges[gra->E-1][1] = to;

  return gra->E-1;
}

unsigned RDL_addUEdge_g(RDL_graph *gra, unsigned from, unsigned to)
{
  unsigned edge_id, i;

  if(from >= gra->V || to >= gra->V) {
    RDL_outputFunc(RDL_ERROR, "Tried to add an edge with atoms not in range.\n");
    RDL_outputFunc(RDL_ERROR,  "edge (%u,%u) can not be added to graph with %u atoms.\n",
        from, to, gra->V);
    return RDL_INVALID_RESULT;
  }

  if (from == to) {
    RDL_outputFunc(RDL_WARNING, "Adding a loop is not allowed, node %u\n", from);
    return RDL_INVALID_RESULT;
  }

  for(i=0; i<gra->degree[from]; ++i) {
    if(gra->adjList[from][i][0] == to) {
      /*edge already exists*/
      return RDL_DUPLICATE_EDGE;
    }
  }
  RDL_addEdge(gra, from, to);
  RDL_addEdge(gra, to, from);
  --gra->E; /*was incremented twice*/

  edge_id = RDL_addToEdgeArray(gra, from, to);

  gra->adjList[from][gra->degree[from]-1][1] = edge_id;
  gra->adjList[to][gra->degree[to]-1][1] = edge_id;

  return edge_id;
}

unsigned RDL_edgeId(const RDL_graph *gra, unsigned from, unsigned to)
{
  unsigned j, edge;

  if(from > to) {
    /*swap order to make from < to*/
    edge = to;
    to = from;
    from = edge;
  }

  edge = RDL_INVALID_RESULT;

  for(j=0; j<gra->degree[from]; ++j) {
    if(gra->adjList[from][j][0] == to) {
      edge = gra->adjList[from][j][1];
      break;
    }
  }

  return edge;
}
