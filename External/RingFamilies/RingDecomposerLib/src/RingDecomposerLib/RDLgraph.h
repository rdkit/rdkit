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

/** graph data structure that is used for the URF calculation.
The vertices have to be numbered from 0 to |V|-1. Call the following functions:

 - RDL_initNewGraph_g(unsigned V) to initialize a new graph with V vertices
    (returns RDL_graph pointer)

 - RDL_addUEdge_g(RDL_graph *, unsigned from, unsigned to) to add a new (undirected) edge from
    the vertex with index "from" to the vertex with index "to". It is NO problem
    to add an edge twice - even with different order of vertices
    (will be ignored).

now RDL_calculate can be called on it

 - RDL_deleteGraph to free all allocated space */

#ifndef RDL_GRAPH_H
#define RDL_GRAPH_H

#include "RingDecomposerLib.h"

struct RDL_graph {
  unsigned V, E; /*number of vertices/edges*/
  unsigned *degree; /*array storing how many vertices are adjacent to each vertex*/
  RDL_edge **adjList; /*the vertices are stored with their index (0 to |V|-1)*/
  unsigned **edges; /*array containing pairs of vertices as edges*/
  unsigned edgesAlloced; /*space alloced for this many edges in 'edges' array*/
};

/** initializes a new Graph that edges can be added to (allocates space for it) */
RDL_graph *RDL_initNewGraph_g(unsigned V, char owns_edges);

/** adds an undirected edge to the graph. */
unsigned RDL_addUEdge_g(RDL_graph *, unsigned from, unsigned to);

/** frees all allocated space. */
void RDL_deleteGraph(RDL_graph *gra);

/*=============================================================================*/
/* This structure can also be used for directed graphs (which is being done in
the calculation of CFs). It that case the degree is the OUT-degree only and the
array edges is not used. */

/** adds a directed edge to the graph. */
void RDL_addEdge(RDL_graph *, unsigned from, unsigned to);

/** returns 1 if the vertex i is adjacent to j in the graph, 0 otherwise. */
int RDL_isAdj(const RDL_graph *, unsigned i, unsigned j);

/** prints the adjacency lists, |V| and |E| and the edges (edges only if the
    graph is undirected) */
void RDL_printGraph(const RDL_graph *graph);

/* checks if the graph is connected */
char RDL_checkGraphConnected(const RDL_graph *gra);

/** Returns the index of the bond connecting the atoms 'from' and 'to' */
unsigned RDL_edgeId(const RDL_graph *gra, unsigned from, unsigned to);
#endif
