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

#ifndef RDL_APSP_H
#define RDL_APSP_H

#include "RDLforward.h"
#include "RingDecomposerLib.h"

struct RDL_sPathInfo {
  unsigned **pred; /*predecessor matrix*/
  unsigned **dist; /*distance matrix*/
  char **reachable; /*reachable vertices for each vertex via shortest paths
     only using vertices preceding in order pi (1 or 0/reachable or not)*/
  RDL_graph **dPaths; /*for each vertex: a subgraph of the original graph that
     stores the shortest paths back to the vertex (in Vismara: U_r).Useful
     to find paths/edges of a RCF.*/
};

/** Solves the "All Pairs Shortest Paths" on the Graph and returns the result.
Allocates space for dPaths, but does not initialize it. */
RDL_sPathInfo *RDL_AllPairsShortestPaths(RDL_graph *gra);

/** Deletes everything that 'RDL_AllPairsShortestPaths' has created. Needs the
number of vertices to delete dPaths. */
void RDL_deleteAPSP(RDL_sPathInfo *info, unsigned number);

#endif
