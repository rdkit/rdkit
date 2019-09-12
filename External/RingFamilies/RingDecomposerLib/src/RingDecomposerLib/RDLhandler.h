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

#ifndef RDL_HANDLER_H
#define RDL_HANDLER_H

#include "RDLforward.h"
#include "RingDecomposerLib.h"

/** writes into the position i of the array of |V| integers if the vertex i is
contained on any shortest path from vertex a to vertex b.
(1 if yes, 0 otherwise)*/
void RDL_giveVertices(unsigned a, unsigned b, char *array, const RDL_sPathInfo *, char *visited);

/** writes into the position i of the array of |E| integers if the edge i is
contained on any shortest path from vertex a to vertex b.
(1 if yes, 0 otherwise)*/
void RDL_giveEdges(unsigned a, unsigned b, char *array, const RDL_graph *,
    const RDL_sPathInfo *, char* visited);

typedef enum RDL_IteratorType {
  RDL_RCF_IT, RDL_URF_IT, RDL_ALL_IT
} RDL_IteratorType;

/** combines all the paths out of 'paths1' with all out of 'paths2' to form all
possible cycles and writes them into the result array as elements of {0,1}^n.
Afterwards all allocated space in 'paths1' and 'paths2' is deallocated. */
RDL_cycleIterator* RDL_initCycleIterator(
    RDL_IteratorType itype,
    unsigned rcf_index,
    unsigned rcf_index_max,
    unsigned urf_index,
    unsigned urf_index_max,
    unsigned bcc_index,
    unsigned bcc_index_max,
    char mode,
    const RDL_data* data);

/** count the number of shortest paths from a to b */
double RDL_countPaths(unsigned a, unsigned b, unsigned V, const RDL_sPathInfo *spi);
#endif
