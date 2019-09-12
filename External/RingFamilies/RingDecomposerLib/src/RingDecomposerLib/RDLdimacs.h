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

#ifndef RDL_DIMACS_H
#define RDL_DIMACS_H

#include "RingDecomposerLib.h"

/*
 * simple reading of a dimacs graph!
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct RDL_DimacsGraph {
  RDL_edge* edges;
  RDL_cycle** cycles;
  char** sorted_cycles;
  unsigned nof_edges;
  unsigned nof_nodes;
  unsigned nof_cycles;
  char* name;
} RDL_DimacsGraph;

RDL_API
void RDL_dimacsGraph_delete(RDL_DimacsGraph* graph);

RDL_API
RDL_DimacsGraph* RDL_dimacsGraph_read(const char* filename);

#ifdef __cplusplus
}
#endif

#endif
