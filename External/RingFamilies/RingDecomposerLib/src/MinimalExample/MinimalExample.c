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

/**
 * @file MinimalExample.c
 *
 * @brief Minimal example for the RingDecomposerLib library.
 *
 */

#include <assert.h>
#include <stdio.h>

#include "RingDecomposerLib.h"

/** @private */
int main()
{
  RDL_graph* graph;
  RDL_data* data;
  unsigned nof_urf, i;
  double nof_rc;
  const unsigned nof_edges=10, nof_nodes=8;
  unsigned edges[10][2] = {
    {0,1}, {1,2}, {2,3}, {3,4}, {4,5},
    {5,0}, {0,6}, {6,2}, {3,7}, {7,5} };
  /* step 1: prepare RDL_graph structure input graph */
  graph = RDL_initNewGraph(nof_nodes);
  for (i = 0; i < nof_edges; ++i) {
    RDL_addUEdge(graph, edges[i][0], edges[i][1]);
  }
  /* step 2: calculate URFs, RCFs */
  data = RDL_calculate(graph);
  /* check that calculation was successful */
  assert(data != NULL);
  /* step 3: how many URFs and RCs does the graph have? */
  nof_urf = RDL_getNofURF(data);
  nof_rc = RDL_getNofRC(data);

  RDL_cycle *cycle;
  RDL_cycleIterator *it = RDL_getRCyclesIterator(data);
  while(!RDL_cycleIteratorAtEnd(it)) {
    cycle = RDL_cycleIteratorGetCycle(it);
    /* <do something with the cycle> */
    RDL_deleteCycle(cycle);
    RDL_cycleIteratorNext(it);
  }
  RDL_deleteCycleIterator(it);

  /* delete data and graph */
  RDL_deleteData(data);

  printf("URFs: %d\n RCs: %.0f\n", nof_urf, nof_rc);

  return 0;
}
