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

#include <stdio.h>
#include <stdlib.h>

#include "TestDemo.h"

void RDL_demo_output(RDL_data* URFdata)
{
  unsigned URFcount, idx, count, obIdx, edgeCount;
  RDL_edge *bondArray;
  RDL_node *atoms;
  RDL_cycle **cycleArray;
  char **otherCycleArray;

  /* some output */
  URFcount = RDL_getNofURF(URFdata);
  printf("==========================================================URF=\n");
  printf("Number of Unique Ring Families: %d\n\n", URFcount);
  for(idx=0; idx<URFcount; ++idx)
  {
    printf("URF %d has weight %d.\n", idx, RDL_getWeightForURF(URFdata, idx));
  }
  /* some more output which might change when the order of the input is changed*/
  printf("\n===============================================================\n");
  printf("The rest of this output might depend on the order of the input:\n\n");
  for(idx=0; idx<URFcount; ++idx)
  {
      count = RDL_getEdgesForURF(URFdata, idx, &bondArray);
      printf("There are %d bonds in URF %d.\n", count, idx);
      free(bondArray);
  }
  printf("\n");

  for(idx=0; idx<URFcount; ++idx)
  {
      printf("Atoms in URF %d: ",idx);
      count = RDL_getNodesForURF(URFdata, idx, &atoms);
      for(obIdx=0; obIdx<count; ++obIdx)
      {
          printf("%d ",atoms[obIdx]);
      }
      printf("\n");
      free(atoms);
  }
  printf("\n");

  printf("A possible MCB (SSSR) ");
  count = RDL_getSSSR(URFdata, &cycleArray);
  printf("(%d rings):\n",count);
  for(idx=0; idx<count; ++idx)
  {
      printf("ring %d: ",idx);
      for(obIdx=0; obIdx<cycleArray[idx]->weight;  ++obIdx)
      {
          printf("(%d ",cycleArray[idx]->edges[obIdx][0]);
          printf("%d), ",cycleArray[idx]->edges[obIdx][1]);
      }
      printf("\n");
  }
  RDL_deleteCycles(cycleArray, count);
  printf("\n");

  printf("The RC Prototypes with bonds as pairs of atoms ");
  count = RDL_getRCPrototypes(URFdata, &cycleArray);
  printf("(%d rings):\n",count);
  for(idx=0; idx<count; ++idx)
  {
      printf("ring %d: ",idx);
      for(obIdx=0; obIdx<cycleArray[idx]->weight; ++obIdx)
      {
          printf("(%d ",cycleArray[idx]->edges[obIdx][0]);
          printf("%d), ",cycleArray[idx]->edges[obIdx][1]);
      }
      printf("\n");
  }
  printf("\n");

  printf("The RC Prototypes as vectors ");
  count = RDL_translateCycArray(URFdata, cycleArray, count, &otherCycleArray);
  printf("(%d rings):\n",count);
  /* To be able to understand the bitsets better: */
  edgeCount = RDL_getEdgeArray(URFdata, &bondArray);
  printf("Edge from");
  for(idx=0; idx<edgeCount; ++idx)
  {
      printf("%2d",bondArray[idx][0]);
  }
  printf("\n     to  ");
  for(idx=0; idx<edgeCount; ++idx)
  {
      printf("%2d",bondArray[idx][1]);
  }
  printf("\n");
  free(bondArray);
  /* the bitsets: */
  for(idx=0; idx<count; ++idx)
  {
      printf("ring %3d: ",idx);
      for(obIdx=0; obIdx<edgeCount; ++obIdx)
      {
          printf("%d ",otherCycleArray[idx][obIdx]);
      }
      printf("\n");
  }
  RDL_deleteCycles(cycleArray, count);
  RDL_deleteEdgeIdxArray(otherCycleArray, count);
  printf("==========================================================URF=\n");
}
