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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "RDLdataStruct.h"
#include "RDLgraph.h"
#include "RDLtarjan.h"
#include "RDLtesting.h"
#include "TestValidate.h"

void RDL_validate(RDL_data* data, RDL_DimacsGraph* dimacs_graph, int timeout,
                  int* urf_result, int* basis_result)
{
  unsigned URFcount, idx, weight, nof_relevant, curr_bcc,
    running_idx, idx2, nof_edges, nof_basis, edge, max_nof_basis_bcc=0;
  RDL_cycle **cycle_array;
  char **other_cycle_array, **final_cycle_array, **sortable_cycle_array, **basis_cycle_array;
  unsigned *urf_numbers, *bcc, *nof_basis_bcc;

  *urf_result = 0;
  *basis_result = 0;

  nof_relevant = RDL_getRCycles(data, &cycle_array);
  RDL_deleteCycles(cycle_array, nof_relevant);

  final_cycle_array = malloc(nof_relevant * sizeof(*final_cycle_array));
  urf_numbers = malloc(nof_relevant * sizeof(*urf_numbers));

  URFcount = RDL_getNofURF(data);

  nof_edges = RDL_getNofEdges(data);

  for (idx = 0, running_idx = 0; idx < URFcount; ++idx) {
    weight = RDL_getRCyclesForURF(data, idx,  &cycle_array);
    weight = RDL_translateCycArray(data, cycle_array, weight, &other_cycle_array);
    RDL_deleteCycles(cycle_array, weight);

    for (idx2 = 0; idx2 < weight; ++idx2, ++running_idx) {
      final_cycle_array[running_idx] = other_cycle_array[idx2];
      urf_numbers[running_idx] = idx;

    }

    free(other_cycle_array);
  }


  sortable_cycle_array = malloc(nof_relevant * sizeof(*sortable_cycle_array));
  for (idx = 0; idx < nof_relevant; ++idx) {
    sortable_cycle_array[idx] = malloc((nof_edges + 1) * sizeof(**sortable_cycle_array));
    memcpy(sortable_cycle_array[idx], final_cycle_array[idx], nof_edges * sizeof(**sortable_cycle_array));
    sortable_cycle_array[idx][nof_edges] = 2;
  }

  qsort(sortable_cycle_array, nof_relevant, sizeof(*sortable_cycle_array), RDL_cmp_cycles);

  if (dimacs_graph->nof_cycles == UINT_MAX) {
    fprintf(stderr, "skipping RC test, no RCs in input file...\n");
  }
  else if (nof_relevant != dimacs_graph->nof_cycles) {
    fprintf(stderr, "different number of relevant cycles %u vs %u!\n",
        nof_relevant, dimacs_graph->nof_cycles);
    *urf_result = 1;
  }
  else {
    for (idx = 0; idx < nof_relevant; ++idx) {
      if (memcmp(sortable_cycle_array[idx], dimacs_graph->sorted_cycles[idx],
          nof_edges * sizeof(**sortable_cycle_array))) {
        fprintf(stderr, "relevant cycle mismatch found!\n");
        *urf_result = 1;
        break;
      }
    }
  }

  if (!*urf_result) {
    RDL_timeout_stop_fun(&timeout);

    *urf_result = RDL_validateRingFamilies((const char**)final_cycle_array, urf_numbers,
        nof_relevant, nof_edges, RDL_timeout_stop_fun);
  }

  nof_basis = RDL_getSSSR(data, &cycle_array);

  bcc = malloc(nof_basis * sizeof(*bcc));
  nof_basis_bcc = malloc(data->bccGraphs->nof_bcc * sizeof(*nof_basis_bcc));
  for (idx = 0; idx < data->bccGraphs->nof_bcc; ++idx) {
    nof_basis_bcc[idx] = 0;
  }
  for (idx = 0; idx < nof_basis; ++idx) {
    edge = RDL_getEdgeId(data, cycle_array[idx]->edges[0][0],
        cycle_array[idx]->edges[0][1]);
    curr_bcc = data->bccGraphs->edge_to_bcc_mapping[edge][0];
    bcc[idx] = curr_bcc;
    ++nof_basis_bcc[curr_bcc];
    max_nof_basis_bcc = nof_basis_bcc[curr_bcc] > max_nof_basis_bcc ? nof_basis_bcc[curr_bcc] : max_nof_basis_bcc;
  }

  RDL_translateCycArray(data, cycle_array, nof_basis, &basis_cycle_array);
  RDL_deleteCycles(cycle_array, nof_basis);

  RDL_timeout_stop_fun(&timeout);
  *basis_result = RDL_validateBasis(
      (const char**)final_cycle_array, nof_relevant,
      (const char**)basis_cycle_array, nof_basis,
      nof_edges, bcc, RDL_timeout_stop_fun);
  free(bcc);
  free(urf_numbers);
  free(nof_basis_bcc);

  RDL_deleteEdgeIdxArray(final_cycle_array, nof_relevant);
  RDL_deleteEdgeIdxArray(basis_cycle_array, nof_basis);
  RDL_deleteEdgeIdxArray(sortable_cycle_array, nof_relevant);
}
