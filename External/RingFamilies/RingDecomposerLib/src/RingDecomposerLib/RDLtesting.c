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

#include "RDLtesting.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>

#include "RDLdataStruct.h"
#include "RDLgraph.h"
#include "RDLstack.h"
#include "RDLutility.h"
#include "RingDecomposerLib.h"

int RDL_no_stop_fun(const void* nothing)
{
  return 0;
}

int RDL_timeout_stop_fun(const void* p_timeout)
{
  static int timeout = -1;
  static time_t start = -1;
  time_t end;

  if (p_timeout) {
    timeout = *(int*)p_timeout;
    start = time(NULL);
  }

  if (timeout < 0 || start < 0) {
    return 0;
  }
  else {
    end = time(NULL);
    if (end - start > timeout) {
      return 1;
    }
    else {
      return 0;
    }
  }
}

int RDL_cmp_cycles(const void* vc1, const void* vc2)
{
  unsigned first_difference = 0;
  unsigned len1 = 0;
  unsigned len2 = 0;
  unsigned i = 0;

  const char* c1 = *((const char**) vc1);
  const char* c2 = *((const char**) vc2);

  while (c1[i] != 2) {
    if (c1[i]) {
      ++len1;
      if (!c2[i] && !first_difference) {
        first_difference = 1;
      }
    }
    if (c2[i]) {
      ++len2;
      if (!c1[i] && !first_difference) {
        first_difference = 2;
      }
    }
    ++i;
  }

  if (len1 < len2) {
    return -1;
  }
  else if (len1 == len2) {
    if (first_difference == 1) {
      return -1;
    }
    else if (first_difference == 0) {
      return 0;
    }
  }

  return 1;
}


typedef struct RDL_stackElement {
    unsigned current_ring_idx;
    unsigned stage;
    const char* current_bonds;
    unsigned owns_bonds;
    unsigned curr_bcc;
} RDL_stackElement;

int RDL_check_property3(const char* start_bonds,
                        unsigned start_bonds_size,
                        const char* ring2_bonds,
                        const char** current_rings,
                        unsigned current_rings_size,
                        unsigned start_ring_idx,
                        int (*stop) (const void*))
{
  const char *current_bonds, *ring;
  char *new_bonds;
  unsigned current_ring_idx, i;
  int family_related;
  RDL_stack* stack;
  RDL_stackElement *top, *new_element;

  stack = RDL_stack_new();

  top = malloc(sizeof(*top));
  top->current_ring_idx = start_ring_idx;
  top->current_bonds = start_bonds;
  top->stage = 0;
  top->owns_bonds = 0;
  RDL_stack_push(stack, top);

  family_related = 0;

  while (!RDL_stack_empty(stack) && !family_related) {
    top = RDL_stack_top(stack);
    current_bonds = top->current_bonds;
    current_ring_idx = top->current_ring_idx;

    family_related = 1;
    for (i = 0; i < start_bonds_size; ++i) {
      if (current_bonds[i] != ring2_bonds[i]) {
        family_related = 0;
      }
    }

    if (!family_related && top->stage < 2 && current_ring_idx < current_rings_size) {
      /* looking at this ring the first time, do not add it and recurse */
      if (top->stage == 0) {
        top->stage = 1;

        new_element = malloc(sizeof(*new_element));
        new_element->stage = 0;
        new_element->current_bonds = current_bonds;
        new_element->current_ring_idx = current_ring_idx + 1;
        new_element->owns_bonds = 0;
        RDL_stack_push(stack, new_element);
      }
      /* looking at this ring the second time, add it and recurse */
      else if (top->stage == 1) {
        top->stage = 2;

        ring = current_rings[current_ring_idx];
        new_bonds = malloc(start_bonds_size * sizeof(*new_bonds));
        for (i = 0; i < start_bonds_size; ++i) {
          new_bonds[i] = !ring[i] != !current_bonds[i] ? 1 : 0; /* logical XOR */
        }

        new_element = malloc(sizeof(*new_element));
        new_element->stage = 0;
        new_element->current_bonds = new_bonds;
        new_element->current_ring_idx = current_ring_idx + 1;
        new_element->owns_bonds = 1;
        RDL_stack_push(stack, new_element);
      }
    }
    /* looking at this ring for the third time or at an non existing ring, end recursion */
    else {
      if (top->owns_bonds) {
        free((char*)top->current_bonds);
      }
      free(top);
      RDL_stack_pop(stack);
    }
    if (stop(NULL)) {
      family_related = -1;
      break;
    }
  }

  while (!RDL_stack_empty(stack)) {
    top = RDL_stack_top(stack);
    if (top->owns_bonds) {
      free((char*)top->current_bonds);
    }
    free(top);
    RDL_stack_pop(stack);
  }

  RDL_stack_delete(stack);

  return family_related;
}

unsigned RDL_count_bonds(const char* bitset, unsigned nofBonds)
{
  unsigned count;
  unsigned idx;

  for (idx = 0, count = 0; idx < nofBonds; ++idx) {
    if (bitset[idx]) {
      ++count;
    }
  }

  return count;
}

char** RDL_calculate_pairwise_family(const char** cycle_array,
                                  unsigned nof_rings,
                                  unsigned nof_bonds,
                                  int (*stop) (const void*))
{
  char* family_relation_data;
  char **family_relation, **ring_sets;
  const char **rings;
  const char *ring1, *ring2, *ring;
  unsigned i, j, k, ring1_count, ring2_count, intersection_count,
    ringcounter, ridx, ring_count;

  size_t matrix_size = nof_rings;

  if (nof_rings >= sqrt(UINT_MAX) - 1) {
    return NULL;
  }

  matrix_size *= nof_rings;

  /* first establish pairwise URF-relations */
  family_relation_data = malloc(matrix_size * sizeof(*family_relation_data));

  /* too much memory requested */
  if (!family_relation_data) {
    return NULL;
  }

  memset(family_relation_data, 0, matrix_size * sizeof(*family_relation_data));
  family_relation = malloc(nof_rings * sizeof(*family_relation));

  for (i = 0; i < nof_rings; ++i) {
    family_relation[i] = family_relation_data + i * nof_rings;
  }

  ring_sets = malloc(nof_rings * sizeof(*ring_sets));

  for (i = 0; i < nof_rings; ++i) {
    ring_sets[i] = malloc(nof_bonds * sizeof(**ring_sets));
    memcpy(ring_sets[i], cycle_array[i], nof_bonds);
  }

  rings = malloc(nof_rings * sizeof(*rings));

  for (i = 0; i < nof_rings && !stop(NULL); ++i) {
    for (j = i+1; j < nof_rings && !stop(NULL); ++j) {
      ring1 = ring_sets[i];
      ring2 = ring_sets[j];
      /* check property 1: |C1| == |C2| */
      ring1_count = RDL_count_bonds(ring1, nof_bonds);
      ring2_count = RDL_count_bonds(ring2, nof_bonds);
      if (ring1_count != ring2_count) {
        continue;
      }
      /* check property 2: C1 ^ C2 != o */
      intersection_count = 0;
      for (k = 0; k < nof_bonds; ++k) {
        if (ring1[k] && ring2[k]) {
          ++intersection_count;
        }
      }
      if (intersection_count == 0) {
        continue;
      }

      /*
       * check property 3: exists any subset of strictly smaller rings,
       * so the symmetric difference of ring 1 and a subset of
       * those smaller rings equals ring 2
       */
      ringcounter = 0;
      for(ridx = 0; ridx < nof_rings; ++ridx) {
        ring = ring_sets[ridx];
        ring_count = RDL_count_bonds(ring, nof_bonds);
        if (ring_count < ring1_count) {
          rings[ringcounter] = ring;
          ++ringcounter;
        }
      }

      family_relation[i][j] = RDL_check_property3(
          ring1, nof_bonds, ring2, rings, ringcounter, 0, stop);
      family_relation[j][i] = family_relation[i][j];
    }
  }

  for (i = 0; i < nof_rings; ++i) {
    free(ring_sets[i]);
  }

  free(ring_sets);
  free((void*)rings);

  return family_relation;
}

/*
 * calculate transitive closure of adjacency matrix
 * with DFS frome each node (O(n^2))
 */
void RDL_calculate_transitive_closure(char** graph, unsigned n)
{
  RDL_stack* dfs_stack;
  char* visited = (char*)malloc(n * sizeof(*visited));
  unsigned i, j, current;
  unsigned *val, *current_element, *new_element;

  for (i = 0; i < n; ++i) {
    dfs_stack = RDL_stack_new();
    memset(visited, 0, n * sizeof(*visited));

    val = malloc(sizeof(*val));
    *val = i;

    RDL_stack_push(dfs_stack, val);

    while (!RDL_stack_empty(dfs_stack)) {
      current_element = RDL_stack_top(dfs_stack);
      current = *current_element;
      RDL_stack_pop(dfs_stack);
      free(current_element);

      visited[current] = 1;
      for (j = 0; j < n; ++j) {
        if (graph[current][j]) {
          graph[i][j] = 1;
          graph[j][i] = 1;
          if (!visited[j]) {
            new_element = malloc(sizeof(*new_element));
            *new_element = j;
            RDL_stack_push(dfs_stack, new_element);
          }
        }
      }
    }

    RDL_stack_delete(dfs_stack);
  }

  free(visited);
}

/*
 * validate the URFs of given ringsystem,
 * StopFunc specifies function to interrupt the process
 * return 0 if successful, 1 on difference, -1 on timeout
 */
int RDL_validateRingFamilies(const char** cycle_array,
                             unsigned* family_numbers,
                             unsigned nof_rings,
                             unsigned nof_bonds,
                             int (*stop) (const void*))
{
  char** family_relation;
  unsigned i, j;
  int return_value = 0;

  /* if there are no rings, don't bother with checking */
  if (!nof_rings) {
    return 0;
  }

  /* calculate pairwise URF relation */
  family_relation = RDL_calculate_pairwise_family(
      cycle_array, nof_rings, nof_bonds, stop);

  if (!family_relation) {
    return -2;
  }

  if (stop(NULL)) {
    return_value = -1;
  }
  else {
    /* calculate transitive closure of pairwise URF relation */
    RDL_calculate_transitive_closure(family_relation, nof_rings);

    /* compare URF relations */
    for (i = 0; i < nof_rings; ++i) {
      for (j = i+1; j < nof_rings; ++j) {
        unsigned family1 = family_numbers[i];
        unsigned family2 = family_numbers[j];

        if (family1 == family2 && !family_relation[i][j]) {
          return_value = 1;
        }
        else if (family1 != family2 && family_relation[i][j]){
          return_value = 1;

        }
      }
    }
  }

  free(*family_relation);
  free(family_relation);

  return return_value;
}

int cmp_unsigned(const void* a, const void* b)
{
  unsigned va, vb;

  va = *((unsigned*)a);
  vb = *((unsigned*)b);

  if (va < vb) {
    return -1;
  }
  else if (va == vb) {
    return 0;
  }
  else {
    return 1;
  }
}

/* check that we can handle invalid input */
int RDL_checkInvalidConsistency(const RDL_data* data)
{
  RDL_edge *edges;
  RDL_node *nodes;

  RDL_cycle **cycle_array;
  unsigned *families;

  unsigned result = 0;
  unsigned result_num;
  double result_num_double;
  const unsigned VeryHighIndex = UINT_MAX/2;
  const unsigned VeryHighIndexPlus1 = VeryHighIndex + 1;

  RDL_graph *empty_graph;
  RDL_data *empty_result;
  RDL_cycleIterator *it;

  /* the next errors on the RDL error function output are expected! */
  RDL_outputFunc(RDL_DEBUG, "The following errors are expected {\n");

  empty_graph = RDL_initNewGraph(0);

  /* test RDL_calculate */
  empty_result = RDL_calculate(NULL);
  if (empty_result) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_calculate(NULL)\n");
    result = 1;
  }

  empty_result = RDL_calculate(empty_graph);
  if (empty_result) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_calculate(empty_graph)\n");
    result = 1;
  }

  RDL_deleteGraph(empty_graph);

  /* test RDL_addUEdge */
  empty_graph = RDL_initNewGraph(2);
  if (RDL_addUEdge(empty_graph, 1, 2) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_addUEdge(empty_graph, 1, 2)\n");
    result = 1;
  }

  RDL_addUEdge(empty_graph, 0, 1);
  if (RDL_addUEdge(empty_graph, 1, 0) != RDL_DUPLICATE_EDGE) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_addUEdge(empty_graph, 0, 1)\n");
    result = 1;
  }

  if (RDL_addUEdge(empty_graph, 1, 1) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_addUEdge(empty_graph, 1, 1)\n");
    result = 1;
  }

  RDL_deleteGraph(empty_graph);

  /* test RDL_getNofURF */
  if (RDL_getNofURF(NULL) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofURF(NULL)\n");
    result = 1;
  }

  /* test RDL_getNofRCF */
  if (RDL_getNofRCF(NULL) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCF(NULL)\n");
    result = 1;
  }

  /* test RDL_getNofRCForURF */
  result_num_double = RDL_getNofRCForURF(data, VeryHighIndex);
  if (result_num_double != RDL_INVALID_RC_COUNT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCForURF(data, VeryHighIndex)\n");
    result = 1;
  }

  result_num_double = RDL_getNofRCForURF(NULL, 0);
  if (result_num_double != RDL_INVALID_RC_COUNT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCForURF(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_edgeId */
  if (RDL_getEdgeId(data, 1, 1) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgeId(data, 1, 1)\n");
    result = 1;
  }

  if (RDL_getEdgeId(NULL, 1, 1) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgeId(data, 1, 1)\n");
    result = 1;
  }

  if (RDL_getEdgeId(data, VeryHighIndex, 1) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgeId(data, VeryHighIndex, 1)\n");
    result = 1;
  }

  if (RDL_getEdgeId(data, 1, VeryHighIndex) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgeId(data, 1, VeryHighIndex)\n");
    result = 1;
  }

  /* test RDL_getNofRCForRCF */
  result_num_double = RDL_getNofRCForRCF(data, VeryHighIndex);
  if (result_num_double != RDL_INVALID_RC_COUNT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCForRCF(data, VeryHighIndex)\n");
    result = 1;
  }

  result_num_double = RDL_getNofRCForRCF(NULL, 0);
  if (result_num_double != RDL_INVALID_RC_COUNT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCForRCF(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getNofRC */
  result_num_double = RDL_getNofRC(NULL);
  if (result_num_double != RDL_INVALID_RC_COUNT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRC(NULL)\n");
    result = 1;
  }

  /* test RDL_getRCyclesForURF */
  result_num = RDL_getRCyclesForURF(data, VeryHighIndex, &cycle_array);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForURF(data, VeryHighIndex, &cycle_array1)\n");
    result = 1;
  }
  RDL_deleteCycles(cycle_array, 0);

  /* test RDL_getRCyclesForURF */
  result_num = RDL_getRCyclesForURF(NULL, 0, &cycle_array);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForURF(NULL, 0, &cycle_array1)\n");
    result = 1;
  }
  RDL_deleteCycles(cycle_array, 0);

  /* test RDL_getRCyclesForURFIterator */
  it = RDL_getRCyclesForURFIterator(data, VeryHighIndex);
  if (it != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForURFIterator(data, VeryHighIndex)\n");
    result = 1;
  }

  /* test RDL_getRCyclesForURFIterator */
  it = RDL_getRCyclesForURFIterator(NULL, 0);
  if (it != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForURFIterator(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getRCyclesForRCF */
  result_num = RDL_getRCyclesForRCF(data, VeryHighIndex, &cycle_array);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForRCF(data, VeryHighIndex, &cycle_array1)\n");
    result = 1;
  }
  RDL_deleteCycles(cycle_array, 0);

  /* test RDL_getRCyclesForRCF */
  result_num = RDL_getRCyclesForRCF(NULL, 0, &cycle_array);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForRCF(NULL, 0, &cycle_array1)\n");
    result = 1;
  }
  RDL_deleteCycles(cycle_array, 0);

  /* test RDL_getRCyclesForRCFIterator */
  it = RDL_getRCyclesForRCFIterator(data, VeryHighIndex);
  if (it != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForRCFIterator(data, VeryHighIndex)\n");
    result = 1;
  }

  /* test RDL_getRCyclesForRCFIterator */
  it = RDL_getRCyclesForRCFIterator(NULL, 0);
  if (it != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesForRCFIterator(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getRCycles */
  result_num = RDL_getRCycles(NULL, &cycle_array);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCycles(NULL, &cycle_array1)\n");
    result = 1;
  }
  RDL_deleteCycles(cycle_array, 0);

  /* test RDL_getRCyclesIterator */
  it = RDL_getRCyclesIterator(NULL);
  if (it != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCyclesIterator(NULL)\n");
    result = 1;
  }

  /* test RDL_cycleIteratorAtEnd */
  if (!RDL_cycleIteratorAtEnd(NULL)) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_cycleIteratorAtEnd(NULL)\n");
    result = 1;
  }

  /* RDL_cycleIteratorNext */
  it = RDL_getRCyclesIterator(data);
  while (!RDL_cycleIteratorAtEnd(it)) {
    RDL_cycleIteratorNext(it);
  }
  if (!it || RDL_cycleIteratorNext(it) != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_cycleIteratorNext(it)\n");
    result = 1;
  }

  if (RDL_cycleIteratorNext(NULL) != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_cycleIteratorNext(NULL)\n");
    result = 1;
  }

  /* RDL_cycleIteratorGetCycle */
  if (RDL_cycleIteratorGetCycle(NULL) != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_cycleIteratorGetCycle(NULL)\n");
    result = 1;
  }

  if (!it || RDL_cycleIteratorGetCycle(it) != NULL) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_cycleIteratorGetCycle(it)\n");
    result = 1;
  }

  RDL_deleteCycleIterator(it);

  /* test RDL_getRCPrototypes*/
  result_num = RDL_getRCPrototypes(NULL, &cycle_array);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCPrototypes(NULL, &cycle_array1)\n");
    result = 1;
  }
  RDL_deleteCycles(cycle_array, 0);

  /* test RDL_getSSSR */
  result_num = RDL_getSSSR(NULL, &cycle_array);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getSSSR(NULL, &cycle_array1)\n");
    result = 1;
  }
  RDL_deleteCycles(cycle_array, 0);

  /* test RDL_getWeightForURF */
  result_num = RDL_getWeightForURF(data, VeryHighIndex);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getWeightForURF(data, VeryHighIndex)\n");
    result = 1;
  }

  result_num = RDL_getWeightForURF(NULL, 0);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getWeightForURF(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getWeightForRCF */
  result_num = RDL_getWeightForRCF(data, VeryHighIndex);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getWeightForRCF(data, VeryHighIndex)\n");
    result = 1;
  }

  result_num = RDL_getWeightForRCF(NULL, 0);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getWeightForRCF(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getNodesForURF */
  result_num = RDL_getNodesForURF(data, VeryHighIndex, &nodes);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNodesForURF(data, VeryHighIndex, &family_nodes1)\n");
    result = 1;
  }
  free(nodes);

  result_num = RDL_getNodesForURF(NULL, 0, &nodes);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNodesForURF(NULL, 0, &family_nodes1)\n");
    result = 1;
  }
  free(nodes);

  /* test RDL_getEdgesForURF */
  result_num = RDL_getEdgesForURF(data, VeryHighIndex, &edges);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgesForURF(data, VeryHighIndex, &family_edges1)\n");
    result = 1;
  }
  free(edges);

  result_num = RDL_getEdgesForURF(NULL, 0, &edges);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,"failed to detect invalid result: RDL_getEdgesForURF(NULL, 0, &family_edges1)\n");
    result = 1;
  }
  free(edges);

  /* test RDL_getNodesForRCF */
  result_num = RDL_getNodesForRCF(data, VeryHighIndex, &nodes);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNodesForRCF(data, VeryHighIndex, &rcf_nodes1)\n");
    result = 1;
  }
  free(nodes);

  result_num = RDL_getNodesForRCF(NULL, 0, &nodes);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNodesForURF(NULL, 0, &family_nodes1)\n");
    result = 1;
  }
  free(nodes);

  /* test RDL_getEdgesForRCF */
  result_num = RDL_getEdgesForRCF(data, VeryHighIndex, &edges);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgesForRCF(data, VeryHighIndex, &rcf_edges1)\n");
    result = 1;
  }
  free(edges);

  result_num = RDL_getEdgesForRCF(NULL, 0, &edges);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgesForRCF(NULL, 0q, &rcf_edges1)\n");
    result = 1;
  }
  free(edges);

  /* test RDL_getURFsContainingNode */
  result_num = RDL_getURFsContainingNode(data, VeryHighIndex, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getURFsContainingNode(data, VeryHighIndex, &families_for_node2)\n");
    result = 1;
  }
  free(families);

  result_num = RDL_getURFsContainingNode(NULL, 0, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getURFsContainingNode(NULL, 0, &families_for_node2)\n");
    result = 1;
  }
  free(families);

  /* test RDL_getNofURFContainingNode */
  result_num = RDL_getNofURFContainingNode(data, VeryHighIndex);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofURFContainingNode(data, VeryHighIndex)\n");
    result = 1;
  }

  result_num = RDL_getNofURFContainingNode(NULL, 0);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofURFContainingNode(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getURFsContainingEdge */
  result_num = RDL_getURFsContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getURFsContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1, &families_for_edge2)\n");
    result = 1;
  }
  free(families);

  result_num = RDL_getURFsContainingEdge(data, 0, 0, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getURFsContainingEdge(data, 0, 0, &families_for_edge2)\n");
    result = 1;
  }
  free(families);

  result_num = RDL_getURFsContainingEdge(NULL, 0, 0, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getURFsContainingEdge(NULL, 0, 0, &families_for_edge2)\n");
    result = 1;
  }
  free(families);

  /* test RDL_getNofURFContainingEdge */
  result_num = RDL_getNofURFContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofURFContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1)\n");
    result = 1;
  }

  result_num = RDL_getNofURFContainingEdge(NULL, 0, 1);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofURFContainingEdge(NULL, 0, 1)\n");
    result = 1;
  }

  result_num = RDL_getNofURFContainingEdge(data, 1, 1);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofURFContainingEdge(data, 1, 1)\n");
    result = 1;
  }

  /* test RDL_getRCFsContainingNode */
  result_num = RDL_getRCFsContainingNode(data, VeryHighIndex, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCFsContainingNode(data, VeryHighIndex, &families_for_node2)\n");
    result = 1;
  }
  free(families);

  result_num = RDL_getRCFsContainingNode(NULL, 0, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCFsContainingNode(NULL, 0, &families_for_node2)\n");
    result = 1;
  }
  free(families);

  /* test RDL_getNofRCFContainingNode */
  result_num = RDL_getNofRCFContainingNode(data, VeryHighIndex);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCFContainingNode(data, VeryHighIndex)\n");
    result = 1;
  }

  result_num = RDL_getNofRCFContainingNode(NULL, 0);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCFContainingNode(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getRCFsContainingEdge */
  result_num = RDL_getRCFsContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCFsContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1, &families_for_edge2)\n");
    result = 1;
  }
  free(families);

  result_num = RDL_getRCFsContainingEdge(data, 0, 0, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCFsContainingEdge(data, 0, 0, &families_for_edge2)\n");
    result = 1;
  }
  free(families);

  result_num = RDL_getRCFsContainingEdge(NULL, 0, 0, &families);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRCFsContainingEdge(NULL, 0, 0, &families_for_edge2)\n");
    result = 1;
  }
  free(families);

  /* test RDL_getNofRCFContainingEdge */
  result_num = RDL_getNofRCFContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCFContainingEdge(data, VeryHighIndex, VeryHighIndexPlus1)\n");
    result = 1;
  }

  result_num = RDL_getNofRCFContainingEdge(NULL, 0, 1);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCFContainingEdge(NULL, 0, 1)\n");
    result = 1;
  }

  result_num = RDL_getNofRCFContainingEdge(data, 1, 1);
  if (result_num != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRCFContainingEdge(data, 1, 1)\n");
    result = 1;
  }

  /* test RDL_getNofRingsystems */
  if (RDL_getNofRingsystems(NULL) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofRingsystems(NULL)\n");
    result = 1;
  }

  /* test RDL_getNofNodesForRingsystem */
  if (RDL_getNofNodesForRingsystem(data, VeryHighIndex) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofNodesForRingsystem(data, VeryHighIndex)\n");
    result = 1;
  }

  if (RDL_getNofNodesForRingsystem(NULL, 0) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofNodesForRingsystem(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getNofedgesForRingsystem */
  if (RDL_getNofEdgesForRingsystem(data, VeryHighIndex) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofEdgesForRingsystem(data, VeryHighIndex)\n");
    result = 1;
  }

  if (RDL_getNofEdgesForRingsystem(NULL, 0) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNofEdgesForRingsystem(NULL, 0)\n");
    result = 1;
  }

  /* test RDL_getNodesForRingsystem */
  if (RDL_getNodesForRingsystem(data, VeryHighIndex, &nodes) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNodesForRingsystem(data, VeryHighIndex, &family_nodes1)\n");
    result = 1;
  }
  free(nodes);

  if (RDL_getNodesForRingsystem(NULL, 0, &nodes) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getNodesForRingsystem(NULL, 0, &family_nodes1)\n");
    result = 1;
  }
  free(nodes);

  /* test RDL_getEdgesForRingsystem */
  if (RDL_getEdgesForRingsystem(data, VeryHighIndex, &edges) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgesForRingsystem(data, VeryHighIndex, &family_edges1)\n");
    result = 1;
  }
  free(edges);

  if (RDL_getEdgesForRingsystem(NULL, 0, &edges) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getEdgesForRingsystem(NULL, 0, &family_edges1)\n");
    result = 1;
  }
  free(edges);

  /* test RDL_getRingsystemForEdge */
  if (RDL_getRingsystemForEdge(data, VeryHighIndex, VeryHighIndexPlus1) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRingsystemForEdge(data, VeryHighIndex, VeryHighIndexPlus1)\n");
    result = 1;
  }

  if (RDL_getRingsystemForEdge(data, 1, 1) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRingsystemForEdge(data, 1, 1)\n");
    result = 1;
  }

  if (RDL_getRingsystemForEdge(NULL, 0, 1) != RDL_INVALID_RESULT) {
    fprintf(stderr,
        "failed to detect invalid result: RDL_getRingsystemForEdge(NULL, 0, 1)\n");
    result = 1;
  }

  RDL_outputFunc(RDL_DEBUG, "} end expected errors\n");

  return result;
}

/* check the consistency of URFs/RCFs */
int RDL_checkFamilyConsistency(const RDL_data* data,
    unsigned (*getNumberPtr) (const RDL_data*),
    double (*getNofRCPtr) (const RDL_data*, unsigned),
    unsigned (*getRCyclesPtr) (const RDL_data*, unsigned, RDL_cycle***),
    unsigned (*getWeightPtr) (const RDL_data*, unsigned),
    unsigned (*getEdgesPtr) (const RDL_data*, unsigned, RDL_edge**),
    unsigned (*getNodesPtr) (const RDL_data*, unsigned, RDL_node**),
    unsigned (*getNofFamiliesContainingNodePtr) (const RDL_data*, RDL_node),
    unsigned (*getFamiliesContainingNodePtr) (const RDL_data*, RDL_node, unsigned**),
    unsigned (*getNofFamiliesContainingEdgePtr) (const RDL_data*, RDL_node, RDL_node),
    unsigned (*getFamiliesContainingEdgePtr) (const RDL_data*, RDL_node, RDL_node, unsigned**),
    unsigned (*getFamilyFromCycle)(const RDL_cycle*))
{
  unsigned idx1, idx2, idx3, idx4, weight,
    nof_edges, nof_nodes, node;
  RDL_edge *family_edges1, *all_edges;
  RDL_node *family_nodes1;
  unsigned nof_families;
  char **char_array1;
  RDL_cycle **cycle_array1;
  char *family_edges2, *family_nodes2;
  double total_rc = 0.0, calc_nof_cycles;

  unsigned **families_for_node1, **families_for_edge1;
  unsigned *families_for_node2, *families_for_edge2;

  unsigned nof_edges_family1, nof_edges_family2, nof_edges_family,
    nof_nodes_family1, nof_nodes_family2, nof_nodes_family;

  unsigned result = 0;

  const RDL_graph* graph = data->graph;

  RDL_getEdgeArray(data, &all_edges);

  nof_families = getNumberPtr(data);

  nof_edges = RDL_getNofEdges(data);
  nof_nodes = graph->V;

  family_edges2 = malloc(nof_edges * sizeof(*family_edges2));
  family_nodes2 = malloc(nof_nodes * sizeof(*family_nodes2));

  /* init families for nodes */
  families_for_node1 = malloc(nof_nodes * sizeof(*families_for_node1));
  for (idx1 = 0; idx1 < nof_nodes; ++idx1) {
    families_for_node1[idx1] = malloc(nof_families * sizeof(*families_for_node1[idx1]));
    for (idx2 = 0; idx2 < nof_families; ++idx2) {
      families_for_node1[idx1][idx2] = UINT_MAX;
    }
  }

  /* init families for edge */
  families_for_edge1 = malloc(nof_edges * sizeof(*families_for_edge1));
  for (idx1 = 0; idx1 < nof_edges; ++idx1) {
    families_for_edge1[idx1] = malloc(nof_families * sizeof(*families_for_edge1[idx1]));
    for (idx2 = 0; idx2 < nof_families; ++idx2) {
      families_for_edge1[idx1][idx2] = UINT_MAX;
    }
  }

  /*
   * check families
   *
   * RDL_getNofRCForURF / RDL_getNofRCForRCF
   * RDL_getRCyclesForURF / RDL_getRCyclesForRCF
   * RDL_getWeightForURF / RDL_getWeightForRCF
   * RDL_getEdgesForURF / RDL_getEdgesForRCF
   * RDL_getNodesForURF / RDL_getNodesForRCF
   *
   */
  for (idx1 = 0; idx1 < nof_families; ++idx1) {
    weight = getRCyclesPtr(data, idx1, &cycle_array1);
    weight = RDL_translateCycArray(data, cycle_array1, weight, &char_array1);

    calc_nof_cycles = getNofRCPtr(data, idx1);
    total_rc += weight;

    /*
     * validate RDL_getNofRCForURF / RDL_getNofRCForRCF consistency
     * (number of cycles calculated and enumerated)
     */
    if (weight != calc_nof_cycles) {
      fprintf(stderr, "mismatch between nof RCs calculated/retrieved: %d vs. %.0f\n", weight, calc_nof_cycles);
      result = 1;
    }

    /* there must be cycles! */
    if (weight == 0) {
      fprintf(stderr, "empty cycle set found!\n");
      result = 1;
    }

    memset(family_edges2, 0, nof_edges * sizeof(*family_edges2));
    memset(family_nodes2, 0, nof_nodes * sizeof(*family_nodes2));

    for (idx2 = 0; idx2 < weight; ++idx2) {
      /* check RDL_cycle family consistency */
      if (getFamilyFromCycle(cycle_array1[idx2]) != idx1) {
        fprintf(stderr, "family number mismatch: %d vs. %d\n", idx1,
            getFamilyFromCycle(cycle_array1[idx2]));
        result = 1;
      }

      /* check that no cycle is empty */
      if (cycle_array1[idx2]->weight == 0) {
        fprintf(stderr, "empty cycle found!\n");
        result = 1;
      }

      /*
       * check consistency of RDL_cycle.weight
       * and RDL_getWeightForURF / RDL_getWeightForRCF
       */
      if (cycle_array1[idx2]->weight != getWeightPtr(data, idx1)) {
        fprintf(stderr, "Family weight mismatch!\n");
        result = 1;
      }

      for (idx3 = 0; idx3 < cycle_array1[idx2]->weight; ++idx3) {
        /* check EXPLICITLY which nodes and edges are in the current cycle */
        node = cycle_array1[idx2]->edges[idx3][0];
        family_nodes2[node] |= 1;
        for (idx4 = 0; idx4 < nof_families; ++idx4) {
          if (families_for_node1[node][idx4] == idx1) {
            break;
          }
          else if (families_for_node1[node][idx4] == UINT_MAX) {
            families_for_node1[node][idx4] = idx1;
            break;
          }
        }

        node = cycle_array1[idx2]->edges[idx3][1];
        family_nodes2[node] |= 1;
        for (idx4 = 0; idx4 < nof_families; ++idx4) {
          if (families_for_node1[node][idx4] == idx1) {
            break;
          }
          else if (families_for_node1[node][idx4] == UINT_MAX) {
            families_for_node1[node][idx4] = idx1;
            break;
          }
        }
      }

      for (idx3 = 0; idx3 < nof_edges; ++idx3) {
        family_edges2[idx3] |= char_array1[idx2][idx3];
        if (char_array1[idx2][idx3]) {
          for (idx4 = 0; idx4 < nof_families; ++idx4) {
            if (families_for_edge1[idx3][idx4] == idx1) {
              break;
            }
            else if (families_for_edge1[idx3][idx4] == UINT_MAX) {
              families_for_edge1[idx3][idx4] = idx1;
              break;
            }
          }
        }

      }
    }

    RDL_deleteCycles(cycle_array1, weight);
    RDL_delete2DCharArray(char_array1, weight);

    /* now enumerate the edges for the current family */
    nof_edges_family1 = getEdgesPtr(data, idx1, &family_edges1);
    nof_edges_family2 = 0;
    for (idx3 = 0; idx3 < nof_edges; ++idx3) {
      if (family_edges2[idx3]) {
        ++nof_edges_family2;
      }
    }

    if (nof_edges_family1 == 0) {
      fprintf(stderr, "family seems to be empty!\n");
      result = 1;
    }

    /* check that the number of edges in the current family is the same */
    if (nof_edges_family1 != nof_edges_family2) {
      fprintf(stderr, "mismatch of number of edges in family %d: %d vs. %d\n",
          idx1, nof_edges_family1, nof_edges_family2);
      result = 1;
    }

    /*
     * and finally verify that the edges are the same
     * (it is sufficient to compare one direction and the number!)
     */
    for (idx3 = 0; idx3 < nof_edges_family1; ++idx3) {
      if (!family_edges2[RDL_edgeId(graph, family_edges1[idx3][0], family_edges1[idx3][1])]) {
        fprintf(stderr, "mismatch of edges in family %d: %d missing\n",
            idx1, RDL_edgeId(graph, family_edges1[idx3][0], family_edges1[idx3][1]));
        result = 1;
        break;
      }
    }

    nof_nodes_family1 = getNodesPtr(data, idx1, &family_nodes1);

    if (nof_nodes_family1 == 0) {
      fprintf(stderr, "family seems to be empty!\n");
      result = 1;
    }

    nof_nodes_family2 = 0;
    for (idx3 = 0; idx3 < nof_nodes; ++idx3) {
      if (family_nodes2[idx3]) {
        ++nof_nodes_family2;
      }
    }

    /* check that the number of nodes in the current family is the same */
    if (nof_nodes_family1 != nof_nodes_family2) {
      fprintf(stderr, "mismatch of number of nodes in family %d: %d vs. %d\n",
          idx1, nof_nodes_family1, nof_nodes_family2);
      result = 1;
    }

    /*
     * and finally verify that the nodes are the same
     * (it is sufficient to compare one direction and the number!)
     */
    for (idx3 = 0; idx3 < nof_nodes_family1; ++idx3) {
      if (!family_nodes2[family_nodes1[idx3]]) {
        fprintf(stderr, "mismatch of edges in family %d: %d missing\n",
            idx1, family_nodes1[idx3]);
        result = 1;
        break;
      }
    }

    free(family_edges1);
    free(family_nodes1);
  }

  free(family_edges2);
  free(family_nodes2);

  if (total_rc != RDL_getNofRC(data)) {
    fprintf(stderr, "inconsistency in the number of RCs!\n");
    result = 1;
  }

  /*
   * now verify the inverse direction for families
   *
   * RDL_getURFsContainingNode / RDL_getRCFsContainingNode
   * RDL_getNofURFsContainingNode / RDL_getNofRCFsContainingNode
   */
  for (idx3 = 0; idx3 < nof_nodes; ++idx3) {
    for (idx4 = 0; idx4 < nof_families; ++idx4) {
      if (families_for_node1[idx3][idx4] == UINT_MAX) {
        break;
      }
    }
    nof_nodes_family = getFamiliesContainingNodePtr(data, idx3, &families_for_node2);

    if (nof_nodes_family != getNofFamiliesContainingNodePtr(data, idx3)) {
      fprintf(stderr, "mismatch between the numbers of families for node!\n");
      result = 1;
    }

    /* check the number of families for node is equal */
    if (nof_nodes_family != idx4) {
      fprintf(stderr, "mismatch between number of families for node %d: %d vs. %d\n",
          idx3, nof_nodes_family, idx4);
      result = 1;
    }
    else {
      /* sort and compare one by one */
      qsort(families_for_node2, nof_nodes_family, sizeof(*families_for_node2), cmp_unsigned);
      qsort(families_for_node1[idx3], nof_nodes_family, sizeof(*families_for_node1[idx3]), cmp_unsigned);

      for (idx4 = 0; idx4 < nof_nodes_family; ++idx4) {
        if (families_for_node1[idx3][idx4] != families_for_node2[idx4]) {
          fprintf(stderr, "mismatch between families for node %d: %d vs %d\n",
              idx3, families_for_node1[idx3][idx4], families_for_node2[idx4]);
          result = 1;
        }
      }
    }

    free(families_for_node2);
    free(families_for_node1[idx3]);
  }

  free(families_for_node1);

  /*
   * now verify the inverse direction for families
   *
   * RDL_getURFsContainingEdge / RDL_getRCFsContainingEdge
   * RDL_getNofURFsContainingEdge / RDL_getNofRCFsContainingEdge
   */
  for (idx3 = 0; idx3 < nof_edges; ++idx3) {
    for (idx4 = 0; idx4 < nof_families; ++idx4) {
      if (families_for_edge1[idx3][idx4] == UINT_MAX) {
        break;
      }
    }

    nof_edges_family = getFamiliesContainingEdgePtr(data, all_edges[idx3][0],
        all_edges[idx3][1], &families_for_edge2);

    if (nof_edges_family != getNofFamiliesContainingEdgePtr(data, all_edges[idx3][0],
        all_edges[idx3][1])) {
      fprintf(stderr, "mismatch between the numbers of families for edge!\n");
      result = 1;
    }

    /* check the number of families for edge is equal */
    if (nof_edges_family != idx4) {
      fprintf(stderr, "mismatch between number of families for edge %d: %d vs. %d\n",
          idx3, nof_edges_family, idx4);
      result = 1;
    }
    else {
      /* sort and compare one by one */
      qsort(families_for_edge2, nof_edges_family, sizeof(*families_for_edge2), cmp_unsigned);
      qsort(families_for_edge1[idx3], nof_edges_family, sizeof(*families_for_edge1[idx3]), cmp_unsigned);

      for (idx4 = 0; idx4 < nof_edges_family; ++idx4) {
        if (families_for_edge1[idx3][idx4] != families_for_edge2[idx4]) {
          fprintf(stderr, "mismatch between families for edge %d: %d vs %d\n",
              idx3, families_for_edge1[idx3][idx4], families_for_edge2[idx4]);
          result = 1;
        }
      }
    }

    free(families_for_edge2);
    free(families_for_edge1[idx3]);
  }

  free(families_for_edge1);
  free(all_edges);

  return result;
}

int RDL_checkRingsystemConsistency(const RDL_data* data, char graph_connected, unsigned z)
{
  unsigned idx1, idx2, idx3, nof_rs_edges;
  unsigned result = 0;

  RDL_edge *rs_edges;
  unsigned rs_edge;
  unsigned *rs_edge_mapping;

  RDL_node *rs_nodes;
  unsigned nof_rs_nodes1, nof_rs_edges2;

  /* verify size of ringsystems RDL_getNofRingsystems for connected graphs */
  if (graph_connected && z && !RDL_getNofRingsystems(data)) {
    fprintf(stderr, "graph has z > 0 but no ringsystem!\n");
    result = 1;
  }

  rs_edge_mapping = malloc(data->graph->E * sizeof(*rs_edge_mapping));
  for (idx1 = 0; idx1 < data->graph->E; ++idx1) {
    rs_edge_mapping[idx1] = RDL_NO_RINGSYSTEM;
  }

  /*
   * check ringsystems
   *
   * RDL_getEdgesForRingsystem
   * RDL_getNofNodesForRingsystem
   *
   */
  for (idx1 = 0; idx1 < RDL_getNofRingsystems(data); ++idx1) {
    nof_rs_edges = RDL_getEdgesForRingsystem(data, idx1, &rs_edges);

    if (nof_rs_edges - RDL_getNofNodesForRingsystem(data, idx1) + 1 == 0) {
      fprintf(stderr, "ring system has z = 0!\n");
      result = 1;
    }

    if (!nof_rs_edges) {
      fprintf(stderr, "empty ring system!\n");
      result = 1;
    }

    for (idx2 = 0; idx2 < nof_rs_edges; ++idx2) {
      rs_edge = RDL_edgeId(data->graph, rs_edges[idx2][0], rs_edges[idx2][1]);
      if (rs_edge_mapping[rs_edge] != RDL_NO_RINGSYSTEM) {
        fprintf(stderr, "inconsistent ringsystem results: %d\n", idx1);
        result = 1;
      }
      else {
        rs_edge_mapping[rs_edge] = idx1;
      }
    }

    nof_rs_nodes1 = RDL_getNodesForRingsystem(data, idx1, &rs_nodes);
    nof_rs_edges2 = 0;

    for (idx2 = 0; idx2 < nof_rs_nodes1; ++idx2) {
      for (idx3 = idx2 + 1; idx3 < nof_rs_nodes1; ++idx3) {
        rs_edge = RDL_edgeId(data->graph, rs_nodes[idx2], rs_nodes[idx3]);
        if (rs_edge != RDL_INVALID_RESULT) {
          ++nof_rs_edges2;
          if (rs_edge_mapping[rs_edge] != idx1) {
            fprintf(stderr, "inconsistent ringsystem results: %d\n", idx1);
            result = 1;
          }
        }
      }
    }

    if (nof_rs_edges != nof_rs_edges2) {
      fprintf(stderr, "inconsistent ringsystem %d results: %d vs. %d \n",
          idx1, nof_rs_edges, nof_rs_edges2);
      result = 1;
    }

    free(rs_nodes);

    free(rs_edges);
  }
  for (idx1 = 0; idx1 < data->graph->E; ++idx1) {
    idx2 = RDL_getRingsystemForEdge(data,
        data->graph->edges[idx1][0], data->graph->edges[idx1][1]);
    if (rs_edge_mapping[idx1] != idx2) {
      fprintf(stderr, "inconsistent ringsystem results for edge: %d\n", idx1);
      result = 1;
    }
  }
  free(rs_edge_mapping);

  return result;
}

int RDL_checkRCPConsistency(const RDL_data* data, unsigned z)
{
  unsigned result = 0;
  unsigned nof_cycles1, nof_cycles2, idx1, idx2, idx3, idx4;
  unsigned equal, no_intersect;
  char **char_array1, **char_array2;
  RDL_cycle **cycle_array1, **cycle_array2;

  nof_cycles1 = RDL_getRCPrototypes(data, &cycle_array1);
  nof_cycles1 = RDL_translateCycArray(data, cycle_array1, nof_cycles1, &char_array1);

  if (z && !nof_cycles1) {
    fprintf(stderr, "no RCPs found, but z != 0\n");
    return 1;
  }

  for (idx1 = 0; idx1 < data->nofRCFs; ++idx1) {
    nof_cycles2 = RDL_getRCyclesForRCF(data, idx1, &cycle_array2);
    nof_cycles2 = RDL_translateCycArray(data, cycle_array2, nof_cycles2, &char_array2);

    for (idx2 = 0; idx2 < nof_cycles2; ++idx2) {

      if (cycle_array2[idx2]->rcf != idx1) {
        fprintf(stderr, "RCF number mismatch!");
        result = 1;
      }

      for (idx3 = 0; idx3 < nof_cycles1; ++idx3) {
        equal = 1;
        no_intersect = 1;

        for (idx4 = 0; idx4 < data->graph->E; ++idx4) {
          if (char_array2[idx2][idx4] != char_array1[idx3][idx4]) {
            equal = 0;
          }
          else {
            no_intersect = 0;
          }
        }
        /* verify that the RCP has the correct rcf ID */
        if (idx3 != cycle_array2[idx2]->rcf && equal) {
          fprintf(stderr, "different RCF number but same ring!\n");
          result = 1;
        }
        if (idx3 == cycle_array2[idx2]->rcf && no_intersect) {
          fprintf(stderr, "same RCF number but no shared edge!\n");
          result = 1;
        }
      }
    }

    RDL_delete2DCharArray(char_array2, nof_cycles2);
    RDL_deleteCycles(cycle_array2, nof_cycles2);
  }

  RDL_deleteCycles(cycle_array1, nof_cycles1);
  RDL_delete2DCharArray(char_array1, nof_cycles1);

  return result;
}

static unsigned RDL_getURFFromCycle(const RDL_cycle *cycle)
{
  return cycle->urf;
}

static unsigned RDL_getRCFFromCycle(const RDL_cycle *cycle)
{
  return cycle->rcf;
}

int RDL_checkConsistency(const RDL_data* data)
{
  unsigned nof_urf, idx1, nof_edges, nof_nodes,
    nof_cycles1, result = 0, z;
  RDL_cycle **cycle_array1;

  const RDL_graph* graph = data->graph;
  char graph_connected;
  RDL_cycleIterator *it;


  nof_urf = RDL_getNofURF(data);
  nof_edges = RDL_getNofEdges(data);

  if (nof_edges != graph->E) {
    fprintf(stderr, "different number of edges!");
    return 1;
  }

  nof_nodes = graph->V;
  z = nof_edges - nof_nodes + 1;
  graph_connected = RDL_checkGraphConnected(data->graph);

  /* there must be an URF if cyclomatic number != 0 */
  if (z && !nof_urf) {
    fprintf(stderr, "no URFs but z != 0!\n");
    return 1;
  }

  if (RDL_checkInvalidConsistency(data)) {
    fprintf(stderr, "Handling of invalid values failed!\n");
    result = 1;
  }
  else {
    fprintf(stderr, "Handling of invalid values successful!\n");
  }

  if (RDL_checkFamilyConsistency(
      data, RDL_getNofURF, RDL_getNofRCForURF,
      RDL_getRCyclesForURF, RDL_getWeightForURF,
      RDL_getEdgesForURF, RDL_getNodesForURF,
      RDL_getNofURFContainingNode, RDL_getURFsContainingNode,
      RDL_getNofURFContainingEdge, RDL_getURFsContainingEdge,
      RDL_getURFFromCycle)) {
    fprintf(stderr, "URF results inconsistent!\n");
    result = 1;
  }
  else {
    fprintf(stderr, "URF results successful!\n");
  }

  if (RDL_checkFamilyConsistency(
      data, RDL_getNofRCF, RDL_getNofRCForRCF,
      RDL_getRCyclesForRCF, RDL_getWeightForRCF,
      RDL_getEdgesForRCF, RDL_getNodesForRCF,
      RDL_getNofRCFContainingNode, RDL_getRCFsContainingNode,
      RDL_getNofRCFContainingEdge, RDL_getRCFsContainingEdge,
      RDL_getRCFFromCycle)) {
    fprintf(stderr, "RCF results inconsistent!\n");
    result = 1;
  }
  else {
    fprintf(stderr, "RCF results successful!\n");
  }

  /* no actual test, detect memory leaks if premature aborting iteration */
  it = RDL_getRCyclesIterator(data);
  for(idx1 = 0; idx1 < 7 && !RDL_cycleIteratorAtEnd(it); ++idx1, RDL_cycleIteratorNext(it)) {
    ;
  }
  RDL_deleteCycleIterator(it);

  /* check RDL_getRCPrototypes */
  if (RDL_checkRCPConsistency(data, z)) {
    fprintf(stderr, "RCP results inconsistent!\n");
    result = 1;
  }
  else {
    fprintf(stderr, "RCP results successful!\n");
  }

  /* verify RDL_getSSSR for connected graphs (correct size) */
  nof_cycles1 = RDL_getSSSR(data, &cycle_array1);

  if (graph_connected && nof_cycles1 != z) {
    fprintf(stderr, "unexpected SSSR size, z = %u, size = %u!\n",
        data->graph->E - data->graph->V + 1, nof_cycles1);
    result = 1;
  }
  RDL_deleteCycles(cycle_array1, nof_cycles1);

  if (RDL_checkRingsystemConsistency(data, graph_connected, z)) {
    fprintf(stderr, "Ringsystem results inconsistent!\n");
    result = 1;
  }
  else {
    fprintf(stderr, "Ringsystem results successful!\n");
  }

  return result;

}

int RDL_validateBasis(const char** relevant_cycle_array,
                      unsigned nof_rc,
                      const char** basis_cycle_array,
                      unsigned nof_basis,
                      unsigned nof_bonds,
                      unsigned *bcc,
                      int (*stop) (const void*))
{
  /*
   * basically enumerates all combinations of the basis cycles
   * until all relevant cycles are found => because the RCs
   * are a basis (the union of all MCBs), our basis was also one
   */
  const char *current_bonds, *ring;
  char *new_bonds;
  unsigned current_ring_idx, i, nof_found_rc=0;
  int cmp, mid, lower, upper;
  RDL_stack* stack;
  RDL_stackElement *top, *new_element;
  unsigned curr_bcc, found;
  char **sortable_cycle_array;
  char *found_rc;
  int result = 0;

  found_rc = malloc(nof_rc * sizeof(*found_rc));
  memset(found_rc, 0, nof_rc * sizeof(*found_rc));

  sortable_cycle_array = malloc(nof_rc * sizeof(*sortable_cycle_array));
  for (i = 0; i < nof_rc; ++i) {
    sortable_cycle_array[i] = malloc((nof_bonds + 1) * sizeof(**sortable_cycle_array));
    memcpy(sortable_cycle_array[i], relevant_cycle_array[i], nof_bonds * sizeof(**sortable_cycle_array));
    sortable_cycle_array[i][nof_bonds] = 2;
  }

  qsort(sortable_cycle_array, nof_rc, sizeof(*sortable_cycle_array), RDL_cmp_cycles);

  stack = RDL_stack_new();

  top = malloc(sizeof(*top));
  top->current_ring_idx = 0;
  new_bonds = malloc((nof_bonds + 1) * sizeof(*new_bonds));
  memset(new_bonds, 0, (nof_bonds + 1) * sizeof(*new_bonds));
  new_bonds[nof_bonds] = 2;
  top->current_bonds = new_bonds;
  top->stage = 0;
  top->owns_bonds = 1;
  top->curr_bcc = UINT_MAX;
  RDL_stack_push(stack, top);

  while (!RDL_stack_empty(stack) && nof_found_rc < nof_rc) {
    top = RDL_stack_top(stack);
    current_bonds = top->current_bonds;
    current_ring_idx = top->current_ring_idx;
    curr_bcc = top->curr_bcc;

    found = 0;
    /* search ring in relevant cycles using binary search */
    lower = 0;
    upper = nof_rc - 1;
    while (lower <= upper && !found) {
      mid = lower + (upper - lower) / 2;
      cmp = RDL_cmp_cycles(&current_bonds, sortable_cycle_array + mid);
      if (cmp == 0) {
        found = 1;
      }
      else if (cmp == -1) {
        upper = mid - 1;
      }
      else if (cmp == 1) {
        lower = mid + 1;
      }
    }

    if (found && !found_rc[mid]) {
      found_rc[mid] = 1;
      ++nof_found_rc;
    }

    if (top->stage < 2 && current_ring_idx < nof_basis) {
      /* looking at this ring the first time, do not add it and recurse */
      if (top->stage == 0) {
        top->stage = 1;

        new_element = malloc(sizeof(*new_element));
        new_element->stage = 0;
        new_element->current_bonds = current_bonds;
        new_element->current_ring_idx = current_ring_idx + 1;
        new_element->owns_bonds = 0;
        new_element->curr_bcc = curr_bcc;
        RDL_stack_push(stack, new_element);
      }
      /* looking at this ring the second time, add it and recurse */
      else if (top->stage == 1) {
        top->stage = 2;

        if (curr_bcc == UINT_MAX || curr_bcc == bcc[current_ring_idx]) {
          ring = basis_cycle_array[current_ring_idx];
          new_bonds = malloc((nof_bonds + 1) * sizeof(*new_bonds));
          for (i = 0; i < nof_bonds; ++i) {
            new_bonds[i] = !ring[i] != !current_bonds[i] ? 1 : 0; /* logical XOR */
          }
          new_bonds[nof_bonds] = 2;
          new_element = malloc(sizeof(*new_element));
          new_element->stage = 0;
          new_element->current_bonds = new_bonds;
          new_element->current_ring_idx = current_ring_idx + 1;
          new_element->owns_bonds = 1;
          new_element->curr_bcc = bcc[current_ring_idx];
          RDL_stack_push(stack, new_element);
        }
      }
    }
    /* looking at this ring for the third time or at an non existing ring, end recursion */
    else {
      if (top->owns_bonds) {
        free((char*)top->current_bonds);
      }
      free(top);
      RDL_stack_pop(stack);
    }
    if (stop(NULL)) {
      break;
    }
  }

  while (!RDL_stack_empty(stack)) {
    top = RDL_stack_top(stack);
    if (top->owns_bonds) {
      free((char*)top->current_bonds);
    }
    free(top);
    RDL_stack_pop(stack);
  }

  RDL_stack_delete(stack);

  if (nof_found_rc != nof_rc) {
    result = 1;
  }

  if (stop(NULL)) {
    result = -1;
  }

  free(found_rc);
  for (i = 0; i < nof_rc; ++i) {
    free(sortable_cycle_array[i]);
  }
  free(sortable_cycle_array);

  return result;
}
