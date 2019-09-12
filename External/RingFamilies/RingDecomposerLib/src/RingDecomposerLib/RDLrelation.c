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

#include "RDLrelation.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "RDLapsp.h"
#include "RDLcycleFams.h"
#include "RDLgraph.h"
#include "RDLhandler.h"
#include "RDLinfo.h"
#include "RDLstack.h"
#include "RDLutility.h"
#include "RDLbitset.h"

/** returns the index in the array of RCFs (fams) that the RCF with the weight
at the index "weight" and position j has */
static unsigned RDL_idxWeights(RDL_URFinfo *uInfo, unsigned weight, unsigned j)
{
  unsigned i, sum = 0;
  for(i=0; i<weight; ++i)
  {
    sum += uInfo->nofProtos[i];
  }
  return sum + j;
}

/** Writes into the URFinfo if two RCFs are potentially URF-related because of
linear dependencies (condition 3 URF-relation). In this step all relevant cycle
families are marked as well.*/
/* 3 sets are maintained (see also Vismara):
  'B' is the current set of basis prototypes (corresponding to a RCF)
  'B<' is the set of cycles in B which are shorter than the cycles we are
  looking at at that moment. RCFs in 'B<' are marked with 1.
  'B=' is the set of cycles in B which have the same length as the cycles we are
  looking at at that moment. RCFs in 'B=' are marked with 2 or 3.
To see whether or not two RCFs are potentially URF-related we look at each RCF
'C' of each weight individually:
  If C is independent of 'B<' we look at it's dependency on all combinations
  of 'B<' with one of 'B='.
  If it is independent of all those combinations we look at if it is independent
  of the combination of 'B<' and all of 'B='.
  If this is the case, it becomes part of the basis B (part of 'B=')
  If we find that it is RDL_dependent on 'B<' and one of 'B=', the element of
  'B=' and C are marked as potentially URF-related.*/
void RDL_checkDependencies(RDL_cfURF *RCFs, RDL_graph *graph, RDL_URFinfo *uInfo)
{
  unsigned i,j,k, index;

  /* this array contains the basis cycles in row echelon form */
  unsigned char **compressedBasisCycles;
  unsigned maxBasisCyclesSize = graph->E - graph->V + 1;
  unsigned currentBasisCyclesSize = 0;
  unsigned currentBasisCyclesSmallerEnd = 0;
  unsigned oldBasisCyclesSize;
  unsigned compressedSize = 0;

  unsigned char **compressedRelevantCycles;
  unsigned currentRelevantCyclesSize = 0;
  unsigned currentRelevantCyclesSmallerEnd = 0;
  unsigned *relevantCyclesMap;

  unsigned char *compressedCycleWithSmallerAdded;
  unsigned char *compressedCycleWithEqualAdded;
  unsigned char **compressedPrototypesForGaussian;

  const unsigned char* empty_cycle;

  if(RCFs->nofFams < 3) {
    /*if only 0, 1 or 2 families exist, they are all relevant and independent*/
    for(i=0; i<uInfo->nofWeights; ++i) {
      for(j=0; j<uInfo->nofProtos[i]; ++j) {
        uInfo->URFrel[i][j][j] = 1; /*URF-related to itself*/
      }
    }
    for(i=0; i<RCFs->nofFams; ++i) {
      RCFs->fams[i]->mark = 1;
    }
    return;
  }

  compressedBasisCycles = malloc(maxBasisCyclesSize * sizeof(*compressedBasisCycles));
  compressedPrototypesForGaussian = malloc(RCFs->nofFams * sizeof(*compressedPrototypesForGaussian));
  compressedRelevantCycles = malloc(RCFs->nofFams * sizeof(*compressedRelevantCycles));
  relevantCyclesMap = malloc(RCFs->nofFams * sizeof(*relevantCyclesMap));

  /* copy the prototypes for gaussian elimination */
  for (i = 0; i < RCFs->nofFams; ++i) {
    compressedSize = RDL_bitset_compressed(
        &(compressedPrototypesForGaussian[i]), RCFs->fams[i]->prototype, graph->E);
  }
  empty_cycle = malloc(compressedSize * sizeof(*empty_cycle));
  memset((unsigned char*)empty_cycle, 0, compressedSize * sizeof(*empty_cycle));
  /* printf("%d %d\n", compressedSize, graph->E); */

  /* iterate over all weights */
  for (i = 0; i < uInfo->nofWeights; ++i) {
    /* reset indices indication the first index of larger cycles */
    currentBasisCyclesSmallerEnd = currentBasisCyclesSize;
    currentRelevantCyclesSmallerEnd = currentRelevantCyclesSize;
    for (j = 0; j < uInfo->nofProtos[i]; ++j) {

      index = RDL_idxWeights(uInfo, i, j);
      /* copy current cycle */
      compressedCycleWithSmallerAdded = malloc(compressedSize * sizeof(*compressedCycleWithSmallerAdded));
      memcpy(compressedCycleWithSmallerAdded, compressedPrototypesForGaussian[index],
          compressedSize * sizeof(*compressedCycleWithSmallerAdded));

      /*
       * add smaller cycles of the basis (set B<)
       *
       * the basisCycles are in row echelon form (see below), so we add each cycle that
       * has a 1 at position k (all smaller positions are 0, row echelon)
       *
       * [this corresponds to the "add row operation" of gaussian elimination]
       */
      for (k = 0; k < currentBasisCyclesSmallerEnd; ++k) {
        if (RDL_bitset_test(compressedCycleWithSmallerAdded, k)) {
          RDL_bitset_xor_inplace(compressedCycleWithSmallerAdded,
              compressedBasisCycles[k], compressedSize);
        }
      }

      /*
       * if the result is empty, then this cycle linearly depends
       * on smaller cycles of the basis B<
       *
       * => trivial
       * <= is also true because adding more smaller cycles can't make
       * the cycle 0 (how could it?)
       */
      if (RDL_bitset_empty(compressedCycleWithSmallerAdded, empty_cycle, compressedSize)) {
        free(compressedCycleWithSmallerAdded);
        continue;
      }

      /* this cycle does not depend on strictly smaller cycles => relevant */
      compressedRelevantCycles[currentRelevantCyclesSize] = compressedCycleWithSmallerAdded;
      relevantCyclesMap[currentRelevantCyclesSize] = j; /* this is no bug, we really need the "internal" index */
      ++currentRelevantCyclesSize;

      RCFs->fams[index]->mark = 1;
      uInfo->URFrel[i][j][j] = 1; /*URF-related to itself*/

      /* make a copy for adding the equal sized cycles */
      compressedCycleWithEqualAdded = malloc(compressedSize * sizeof(*compressedCycleWithEqualAdded));
      memcpy(compressedCycleWithEqualAdded, compressedCycleWithSmallerAdded,
          compressedSize * sizeof(*compressedCycleWithSmallerAdded));

      /*
       * add equal sized cycles of the basis (set B=)
       *
       * the same arguments as for the smaller cycles apply
       */
      for (k = currentBasisCyclesSmallerEnd; k < currentBasisCyclesSize; ++k) {
        if (RDL_bitset_test(compressedCycleWithEqualAdded, k)) {
          RDL_bitset_xor_inplace(compressedCycleWithEqualAdded,
              compressedBasisCycles[k], compressedSize);
        }
      }

      if (RDL_bitset_empty(compressedCycleWithEqualAdded, empty_cycle, compressedSize)) {
        /*
         * if the cycle depends linearly on any equal sized cycle of the
         * basis (B=), we check all RELEVANT cycles we've seen so far for
         * linear dependence
         *
         * (-1 because we don't have to add the ring to itself)
         */

        for (k = currentRelevantCyclesSmallerEnd; k < currentRelevantCyclesSize-1; ++k) {
          memcpy(compressedCycleWithEqualAdded, compressedCycleWithSmallerAdded,
              compressedSize * sizeof(*compressedCycleWithSmallerAdded));
          RDL_bitset_xor_inplace(compressedCycleWithEqualAdded,
              compressedRelevantCycles[k], compressedSize);
          if (RDL_bitset_empty(compressedCycleWithEqualAdded, empty_cycle, compressedSize)) {
            /*
             * if there is any equal sized cycle that is in fact equal to the current cycle
             * (after adding smaller cycles) => those cycles fulfill URF condition 1 (equal size)
             * and also condition 3 (depended on each other and smaller cycles)
             */
            uInfo->URFrel[i][j][relevantCyclesMap[k]] = 1;
            uInfo->URFrel[i][relevantCyclesMap[k]][j] = 1;
          }
        }

        free(compressedCycleWithEqualAdded);
      }
      else {
        /*
         * the current cycle is also independent of equal sized cycles,
         * so it's part of the basis, add it to B=
         *
         * (so it can't be depedent of any equal sized ring seen so far,
         * so URF relation is out of discussion)
         */
        oldBasisCyclesSize = currentBasisCyclesSize;
        compressedBasisCycles[currentBasisCyclesSize] = compressedCycleWithEqualAdded;
        ++currentBasisCyclesSize;

        /*
         * this is where the magic is happening: ensure that the basis is
         * indeed in row echelon form (we use the property throughout the algorithm)
         *
         * swap columns such that indeed the current ring with index currentBasisCycleSize-1
         * has a 1 at column currentBasisCycleSize-1
         * => ROW ECHELON FORM
         *
         * [this corresponds to the "swap columns operation" of gaussian elimination]
         */
        if (!RDL_bitset_test(compressedCycleWithEqualAdded, oldBasisCyclesSize)) {
          for (k = oldBasisCyclesSize + 1; k < graph->E; ++k) {
            if (RDL_bitset_test(compressedCycleWithEqualAdded, k)) {
              RDL_swap_columns(compressedBasisCycles, currentBasisCyclesSize, oldBasisCyclesSize, k);
              RDL_swap_columns(compressedRelevantCycles, currentRelevantCyclesSize, oldBasisCyclesSize, k);
              RDL_swap_columns(compressedPrototypesForGaussian, RCFs->nofFams, oldBasisCyclesSize, k);
              break;
            }
          }
        }
      }

    }
  }

  for (i = 0; i < currentBasisCyclesSize; ++i) {
    free(compressedBasisCycles[i]);
  }
  free(compressedBasisCycles);

  for (i = 0; i < currentRelevantCyclesSize; ++i) {
    free(compressedRelevantCycles[i]);
  }
  free(compressedRelevantCycles);
  free(relevantCyclesMap);

  for (i = 0; i < RCFs->nofFams; ++i) {
    free(compressedPrototypesForGaussian[i]);
  }
  free(compressedPrototypesForGaussian);
  free((void*)empty_cycle);
}

/** finds all edges in the RCF and stores the result in the array edges.*/
void RDL_findEdges(char *edges, RDL_cfam *RCF, RDL_graph *gra, RDL_sPathInfo *spi)
{
  char* visited;
  visited = malloc(gra->E * sizeof(*visited));
  memset(visited, 0, gra->E * sizeof(*visited));
  RDL_giveEdges(RCF->r, RCF->p, edges, gra, spi, visited);
  memset(visited, 0, gra->E * sizeof(*visited));
  RDL_giveEdges(RCF->r, RCF->q, edges, gra, spi, visited);
  /* we need to add the edges */
  if(RCF->x < UINT_MAX) /*even family*/
  {
    edges[RDL_edgeId(gra, RCF->x, RCF->p)] = 1;
    edges[RDL_edgeId(gra, RCF->x, RCF->q)] = 1;
  }
  else
  {
    edges[RDL_edgeId(gra, RCF->p, RCF->q)] = 1;
  }
  free(visited);
}

/** Checks if the RCFs with indices idx1 and idx2 share at least one edge.
returns 1 if yes and 0 otherwise. */
char RDL_shareEdges(RDL_cfURF *RCFs, unsigned idx1, unsigned idx2, RDL_graph *graph,
                    RDL_sPathInfo *spi)
{
  char *edges1;
  char *edges2; /*arrays of {0,1}^m that can store a set of edges with entries
  being 1 if the edge is contained and 0 otherwise*/
  unsigned i;
  char result = 0;
  edges1 = malloc(graph->E * sizeof(*edges1));
  edges2 = malloc(graph->E * sizeof(*edges2));

  memset(edges1, 0, graph->E * sizeof(*edges1));
  memset(edges2, 0, graph->E * sizeof(*edges2));

  RDL_findEdges(edges1, RCFs->fams[idx1], graph, spi);
  RDL_findEdges(edges2, RCFs->fams[idx2], graph, spi);

  for(i=0; i<graph->E; ++i)
  {
    if(edges1[i] == 1 && edges2[i] == 1)
    {
      result = 1;
      break;
    }
  }
  free(edges1);
  free(edges2);
  return result;
}

static void make_edge_list(
    unsigned** edge_list,
    unsigned* edge_list_size,
    char* edges, RDL_graph *graph,
    RDL_cfam* rcf, RDL_sPathInfo* spi)
{
  unsigned alloced, i;

  memset(edges, 0, graph->E * sizeof(*edges));
  RDL_findEdges(edges, rcf, graph, spi);

  alloced = 64;
  *edge_list = malloc(alloced * sizeof(**edge_list));

  for (i = 0; i < graph->E; ++i) {
    if (edges[i]) {
      if (*edge_list_size == alloced) {
        alloced *= 2;
        *edge_list = realloc(*edge_list, alloced * sizeof(**edge_list));
      }
      (*edge_list)[*edge_list_size] = i;
      ++(*edge_list_size);
    }
  }
  *edge_list = realloc(*edge_list, *edge_list_size * sizeof(**edge_list));
}

/** checks the previously 'marked as potentially URF-related' RCFs for edges
that they have in common which proves that they are URF-related. */
void RDL_checkEdges(RDL_cfURF *RCFs, RDL_graph *graph, RDL_URFinfo *uInfo, RDL_sPathInfo* spi)
{
  unsigned i,j,k,shared, list_idx1, list_idx2;
  char *edges; /*arrays of {0,1}^m that can store a set of edges with entries
  being 1 if the edge is contained and 0 otherwise*/

  unsigned **edge_list;
  unsigned *edge_list_size;

  edges = malloc(graph->E * sizeof(*edges));

  for(i=0; i<uInfo->nofWeights; ++i) /*go through all matrices*/
  {
    edge_list = malloc(uInfo->nofProtos[i] * sizeof(*edge_list));
    edge_list_size = malloc(uInfo->nofProtos[i] * sizeof(*edge_list_size));
    for(j=0; j<uInfo->nofProtos[i]; ++j) {
      edge_list[j] = NULL;
      edge_list_size[j] = 0;
    }
    for(j=0; j<uInfo->nofProtos[i]; ++j) {
      for(k=j+1; k<uInfo->nofProtos[i]; ++k) /*since matrix is symmetric: only
                                              look into top right half*/
      {
        if(uInfo->URFrel[i][j][k] == 1) /*find entries which indicate potential
                                        URF-relations*/
        {
          /* calculate the edge lists only if necessary */
          if (edge_list[j] == NULL) {
            make_edge_list(&edge_list[j], &edge_list_size[j], edges,
                graph, RCFs->fams[RDL_idxWeights(uInfo,i,j)], spi);
          }
          if (edge_list[k] == NULL) {
            make_edge_list(&edge_list[k], &edge_list_size[k], edges,
                graph, RCFs->fams[RDL_idxWeights(uInfo,i,k)], spi);
          }

          shared = 0;
          list_idx1 = list_idx2 = 0;
          while (list_idx1 < edge_list_size[j] && list_idx2 < edge_list_size[k] && !shared) {
            if (edge_list[j][list_idx1] < edge_list[k][list_idx2]) {
              ++list_idx1;
            }
            else if (edge_list[j][list_idx1] > edge_list[k][list_idx2]) {
              ++list_idx2;
            }
            else {
             shared = 1;
            }
          }

          /* if shared edge, relation is confirmed */
          if (shared) {
            uInfo->URFrel[i][j][k] = 1;
            uInfo->URFrel[i][k][j] = 1;
          }
          else {
            uInfo->URFrel[i][j][k] = 0;
            uInfo->URFrel[i][k][j] = 0;
          }
        }
      }
    }

    for(j=0; j<uInfo->nofProtos[i]; ++j) {
      if(edge_list[j]) {
        free(edge_list[j]);
      }
    }
    free(edge_list);
    free(edge_list_size);
  }

  free(edges);
}

/** Takes the matrices in "URFrel" which contain the URF-pair-relation and finds
the transitive closure - the URF-relation. */
void RDL_findTransitiveClosure(RDL_URFinfo *uInfo)
{
  unsigned i,j,k,l;
  RDL_stack* dfs_stack;
  char* visited;
  unsigned cc_id;
  unsigned **cc;
  unsigned *cc_size;
  unsigned u, v;
  unsigned *elements;
  unsigned next_free_element;
  unsigned *curr_elem;

  dfs_stack = RDL_stack_new();

  for(i = 0; i < uInfo->nofWeights; ++i) {
    visited = malloc(uInfo->nofProtos[i] * sizeof(*visited));
    memset(visited, 0, uInfo->nofProtos[i] * sizeof(*visited));

    cc_id = 0;
    cc = malloc(uInfo->nofProtos[i] * sizeof(*cc));
    cc_size = malloc(uInfo->nofProtos[i] * sizeof(*cc_size));
    elements = malloc(uInfo->nofProtos[i] * sizeof(*elements));
    next_free_element = 0;

    /*
     * perform a DFS search on the URF relations
     * and collect the connected components
     */
    for(j = 0; j < uInfo->nofProtos[i]; ++j) {
      if (!visited[j]) {
        cc[cc_id] = malloc(uInfo->nofProtos[i] * sizeof(**cc));
        cc_size[cc_id] = 0;

        elements[next_free_element] = j;
        RDL_stack_push(dfs_stack, &(elements[next_free_element]));
        ++next_free_element;
        visited[j] = 1;
        while (!RDL_stack_empty(dfs_stack)) {
          curr_elem = RDL_stack_top(dfs_stack);
          k = *curr_elem;
          RDL_stack_pop(dfs_stack);

          /* save the connected components */
          cc[cc_id][cc_size[cc_id]] = k;
          ++cc_size[cc_id];

          for (l = 0; l < uInfo->nofProtos[i]; ++l) {
            if (!visited[l] && uInfo->URFrel[i][k][l]) {
              visited[l] = 1;
              elements[next_free_element] = l;
              RDL_stack_push(dfs_stack, &(elements[next_free_element]));
              ++next_free_element;
            }
          }
        }

        ++cc_id;
      }
    }

    /* every pair inside a CC is URF related */
    for(j = 0; j < cc_id; ++j) {
      for (k = 0; k < cc_size[j]; ++k) {
        for (l = k+1; l < cc_size[j]; ++l) {
          u = cc[j][k];
          v = cc[j][l];
          uInfo->URFrel[i][u][v] = 1;
          uInfo->URFrel[i][v][u] = 1;
        }
      }
    }

    for(j = 0; j < cc_id; ++j) {
      free(cc[j]);
    }
    free(cc);
    free(cc_size);
    free(visited);
    free(elements);
  }

  RDL_stack_delete(dfs_stack);
}

void RDL_findRelations(RDL_cfURF *RCFs, RDL_graph *graph,
    RDL_URFinfo *uInfo, RDL_sPathInfo* spi)
{
  RDL_checkDependencies(RCFs, graph, uInfo);
  RDL_checkEdges(RCFs, graph, uInfo, spi);
  RDL_findTransitiveClosure(uInfo);
}

