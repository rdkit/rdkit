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

#include "RDLdimacs.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "RDLtesting.h"

#define RDL_MAXSIZE 256

void RDL_dimacsGraph_delete(RDL_DimacsGraph* graph)
{
  if (graph->cycles) {
    RDL_deleteCycles(graph->cycles, graph->nof_cycles);
  }
  if (graph->sorted_cycles) {
    RDL_deleteEdgeIdxArray(graph->sorted_cycles, graph->nof_cycles);
  }
  if (graph->edges) {
    free(graph->edges);
  }
  if (graph->name) {
    free(graph->name);
  }
  free(graph);
}

RDL_DimacsGraph* RDL_dimacsGraph_read(const char* filename)
{
  FILE* infile;
  RDL_DimacsGraph* graph = NULL;
  char line[RDL_MAXSIZE + 1];
  unsigned node1, node2, ring;
  int error = 0;
  int stage = 0;
  unsigned edge_counter = 0;
  unsigned ring_counter = 0;
  unsigned ring_position_counter = 0;
  unsigned last_ring = 0;
  unsigned ring_size = 0;
  unsigned ring_idx, edge_idx;
  unsigned eidx;
  int found_eidx = -1;
  int read_initializer = -1;
  unsigned len;
  RDL_edge edge;

  infile = fopen(filename, "r");

  if (!infile) {
    return graph;
  }

  graph = malloc(sizeof(*graph));
  graph->nof_cycles = 0;
  graph->edges = NULL;
  graph->cycles = NULL;
  graph->sorted_cycles = NULL;
  graph->name = NULL;

  while (fgets(line, RDL_MAXSIZE, infile) && !error) {
    if (line[0] == 'p') {
      if (!stage) {
        stage = 1;

        read_initializer = sscanf(line, "p edge %u %u %u", &graph->nof_nodes, &graph->nof_edges, &graph->nof_cycles);
        if (read_initializer == 3) {
          graph->edges = malloc(graph->nof_edges * sizeof(*graph->edges));
          graph->cycles = malloc(graph->nof_cycles * sizeof(*graph->cycles));
          memset(graph->cycles, 0, graph->nof_cycles * sizeof(*graph->cycles));
        }
        else if (read_initializer == 2) {
          graph->edges = malloc(graph->nof_edges * sizeof(*graph->edges));
          graph->nof_cycles = UINT_MAX;
          graph->cycles = NULL;
        }
        else {
          RDL_outputFunc(RDL_ERROR, "invalid line starting with 'p'\n");
          error = 1;
        }
      }
      else {
        RDL_outputFunc(RDL_ERROR, "duplicate 'p' line found!\n");
        error = 1;
      }
    }
    else if (line[0] == 'e') {
      if (stage != 1) {
        RDL_outputFunc(RDL_ERROR, "edge 'e' found without 'p' definition!\n");
        error = 1;
      }
      else {
        if (!sscanf(line, "e %u %u", &node1, &node2)) {
          RDL_outputFunc(RDL_ERROR, "invalid edge specification line\n");
          error = 1;
        }
        else if (edge_counter == graph->nof_edges) {
          RDL_outputFunc(RDL_ERROR, "too many edges found!\n");
          error = 1;
        }
        else if (node1 < 1 || node1 > graph->nof_nodes || node2 < 1 || node2 > graph->nof_nodes) {
          RDL_outputFunc(RDL_ERROR, "illegal edge specification!\n");
          error = 1;
        }
        else {
          graph->edges[edge_counter][0] = node1;
          graph->edges[edge_counter][1] = node2;
          ++edge_counter;
        }
      }
    }
    else if (line[0] == 'r') {
      if (!stage) {
        RDL_outputFunc(RDL_ERROR, "ring 'r' found without 'p' definition!\n");
        error = 1;
      }
      else if (!graph->cycles) {
        RDL_outputFunc(RDL_ERROR, "ring 'r' found without 'p' definition for rings!\n");
        error = 1;
      }
      else {
        stage = 2;
        if (sscanf(line, "r %u %u %u %u", &ring, &ring_size, &node1, &node2) != 4) {
          RDL_outputFunc(RDL_ERROR, "invalid ring specification line\n");
          error = 1;
        }
        else if (ring < 1|| ring > graph->nof_cycles) {
          RDL_outputFunc(RDL_ERROR, "invalid ring counter '%u'!\n", ring);
          error = 1;
        }
        else if (node1 < 1 || node1 > graph->nof_nodes || node2 < 1 || node2 > graph->nof_nodes) {
          RDL_outputFunc(RDL_ERROR, "illegal edge specification!\n");
          error = 1;
        }
        else {
          if (ring != last_ring) {
            if (last_ring) {
              if (ring_position_counter != graph->cycles[last_ring-1]->weight) {
                RDL_outputFunc(RDL_ERROR, "last ring %u wasn't filled up!\n", last_ring);
                error = 1;
              }
              else if (ring != last_ring + 1) {
                RDL_outputFunc(RDL_ERROR, "rings have to be continuous, '%u' follows on '%u'", ring, last_ring);
                error = 1;
              }
            }

            ring_position_counter = 0;
            graph->cycles[ring-1] = malloc(sizeof(**graph->cycles));
            graph->cycles[ring-1]->weight = ring_size;
            graph->cycles[ring-1]->edges = malloc(ring_size * sizeof(*graph->cycles[ring-1]->edges));
            last_ring = ring;
            ++ring_counter;
          }

          if (ring_position_counter == graph->cycles[ring-1]->weight) {
            RDL_outputFunc(RDL_ERROR, "too many edges for ring!\n");
            error = 1;
          }
          else {
            graph->cycles[ring-1]->edges[ring_position_counter][0] = node1;
            graph->cycles[ring-1]->edges[ring_position_counter][1] = node2;
            ++ring_position_counter;
          }
        }
      }
    }
    else if (line[0] == 'c') {
      len = strlen(line);
      if (!graph->name && len > 3) {
        /* skip first two characters (but reserve +1 for \0 and -1 for the newline) */
        graph->name = malloc((len - 2) * sizeof(*graph->name));
        strncpy(graph->name, line + 2, len - 3);
        graph->name[len - 3] = '\0';
      }
    }
    else {
      RDL_outputFunc(RDL_ERROR, "invalid line!\n");
      error = 1;
    }
  }

  if (!error && graph->nof_cycles != UINT_MAX) {
    if (ring_counter != graph->nof_cycles) {
      RDL_outputFunc(RDL_ERROR, "mismatch between number of cycles!\n");
      error = 1;
    }

    if (last_ring) {
      if (ring_position_counter != graph->cycles[last_ring-1]->weight) {
        RDL_outputFunc(RDL_ERROR, "last ring %u wasn't filled up!\n", last_ring);
        error = 1;
      }
    }
  }

  if (!error && edge_counter != graph->nof_edges) {
    RDL_outputFunc(RDL_ERROR, "mismatch between number of edges!\n");
    error = 1;
  }

  if (error && graph) {
    RDL_dimacsGraph_delete(graph);
    graph = NULL;
  }

  fclose(infile);

  if (!error) {
    if (graph->cycles) {
      graph->sorted_cycles = malloc(graph->nof_cycles * sizeof(*graph->sorted_cycles));
      for (ring_idx = 0; ring_idx < graph->nof_cycles; ++ring_idx) {
        graph->sorted_cycles[ring_idx] = malloc((graph->nof_edges + 1) * sizeof(**graph->sorted_cycles));
        memset(graph->sorted_cycles[ring_idx], 0, graph->nof_edges * sizeof(**graph->sorted_cycles));
        graph->sorted_cycles[ring_idx][graph->nof_edges] = 2;
        for (edge_idx = 0; edge_idx < graph->cycles[ring_idx]->weight; ++edge_idx) {
          edge[0] = graph->cycles[ring_idx]->edges[edge_idx][0];
          edge[1] = graph->cycles[ring_idx]->edges[edge_idx][1];
          for (eidx = 0; eidx < graph->nof_edges; ++eidx) {
            if ((graph->edges[eidx][0] == edge[0]&& graph->edges[eidx][1] == edge[1])
                || (graph->edges[eidx][0] == edge[1] && graph->edges[eidx][1] == edge[0])) {
              found_eidx = eidx;
            }
          }
          graph->sorted_cycles[ring_idx][found_eidx] = 1;
        }
      }

      qsort(graph->sorted_cycles, graph->nof_cycles, sizeof(*graph->sorted_cycles), RDL_cmp_cycles);
    }
    else {
      graph->sorted_cycles = NULL;
    }
  }

  return graph;
}
