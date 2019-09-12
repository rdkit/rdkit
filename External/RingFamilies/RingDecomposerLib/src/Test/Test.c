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
 * @file
 *
 * @brief Demo and validation tool for the RingDecomposerLib library
 *
 * This tool can be used for demo output and for testing/validation.
 *
 * Demo output
 * -----------
 * Start the tool with
 *
 *     Test demo <filename>
 *
 * The program will output the ring topology calculated by
 * the library. The input file must be in DIMACS format
 * as described below.
 *
 *
 * Validation
 * ----------
 * Start the tool with
 *
 *     Test validate <filename> [<timeout>]
 *
 * The program will compare the relevant cycles to the cycles
 * present in the input file.
 * Furthermore the URFs will be validated, i.e.
 * calculated with an independent exponential algorithm
 * (definition of the URFs). The same procedure is applied
 * to the SSSR (verify it's a cycle base). A number of consistency
 * tests is executed to ensure integrity and robustness.
 * See the file `test/README` for a descriptions of the test files.
 *
 * The input file must be in DIMACS format
 * as described below.
 * A number of interesting example graphs are in the folder
 * `test`.
 *
 * The parameter `timeout` is optional and specifies a timeout
 * in seconds for the exponential URF validation and the SSSR validation algorithm.
 * (Timeouts are failures!)
 *
 * File format
 * ------------
 *
 * The input file format is a (modified) DIMACS graph format.
 * See folder `test` for example files.
 * Graph format specification:
 *
 * First line:
 *
 *     p <number of nodes> <number of edges> [<nof cycles>]
 *
 * with `<nof cycles>` optional.
 *
 * Then the edges of the graph follow
 *
 *     e <node1> <node2>
 *
 * If `<nof cycles>` was specified, then the relevant cycles
 * of the graph are listed:
 *
 *     r <ring id> <ring size> <node1> <node2>
 *
 * `<ring id>` and `<ring size>` are repeated for each edge
 * of the graph as specified by `<node1>` and `<node2>`.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "RDLtesting.h"
#include "TestDemo.h"
#include "TestValidate.h"

/** @private */
int main(int argc, char **argv)
{
  const char* filename;
  RDL_graph *graph;
  RDL_data *data;
  const char* mode;
  unsigned idx;
  int timeout = -1, urf_validation_result, basis_validation_result, consistency_result;
  int return_code;

  if (argc < 2) {
    fprintf(stderr, "usage: %s mode filename [timeout=-1]\nwith mode 'demo' "
        "for demo output of URF calculations\n"
        "or\nwith mode 'validate' for testing/validation\n", argv[0]);
    return EXIT_FAILURE;
  }

  mode = argv[1];

  if (!strcmp(mode, "demo")) {
    if (argc < 3) {
      fprintf(stderr, "usage: %s demo filename\n", argv[0]);
      return EXIT_FAILURE;
    }
  }
  else if (!strcmp(mode, "validate")) {
    if (argc < 3) {
      fprintf(stderr, "usage: %s demo filename [timeout=-1]\n"
          "with timeout is validation timeout in seconds, default -1 (no timeout)\n", argv[0]);
      return EXIT_FAILURE;
    }
    if (argc > 3) {
      if (!sscanf(argv[3], "%d", &timeout)) {
        fprintf(stderr, "invalid timeout in seconds: '%s'\n", argv[3]);
        return EXIT_FAILURE;
      }
    }
  }
  else {
    fprintf(stderr, "invalid mode '%s', allowed are 'demo' and 'validate'\n", mode);
    return EXIT_FAILURE;
  }

  filename = argv[2];

  RDL_DimacsGraph* dimacs_graph = RDL_dimacsGraph_read(filename);

  if (!dimacs_graph) {
    fprintf(stderr, "unable to read file '%s'\n", filename);
    return EXIT_FAILURE;
  }

  graph = RDL_initNewGraph(dimacs_graph->nof_nodes);

  for (idx = 0; idx < dimacs_graph->nof_edges; ++idx) {
    RDL_addUEdge(graph, dimacs_graph->edges[idx][0]-1, dimacs_graph->edges[idx][1]-1);
  }

  RDL_setOutputFunction(RDL_writeToStderr);
  /* calculate Unique Ring Families */
  data = RDL_calculate(graph);
  if (!data) {
    fprintf(stderr, "Calculation failed!\n");
    RDL_deleteGraph(graph);
    RDL_dimacsGraph_delete(dimacs_graph);
    return EXIT_FAILURE;
  }

  if (!strcmp(mode, "demo")) {
    RDL_demo_output(data);
    return_code = EXIT_SUCCESS;
  }
  else if (!strcmp(mode, "validate")) {
    RDL_validate(data, dimacs_graph, timeout,
        &urf_validation_result, &basis_validation_result);
    if (urf_validation_result == -1) {
      fprintf(stdout, "URF validation timed out!\n");
    }
    else if (urf_validation_result == -2) {
      fprintf(stdout, "URF validation ran out of memory!\n");
    }
    else if (urf_validation_result == 0) {
      fprintf(stdout, "URF validation successful!\n");
    }
    else {
      fprintf(stderr, "URF validation failed!\n");
    }

    if (basis_validation_result == -1) {
      fprintf(stdout, "MCB validation timed out!\n");
    }
    else if (basis_validation_result == 0) {
      fprintf(stdout, "MCB validation successful!\n");
    }
    else {
      fprintf(stderr, "MCB validation failed!\n");
    }

    consistency_result = RDL_checkConsistency(data);
    if (consistency_result == 0) {
      fprintf(stdout, "consistency validation successful!\n");
    }
    else {
      fprintf(stderr, "consistency validation failed!\n");
    }

    if (urf_validation_result || consistency_result
        || basis_validation_result) {
      return_code = EXIT_FAILURE;
    }
    else {
      return_code = EXIT_SUCCESS;
    }
  }
  else {
    return_code = EXIT_FAILURE;
  }

  /* delete URFdata and the graph */
  RDL_deleteData(data);
  RDL_dimacsGraph_delete(dimacs_graph);

  return return_code;
}
