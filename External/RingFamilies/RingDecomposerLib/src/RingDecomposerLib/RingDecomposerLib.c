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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "RingDecomposerLib.h"

#include "RDLapsp.h"
#include "RDLbitset.h"
#include "RDLcycleFams.h"
#include "RDLdataStruct.h"
#include "RDLgraph.h"
#include "RDLhandler.h"
#include "RDLinfo.h"
#include "RDLrelation.h"
#include "RDLtarjan.h"
#include "RDLutility.h"

const unsigned RDL_RESERVED_START = 64;

/* define some constants */
const unsigned RDL_INVALID_RESULT = UINT_MAX;
const unsigned RDL_DUPLICATE_EDGE = UINT_MAX - 1;
const unsigned RDL_NO_RINGSYSTEM = UINT_MAX - 2;
const double RDL_INVALID_RC_COUNT = DBL_MAX;

/* public, set the global output function */
void RDL_setOutputFunction(RDL_outputFunction func)
{
  RDL_outputFunc = func;
}

/* public */
RDL_graph *RDL_initNewGraph(unsigned V)
{
  return RDL_initNewGraph_g(V, 1);
}

/* public, add an edge with range check! */
unsigned RDL_addUEdge(RDL_graph *gra, RDL_node from, RDL_node to)
{
  return RDL_addUEdge_g(gra, from, to);
}

/* public, the calculation routine */
RDL_data *RDL_calculate(RDL_graph *gra)
{
  RDL_data *data;
  unsigned i, j, k, urf_index, rcf_index;
  unsigned nof_relevant_fams, nof_relevant_fams_sum;

  if (!gra) {
    RDL_outputFunc(RDL_ERROR, "The graph is NULL.\n");
    return NULL;
  }

  /* we can't calculate anything if there isn't at least one node */
  if(!gra->V) {
    RDL_outputFunc(RDL_ERROR, "The graph has no nodes.\n");
    return NULL;
  }

  data = malloc(sizeof(*data));

  /* FIRST STEP: TARJAN */
  data->bccGraphs = RDL_tarjanBCC(gra);
  data->nofURFs = 0;
  data->nofRCFs = 0;

  /* allocate result structures */
  data->spiPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->spiPerBCC));
  data->CFsPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->CFsPerBCC));
  data->urfInfoPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->urfInfoPerBCC));
  data->nofURFsPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->nofURFsPerBCC));
  data->nofRCFsPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->nofRCFsPerBCC));

  /* the followeing steps are performed for each BCC */
  for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
    /* solve APSP problem */
    data->spiPerBCC[i] = RDL_AllPairsShortestPaths(data->bccGraphs->bcc_graphs[i]);
    /* calculate RCFs with Vismara's algortihm */
    data->CFsPerBCC[i] = RDL_findCycleFams(data->bccGraphs->bcc_graphs[i], data->spiPerBCC[i]);
    /* this can fail if the BCC is to large */
    if (!data->CFsPerBCC[i]) {
      /* delete families and APSP information UPTO current index */
      for (j = 0; j < i; ++j) {
        RDL_deleteAPSP(data->spiPerBCC[j], data->bccGraphs->bcc_graphs[j]->V);
        RDL_deleteCycleFams(data->CFsPerBCC[j]);
        if (data->nofURFsPerBCC[j] > 0) {
          RDL_deleteURFInfo(data->urfInfoPerBCC[j]);
        }
      }
      /* delete APSP for current index (there is no RCF info) */
      RDL_deleteAPSP(data->spiPerBCC[i], data->bccGraphs->bcc_graphs[i]->V);

      /* free alloced per BCC structures */
      free(data->spiPerBCC);
      free(data->CFsPerBCC);
      free(data->nofURFsPerBCC);
      free(data->nofRCFsPerBCC);
      free(data->urfInfoPerBCC);
      /* and of course the BCC graph */
      RDL_deleteBCCGraph(data->bccGraphs);
      /* free resulting struct */
      free(data);

      return NULL;
    }
    if(data->CFsPerBCC[i]->nofFams > 0) {
      /* if there is at least one RCF, check URF relation */
      data->urfInfoPerBCC[i] = RDL_checkURFRelation(data->CFsPerBCC[i],
          data->bccGraphs->bcc_graphs[i], data->spiPerBCC[i]);
      data->nofURFsPerBCC[i] = data->urfInfoPerBCC[i]->nofURFs;

      nof_relevant_fams = 0;
      /* count RCFs */
      for (j = 0; j < data->CFsPerBCC[i]->nofFams; ++j) {
        if (data->CFsPerBCC[i]->fams[j]->mark) {
          ++nof_relevant_fams;
        }
      }

      nof_relevant_fams_sum = 0;
      for (j = 0; j < data->nofURFsPerBCC[i]; ++j) {
        nof_relevant_fams_sum += data->urfInfoPerBCC[i]->nofCFsPerURF[j];
      }

      if (nof_relevant_fams != nof_relevant_fams_sum) {
        RDL_outputFunc(RDL_ERROR, "different number of relevant families!\n");
        /* internal check, should never happen */
        assert(0);
      }
      data->nofRCFsPerBCC[i] = nof_relevant_fams;
    }
    else {
      data->nofURFsPerBCC[i] = 0;
      data->nofRCFsPerBCC[i] = 0;
    }
    data->nofURFs += data->nofURFsPerBCC[i];
    data->nofRCFs += data->nofRCFsPerBCC[i];
  }

  /* create a mapping from URFs to BCCs */
  data->urf_to_bcc = malloc(data->nofURFs * sizeof(*data->urf_to_bcc));
  /* create a mapping from RCFs to URFs */
  data->rcf_to_urf = malloc(data->nofRCFs * sizeof(*data->rcf_to_urf));
  urf_index = 0;
  rcf_index = 0;
  for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
    for (j = 0; j < data->nofURFsPerBCC[i]; ++j, ++urf_index) {
      data->urf_to_bcc[urf_index][0] = i;
      data->urf_to_bcc[urf_index][1] = j;
      for (k = 0; k < data->urfInfoPerBCC[i]->nofCFsPerURF[j]; ++k, ++rcf_index) {
        data->rcf_to_urf[rcf_index][0] = urf_index;
        data->rcf_to_urf[rcf_index][1] = k;
      }
    }
  }

  data->graph = gra;

  return data;
}

/* public, free memory */
void RDL_deleteData(RDL_data *data)
{
  unsigned i;

  for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
    RDL_deleteAPSP(data->spiPerBCC[i], data->bccGraphs->bcc_graphs[i]->V);
    RDL_deleteCycleFams(data->CFsPerBCC[i]);
    if (data->nofURFsPerBCC[i] > 0) {
      RDL_deleteURFInfo(data->urfInfoPerBCC[i]);
    }
  }
  free(data->spiPerBCC);
  free(data->CFsPerBCC);
  free(data->nofURFsPerBCC);
  free(data->nofRCFsPerBCC);
  free(data->urfInfoPerBCC);
  RDL_deleteBCCGraph(data->bccGraphs);
  free(data->urf_to_bcc);
  free(data->rcf_to_urf);

  RDL_deleteGraph(data->graph);
  free(data);
}

/* public, get the total number of URFs */
unsigned RDL_getNofURF(const RDL_data *data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }
  return data->nofURFs;
}

/* public, get the total number of RCFs */
unsigned RDL_getNofRCF(const RDL_data *data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }
  return data->nofRCFs;
}

/*
 * private
 * return the weight for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static unsigned RDL_getWeight_internal(const RDL_data *data,
    unsigned bcc_index, unsigned urf_internal_index, unsigned rcf_internal_index)
{
  return data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index][rcf_internal_index]->weight;
}

/* public, get weight for the given URF */
unsigned RDL_getWeightForURF(const RDL_data *data, unsigned index)
{
  unsigned bcc_index, internal_index;
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RESULT;
  }
  /* ger internal index */
  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];
  return RDL_getWeight_internal(data, bcc_index, internal_index, 0);
}

/* public, get weight for the given RCF */
unsigned RDL_getWeightForRCF(const RDL_data *data, unsigned index)
{
  unsigned bcc_index, urf_internal_index, rcf_internal_index, urf_index;
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RESULT;
  }

  /* map rcf to urf and urf to bcc */
  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];
  return RDL_getWeight_internal(data, bcc_index,
      urf_internal_index, rcf_internal_index);
}

/*
 * private
 * return the nodes for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static void RDL_getNodes_internal(const RDL_data *data, unsigned bcc_index,
    unsigned urf_internal_index, unsigned rcf_internal_index, char* atoms)
{
  const RDL_cfam **URF;
  char *visited;
  const RDL_graph* graph;

  graph = data->bccGraphs->bcc_graphs[bcc_index];

  URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];

  visited = malloc(graph->V *  sizeof(*visited));

  /* BFS on both shortest path graphs */
  memset(visited, 0, graph->V *  sizeof(*visited));
  RDL_giveVertices(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
      atoms, data->spiPerBCC[bcc_index], visited);

  memset(visited, 0, graph->V *  sizeof(*visited));
  RDL_giveVertices(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
      atoms, data->spiPerBCC[bcc_index], visited);

  if(URF[rcf_internal_index]->x < UINT_MAX) /*even cycle*/ {
    atoms[URF[rcf_internal_index]->x] = 1;
  }
  free(visited);
}


/* private
 * gives an array of indices of nodes that are contained in the URF with the
 * given index. Array is terminated by UINT_MAX
 */
static RDL_node *RDL_getNodesURF(const RDL_data *data, unsigned index)
{
  unsigned i,nofFams,nextfree=0,alloced, bcc_index, internal_index;
  char *atoms;
  RDL_node *result;
  const RDL_graph* graph;

  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  atoms = malloc(graph->V * sizeof(*atoms));
  memset(atoms, 0, graph->V * sizeof(*atoms));

  nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];
  alloced = RDL_RESERVED_START;
  result = malloc(alloced * sizeof(*result));

  /* iterate over all RCFs in this URFs and collect data in one bitset (atoms) */
  for(i=0; i<nofFams; ++i) {
    RDL_getNodes_internal(data, bcc_index, internal_index, i, atoms);
  }

  /* translate into a dynamic table */
  for(i=0; i<graph->V; ++i) {
    if(atoms[i] == 1) {
      if(nextfree == alloced) {/*double the size*/
        alloced *= 2;
        result = realloc(result, alloced*sizeof(*result));
      }
      result[nextfree++] = data->bccGraphs->node_from_bcc_mapping[bcc_index][i];
    }
  }

  result = realloc(result, (nextfree+1)*sizeof(*result));

  result[nextfree] = UINT_MAX;
  free(atoms);

  return result;
}

/* private
 * gives an array of indices of nodes that are contained in the RRF with the
 * given index. Array is terminated by UINT_MAX
 */
static RDL_node *RDL_getNodesRCF(const RDL_data *data, unsigned index)
{
  unsigned i,nextfree=0,alloced, bcc_index,
      urf_internal_index, rcf_internal_index, urf_index;
  char *atoms;
  RDL_node *result;
  const RDL_graph* graph;

  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  atoms = malloc(graph->V * sizeof(*atoms));
  memset(atoms, 0, graph->V * sizeof(*atoms));

  alloced = RDL_RESERVED_START;
  result = malloc(alloced * sizeof(*result));

  /* save the nodes into bitset */
  RDL_getNodes_internal(data, bcc_index, urf_internal_index,
      rcf_internal_index, atoms);

  /* translate into a dynamic table */
  for(i=0; i<graph->V; ++i) {
    if(atoms[i] == 1) {
      if(nextfree == alloced) {/*double the size*/
        alloced *= 2;
        result = realloc(result, alloced*sizeof(*result));
      }
      result[nextfree++] = data->bccGraphs->node_from_bcc_mapping[bcc_index][i];
    }
  }

  result = realloc(result, (nextfree+1)*sizeof(*result));

  result[nextfree] = UINT_MAX;
  free(atoms);

  return result;
}

/* public, calculate the number of RCs */
double RDL_getNofRC(const RDL_data *data)
{
  double result = 0, intermediate;
  unsigned i;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RC_COUNT;
  }

  /* add up results for URFs */
  for (i = 0; i < data->nofURFs; ++i) {
    intermediate = RDL_getNofRCForURF(data, i);
    if (intermediate == RDL_INVALID_RC_COUNT) {
      return RDL_INVALID_RC_COUNT;
    }

    result += intermediate;
  }

  return result;
}

/*
 * private
 * return the # RCs for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static double RDL_getNofRCForRCF_internal(const RDL_data *data,
    unsigned bcc_index, unsigned urf_internal_index,
    unsigned rcf_internal_index)
{
  double nofPaths1, nofPaths2, result;
  const double prod_limit = sqrt(DBL_MAX) - 1.0;
  const RDL_cfam** URF;
  const RDL_graph* graph;

  result = 0.0;

  graph = data->bccGraphs->bcc_graphs[bcc_index];

  URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];

  /* count the number of paths one both sides */
  nofPaths1 = RDL_countPaths(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
      graph->V, data->spiPerBCC[bcc_index]);
  nofPaths2 = RDL_countPaths(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
      graph->V, data->spiPerBCC[bcc_index]);

  result = nofPaths1 * nofPaths2;

  /* check if either of the number paths is larger than what we can multiply */
  if (nofPaths1 >= prod_limit || nofPaths2 >= prod_limit) {
    RDL_outputFunc(RDL_WARNING, "result overflow when counting paths!\n");
    return RDL_INVALID_RC_COUNT;
  }

  return result;
}

/* public, return the # RCs for given RCF */
double RDL_getNofRCForRCF(const RDL_data *data, unsigned index)
{
  unsigned bcc_index, urf_id, urf_internal_index, rcf_internal_index;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RC_COUNT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RC_COUNT;
  }

  /* map indices */
  urf_id = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];

  bcc_index = data->urf_to_bcc[urf_id][0];
  urf_internal_index = data->urf_to_bcc[urf_id][1];

  return RDL_getNofRCForRCF_internal(data, bcc_index,
      urf_internal_index, rcf_internal_index);
}

/* public, return the # RCs for given URF */
double RDL_getNofRCForURF(const RDL_data *data, unsigned index)
{
  unsigned i,nofFams, bcc_index, internal_index;
  double result=0, prod;

  const double sum_limit = DBL_MAX/2.0 - 1.0;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RC_COUNT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RC_COUNT;
  }

  /* map indices */
  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];

  nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];

  /* iterate all RCFs in this URF */
  for(i=0; i<nofFams; ++i) {
    prod = RDL_getNofRCForRCF_internal(data, bcc_index, internal_index, i);

    /* check if any summand is larger than limit */
    if (prod >= sum_limit || result >= sum_limit) {
      RDL_outputFunc(RDL_WARNING, "result overflow when counting paths!\n");
      return RDL_INVALID_RC_COUNT;
    }

    result += prod;
  }

  return result;
}

/* private
 * return the edges for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static void RDL_getEdges_internal(const RDL_data *data,
    unsigned bcc_index, unsigned urf_internal_index,
    unsigned rcf_internal_index, char* edges)
{
  const RDL_cfam **URF;
  char *visited;
  const RDL_graph* graph;

  graph = data->bccGraphs->bcc_graphs[bcc_index];

  URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];

  /* perform BFS on both paths and collect edges */
  visited = malloc(graph->V *  sizeof(*visited));
  memset(visited, 0, graph->V *  sizeof(*visited));

  RDL_giveEdges(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
      edges, graph, data->spiPerBCC[bcc_index], visited);

  memset(visited, 0, graph->V *  sizeof(*visited));
  RDL_giveEdges(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
      edges, graph, data->spiPerBCC[bcc_index], visited);

  /*
   * in contrast to the nodes, we have to add additional edges, because we're
   * looking at cycles...
   */
  if(URF[rcf_internal_index]->x < UINT_MAX) /*even cycle*/ {
    edges[RDL_edgeId(graph,URF[rcf_internal_index]->q,URF[rcf_internal_index]->x)] = 1;
    edges[RDL_edgeId(graph,URF[rcf_internal_index]->p,URF[rcf_internal_index]->x)] = 1;
  }
  else /*odd cycle*/ {
    edges[RDL_edgeId(graph,URF[rcf_internal_index]->q,URF[rcf_internal_index]->p)] = 1;
  }
  free(visited);
}


/* private
 * gives an array of indices of edges that are contained in the URF with the
 * given index. Array is terminated by UINT_MAX
 */
static unsigned *RDL_getEdgesURF(const RDL_data *data, unsigned index)
{
  unsigned i,nofFams,nextfree=0,alloced, bcc_index, internal_index;
  char *edges;
  unsigned *result;
  const RDL_graph* graph;

  /* map indices */
  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  edges = malloc(graph->E * sizeof(*edges));
  memset(edges, 0, graph->E * sizeof(*edges));

  nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];
  alloced = RDL_RESERVED_START;
  result = malloc(alloced * sizeof(*result));

  /* iterate over all RCFs in this URF and collect results in bitset (edges) */
  for(i=0; i<nofFams; ++i) {
    RDL_getEdges_internal(data, bcc_index, internal_index, i, edges);
  }

  /* translate into dynamic table */
  for(i=0; i<graph->E; ++i) {
    if(edges[i] == 1) {
      if(nextfree == alloced)
      {/*double the size*/
        alloced *= 2;
        result = realloc(result, alloced*sizeof(*result));
      }
      result[nextfree++] = data->bccGraphs->edge_from_bcc_mapping[bcc_index][i];
    }
  }
  result = realloc(result, (nextfree+1)*sizeof(*result));

  result[nextfree] = UINT_MAX;

  free(edges);
  return result;
}

/* private
 * gives an array of indices of edges that are contained in the RCF with the
 * given index. Array is terminated by UINT_MAX
 */
static unsigned *RDL_getEdgesRCF(const RDL_data *data, unsigned index)
{
  unsigned j,nextfree=0,alloced, urf_index, bcc_index,
      urf_internal_index, rcf_internal_index;
  char *edges;
  unsigned *result;
  const RDL_graph* graph;

  /* index mapping */
  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  edges = malloc(graph->E * sizeof(*edges));
  memset(edges, 0, graph->E * sizeof(*edges));

  alloced = RDL_RESERVED_START;
  result = malloc(alloced * sizeof(*result));

  RDL_getEdges_internal(data, bcc_index, urf_internal_index,
      rcf_internal_index, edges);

  /* translate into dynamic table */
  for(j=0; j<graph->E; ++j) {
    if(edges[j] == 1) {
      if(nextfree == alloced)
      {/*double the size*/
        alloced *= 2;
        result = realloc(result, alloced*sizeof(*result));
      }
      result[nextfree++] = data->bccGraphs->edge_from_bcc_mapping[bcc_index][j];
    }
  }

  result = realloc(result, (nextfree+1)*sizeof(*result));

  result[nextfree] = UINT_MAX;

  free(edges);
  return result;
}


/* private, calls getNodesURF() or getEdgesURF() depending on mode 'a' or 'b' */
static unsigned *RDL_giveURF(const RDL_data *data, unsigned index, char mode)
{
  unsigned *result;
  if(mode == 'a')
  {
    result = RDL_getNodesURF(data, index);
  }
  else if(mode == 'b')
  {
    result = RDL_getEdgesURF(data, index);
  }
  else
  {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_giveURF()' with invalid mode '%c'\n", mode);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  return result;
}

/* private, calls getNodesRCF() or getEdgesRCF() depending on mode 'a' or 'b' */
static unsigned *RDL_giveRCF(const RDL_data *data, unsigned index, char mode)
{
  unsigned *result;
  if(mode == 'a')
  {
    result = RDL_getNodesRCF(data, index);
  }
  else if(mode == 'b')
  {
    result = RDL_getEdgesRCF(data, index);
  }
  else
  {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_giveRCF()' with invalid mode '%c'\n", mode);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  return result;
}

/* public, function for retrieving nodes in an URF */
unsigned RDL_getNodesForURF(const RDL_data *data, unsigned index, RDL_node **ptr)
{
  unsigned i;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  (*ptr) = RDL_getNodesURF(data, index);
  for(i=0; (*ptr)[i]<UINT_MAX; ++i); /*counts the number of atoms*/
  return i;
}

/* public, function for retrieving nodes in an RCF */
unsigned RDL_getNodesForRCF(const RDL_data *data, unsigned index, RDL_node **ptr)
{
  unsigned i;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  (*ptr) = RDL_getNodesRCF(data, index);
  for(i=0; (*ptr)[i]<UINT_MAX; ++i); /*counts the number of atoms*/
  return i;
}

/* public, function for retrieving edges of URF */
unsigned RDL_getEdgesForURF(const RDL_data *data, unsigned index, RDL_edge **ptr)
{
  unsigned nextfree, alloced;
  RDL_edge *result;
  unsigned *edgeIndices;
  unsigned i;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  nextfree = 0;
  alloced = RDL_RESERVED_START;
  result = malloc(alloced * sizeof(*result));
  edgeIndices = RDL_getEdgesURF(data, index);
  /* translate result from edge indices to edges */
  for(i=0; edgeIndices[i]<UINT_MAX; ++i)
  {
    if(nextfree == alloced)/*more space needed in result*/
    {
      alloced *= 2; /* double the space */
      result = realloc(result, alloced * sizeof(*result));
    }
    result[nextfree][0] = data->graph->edges[edgeIndices[i]][0];
    result[nextfree][1] = data->graph->edges[edgeIndices[i]][1];
    ++nextfree;
  }
  result = realloc(result, nextfree * sizeof(*result));
  free(edgeIndices);
  (*ptr) = result;
  return nextfree;
}

/* public, function for retrieving edges of RCF */
unsigned RDL_getEdgesForRCF(const RDL_data *data, unsigned index, RDL_edge **ptr)
{
  unsigned nextfree, alloced;
  RDL_edge *result;
  unsigned *edgeIndices;
  unsigned i;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  nextfree = 0;
  alloced = RDL_RESERVED_START;
  result = malloc(alloced * sizeof(*result));
  edgeIndices = RDL_getEdgesRCF(data, index);
  /* translate result from edge indices to edges */
  for(i=0; edgeIndices[i]<UINT_MAX; ++i)
  {
    if(nextfree == alloced)/*more space needed in result*/
    {
      alloced *= 2; /* double the space */
      result = realloc(result, alloced * sizeof(*result));
    }
    result[nextfree][0] = data->graph->edges[edgeIndices[i]][0];
    result[nextfree][1] = data->graph->edges[edgeIndices[i]][1];
    ++nextfree;
  }
  result = realloc(result, nextfree * sizeof(*result));
  free(edgeIndices);
  (*ptr) = result;
  return nextfree;
}

/* public, to delete result of RDL_getCyclesChar() */
void RDL_deleteCyclesChar(char **cycles)
{
  unsigned i;
  for(i=0; cycles[i]!=NULL; ++i)
  {
    free(cycles[i]);
  }
  free(cycles);
}

/* public, delete cycle */
void RDL_deleteCycle(RDL_cycle *cycle)
{
  free(cycle->edges);
  free(cycle);
}

/* public, delete cycles */
void RDL_deleteCycles(RDL_cycle **cycles, unsigned number)
{
  unsigned i;
  for(i=0; i<number; ++i)
  {
    RDL_deleteCycle(cycles[i]);
  }
  free(cycles);
}

/*
 * private
 * Gives all families containing the object, which can be an atom or a bond.
mode:
  - 'a': nodes
  - 'b': edges
returns an array of integers containing all indices of families containing the
object. The array ends with a terminating UINT_MAX on its last position and has
to be deallocated with 'free()'
*/
static unsigned *RDL_listFamilies(const RDL_data *data, unsigned object,
    char mode, char family)
{
  unsigned *result;
  unsigned nof_families;
  unsigned nextfree=0, alloced=RDL_RESERVED_START;
  unsigned *objects;
  unsigned i,j, bcc_index, target_bcc, urf_index;
  char contained, can_be_contained;

  if(!(mode == 'a' || mode == 'b'))
  {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_listFamilies()' with invalid mode '%c'\n",mode);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  unsigned* (*giveFamily) (const RDL_data*, unsigned, char) = NULL;

  if (family == 'u') {
    nof_families = data->nofURFs;
    giveFamily = RDL_giveURF;
  }
  else if (family == 'r') {
    nof_families = data->nofRCFs;
    giveFamily = RDL_giveRCF;
  }
  else {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_listFamilies()' with invalid family '%c'\n",family);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  result=malloc(alloced * sizeof(*result));
  for(i = 0; i < nof_families; ++i)
  {
    if (family == 'u') {
      urf_index = i;
    }
    else if (family == 'r') {
      urf_index = data->rcf_to_urf[i][0];
    }

    bcc_index = data->urf_to_bcc[urf_index][0];

    /*
     * check if object belongs to this BCC at all,
     * if not it can't be in the URF either
     */
    can_be_contained = 0;
    if (mode == 'a') {
      for (j = 0; j < data->bccGraphs->nof_bcc_per_node[object]; ++j) {
        target_bcc = data->bccGraphs->node_to_bcc_mapping[object][2*j];
        if (target_bcc == bcc_index) {
          can_be_contained = 1;
          break;
        }
      }
    }
    else {
      target_bcc = data->bccGraphs->edge_to_bcc_mapping[object][0];
      if (target_bcc == bcc_index) {
        can_be_contained = 1;
      }
    }

    contained = 0;
    if (can_be_contained) {
      objects = giveFamily(data, i, mode);
      for(j=0; objects[j]<UINT_MAX; ++j)
      {
        if(objects[j] == object)
        {
          contained = 1;
          break;
        }
      }
      free(objects);
    }
    if (contained) {
      if(nextfree == alloced) {
        alloced *= 2;
        result = realloc(result, alloced*sizeof(*result));
      }
      result[nextfree++] = i;
    }
  }

  result = realloc(result, (nextfree+1)*sizeof(*result));
  result[nextfree] = UINT_MAX;
  return result;
}

/* private, get families containing a node */
static unsigned RDL_getFamiliesContainingNode_internal(
    const RDL_data *data, RDL_node object, unsigned **ptr,
    char family)
{
  unsigned i;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (object >= data->graph->V) {
    RDL_outputFunc(RDL_ERROR, "invalid node: %u\n", object);
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = malloc(sizeof(**ptr));
    return 0;
  }

  *ptr = RDL_listFamilies(data, object, 'a', family);
  for(i=0; (*ptr)[i]<UINT_MAX; ++i);
  return i;
}

/* public, get the RCF indices with this node */
unsigned RDL_getRCFsContainingNode(const RDL_data *data, RDL_node object,
                                   unsigned **ptr)
{
  return RDL_getFamiliesContainingNode_internal(data, object, ptr, 'r');
}

/* public, get the URF indices with this node */
unsigned RDL_getURFsContainingNode(const RDL_data *data, RDL_node object,
                                   unsigned **ptr)
{
  return RDL_getFamiliesContainingNode_internal(data, object, ptr, 'u');
}

/* private, get the family indices with this edge */
static unsigned RDL_getFamiliesContainingEdge_internal(
    const RDL_data *data, RDL_node node1, RDL_node node2, unsigned **ptr,
    char family)
{
  unsigned i, edge;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (node1 >= data->graph->V || node2 >= data->graph->V) {
    RDL_outputFunc(RDL_ERROR, "invalid edge: %u %u\n", node1, node2);
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  edge = RDL_edgeId(data->graph, node1, node2);

  /* check if this edge exists at all */
  if (edge == RDL_INVALID_RESULT) {
    RDL_outputFunc(RDL_ERROR, "invalid edge: %u %u\n", node1, node2);
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = malloc(sizeof(**ptr));
    return 0;
  }

  *ptr = RDL_listFamilies(data, edge, 'b', family);
  for(i=0; (*ptr)[i]<UINT_MAX; ++i);
  return i;
}


/* public, get the RCF indices with this edge */
unsigned RDL_getRCFsContainingEdge(const RDL_data *data, RDL_node node1,
                                   RDL_node node2, unsigned **ptr)
{
  return RDL_getFamiliesContainingEdge_internal(data, node1, node2, ptr, 'r');
}

/* public, get the URF indices with this edge */
unsigned RDL_getURFsContainingEdge(const RDL_data *data, RDL_node node1,
                                   RDL_node node2, unsigned **ptr)
{
  return RDL_getFamiliesContainingEdge_internal(data, node1, node2, ptr, 'u');
}

/* public, get the number of URFs that this node is part of */
unsigned RDL_getNofURFContainingNode(const RDL_data *data, RDL_node node)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getURFsContainingNode(data, node, &arr);
  free(arr);
  return number;
}

/* public, get the number of URFs that this edge is part of */
unsigned RDL_getNofURFContainingEdge(const RDL_data *data, RDL_node node1, RDL_node node2)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getURFsContainingEdge(data, node1, node2, &arr);
  free(arr);
  return number;
}

/* public, get the number of RCFs that this node is part of */
unsigned RDL_getNofRCFContainingNode(const RDL_data *data, RDL_node node)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getRCFsContainingNode(data, node, &arr);
  free(arr);
  return number;
}

/* public, get the number of RCFs that this edge is part of */
unsigned RDL_getNofRCFContainingEdge(const RDL_data *data, RDL_node node1, RDL_node node2)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getRCFsContainingEdge(data, node1, node2, &arr);
  free(arr);
  return number;
}

/* construct a SSSR */
unsigned RDL_getSSSR(const RDL_data *data, RDL_cycle ***ptr)
{
  unsigned j, k, rcf=0, urf=0;
  RDL_cycle **result;
  unsigned size, added, currBond;
  unsigned bcc_index, internal_index;
  char* prototype;
  const RDL_graph* curr_graph;
  unsigned total_counter, alloced;
  unsigned edge;
  unsigned char *currentCycle;
  unsigned char **basisCycles;
  unsigned char **prototypesForGaussian;
  unsigned currentBasisCyclesSize;
  unsigned oldBasisCyclesSize;
  unsigned compressedSize = 0;
  const unsigned char* empty_cycle;


  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = malloc(sizeof(**ptr));
    return 0;
  }

  total_counter = 0;
  if (data->graph->E < data->graph->V) {
    /*
     * if E - V + 1 < 1  <=> E - V < 0 <=> E < V
     *
     * graph is not connected and our heuristic for reservation
     * fails => allocate something ;)
     * */
    alloced = RDL_RESERVED_START;
  }
  else {
    alloced = data->graph->E - data->graph->V + 1;
  }
  result = malloc(alloced * sizeof(*result));

  /* calculate basis for each BCC (and therefore for each CC) */

  for (bcc_index = 0; bcc_index < data->bccGraphs->nof_bcc; ++bcc_index) {
    curr_graph = data->bccGraphs->bcc_graphs[bcc_index];

    size=curr_graph->E - curr_graph->V + 1;

    basisCycles = malloc(size * sizeof(*basisCycles));
    prototypesForGaussian =
        malloc(data->nofURFsPerBCC[bcc_index] * sizeof(*prototypesForGaussian));

    /* copy the prototypes for gaussian elimination */
    /*
     * it is sufficient to take the first prototype for each URF
     * (the other's are linearly dependent anyway
     */
    for (internal_index = 0;
        internal_index < data->nofURFsPerBCC[bcc_index]; ++internal_index) {
      compressedSize = RDL_bitset_compressed(&(prototypesForGaussian[internal_index]),
          data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->prototype,
          curr_graph->E);
    }

    empty_cycle = malloc(compressedSize * sizeof(*empty_cycle));
    memset((unsigned char*)empty_cycle, 0, compressedSize * sizeof(*empty_cycle));

    added = 0;
    currentBasisCyclesSize = 0;

    for (internal_index = 0;
        internal_index < data->nofURFsPerBCC[bcc_index]; ++internal_index, ++urf) {

      currentCycle = malloc(compressedSize * sizeof(*currentCycle));
      memcpy(currentCycle, prototypesForGaussian[internal_index],
          compressedSize * sizeof(*currentCycle));

      /*
       * basisCycles is in row echelon form!result = realloc(result, nextfree * sizeof(*result));
       * [this corresponds to "add row operation"]
       */
      for (k = 0; k < currentBasisCyclesSize; ++k) {
        if (RDL_bitset_test(currentCycle, k)) {
          RDL_bitset_xor_inplace(currentCycle, basisCycles[k], compressedSize);
        }
      }

      if (RDL_bitset_empty(currentCycle, empty_cycle, compressedSize)) {
        free(currentCycle);
        continue;
      }
      /*
       * the current cycle is indepedent of equal and smaller
       * sized cycles => part of the basis
       */
      oldBasisCyclesSize = currentBasisCyclesSize;
      basisCycles[currentBasisCyclesSize] = currentCycle;
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
      if (!RDL_bitset_test(currentCycle, oldBasisCyclesSize)) {
        for (k = oldBasisCyclesSize + 1; k < curr_graph->E; ++k) {
          if (RDL_bitset_test(currentCycle, k)) {
            RDL_swap_columns(basisCycles, currentBasisCyclesSize, oldBasisCyclesSize, k);
            RDL_swap_columns(prototypesForGaussian, data->nofURFsPerBCC[bcc_index], oldBasisCyclesSize, k);
            break;
          }
        }
      }

      prototype = data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->prototype;

      /* can happen for unconnected graphs */
      if (total_counter >= alloced) {
        alloced *= 2;
        result = realloc(result, (alloced) * sizeof(*result));
      }

      currBond = 0;
      result[total_counter] = malloc(sizeof(**result));
      result[total_counter]->edges = malloc(data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->weight *
                                   sizeof(*(result[total_counter]->edges)));
      result[total_counter]->weight = data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->weight;
      result[total_counter]->urf = urf;
      result[total_counter]->rcf = rcf;
      for(j=0; j < curr_graph->E; ++j) {
        if(prototype[j] == 1) {
          edge = data->bccGraphs->edge_from_bcc_mapping[bcc_index][j];
          result[total_counter]->edges[currBond][0] = data->graph->edges[edge][0];
          result[total_counter]->edges[currBond][1] = data->graph->edges[edge][1];
          ++currBond;
        }
      }

      rcf += data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];

      /*if enough inRDL_dependent cycles were found, break out of the loop*/
      ++total_counter;
      if(++added == size) {
        break;
      }
    }

    for (k = 0; k < currentBasisCyclesSize; ++k) {
      free(basisCycles[k]);
    }
    free(basisCycles);

    for (k = 0; k < data->nofURFsPerBCC[bcc_index]; ++k) {
      free(prototypesForGaussian[k]);
    }
    free(prototypesForGaussian);
    free((void*)empty_cycle);
  }

  result = realloc(result, total_counter * sizeof(*result));

  (*ptr) = result;
  return total_counter;
}

/* public, get the RCPs */
unsigned RDL_getRCPrototypes(const RDL_data *data, RDL_cycle ***ptr)
{
  RDL_cycle **result;
  unsigned nof_urf;
  unsigned nofRel=0, rcf=0;
  unsigned currFam=0, currEdge;
  unsigned i,j,k;
  const char* prototype;
  unsigned bcc_index, mapped_edge;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = malloc(sizeof(**ptr));
    return 0;
  }

  for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
    for (j = 0; j < data->nofURFsPerBCC[i]; ++j) {
      nofRel += data->urfInfoPerBCC[i]->nofCFsPerURF[j];
    }
  }

  /*allocate space*/
  result = malloc(nofRel * sizeof(*result));

  for (bcc_index = 0; bcc_index < data->bccGraphs->nof_bcc; ++bcc_index) {
    nof_urf = data->nofURFsPerBCC[bcc_index];
    for (i = 0; i < nof_urf; ++i) {
      for (j = 0; j < data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[i]; ++j, ++rcf) {
        prototype = data->urfInfoPerBCC[bcc_index]->URFs[i][j]->prototype;
        result[currFam] = malloc(sizeof(**result));
        result[currFam]->edges = malloc(data->urfInfoPerBCC[bcc_index]->URFs[i][j]->weight *
                                  sizeof(*result[currFam]->edges));
        result[currFam]->weight = data->urfInfoPerBCC[bcc_index]->URFs[i][j]->weight;
        result[currFam]->urf = i;
        result[currFam]->rcf = rcf;
        currEdge = 0;
        for (k = 0; k < data->bccGraphs->bcc_graphs[bcc_index]->E; ++k) {
          if (prototype[k] == 1) {
            mapped_edge = data->bccGraphs->edge_from_bcc_mapping[bcc_index][k];
            result[currFam]->edges[currEdge][0] = data->graph->edges[mapped_edge][0];
            result[currFam]->edges[currEdge][1] = data->graph->edges[mapped_edge][1];
            ++currEdge;
          }
        }
        ++currFam;
      }
    }
  }

  (*ptr) = result;
  return nofRel;
}

/* public, get an RDL_cycleIterator for RCFs */
RDL_cycleIterator* RDL_getRCyclesForRCFIterator(const RDL_data *data, unsigned index)
{
  RDL_cycleIterator* it;
  unsigned rcf_internal_index, urf_internal_index, urf_index;
  unsigned bcc_index;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return NULL;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return NULL;
  }

  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];

  it = RDL_initCycleIterator(RDL_RCF_IT,
      rcf_internal_index, rcf_internal_index,
      urf_internal_index, urf_internal_index,
      bcc_index, bcc_index,
      'b', data);

  return it;
}

/* public, get an RDL_cycleIterator for URFs */
RDL_cycleIterator* RDL_getRCyclesForURFIterator(const RDL_data *data, unsigned index)
{
  RDL_cycleIterator* it;
  unsigned urf_internal_index, bcc_index;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return NULL;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return NULL;
  }

  bcc_index = data->urf_to_bcc[index][0];
  urf_internal_index = data->urf_to_bcc[index][1];

  it = RDL_initCycleIterator(
      RDL_URF_IT,
      0, 0,
      urf_internal_index, urf_internal_index,
      bcc_index, bcc_index,
      'b', data);

  return it;
}

/* public, get an RDL_cycleIterator for all RCs */
RDL_cycleIterator* RDL_getRCyclesIterator(const RDL_data *data)
{
  RDL_cycleIterator* it;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return NULL;
  }

  it = RDL_initCycleIterator(
      RDL_ALL_IT,
      0, 0,
      0, 0,
      0, data->bccGraphs->nof_bcc-1,
      'b', data);

  return it;
}

/* private, get RCs from an RDL_cycleIterator */
static unsigned RDL_getRCycles_internal(
    const RDL_data *data, RDL_cycle ***ptr, RDL_cycleIterator* it)
{
  RDL_cycle **result;
  unsigned alloced, nextfree;

  alloced = RDL_RESERVED_START;
  nextfree = 0;
  result = malloc(alloced * sizeof(*result));

  while (!RDL_cycleIteratorAtEnd(it)) {
    if(nextfree == alloced) {
      alloced *= 2;
      result = realloc(result, alloced * sizeof(*result));
    }
    result[nextfree] = RDL_cycleIteratorGetCycle(it);
    ++nextfree;

    RDL_cycleIteratorNext(it);
  }
  RDL_deleteCycleIterator(it);

  result = realloc(result, nextfree * sizeof(*result));

  *ptr = result;
  return nextfree;
}

/* public, get all RCs */
unsigned RDL_getRCycles(const RDL_data *data, RDL_cycle ***ptr)
{
  RDL_cycleIterator* it;

  it = RDL_getRCyclesIterator(data);

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  return RDL_getRCycles_internal(data, ptr, it);
}

/* public, get all RCs for given URF */
unsigned RDL_getRCyclesForURF(const RDL_data *data, unsigned index, RDL_cycle ***ptr)
{
  RDL_cycleIterator* it;

  it = RDL_getRCyclesForURFIterator(data, index);

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  return RDL_getRCycles_internal(data, ptr, it);
}

/* public, get all RCs for given RCF */
unsigned RDL_getRCyclesForRCF(const RDL_data *data, unsigned index, RDL_cycle ***ptr)
{
  RDL_cycleIterator* it;

  it = RDL_getRCyclesForRCFIterator(data, index);

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  return RDL_getRCycles_internal(data, ptr, it);
}

/* public, utils function for converting cycle to bitset representation */
unsigned RDL_translateCycArray(const RDL_data *data, RDL_cycle **array,
                               unsigned number, char ***ptr)
{
  unsigned i,j,RDL_edgeIdx;
  char **result;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = RDL_alloc2DCharArray(0, 0);
    return RDL_INVALID_RESULT;
  }

  if (number < 1) {
    (*ptr) = RDL_alloc2DCharArray(0, 0);
    return 0;
  }

  result = RDL_alloc2DCharArray(number, data->graph->E);
  for(i=0; i<number;  ++i)/*initialize to 0s*/
  {
    for(j=0; j<data->graph->E; ++j)
    {
      result[i][j] = 0;
    }
  }

  for(i=0; i<number; ++i)
  {
    for(j=0; j<array[i]->weight; ++j)
    {
      RDL_edgeIdx = RDL_edgeId(data->graph, array[i]->edges[j][0],
                               array[i]->edges[j][1]);
      result[i][RDL_edgeIdx] = 1;
    }
  }
  (*ptr) = result;
  return number;
}

/* public, delete an edge array */
void RDL_deleteEdgeIdxArray(char **cycles, unsigned num)
{
  RDL_delete2DCharArray(cycles, num);
}

/* public, get the edges of the graph */
unsigned RDL_getEdgeArray(const RDL_data *data, RDL_edge **ptr)
{
  RDL_edge *result;
  unsigned idx;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = malloc(sizeof(**ptr));
    return RDL_INVALID_RESULT;
  }

  result = malloc(data->graph->E * sizeof(*result));
  for(idx=0; idx<data->graph->E; ++idx)
  {
    result[idx][0] = data->graph->edges[idx][0];
    result[idx][1] = data->graph->edges[idx][1];
  }
  (*ptr) = result;
  return data->graph->E;
}

/* public, get the number of edges in the graph */
unsigned RDL_getNofEdges(const RDL_data *data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }
  return data->graph->E;
}

/* public, get the number of ringsystems calculated */
unsigned RDL_getNofRingsystems(const RDL_data* data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->nof_bcc;
}

/* public, get the number of nodes for the ring sytem */
unsigned RDL_getNofNodesForRingsystem(const RDL_data *data, unsigned idx)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->bcc_graphs[idx]->V;
}

/* public, get the number of edges for the ring sytem */
unsigned RDL_getNofEdgesForRingsystem(const RDL_data *data, unsigned idx)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->bcc_graphs[idx]->E;

}

/* public, get the edge for the ring sytem */
unsigned RDL_getEdgesForRingsystem(
    const RDL_data *data, unsigned idx, RDL_edge** edges)
{
  unsigned i, edge;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*edges) = malloc(sizeof(**edges));
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    (*edges) = malloc(sizeof(**edges));
    return RDL_INVALID_RESULT;
  }

  (*edges) = malloc(data->bccGraphs->bcc_graphs[idx]->E * sizeof(**edges));

  for (i = 0; i < data->bccGraphs->bcc_graphs[idx]->E; ++i) {
    edge = data->bccGraphs->edge_from_bcc_mapping[idx][i];
    (*edges)[i][0] = data->graph->edges[edge][0];
    (*edges)[i][1] = data->graph->edges[edge][1];
  }

  return data->bccGraphs->bcc_graphs[idx]->E;
}

/* public, get the nodes for the ring sytem */
unsigned RDL_getNodesForRingsystem(
    const RDL_data *data, unsigned idx, RDL_node** nodes)
{
  unsigned i, node;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*nodes) = malloc(sizeof(**nodes));
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    (*nodes) = malloc(sizeof(**nodes));
    return RDL_INVALID_RESULT;
  }

  (*nodes) = malloc(data->bccGraphs->bcc_graphs[idx]->V * sizeof(*nodes));

  for (i = 0; i < data->bccGraphs->bcc_graphs[idx]->V; ++i) {
    node = data->bccGraphs->node_from_bcc_mapping[idx][i];
    (*nodes)[i] = node;
  }

  return data->bccGraphs->bcc_graphs[idx]->V;
}

/* public, get the ringsystem an edge belongs to */
unsigned RDL_getRingsystemForEdge(
    const RDL_data* data, unsigned from, unsigned to)
{
  unsigned edge;

  edge = RDL_getEdgeId(data, from, to);

  if (edge == RDL_INVALID_RESULT) {
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->edge_to_bcc_mapping[edge][0];
}

/* public, get the ID of an edge */
unsigned RDL_getEdgeId(const RDL_data *data, unsigned from, unsigned to)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (!data->graph) {
    RDL_outputFunc(RDL_ERROR, "RDL_graph is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (from >= data->graph->V || to >= data->graph->V) {
    RDL_outputFunc(RDL_ERROR, "invalid edge %u %u\n", from, to);
    return RDL_INVALID_RESULT;
  }

  return RDL_edgeId(data->graph, from, to);
}
