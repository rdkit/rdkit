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

#ifndef RDL_TARJAN_H
#define RDL_TARJAN_H

#include "RDLforward.h"
#include "RingDecomposerLib.h"

struct RDL_BCCGraph {
    unsigned nof_bcc;

    RDL_graph** bcc_graphs;

    unsigned** edge_to_bcc_mapping;
    unsigned** node_to_bcc_mapping;
    unsigned* nof_bcc_per_node;

    unsigned** edge_from_bcc_mapping;
    unsigned** node_from_bcc_mapping;
    unsigned* nof_nodes_per_bcc;
    unsigned* nof_edges_per_bcc;

    const RDL_graph* complete_graph;
};

void RDL_deleteBCCGraph(RDL_BCCGraph* graph);

RDL_BCCGraph* RDL_tarjanBCC(const RDL_graph* graph);

#endif
