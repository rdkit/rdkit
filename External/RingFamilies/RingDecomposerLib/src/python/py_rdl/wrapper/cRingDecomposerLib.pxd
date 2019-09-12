# This file is part of the RingDecomposerLib, licensed
# under BSD New license (see LICENSE in the root directory).
# Copyright (c) 2016
# University of Hamburg, ZBH - Center for Bioinformatics
# Niek Andresen, Florian Flachsenberg, Matthias Rarey
# 
# Please cite:
# 
# Kolodzik, A.; Urbaczek, S.; Rarey, M.
# Unique Ring Families: A Chemically Meaningful Description
# of Molecular Ring Topologies.
# J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
# 
# Flachsenberg, F.; Andresen, N.; Rarey, M.
# RingDecomposerLib: An Open-Source Implementation of
# Unique Ring Families and Other Cycle Bases.
# J. Chem. Inf. Model., 2017, 57 (2), pp 122-126

# import all the functions we need
cdef extern from "RingDecomposerLib.h":
    ctypedef struct RDL_graph
    ctypedef struct RDL_data

    ctypedef unsigned int RDL_node
    ctypedef RDL_node RDL_edge[2]
    ctypedef struct RDL_cycle:
        RDL_edge* edges
        unsigned weight
        unsigned urf
        unsigned rcf
    ctypedef struct RDL_cycleIterator

    const unsigned RDL_INVALID_RESULT
    const unsigned RDL_DUPLICATE_EDGE

    RDL_data *RDL_calculate(RDL_graph *)
    void RDL_deleteData(RDL_data *)
    void RDL_deleteGraph(RDL_graph *)
    unsigned RDL_getNofURF(const RDL_data *)
    unsigned RDL_getNofRCF(const RDL_data *)
    unsigned RDL_getWeightForURF(const RDL_data *, unsigned)
    unsigned RDL_getWeightForRCF(const RDL_data *, unsigned)
    unsigned RDL_getNodesForURF(const RDL_data *, unsigned, RDL_node **)
    unsigned RDL_getEdgesForURF(const RDL_data *, unsigned, RDL_edge **)
    unsigned RDL_getNodesForRCF(const RDL_data *, unsigned, RDL_node **)
    unsigned RDL_getEdgesForRCF(const RDL_data *, unsigned, RDL_edge **)
        
    unsigned RDL_getURFsContainingNode(const RDL_data *, RDL_node, unsigned **)
    unsigned RDL_getURFsContainingEdge(const RDL_data *, RDL_node, RDL_node, unsigned **)
    unsigned RDL_getRCFsContainingNode(const RDL_data *, RDL_node, unsigned **)
    unsigned RDL_getRCFsContainingEdge(const RDL_data *, RDL_node, RDL_node, unsigned **)
    unsigned RDL_getSSSR(const RDL_data *, RDL_cycle ***)
    unsigned RDL_getRCPrototypes(const RDL_data *, RDL_cycle ***)
    
    double RDL_getNofRC(const RDL_data *)
    double RDL_getNofRCForURF(const RDL_data *, unsigned)

    unsigned RDL_getNofRingsystems(const RDL_data *)
    unsigned RDL_getNodesForRingsystem(const RDL_data *, unsigned, RDL_node **)
    unsigned RDL_getEdgesForRingsystem(const RDL_data *, unsigned, RDL_edge **)

    void RDL_cycleIteratorNext(RDL_cycleIterator *)
    RDL_cycle* RDL_cycleIteratorGetCycle(RDL_cycleIterator *)
    int RDL_cycleIteratorAtEnd(RDL_cycleIterator *)
    void RDL_deleteCycleIterator(RDL_cycleIterator *)

    RDL_cycleIterator* RDL_getRCyclesForRCFIterator(const RDL_data *, unsigned)
    RDL_cycleIterator* RDL_getRCyclesForURFIterator(const RDL_data *, unsigned)
    RDL_cycleIterator* RDL_getRCyclesIterator(const RDL_data *)
    void RDL_deleteCycle(RDL_cycle *)
    void RDL_deleteCycles(RDL_cycle **, unsigned)

    RDL_graph* RDL_initNewGraph(unsigned)
    unsigned RDL_addUEdge(RDL_graph *, RDL_node, RDL_node)