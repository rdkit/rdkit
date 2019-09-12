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

from cRingDecomposerLib cimport *
from libc.stdlib cimport free
from .RDLError import RDLError
from .Cycle import Cycle


cdef convert_cycle(RDL_cycle* cycle):
    """
    Convert the internal cycle representation Cycle python object
    """
    nof_edges = cycle[0].weight
    edges = [(int(cycle[0].edges[i][0]), int(cycle[0].edges[i][1])) for i in range(nof_edges)]
    nodes = set()
    for e1, e2 in edges:
        nodes.add(e1)
        nodes.add(e2)

    return Cycle(edges, nodes, int(cycle[0].urf), int(cycle[0].rcf))


cdef class CycleIteratorInternal:
    cdef RDL_cycleIterator* _c_iterator
    
    def __cinit__(self):
        self._c_iterator = NULL
        
    def __dealloc__(self):
        if self._c_iterator is not NULL:
            RDL_deleteCycleIterator(self._c_iterator)
    
    # we cannot set the iterator in the constructor, because
    # it's a C struct
    cdef set_iterator(self, RDL_cycleIterator *it):
        self._c_iterator = it
    
    def __call__(self):
        if self._c_iterator is NULL:
            return None
        if RDL_cycleIteratorAtEnd(self._c_iterator):
            return None
        
        cdef RDL_cycle *cycle = RDL_cycleIteratorGetCycle(self._c_iterator)
        conv_cycle = convert_cycle(cycle)
        RDL_deleteCycle(cycle)

        RDL_cycleIteratorNext(self._c_iterator)

        return conv_cycle


cdef class DataInternal:
    cdef RDL_data* _c_rdl_data
    cdef RDL_graph* _c_rdl_graph
    cdef object _rdl_graph

    def __cinit__(self, nof_nodes, edges):
        self._c_rdl_data = NULL
        self._c_rdl_graph = NULL
        self.set_graph(nof_nodes, edges)

    def __dealloc__(self):
        self._delete_c_structs()

    def _delete_c_structs(self):
        if self._c_rdl_data is not NULL:
            RDL_deleteData(self._c_rdl_data)
            self._c_rdl_graph = NULL
            self._c_rdl_data = NULL
        if self._c_rdl_graph is not NULL:
            RDL_deleteGraph(self._c_rdl_graph)
            self._c_rdl_graph = NULL
            self._c_rdl_data = NULL

    cpdef set_graph(self, nof_nodes, edges):
        """
        Set the and trigger new calculation.
        """
        self._delete_c_structs()
        # if there are no nodes, add one that calculation is possible
        if nof_nodes == 0:
            nof_nodes = 1
        self._c_rdl_graph = RDL_initNewGraph(nof_nodes)
        for index1, index2 in edges:
            index = RDL_addUEdge(self._c_rdl_graph, index1, index2)
            if index == RDL_INVALID_RESULT:
                raise RDLError('Internal error, could not add edge!')
            elif index == RDL_DUPLICATE_EDGE:
                # this is an error because wanna catch this before
                raise RDLError('Internal error, tried to add edge twice!')

    cpdef calculate(self):
        """
        Actual calculation.
        """
        if self._c_rdl_graph is NULL:
            raise RDLError('No graph specified!')
        self._c_rdl_data = RDL_calculate(self._c_rdl_graph)
        if not self.is_calculated():
            raise RDLError('Calculation failed!')

    cpdef is_calculated(self):
        """
        Check if calculation was successful (i.e. not NULL)
        """
        return self._c_rdl_data is not NULL

    cpdef _assert_calculated(self):
        if not self.is_calculated():
            raise RDLError('No data calculated! Calculate it first!')

    cpdef _assert_success(self, retval):
        """
        Assert that return value was not INVALID
        """
        if retval == RDL_INVALID_RESULT:
            raise RDLError('Internal error, could not retrieve data!')

    cpdef get_nof_urf(self):
        self._assert_calculated()
        number = RDL_getNofURF(self._c_rdl_data)
        self._assert_success(number)

        return int(number)

    cpdef get_nof_ringsystems(self):
        self._assert_calculated()
        number = RDL_getNofRingsystems(self._c_rdl_data)
        self._assert_success(number)

        return int(number)

    cpdef get_nof_rcf(self):
        self._assert_calculated()
        number = RDL_getNofRCF(self._c_rdl_data)
        self._assert_success(number)

        return int(number)

    cpdef get_nof_relevant_cycles(self):
        self._assert_calculated()
        number = RDL_getNofRC(self._c_rdl_data)
        self._assert_success(number)

        return float(number)

    cpdef get_nof_relevant_cycles_for_urf(self, unsigned urf_index):
        self._assert_calculated()
        self._check_urf_index(urf_index)
        number = RDL_getNofRCForURF(self._c_rdl_data, urf_index)
        self._assert_success(number)

        return float(number)

    cpdef _check_urf_index(self, unsigned index):
        if index >= self.get_nof_urf():
            raise IndexError("invalid index!")

    cpdef _check_rcf_index(self, unsigned index):
        if index >= self.get_nof_rcf():
            raise IndexError("invalid index!")

    cpdef _check_rs_index(self, unsigned index):
        if index >= self.get_nof_ringsystems():
            raise IndexError("invalid index!")

    cdef _convert_cycles(self, RDL_cycle** cycles, unsigned nof_cycles):
        """
        Convert the internal cycle representation Cycle python object
        """
        result = []
        for cycle_idx in range(nof_cycles):
            cycle = cycles[cycle_idx]
            result.append(convert_cycle(cycle))

        return result

    # from here on the wrapping of methods follows

    cpdef get_weight_for_urf(self, unsigned urf_index):
        self._assert_calculated()
        self._check_urf_index(urf_index)

        cdef unsigned weight
        weight = RDL_getWeightForURF(self._c_rdl_data, urf_index)
        self._assert_success(weight)
        return int(weight)

    cpdef get_weight_for_rcf(self, unsigned rcf_index):
        self._assert_calculated()
        self._check_rcf_index(rcf_index)

        cdef unsigned weight
        weight = RDL_getWeightForRCF(self._c_rdl_data, rcf_index)
        self._assert_success(weight)
        return int(weight)

    cpdef get_node_indices_for_urf(self, unsigned urf_index):
        self._assert_calculated()
        self._check_urf_index(urf_index)

        cdef RDL_node* nodes
        nof_nodes = RDL_getNodesForURF(self._c_rdl_data, urf_index, &nodes)
        self._assert_success(nof_nodes)
        result = [int(nodes[i]) for i in range(nof_nodes)]
        free(nodes)
        return result

    cpdef get_node_indices_for_rcf(self, unsigned rcf_index):
        self._assert_calculated()
        self._check_rcf_index(rcf_index)

        cdef RDL_node* nodes
        nof_nodes = RDL_getNodesForRCF(self._c_rdl_data, rcf_index, &nodes)
        self._assert_success(nof_nodes)
        result = [int(nodes[i]) for i in range(nof_nodes)]
        free(nodes)
        return result

    cpdef get_node_indices_for_ringsystem(self, unsigned rs_index):
        self._assert_calculated()
        self._check_rs_index(rs_index)

        cdef RDL_node* nodes
        nof_nodes = RDL_getNodesForRingsystem(self._c_rdl_data, rs_index, &nodes)
        self._assert_success(nof_nodes)
        result = [int(nodes[i]) for i in range(nof_nodes)]
        free(nodes)
        return result

    cpdef get_urf_indices_for_node(self, RDL_node node):
        self._assert_calculated()

        cdef unsigned* urfs
        nof_urfs = RDL_getURFsContainingNode(self._c_rdl_data, node, &urfs)
        self._assert_success(nof_urfs)

        result = [int(urfs[i]) for i in range(nof_urfs)]
        free(urfs)
        return result

    cpdef get_rcf_indices_for_node(self, RDL_node node):
        self._assert_calculated()

        cdef unsigned* rcfs
        nof_rcfs = RDL_getRCFsContainingNode(self._c_rdl_data, node, &rcfs)
        self._assert_success(nof_rcfs)

        result = [int(rcfs[i]) for i in range(nof_rcfs)]
        free(rcfs)
        return result

    cpdef get_edge_indices_for_urf(self, unsigned urf_index):
        self._assert_calculated()
        self._check_urf_index(urf_index)

        cdef RDL_edge* edges
        nof_edges = RDL_getEdgesForURF(self._c_rdl_data, urf_index, &edges)
        self._assert_success(nof_edges)

        result = [(int(edges[i][0]), int(edges[i][1])) for i in range(nof_edges)]
        free(edges)
        return result

    cpdef get_edge_indices_for_rcf(self, unsigned rcf_index):
        self._assert_calculated()
        self._check_rcf_index(rcf_index)

        cdef RDL_edge* edges
        nof_edges = RDL_getEdgesForRCF(self._c_rdl_data, rcf_index, &edges)
        self._assert_success(nof_edges)

        result = [(int(edges[i][0]), int(edges[i][1])) for i in range(nof_edges)]
        free(edges)
        return result

    cpdef get_edge_indices_for_ringsystem(self, unsigned rs_index):
        self._assert_calculated()
        self._check_rs_index(rs_index)

        cdef RDL_edge* edges
        nof_edges = RDL_getEdgesForRingsystem(self._c_rdl_data, rs_index, &edges)
        self._assert_success(nof_edges)

        result = [(int(edges[i][0]), int(edges[i][1])) for i in range(nof_edges)]
        free(edges)
        return result

    cpdef get_urf_indices_for_edge(self, RDL_node node1, RDL_node node2):
        self._assert_calculated()

        cdef unsigned* urfs
        nof_urfs = RDL_getURFsContainingEdge(self._c_rdl_data, node1, node2, &urfs)
        self._assert_success(nof_urfs)

        result = [int(urfs[i]) for i in range(nof_urfs)]
        free(urfs)
        return result

    cpdef get_rcf_indices_for_edge(self, RDL_node node1, RDL_node node2):
        self._assert_calculated()

        cdef unsigned* rcfs
        nof_rcfs = RDL_getRCFsContainingEdge(self._c_rdl_data, node1, node2, &rcfs)
        self._assert_success(nof_rcfs)

        result = [int(rcfs[i]) for i in range(nof_rcfs)]
        free(rcfs)
        return result

    cpdef get_relevant_cycles_with_indices_for_urf(self, unsigned urf_index):
        self._assert_calculated()
        self._check_urf_index(urf_index)

        cdef RDL_cycleIterator* it = RDL_getRCyclesForURFIterator(self._c_rdl_data, urf_index)
        py_iterator = CycleIteratorInternal()
        py_iterator.set_iterator(it)
        
        return iter(py_iterator, None)

    cpdef get_relevant_cycles_with_indices_for_rcf(self, unsigned rcf_index):
        self._assert_calculated()
        self._check_rcf_index(rcf_index)

        cdef RDL_cycleIterator* it = RDL_getRCyclesForRCFIterator(self._c_rdl_data, rcf_index)
        py_iterator = CycleIteratorInternal()
        py_iterator.set_iterator(it)
        
        return iter(py_iterator, None)

    cpdef get_sssr(self):
        self._assert_calculated()

        cdef RDL_cycle** cycles
        nof_cycles = RDL_getSSSR(self._c_rdl_data, &cycles)
        self._assert_success(nof_cycles)

        result = self._convert_cycles(cycles, nof_cycles)
        RDL_deleteCycles(cycles, nof_cycles)
        return result

    cpdef get_rcs(self,):
        self._assert_calculated()

        cdef RDL_cycleIterator* it = RDL_getRCyclesIterator(self._c_rdl_data)
        py_iterator = CycleIteratorInternal()
        py_iterator.set_iterator(it)
        
        return iter(py_iterator, None)

    cpdef get_rcps(self):
        self._assert_calculated()

        cdef RDL_cycle** cycles
        nof_cycles = RDL_getRCPrototypes(self._c_rdl_data, &cycles)
        self._assert_success(nof_cycles)

        result = self._convert_cycles(cycles, nof_cycles)
        RDL_deleteCycles(cycles, nof_cycles)
        return result
