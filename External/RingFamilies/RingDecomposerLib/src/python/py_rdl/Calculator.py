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

import operator

# force new style classes for python3 compatibility
__metaclass__ = type

from .CycleFamily import URF, RCF, CycleFamily, Ringsystem
from .Graph import Graph
from .wrapper import DataInternal
from .wrapper import RDLError
from .wrapper import Cycle


class Calculator():
    """
    Object for representing a calculation result.
    """
    def __init__(self, graph):
        """
        Initialize a new Calculator object with a Graph.

        Use @classmethod provided instead
        (get_calculated_result_for_graph, get_calculated_result).
        
        You can access URFs, RCFs and ring systems with the
        public attributes `urfs`, `rcfs` and `ringsystems`.
        The objects can be uses as indices for all functions below,
        that operate on the respective ring descriptions.
        """
        ## @brief List of calculated URFs.
        ## You can access individual URFs using this list and
        ## use the URF objects as index to all functions below,
        ## that take an URF as argument.
        self.urfs = None
        ## @brief List of calculated RCFs.
        ## You can access individual RCFs using this list and
        ## use the RCF objects as index to all functions below,
        ## that take an RCF as argument.
        self.rcfs = None
        ## @brief List of calculated ring systems.
        ## You can access individual ring systems using this list and
        ## use the Ringsystem objects as index to all functions below,
        ## that take an ring system as argument.
        self.ringsystems = None
        self._graph = None
        self._data_internal = None
        self.set_graph(graph)

    @classmethod
    def get_calculated_result_for_graph(cls, graph):
        """
        Calculate the ring topologies for a Graph
        and return the results.
        """
        urf_data = cls(graph)
        urf_data.calculate()
        return urf_data

    @classmethod
    def get_calculated_result(cls,
                              edge_iterable,
                              get_node_1=operator.itemgetter(0),
                              get_node_2=operator.itemgetter(1),
                              get_node_id=lambda x: x,
                              get_edge_id=lambda x: x):
        """
        Calculate the ring topologies and return the results.

        edge_iterable: an iterable containing the edges
        get_node_1: function retrieving the first node of an edge
        get_node_2: function retrieving the second node of an edge
        get_node_id: function for retrieving node identifier
        get_edge_id: function for retrieving edge identifier
        returns the calculated Calculator

        This @classmethod creates a new Graph.
        The edge_iterable can be any iterable (list, generator etc.).
        An edge can be as simple as a pair of nodes.
        The get_node function must return the respective nodes adjacent
        to an edge.
        The node must be hashable. If it's not, use some unique ID
        instead.
        The edge must be hashable. If it's not, use some unique ID
        instead.
        """
        graph = Graph.from_edges(edge_iterable, get_node_1,
            get_node_2, get_node_id, get_edge_id)
        return cls.get_calculated_result_for_graph(graph)

    def set_graph(self, graph):
        """
        Set the graph data structure.

        This has to be a Graph object. See Graph
        for a generic interface.
        """
        self._graph = graph
        try:
            keys = graph.get_edges().iterkeys()
        except AttributeError:
            keys = graph.get_edges().keys()

        self._data_internal = DataInternal(graph.get_nof_nodes(), keys)

    def _convert_urfs(self):
        """
        Convert the internal URF representation.
        """
        self.urfs = []
        for i in range(self.get_nof_urf()):
            edges = self.get_edges_for_urf(i)
            nodes = self.get_nodes_for_urf(i)
            weight = self.get_weight_for_urf(i)
            urf = URF(i, nodes, edges, weight)
            self.urfs.append(urf)

    def _convert_rcfs(self):
        """
        Convert the internal RCF representation.
        """
        self.rcfs = []
        for i in range(self.get_nof_rcf()):
            edges = self.get_edges_for_rcf(i)
            nodes = self.get_nodes_for_rcf(i)
            weight = self.get_weight_for_rcf(i)
            rcf = RCF(i, nodes, edges, weight)
            self.rcfs.append(rcf)

    def _convert_rs(self):
        """
        Convert ringsystems
        """
        self.ringsystems = []
        for i in range(self.get_nof_ringsystems()):
            edges = self.get_edges_for_ringsystem(i)
            nodes = self.get_nodes_for_ringsystem(i)
            rs = Ringsystem(i, nodes, edges)
            self.ringsystems.append(rs)

    def calculate(self):
        """
        Calculate results for the graph datastructure.
        raises RDLError, if the calculation fails

        When creating a Calculator object with the
        @classmethod provided
        (get_calculated_result_for_graph, get_calculated_result)
        calling this function
        is not necessary.
        """
        self._data_internal.calculate()
        self._convert_urfs()
        self._convert_rcfs()
        self._convert_rs()

    def is_calculated(self):
        """
        Check calculation status.
        returns True if calculation was successful, False otherwise.
        """
        return self._data_internal.is_calculated()

    def get_nof_urf(self):
        """
        Get the number of URFs.
        raises RDLError if calculation wasn't successful
        """
        return self._data_internal.get_nof_urf()

    def get_nof_rcf(self):
        """
        Get the number of RCFs.
        raises RDLError if calculation wasn't successful
        """
        return self._data_internal.get_nof_rcf()

    def get_nof_ringsystems(self):
        """
        Get the number of Ringsystems.
        A ring system is a 2-connected component of the graph.
        raises RDLError if calculation wasn't successful
        """
        return self._data_internal.get_nof_ringsystems()

    def __len__(self):
        """
        Builtin for getting the length (number of URFs).
        """
        return self.get_nof_urf()

    def _convert_indices_to_nodes(self, node_indices):
        """
        Converts the internal indices back to the node identifiers.
        """
        return [self._graph.get_node_for_index(index) for index in node_indices]

    def _convert_indices_to_edges(self, edge_indices):
        """
        Converts the internal edge indicies back to edge identifier.
        """
        return [self._graph.get_edge_for_indices(index1, index2)
                for index1, index2 in edge_indices]

    def _convert_cycles(self, edge_indices):
        """
        Convert the internal cycles.
        """
        return Cycle(self._convert_indices_to_edges(edge_indices.edges),
                     self._convert_indices_to_nodes(edge_indices.nodes),
                     self.urfs[edge_indices.urf], self.rcfs[edge_indices.rcf])

    def get_nodes_for_urf(self, urf_index):
        """
        Get the nodes in this URF.
        raises RDLError if calculation wasn't successful
        returns list of nodes in this URF
        """
        if isinstance(urf_index, URF):
            urf_index = urf_index.index
        node_indices = self._data_internal.get_node_indices_for_urf(urf_index)
        return self._convert_indices_to_nodes(node_indices)

    def get_nodes_for_rcf(self, rcf_index):
        """
        Get the nodes in this RCF.
        raises RDLError if calculation wasn't successful
        returns list of nodes in this RCF
        """
        if isinstance(rcf_index, RCF):
            rcf_index = rcf_index.index
        node_indices = self._data_internal.get_node_indices_for_rcf(rcf_index)
        return self._convert_indices_to_nodes(node_indices)

    def get_nodes_for_ringsystem(self, rs_index):
        """
        Get the nodes in this Ringsystem.
        A ring system is a 2-connected component of the graph.
        raises RDLError if calculation wasn't successful
        returns list of nodes in this Ringsystem
        """
        if isinstance(rs_index, Ringsystem):
            rs_index = rs_index.index
        node_indices = self._data_internal.get_node_indices_for_ringsystem(rs_index)
        return self._convert_indices_to_nodes(node_indices)

    def get_edges_for_urf(self, urf_index):
        """
        Get the edges in this URF.
        raises RDLError if calculation wasn't successful
        returns list of edges in this URF
        """
        if isinstance(urf_index, URF):
            urf_index = urf_index.index
        edge_indices = self._data_internal.get_edge_indices_for_urf(urf_index)
        return self._convert_indices_to_edges(edge_indices)

    def get_edges_for_rcf(self, rcf_index):
        """
        Get the edges in this RCF.
        raises RDLError if calculation wasn't successful
        returns list of edges in this RCF
        """
        if isinstance(rcf_index, RCF):
            rcf_index = rcf_index.index
        edge_indices = self._data_internal.get_edge_indices_for_rcf(rcf_index)
        return self._convert_indices_to_edges(edge_indices)

    def get_edges_for_ringsystem(self, rs_index):
        """
        Get the edges in this Ringsystem.
        A ring system is a 2-connected component of the graph.
        raises RDLError if calculation wasn't successful
        returns list of edges in this Ringsystem
        """
        if isinstance(rs_index, Ringsystem):
            rs_index = rs_index.index
        edge_indices = self._data_internal.get_edge_indices_for_ringsystem(rs_index)
        return self._convert_indices_to_edges(edge_indices)

    def get_weight_for_urf(self, urf_index):
        """
        Get weight of this URF.
        raises RDLError if calculation wasn't successful
        return weight for given URF
        """
        if isinstance(urf_index, URF):
            urf_index = urf_index.index
        return self._data_internal.get_weight_for_urf(urf_index)

    def get_weight_for_rcf(self, rcf_index):
        """
        Get weight of this RCF.
        raises RDLError if calculation wasn't successful
        return weight for given RCF
        """
        if isinstance(rcf_index, RCF):
            rcf_index = rcf_index.index
        return self._data_internal.get_weight_for_rcf(rcf_index)

    def get_relevant_cycles_for_urf(self, urf_index):
        """
        Get the cycles in this of URF.
        raises RDLError if calculation wasn't successful
        returns generator for enumerating relevant cycles in this URF
        """
        if isinstance(urf_index, URF):
            urf_index = urf_index.index
        for ring in self._data_internal.get_relevant_cycles_with_indices_for_urf(urf_index):
            yield self._convert_cycles(ring)

    def get_relevant_cycles_for_rcf(self, rcf_index):
        """
        Get the cycles in this of RCF.
        raises RDLError if calculation wasn't successful
        returns generator for enumerating relevant cycles in this RCF
        """
        if isinstance(rcf_index, RCF):
            rcf_index = rcf_index.index
        for ring in self._data_internal.get_relevant_cycles_with_indices_for_rcf(rcf_index):
            yield self._convert_cycles(ring)

    def get_urfs_for_node(self, node):
        """
        Get the URFs this node is part of.
        raises RDLError if calculation wasn't successful
        returns a list of URFs
        """
        try:
            node_index = self._graph.get_index_for_node(node)
            indices = self._data_internal.get_urf_indices_for_node(node_index)
        # can happen if graph is not connected and node was never added
        except KeyError:
            return []
        return [self.__getitem__(i) for i in indices]

    def get_urfs_for_edge(self, edge):
        """
        Get the URFs this edge is part of.
        raises RDLError if calculation wasn't successful
        returns list of URFs
        """
        index1, index2 = self._graph.get_indices_for_edge(edge)
        indices = self._data_internal.get_urf_indices_for_edge(index1, index2)
        return [self.__getitem__(i) for i in indices]

    def get_rcfs_for_node(self, node):
        """
        Get the RCFs this node is part of.
        raises RDLError if calculation wasn't successful
        returns a list of RCFs
        """
        try:
            node_index = self._graph.get_index_for_node(node)
            indices = self._data_internal.get_rcf_indices_for_node(node_index)
        # can happen if graph is not connected and node was never added
        except KeyError:
            return []
        return [self.rcfs[i] for i in indices]

    def get_rcfs_for_edge(self, edge):
        """
        Get the RCFs this edge is part of.
        raises RDLError if calculation wasn't successful
        returns list of RCFs
        """
        index1, index2 = self._graph.get_indices_for_edge(edge)
        indices = self._data_internal.get_rcf_indices_for_edge(index1, index2)
        return [self.rcfs[i] for i in indices]

    def get_sssr(self):
        """
        Get a minimal cycle base.
        raises RDLError if calculation wasn't successful
        returns list of cycles
        """
        indices_rings = self._data_internal.get_sssr()
        return [self._convert_cycles(ring) for ring in indices_rings]

    def get_relevant_cycles(self):
        """
        Get relevant cycles.
        raises RDLError if calculation wasn't successful
        returns generator for enumerating relevant cycles
        """
        for ring in self._data_internal.get_rcs():
            yield self._convert_cycles(ring)

    def get_nof_relevant_cycles(self):
        """
        Get the number of relevant cycles.
        raises RDLError if calculation wasn't successful
        returns number of RCs
        """
        return self._data_internal.get_nof_relevant_cycles()

    def get_nof_relevant_cycles_for_urf(self, urf_index):
        """
        Get the number of relevant cycles for given URF.
        raises RDLError if calculation wasn't successful
        returns number of RCs
        """
        if isinstance(urf_index, URF):
            urf_index = urf_index.index
        return self._data_internal.get_nof_relevant_cycles_for_urf(urf_index)

    def get_relevant_cycle_prototypes(self):
        """
        Get relevant cycle prototypes (one for each RCF).
        raises RDLError if calculation wasn't successful
        returns list of relevant cycle prototypes
        """
        indices_rings = self._data_internal.get_rcps()
        return [self._convert_cycles(ring) for ring in indices_rings]

    def __getitem__(self, item):
        """
        Get the item-th URF.
        item: The index of the URF
        raises RDLError if calculation wasn't successful
        raises IndexError if the index is out of range
        returns the item-th URF.
        """
        if not self.is_calculated():
            raise RDLError('No calculated result! Calculate first!')
        return self.urfs.__getitem__(item)

