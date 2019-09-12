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
import warnings

# force new style classes for python3 compatibility
__metaclass__ = type

class Graph():
    """
    Object for holding a graph for calculations.
    """
    def __init__(self):
        """
        Create a new (empty) graph.

        Use this method and add individual edges or
        use from_edges instead.
        """
        # mapping of nodes to indices
        self._node_map = {}
        # mapping of indices to nodes
        self._nodes = []
        # mapping of indices to edges
        self._edges = {}
        # mapping of edges to indices
        self._edge_map = {}

    @classmethod
    def from_edges(cls, edge_iterable,
                   get_node_1=operator.itemgetter(0),
                   get_node_2=operator.itemgetter(1),
                   get_node_id=lambda x: x,
                   get_edge_id=lambda x: x):
        """
        Create a new Graph datastructure
        edge_iterable: an iterable containing the edges
        get_node_1: function retrieving the first node of an edge
        get_node_2: function retrieving the second node of an edge
        get_node_id: function for retrieven node identifier
        get_edge_id: function for retrieving edge identifier

        This classmethod creates a new Graph.
        The edge_iterable can be any iterable (list, generator etc.).
        An edge can be as simple as a pair of nodes.
        The get_node_ function must return the respective nodes adjacent
        to an edge (default: first and second element of the edge).
        The node must be hashable. If it's not, use some unique ID
        instead.
        The edge must be hashable. If it's not, use some unique ID
        instead (provided by get_edge_id).
        """
        graph = cls()
        for edge in edge_iterable:
            graph.add_edge(edge, get_node_1, get_node_2, get_node_id, get_edge_id)

        return graph

    def _add_or_get_node_index(self, node):
        if node not in self._node_map:
            self._node_map[node] = len(self._node_map)
            self._nodes.append(node)
        return self._node_map[node]

    def add_edge(self, edge,
                 get_node_1=operator.itemgetter(0),
                 get_node_2=operator.itemgetter(1),
                 get_node_id=lambda x: x,
                 get_edge_id=lambda x: x):
        """
        Add an edge to this graph.
        edge: edge to be added to the graph
        get_node_1: function retrieving the first node of an edge
        get_node_2: function retrieving the second node of an edge
        get_node_id: fucntion for retrieving node identifier
        get_edge_id: function for retrieving edge identifier
        returns the edge_id if successfull, None if edge already present
        (and warns in that case)

        This function adds an edge to the graph.
        DO NOT USE AFTER CALCULATION
        """
        node1 = get_node_1(edge)
        node2 = get_node_2(edge)
        node1_id = get_node_id(node1)
        node2_id = get_node_id(node2)

        edge_id = get_edge_id(edge)

        index1 = self._add_or_get_node_index(node1_id)
        index2 = self._add_or_get_node_index(node2_id)

        index_pair = (min(index1, index2), max(index1, index2))

        if index_pair in self._edges:
            warnings.warn("trying to add edge that's already present!")
            return None

        self._edges[index_pair] = edge_id

        self._edge_map[edge_id] = index_pair

        return edge_id

    def get_edge_for_indices(self, node_index1, node_index2):
        """
        Returns the edge which is formed by the two nodes
        (as indices).
        node_index1: index of the first node
        node_index2: index of the second node
        returns edge adjacent by two nodes
        raises KeyError if the edge does not exist

        This function is used internally by Calculator, you
        probably won't need this function.
        """
        index_pair = (min(node_index1, node_index2),
                      max(node_index1, node_index2))
        return self._edges[index_pair]

    def get_indices_for_edge(self, edge):
        """
        Returns the edge which is formed by the two nodes
        (as indices).
        edge: edge
        returns edge adjacent by two nodes
        raises KeyError if the edge does not exist

        This function is used internally by Calculator, you
        probably won't need this function.
        """
        return self._edge_map[edge]

    def get_node_for_index(self, node_index):
        """
        Returns the node with the internal index.
        node_index: index of the node
        returns the node object with this index
        raises IndexError if node does not exist

        This function is used internally by Calculator, you
        probably won't need this function.
        """
        return self._nodes[node_index]

    def get_index_for_node(self, node):
        """
        Returns the internal index for the node.
        node: node object
        returns the index for the node object
        raises KeyError if node does not exist

        This function is used internally by Calculator, you
        probably won't need this function.
        """
        return self._node_map[node]

    def get_nof_nodes(self):
        """
        Get the number of nodes in the graph.
        """
        return len(self._nodes)

    def get_edges(self):
        """
        Get the edges in the graph.

        This function is used internally by Calculator, you
        won't need this function.
        """
        return self._edges
