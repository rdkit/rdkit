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

# force new style classes for python3 compatibility
__metaclass__ = type


class ResultCollection:
    """
    Object for representing a result.
    """
    def __init__(self, index, nodes, edges):
        """
        index: index of result
        nodes: list or set of nodes
        edges: list or set of edges
        """
        ## index of the result
        self.index = index
        ## set of nodes of the result
        self.nodes = set(nodes)
        ## set of edges of the result
        self.edges = set(edges)

    def __str__(self):
        """
        to string.
        """
        return "UNKNOWN {:d}".format(self.index)

    def __repr__(self):
        """
        Representation
        """
        return self.__str__()


class Ringsystem(ResultCollection):
    """
    Object for representing a ringsystem.
    A ring system is a 2-connected component of the graph.
    """
    def __init__(self, index, nodes, edges):
        """
        index: index of result
        nodes: list or set of nodes
        edges: list or set of edges
        """
        super(Ringsystem, self).__init__(index, nodes, edges)


class CycleFamily(ResultCollection):
    """
    Object for representing a cycle family.
    """
    def __init__(self, index, nodes, edges, weight):
        """
        index: index of result
        nodes: list or set of nodes
        edges: list or set of edges
        weight: weight of the cycles in the family
        """
        super(CycleFamily, self).__init__(index, nodes, edges)
        ## weight of the cycle family
        self.weight = weight


class URF(CycleFamily):
    """
    Object for representing a unique ring family.
    """
    def __str__(self):
        """
        to string.
        """
        return "URF {:d}".format(self.index)

    def __init__(self, index, nodes, edges, weight):
        """
        index: index of result
        nodes: list or set of nodes
        edges: list or set of edges
        weight: weight of the cycles in the family
        """
        super(URF, self).__init__(index, nodes, edges, weight)


class RCF(CycleFamily):
    """
    Object for representing a relevant cycle family.
    """

    def __str__(self):
        """
        to string.
        """
        return "RCF {:d}".format(self.index)

    def __init__(self, index, nodes, edges, weight):
        """
        index: index of result
        nodes: list or set of nodes
        edges: list or set of edges
        weight: weight of the cycles in the family
        """
        super(RCF, self).__init__(index, nodes, edges, weight)

