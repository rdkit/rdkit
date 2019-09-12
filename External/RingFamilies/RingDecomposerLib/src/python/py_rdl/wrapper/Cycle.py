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

class Cycle:
    """
    This class represents an cycle by it's edges,
    nodes, and the urf and rcf it is part of.
    """
    def __init__(self, edges, nodes, urf, rcf):
        """
        Initialize a cycle with edges, nodes and
        ring families.
        """
        ## set of edges
        self.edges = set(edges)
        ## set of nodes
        self.nodes = set(nodes)
        ## URF this cycle belongs to
        self.urf = urf
        ## RCF this cycle belongs to
        self.rcf = rcf

    @property
    ## Get the weight of the cycle.
    def weight(self):
        return len(self.edges)

    def __str__(self):
        string = "Cycle [WEIGHT {:d}, URF {:d}, RCF {:d}]".format(
            self.weight, self.urf.index, self.rcf.index)
        return string

    def __repr__(self):
      return self.__str__()
