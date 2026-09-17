"""Module containing the diversity and similarity pickers"""

from collections.abc import Callable, Sequence
import enum
from typing import Annotated

import numpy
from numpy.typing import NDArray

import rdkit.DataStructs.cDataStructs


class MaxMinPicker:
    """A class for diversity picking of items using the MaxMin Algorithm"""

    def __init__(self) -> None: ...

    def Pick(self, distMat: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C')], poolSize: int, pickSize: int, firstPicks: Sequence[int] = (), seed: int = -1) -> list[int]:
        """
        Pick a subset of items from a pool of items using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604

        ARGUMENTS:
          - distMat: 1D distance matrix (only the lower triangle elements)
          - poolSize: number of items in the pool
          - pickSize: number of items to pick from the pool
          - firstPicks: (optional) the first items to be picked (seeds the list)
          - seed: (optional) seed for the random number generator
        """

    def LazyPick(self, distFunc: Callable[[int, int], float], poolSize: int, pickSize: int, firstPicks: Sequence[int] = (), seed: int = -1, useCache: object | None = None) -> list[int]:
        """
        Pick a subset of items from a pool of items using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

          - distFunc: a function that should take two indices and return the
                      distance between those two points.
                      NOTE: the implementation caches distance values, so the
                      client code does not need to do so; indeed, it should not.
          - poolSize: number of items in the pool
          - pickSize: number of items to pick from the pool
          - firstPicks: (optional) the first items to be picked (seeds the list)
          - seed: (optional) seed for the random number generator
          - useCache: IGNORED
        """

    def LazyBitVectorPick(self, objects: Sequence[rdkit.DataStructs.cDataStructs.ExplicitBitVect], poolSize: int, pickSize: int, firstPicks: Sequence[int] = (), seed: int = -1, useCache: object | None = None) -> list[int]:
        """
        Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

          - vectors: a sequence of the bit vectors that should be picked from.
          - poolSize: number of items in the pool
          - pickSize: number of items to pick from the pool
          - firstPicks: (optional) the first items to be picked (seeds the list)
          - seed: (optional) seed for the random number generator
          - useCache: IGNORED.
        """

    def LazyPickWithThreshold(self, distFunc: Callable[[int, int], float], poolSize: int, pickSize: int, threshold: float, firstPicks: Sequence[int] = (), seed: int = -1) -> tuple[list[int], float]:
        """
        Pick a subset of items from a pool of items using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

          - distFunc: a function that should take two indices and return the
                      distance between those two points.
                      NOTE: the implementation caches distance values, so the
                      client code does not need to do so; indeed, it should not.
          - poolSize: number of items in the pool
          - pickSize: number of items to pick from the pool
          - threshold: stop picking when the distance goes below this value
          - firstPicks: (optional) the first items to be picked (seeds the list)
          - seed: (optional) seed for the random number generator
        """

    def LazyBitVectorPickWithThreshold(self, objects: Sequence[rdkit.DataStructs.cDataStructs.ExplicitBitVect], poolSize: int, pickSize: int, threshold: float, firstPicks: Sequence[int] = (), seed: int = -1) -> tuple[list[int], float]:
        """
        Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

          - vectors: a sequence of the bit vectors that should be picked from.
          - poolSize: number of items in the pool
          - pickSize: number of items to pick from the pool
          - threshold: stop picking when the distance goes below this value
          - firstPicks: (optional) the first items to be picked (seeds the list)
          - seed: (optional) seed for the random number generator
        """

class LeaderPicker:
    """
    A class for diversity picking of items using Roger Sayle's Leader
    algorithm (analogous to sphere exclusion). The algorithm is
    currently unpublished, but a description is available in this
    presentation from the 2019 RDKit UGM:
    https://github.com/rdkit/UGM_2019/raw/master/Presentations/Sayle_Clustering.pdf
    """

    def __init__(self) -> None: ...

    def LazyBitVectorPick(self, objects: Sequence[rdkit.DataStructs.cDataStructs.ExplicitBitVect], poolSize: int, threshold: float, pickSize: int = 0, firstPicks: Sequence[int] = (), numThreads: int = 1) -> list[int]:
        """
        Pick a subset of items from a collection of bit vectors using Tanimoto distance. The threshold value is a *distance* (i.e. 1-similarity). Note that the numThreads argument is currently ignored.
        """

    def LazyPick(self, distFunc: Callable[[int, int], float], poolSize: int, threshold: float, pickSize: int = 0, firstPicks: Sequence[int] = (), numThreads: int = 1) -> list[int]:
        """
        Pick a subset of items from a pool of items using the user-provided function to determine distances. Note that the numThreads argument is currently ignored.
        """

class HierarchicalClusterPicker:
    """A class for diversity picking of items using Hierarchical Clustering"""

    def __init__(self, clusterMethod: ClusterMethod) -> None: ...

    def Pick(self, distMat: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)], poolSize: int, pickSize: int) -> list[int]:
        """
        Pick a diverse subset of items from a pool of items using hierarchical clustering

        ARGUMENTS:
          - distMat: 1D distance matrix (only the lower triangle elements)
          - poolSize: number of items in the pool
          - pickSize: number of items to pick from the pool
        """

    def Cluster(self, distMat: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)], poolSize: int, pickSize: int) -> list[list[int]]:
        """
        Return a list of clusters of item from the pool using hierarchical clustering

        ARGUMENTS:
          - distMat: 1D distance matrix (only the lower triangle elements)
          - poolSize: number of items in the pool
          - pickSize: number of items to pick from the pool
        """

class ClusterMethod(enum.Enum):
    WARD = 1

    SLINK = 2

    CLINK = 3

    UPGMA = 4

    MCQUITTY = 5

    GOWER = 6

    CENTROID = 7

WARD: ClusterMethod = ClusterMethod.WARD

SLINK: ClusterMethod = ClusterMethod.SLINK

CLINK: ClusterMethod = ClusterMethod.CLINK

UPGMA: ClusterMethod = ClusterMethod.UPGMA

MCQUITTY: ClusterMethod = ClusterMethod.MCQUITTY

GOWER: ClusterMethod = ClusterMethod.GOWER

CENTROID: ClusterMethod = ClusterMethod.CENTROID
