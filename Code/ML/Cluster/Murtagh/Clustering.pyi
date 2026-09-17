from typing import Annotated

import numpy
from numpy.typing import NDArray


def MurtaghCluster(data: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)], nPts: int, sz: int, option: int) -> tuple:
    """
    Cluster points using Murtagh's hierarchical clustering algorithm.

    ARGUMENTS:
      - data: 2D NumPy array of point coordinates
      - nPts: number of points in the array
      - sz: number of coordinate values per point
      - option: clustering option passed to the underlying driver

    RETURNS:
      - tuple of three 1D NumPy arrays containing the cluster results
    """

def MurtaghDistCluster(data: Annotated[NDArray[numpy.float64], dict(shape=(None,), order='C', writable=False)], nPts: int, option: int) -> tuple:
    """
    Cluster using a precomputed condensed distance matrix.

    ARGUMENTS:
      - data: 1D NumPy array containing the lower triangle of the distance matrix
      - nPts: number of points in the array
      - option: clustering option passed to the underlying driver

    RETURNS:
      - tuple of three 1D NumPy arrays containing the cluster results
    """
