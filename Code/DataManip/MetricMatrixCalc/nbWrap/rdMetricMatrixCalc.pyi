"""
Module containing the calculator for metric matrix calculation,
e.g. similarity and distance matrices
"""

from typing import Annotated

import numpy
from numpy.typing import NDArray


def GetEuclideanDistMat(descripMat: object) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
    """
    Compute the distance matrix from a descriptor matrix using the Euclidean distance metric

    ARGUMENTS:

      descripMat - A python object of any one of the following types
                   1. A numeric array of dimensions n by m where n is the number of items in the data set
                      and m is the number of descriptors
                   2. A list of Numeric Vectors (or 1D arrays), each entry in the list corresponds
                      to descriptor vector for one item
                   3. A list (or tuple) of lists (or tuples) of values, where the values can be extracted to
                      double.

    RETURNS:
      A numeric one-dimensional array containing the lower triangle elements of the symmetric distance matrix
    """

def GetTanimotoDistMat(bitVectList: object) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
    """
    Compute the distance matrix from a list of BitVects using the Tanimoto distance metric

    ARGUMENTS:

      bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors,
                    needs to be expanded to support a list of SparseBitVects

    RETURNS:
      A numeric 1 dimensional array containing the lower triangle elements of the
      symmetric distance matrix
    """

def GetTanimotoSimMat(bitVectList: object) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
    """
    Compute the similarity matrix from a list of BitVects

    ARGUMENTS:

      bitVectList - a list of bit vectors. Currently this works only for a list of explicit bit vectors,
                    needs to be expanded to support a list of SparseBitVects

    RETURNS:
      A numeric 1 dimensional array containing the lower triangle elements of the symmetric similarity matrix
    """
