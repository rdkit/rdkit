"""Module containing functions for basic distance geometry operations"""

from typing import Annotated

import numpy
from numpy.typing import NDArray


def DoTriangleSmoothing(boundsMatrix: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C')], tol: float = 0.0) -> bool:
    """
    Do triangle smoothing on a bounds matrix

    ARGUMENTS:

       - mat: a square Numeric array of doubles containing the bounds matrix, this matrix
              *is* modified by the smoothing

    RETURNS:

       a boolean indicating whether or not the smoothing worked.
    """

def EmbedBoundsMatrix(boundsMatrix: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C')], maxIters: int = 10, randomizeOnFailure: bool = False, numZeroFail: int = 2, weights: list = [], randomSeed: int = -1) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
    """
    Embed a bounds matrix and return the coordinates

    ARGUMENTS:

       - boundsMatrix: a square Numeric array of doubles containing the bounds matrix, this matrix
              should already be smoothed
       - maxIters: (optional) the maximum number of random distance matrices to try
       - randomizeOnFailure: (optional) toggles using random coords if a matrix fails to embed
       - numZeroFail: (optional) sets the number of zero eigenvalues to be considered a failure
       - weights: (optional) a sequence of 3 sequences (i,j,weight) indicating elements of
          the bounds matrix whose weights should be adjusted
       - randomSeed: (optional) sets the random number seed used for embedding

    RETURNS:

       a Numeric array of doubles with the coordinates
    """
