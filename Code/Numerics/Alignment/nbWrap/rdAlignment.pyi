"""Module containing functions to align pairs of points in 3D"""

from collections.abc import Sequence


def GetAlignmentTransform(refPoints: object, probePoints: object, weights: Sequence[float] | None = None, reflect: bool = False, maxIterations: int = 50) -> tuple:
    """
    Compute the optimal alignment (minimum RMSD) between two set of points using the quaternion algorithm

    ARGUMENTS:

      - refPoints : reference points specified as a sequence of 3-sequences or sequence of Point3Ds
      - probePoints : probe points to align to reference points - same format
        restrictions as reference points apply here
      - weights : optional sequence of weights to associate to each pair of points
      - reflect : reflect the probe points before attempting alignment
      - maxIterations : maximum number of iterations for the eigen solver

    RETURNS:

      a 2-tuple:
        - SSD value for the alignment
        - the 4x4 transform matrix, as a list of lists
    """
