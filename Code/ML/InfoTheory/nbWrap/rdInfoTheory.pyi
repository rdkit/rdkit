"""
Module containing bunch of functions for information metrics and a ranker to rank bits
"""

from collections.abc import Iterable
import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray


class InfoBitRanker:
    """
    A class to rank the bits from a series of labelled fingerprints
    A simple demonstration may help clarify what this class does.
    Here's a small set of vectors:

    >>> for i,bv in enumerate(bvs): print(bv.ToBitString(),acts[i])
    ...
    0001 0
    0101 0
    0010 1
    1110 1

    Default ranker, using infogain:

    >>> ranker = InfoBitRanker(4,2)
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ...
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ...
    3 1.000 2 0
    2 1.000 0 2
    0 0.311 0 1

    Using the biased infogain:

    >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASENTROPY)
    >>> ranker.SetBiasList((1,))
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ...
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ...
    2 1.000 0 2
    0 0.311 0 1
    1 0.000 1 1

    A chi squared ranker is also available:

    >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.CHISQUARE)
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ...
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ...
    3 4.000 2 0
    2 4.000 0 2
    0 1.333 0 1

    As is a biased chi squared:

    >>> ranker = InfoBitRanker(4,2,InfoTheory.InfoType.BIASCHISQUARE)
    >>> ranker.SetBiasList((1,))
    >>> for i,bv in enumerate(bvs): ranker.AccumulateVotes(bv,acts[i])
    ...
    >>> for bit,gain,n0,n1 in ranker.GetTopN(3): print(int(bit),'%.3f'%gain,int(n0),int(n1))
    ...
    2 4.000 0 2
    0 1.333 0 1
    1 0.000 1 1
    """

    @overload
    def __init__(self, nBits: int, nClasses: int) -> None: ...

    @overload
    def __init__(self, nBits: int, nClasses: int, infoType: InfoType) -> None: ...

    def AccumulateVotes(self, bitVect: object, label: int) -> None:
        """
        Accumulate the votes for all the bits turned on in a bit vector

        ARGUMENTS:

          - bv : bit vector either ExplicitBitVect or SparseBitVect operator
          - label : the class label for the bit vector. It is assumed that 0 <= class < nClasses
        """

    def SetBiasList(self, classList: Iterable[int]) -> None:
        """
        Set the classes to which the entropy calculation should be biased

        This list contains a set of class ids used when in the BIASENTROPY mode of ranking bits.
        In this mode, a bit must be correlated higher with one of the biased classes than all the
        other classes. For example, in a two class problem with actives and inactives, the fraction of
        actives that hit the bit has to be greater than the fraction of inactives that hit the bit

        ARGUMENTS:

          - classList : list of class ids that we want a bias towards
        """

    def SetMaskBits(self, maskBits: Iterable[int]) -> None:
        """
        Set the mask bits for the calculation

        ARGUMENTS:

          - maskBits : list of mask bits to use
        """

    def GetTopN(self, num: int) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
        """
        Returns the top n bits ranked by the information metric
        This is actually the function where most of the work of ranking is happening

        ARGUMENTS:

          - num : the number of top ranked bits that are required
        """

    def WriteTopBitsToFile(self, fileName: str) -> None:
        """Write the bits that have been ranked to a file"""

    def Tester(self, bitVect: object) -> None: ...

class InfoType(enum.Enum):
    ENTROPY = 1

    BIASENTROPY = 2

    CHISQUARE = 3

    BIASCHISQUARE = 4

ENTROPY: InfoType = InfoType.ENTROPY

BIASENTROPY: InfoType = InfoType.BIASENTROPY

CHISQUARE: InfoType = InfoType.CHISQUARE

BIASCHISQUARE: InfoType = InfoType.BIASCHISQUARE

class BitCorrMatGenerator:
    """
    A class to generate a pairwise correlation matrix between a list of bits
    The mode of operation for this class is something like this

       >>> cmg = BitCorrMatGenerator()
       >>> cmg.SetBitList(blist)
       >>> for fp in fpList:
       >>>    cmg.CollectVotes(fp)
       >>> corrMat = cmg.GetCorrMatrix()

       The resulting correlation matrix is a one dimensional nummeric array containing the
       lower triangle elements
    """

    def __init__(self) -> None: ...

    def SetBitList(self, bitList: Iterable[int]) -> None:
        """
        Set the list of bits that need to be correllated

        This may for example be their top ranking ensemble bits

        ARGUMENTS:

          - bitList : an integer list of bit IDs
        """

    def CollectVotes(self, bitVect: object) -> None:
        """
        For each pair of on bits (bi, bj) in fp increase the correlation count for the pair by 1

        ARGUMENTS:

          - fp : a bit vector to collect the fingerprints from
        """

    def GetCorrMatrix(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
        """
        Get the correlation matrix following the collection of votes from a bunch of fingerprints
        """

def InfoEntropy(resArr: Annotated[NDArray, dict(shape=(None,))]) -> float:
    """
    calculates the informational entropy of the values in an array

    ARGUMENTS:

      - resMat: pointer to a long int array containing the data
      - dim: long int containing the length of the _tPtr_ array.

    RETURNS:

      a double
    """

def InfoGain(resArr: Annotated[NDArray, dict(shape=(None, None))]) -> float:
    """
    Calculates the information gain for a variable

    ARGUMENTS:

      - varMat: a Numeric Array object
        varMat is a Numeric array with the number of possible occurrences
          of each result for reach possible value of the given variable.

        So, for a variable which adopts 4 possible values and a result which
          has 3 possible values, varMat would be 4x3

    RETURNS:

      - a Python float object

    NOTES

      - this is a dropin replacement for _PyInfoGain()_ in entropy.py
    """

def ChiSquare(resArr: Annotated[NDArray, dict(shape=(None, None))]) -> float:
    """
    Calculates the chi squared value for a variable

    ARGUMENTS:

      - varMat: a Numeric Array object
        varMat is a Numeric array with the number of possible occurrences
          of each result for reach possible value of the given variable.

        So, for a variable which adopts 4 possible values and a result which
          has 3 possible values, varMat would be 4x3

    RETURNS:

      - a Python float object
    """
