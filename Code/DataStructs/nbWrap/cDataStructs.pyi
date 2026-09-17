"""
Module containing an assortment of functionality for basic data structures.

At the moment the data structures defined are:
  Bit Vector classes (for storing signatures, fingerprints and the like:
    - ExplicitBitVect: class for relatively small (10s of thousands of bits) or
                       dense bit vectors.
    - SparseBitVect:   class for large, sparse bit vectors
  DiscreteValueVect:   class for storing vectors of integers
  SparseIntVect:       class for storing sparse vectors of integers
"""

from collections.abc import Iterable
import enum
from typing import overload

from rdkit.DataStructs import bvtotext as ToBitString


def ConvertToExplicit(sbv: SparseBitVect) -> ExplicitBitVect:
    """
    Converts a SparseBitVector to an ExplicitBitVector and returns the ExplicitBitVector
    """

def CreateFromBitString(bits: str) -> ExplicitBitVect:
    """Creates an ExplicitBitVect from a bit string (string of 0s and 1s)."""

def CreateFromFPSText(fps: str) -> ExplicitBitVect:
    """Creates an ExplicitBitVect from an FPS string."""

@overload
def CreateFromBinaryText(fps: bytes) -> ExplicitBitVect:
    """Creates an ExplicitBitVect from a binary string (byte array)."""

@overload
def CreateFromBinaryText(fps: str) -> ExplicitBitVect:
    """Creates an ExplicitBitVect from a string (byte array)."""

@overload
def InitFromDaylightString(sbv: SparseBitVect, s: str) -> None: ...

@overload
def InitFromDaylightString(sbv: ExplicitBitVect, s: str) -> None:
    """
    Fill a BitVect using an ASCII (Daylight) encoding of a fingerprint.

       **Arguments**
         - bv: either a _SparseBitVect_ or an _ExplicitBitVect_
         - txt: a string with the Daylight encoding (this is the text that
        the Daylight tools put in the FP field of a TDT)
    """

class SparseBitVect:
    """
    A class to store sparse bit vectors.

    This class is most useful for situations where the size of the vector
    is large and relatively few bits are set

    For smaller or denser vectors, the _ExplicitBitVect_ class is much faster.

    As you would expect, _SparseBitVects_ support a set of binary operations
    so you can do things like:
         bv3 = bv1 & bv2  (bitwise and)
         bv3 = bv1 | bv2  (bitwise or)
         bv3 = bv1 ^ bv2  (bitwise xor)
         bv3 = ~bv1       (bitwise negation) NOTE: this operation is likely
                                                      to be VERY slow and inefficient.

    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, size: int) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    def SetBit(self, which: int) -> bool:
        """Turns on a particular bit. Returns the original state of the bit."""

    def SetBitsFromList(self, onBitList: Iterable[int]) -> None:
        """
        Turns on a set of bits. The argument should be a tuple or list of bit ids.
        """

    def UnSetBit(self, which: int) -> bool:
        """Turns off a particular bit. Returns the original state of the bit."""

    def UnSetBitsFromList(self, offBitList: Iterable[int]) -> None:
        """
        Turns off a set of bits. The argument should be a tuple or list of bit ids.
        """

    def GetBit(self, which: int) -> bool:
        """Returns the value of a bit."""

    def GetNumBits(self) -> int:
        """Returns the number of bits in the vector (the vector's size)."""

    def __len__(self) -> int: ...

    def GetNumOnBits(self) -> int:
        """Returns the number of on bits."""

    def GetNumOffBits(self) -> int:
        """Returns the number of off bits."""

    def __getitem__(self, which: int) -> int: ...

    def __setitem__(self, which: int, val: int) -> int: ...

    def GetOnBits(self) -> list[int]:
        """Returns a tuple containing IDs of the on bits."""

    def ToBinary(self) -> bytes:
        """Returns an internal binary representation of the vector."""

    def FromBase64(self, inD: str) -> None:
        """Initializes the vector from a base64 encoded binary string."""

    def ToBase64(self) -> str:
        """
        Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
        """

    def ToList(self) -> list[int]:
        """Return the BitVector as a python list."""

    def __and__(self, arg: SparseBitVect, /) -> SparseBitVect: ...

    def __or__(self, arg: SparseBitVect, /) -> SparseBitVect: ...

    def __xor__(self, arg: SparseBitVect, /) -> SparseBitVect: ...

    def __invert__(self) -> SparseBitVect: ...

    def __eq__(self, arg: SparseBitVect, /) -> bool: ...

    def __ne__(self, arg: SparseBitVect, /) -> bool: ...

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class ExplicitBitVect:
    """
    A class to store explicit bit vectors.

    This class is most useful for situations where the size of the vector
    is relatively small (tens of thousands or smaller).

    For larger vectors, use the _SparseBitVect_ class instead.

    As you would expect, _ExplicitBitVects_ support a set of binary operations
    so you can do things like:
         bv3 = bv1 & bv2  (bitwise and)
         bv3 = bv1 | bv2  (bitwise or)
         bv3 = bv1 ^ bv2  (bitwise xor)
         bv3 = ~bv1       (bitwise negation)

    Bits can be set and read using either the Set/UnsetBit() and GetBit() methods
    or by indexing (i.e. bv[i] = 1 or if bv[i]).
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, size: int) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    @overload
    def __init__(self, size: int, bitsSet: bool) -> None: ...

    def SetBit(self, which: int) -> bool:
        """Turns on a particular bit. Returns the original state of the bit."""

    def SetBitsFromList(self, onBitList: Iterable[int]) -> None:
        """
        Turns on a set of bits. The argument should be a tuple or list of bit ids.
        """

    def UnSetBit(self, which: int) -> bool:
        """Turns off a particular bit. Returns the original state of the bit."""

    def UnSetBitsFromList(self, offBitList: Iterable[int]) -> None:
        """
        Turns off a set of bits. The argument should be a tuple or list of bit ids.
        """

    def GetBit(self, which: int) -> bool:
        """Returns the value of a bit."""

    def GetNumBits(self) -> int:
        """Returns the number of bits in the vector (the vector's size)."""

    def __len__(self) -> int: ...

    def GetNumOnBits(self) -> int:
        """Returns the number of on bits."""

    def GetNumOffBits(self) -> int:
        """Returns the number of off bits."""

    def __getitem__(self, which: int) -> int: ...

    def __setitem__(self, which: int, val: int) -> int: ...

    def GetOnBits(self) -> list[int]:
        """Returns a tuple containing IDs of the on bits."""

    def ToBinary(self) -> bytes:
        """Returns an internal binary representation of the vector."""

    def FromBase64(self, inD: str) -> None:
        """Initializes the vector from a base64 encoded binary string."""

    def ToBase64(self) -> str:
        """
        Converts the vector to a base64 string (the base64 encoded version of the results of ToString()).
        """

    def ToList(self) -> list[int]:
        """Return the Bitvector as a python list (faster than list(vect))"""

    def __and__(self, arg: ExplicitBitVect, /) -> ExplicitBitVect: ...

    def __or__(self, arg: ExplicitBitVect, /) -> ExplicitBitVect: ...

    def __xor__(self, arg: ExplicitBitVect, /) -> ExplicitBitVect: ...

    def __add__(self, arg: ExplicitBitVect, /) -> ExplicitBitVect: ...

    def __invert__(self) -> ExplicitBitVect: ...

    def __eq__(self, arg: ExplicitBitVect, /) -> bool: ...

    def __ne__(self, arg: ExplicitBitVect, /) -> bool: ...

    def __iadd__(self, arg: ExplicitBitVect, /) -> ExplicitBitVect: ...

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

@overload
def TanimotoSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def TanimotoSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def TanimotoSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def TanimotoSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))"""

@overload
def TanimotoSimilarity(siv1: IntSparseIntVect, siv2: IntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def TanimotoSimilarity(siv1: LongSparseIntVect, siv2: LongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def TanimotoSimilarity(siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def TanimotoSimilarity(siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """return the Tanimoto similarity between two vectors"""

@overload
def BulkTanimotoSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTanimotoSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))"""

@overload
def BulkTanimotoSimilarity(v1: IntSparseIntVect, v2: Iterable[IntSparseIntVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTanimotoSimilarity(v1: LongSparseIntVect, v2: Iterable[LongSparseIntVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTanimotoSimilarity(v1: UIntSparseIntVect, v2: Iterable[UIntSparseIntVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTanimotoSimilarity(v1: ULongSparseIntVect, v2: Iterable[ULongSparseIntVect], returnDistance: bool = False) -> list[float]:
    """
    return the Tanimoto similarities between one vector and a sequence of others
    """

def TanimotoSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))"""

def TanimotoSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / (B(bv1) + B(bv2) - B(bv1&bv2))"""

@overload
def CosineSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def CosineSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def CosineSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def CosineSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / sqrt(B(bv1) * B(bv2))"""

@overload
def BulkCosineSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkCosineSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / sqrt(B(bv1) * B(bv2))"""

def CosineSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / sqrt(B(bv1) * B(bv2))"""

def CosineSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / sqrt(B(bv1) * B(bv2))"""

@overload
def KulczynskiSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def KulczynskiSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def KulczynskiSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def KulczynskiSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))"""

@overload
def BulkKulczynskiSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkKulczynskiSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))"""

def KulczynskiSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))"""

def KulczynskiSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2)*(B(bv1) + B(bv2)) / (2 * B(bv1) * B(bv2))"""

@overload
def DiceSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def DiceSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def DiceSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def DiceSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """2*B(bv1&bv2) / (B(bv1) + B(bv2))"""

@overload
def DiceSimilarity(siv1: IntSparseIntVect, siv2: IntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def DiceSimilarity(siv1: LongSparseIntVect, siv2: LongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def DiceSimilarity(siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def DiceSimilarity(siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """return the Dice similarity between two vectors"""

@overload
def BulkDiceSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkDiceSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """2*B(bv1&bv2) / (B(bv1) + B(bv2))"""

@overload
def BulkDiceSimilarity(v1: IntSparseIntVect, v2: Iterable[IntSparseIntVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkDiceSimilarity(v1: LongSparseIntVect, v2: Iterable[LongSparseIntVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkDiceSimilarity(v1: UIntSparseIntVect, v2: Iterable[UIntSparseIntVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkDiceSimilarity(v1: ULongSparseIntVect, v2: Iterable[ULongSparseIntVect], returnDistance: bool = False) -> list[float]:
    """
    return the Dice similarities between one vector and a sequence of others
    """

def DiceSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """2*B(bv1&bv2) / (B(bv1) + B(bv2))"""

def DiceSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """2*B(bv1&bv2) / (B(bv1) + B(bv2))"""

@overload
def SokalSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def SokalSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def SokalSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def SokalSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))"""

@overload
def BulkSokalSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkSokalSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))"""

def SokalSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))"""

def SokalSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / (2*B(bv1) + 2*B(bv2) - 3*B(bv1&bv2))"""

@overload
def McConnaugheySimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def McConnaugheySimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def McConnaugheySimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def McConnaugheySimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))"""

@overload
def BulkMcConnaugheySimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkMcConnaugheySimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))"""

def McConnaugheySimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))"""

def McConnaugheySimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """(B(bv1&bv2) * (B(bv1)+B(bv2)) - B(bv1)*B(bv2)) / (B(bv1) * B(bv2))"""

@overload
def AsymmetricSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def AsymmetricSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def AsymmetricSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def AsymmetricSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / min(B(bv1),B(bv2))"""

@overload
def BulkAsymmetricSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkAsymmetricSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / min(B(bv1),B(bv2))"""

def AsymmetricSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / min(B(bv1),B(bv2))"""

def AsymmetricSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / min(B(bv1),B(bv2))"""

@overload
def BraunBlanquetSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def BraunBlanquetSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def BraunBlanquetSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def BraunBlanquetSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / max(B(bv1),B(bv2))"""

@overload
def BulkBraunBlanquetSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkBraunBlanquetSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / max(B(bv1),B(bv2))"""

def BraunBlanquetSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / max(B(bv1),B(bv2))"""

def BraunBlanquetSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / max(B(bv1),B(bv2))"""

@overload
def RusselSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def RusselSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def RusselSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def RusselSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / B(bv1)"""

@overload
def BulkRusselSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkRusselSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / B(bv1)"""

def RusselSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / B(bv1)"""

def RusselSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / B(bv1)"""

@overload
def RogotGoldbergSimilarity(bv1: SparseBitVect, bv2: SparseBitVect, returnDistance: bool = False) -> float: ...

@overload
def RogotGoldbergSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, returnDistance: bool = False) -> float: ...

@overload
def RogotGoldbergSimilarity(bv1: SparseBitVect, pkl: bytes, returnDistance: bool = False) -> float: ...

@overload
def RogotGoldbergSimilarity(bv1: ExplicitBitVect, pkl: bytes, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / B(bv1)"""

@overload
def BulkRogotGoldbergSimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkRogotGoldbergSimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / B(bv1)"""

def RogotGoldbergSimilarityNeighbors(bvqueries: Iterable[ExplicitBitVect], bvList: Iterable[ExplicitBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / B(bv1)"""

def RogotGoldbergSimilarityNeighbors_sparse(bvqueries: Iterable[SparseBitVect], bvList: Iterable[SparseBitVect]) -> list[tuple[int, float]]:
    """B(bv1&bv2) / B(bv1)"""

@overload
def TverskySimilarity(bv1: SparseBitVect, bv2: SparseBitVect, a: float, b: float, returnDistance: bool = False) -> float: ...

@overload
def TverskySimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect, a: float, b: float, returnDistance: bool = False) -> float: ...

@overload
def TverskySimilarity(bv1: SparseBitVect, pkl: bytes, a: float, b: float, returnDistance: bool = False) -> float: ...

@overload
def TverskySimilarity(bv1: ExplicitBitVect, pkl: bytes, a: float, b: float, returnDistance: bool = False) -> float:
    """B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2))"""

@overload
def TverskySimilarity(siv1: IntSparseIntVect, siv2: IntSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def TverskySimilarity(siv1: LongSparseIntVect, siv2: LongSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def TverskySimilarity(siv1: UIntSparseIntVect, siv2: UIntSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float: ...

@overload
def TverskySimilarity(siv1: ULongSparseIntVect, siv2: ULongSparseIntVect, a: float, b: float, returnDistance: bool = False, bounds: float = 0.0) -> float:
    """return the Tversky similarity between two vectors"""

@overload
def BulkTverskySimilarity(bv1: SparseBitVect, bvList: Iterable[SparseBitVect], a: float, b: float, returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTverskySimilarity(bv1: ExplicitBitVect, bvList: Iterable[ExplicitBitVect], a: float, b: float, returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / (a*B(bv1)+b*B(bv2)+(1-a-b)*B(bv1&bv2))"""

@overload
def BulkTverskySimilarity(v1: IntSparseIntVect, v2: Iterable[IntSparseIntVect], a: float, b: float, returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTverskySimilarity(v1: LongSparseIntVect, v2: Iterable[LongSparseIntVect], a: float, b: float, returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTverskySimilarity(v1: UIntSparseIntVect, v2: Iterable[UIntSparseIntVect], a: float, b: float, returnDistance: bool = False) -> list[float]: ...

@overload
def BulkTverskySimilarity(v1: ULongSparseIntVect, v2: Iterable[ULongSparseIntVect], a: float, b: float, returnDistance: bool = False) -> list[float]:
    """
    return the Tversky similarities between one vector and a sequence of others
    """

@overload
def OnBitSimilarity(v1: SparseBitVect, v2: SparseBitVect) -> float: ...

@overload
def OnBitSimilarity(v1: ExplicitBitVect, v2: ExplicitBitVect) -> float:
    """B(bv1&bv2) / B(bv1|bv2)"""

@overload
def BulkOnBitSimilarity(v1: SparseBitVect, v2: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkOnBitSimilarity(v1: ExplicitBitVect, v2: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """B(bv1&bv2) / B(bv1|bv2)"""

@overload
def AllBitSimilarity(v1: SparseBitVect, v2: SparseBitVect) -> float: ...

@overload
def AllBitSimilarity(v1: ExplicitBitVect, v2: ExplicitBitVect) -> float:
    """(B(bv1) - B(bv1^bv2)) / B(bv1)"""

@overload
def BulkAllBitSimilarity(v1: SparseBitVect, v2: Iterable[SparseBitVect], returnDistance: bool = False) -> list[float]: ...

@overload
def BulkAllBitSimilarity(v1: ExplicitBitVect, v2: Iterable[ExplicitBitVect], returnDistance: bool = False) -> list[float]:
    """(B(bv1) - B(bv1^bv2)) / B(bv1)"""

@overload
def OnBitProjSimilarity(bv1: SparseBitVect, bv2: SparseBitVect) -> list[float]: ...

@overload
def OnBitProjSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> list[float]:
    """Returns a 2-tuple: (B(bv1&bv2) / B(bv1), B(bv1&bv2) / B(bv2))"""

@overload
def OffBitProjSimilarity(bv1: SparseBitVect, bv2: SparseBitVect) -> list[float]: ...

@overload
def OffBitProjSimilarity(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> list[float]: ...

@overload
def NumBitsInCommon(bv1: SparseBitVect, bv2: SparseBitVect) -> int: ...

@overload
def NumBitsInCommon(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> int:
    """Returns the total number of bits in common between the two bit vectors"""

@overload
def OnBitsInCommon(bv1: SparseBitVect, bv2: SparseBitVect) -> list[int]: ...

@overload
def OnBitsInCommon(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> list[int]:
    """Returns the number of on bits in common between the two bit vectors"""

@overload
def OffBitsInCommon(bv1: SparseBitVect, bv2: SparseBitVect) -> list[int]: ...

@overload
def OffBitsInCommon(bv1: ExplicitBitVect, bv2: ExplicitBitVect) -> list[int]:
    """Returns the number of off bits in common between the two bit vectors"""

@overload
def FoldFingerprint(bv: SparseBitVect, foldFactor: int = 2) -> SparseBitVect: ...

@overload
def FoldFingerprint(bv: ExplicitBitVect, foldFactor: int = 2) -> ExplicitBitVect:
    """
    Folds the fingerprint by the provided amount. The default, foldFactor=2, returns a fingerprint that is half the size of the original.
    """

@overload
def AllProbeBitsMatch(probe: SparseBitVect, ref: SparseBitVect) -> bool: ...

@overload
def AllProbeBitsMatch(probe: ExplicitBitVect, ref: ExplicitBitVect) -> bool: ...

@overload
def AllProbeBitsMatch(probe: SparseBitVect, ref: bytes) -> bool: ...

@overload
def AllProbeBitsMatch(probe: ExplicitBitVect, ref: bytes) -> bool:
    """
    Returns True if all bits in the first argument match all bits in the
    vector defined by the pickle in the second argument.
    """

@overload
def BitVectToText(bv1: SparseBitVect) -> str: ...

@overload
def BitVectToText(bv1: ExplicitBitVect) -> str:
    """Returns a string of zeros and ones representing the bit vector."""

@overload
def BitVectToFPSText(bv1: SparseBitVect) -> str: ...

@overload
def BitVectToFPSText(bv1: ExplicitBitVect) -> str:
    """Returns an FPS string representing the bit vector."""

@overload
def BitVectToBinaryText(bv: SparseBitVect) -> bytes: ...

@overload
def BitVectToBinaryText(bv: ExplicitBitVect) -> bytes:
    """Returns a binary string (byte array) representing the bit vector."""

class DiscreteValueType(enum.Enum):
    ONEBITVALUE = 0

    TWOBITVALUE = 1

    FOURBITVALUE = 2

    EIGHTBITVALUE = 3

    SIXTEENBITVALUE = 4

ONEBITVALUE: DiscreteValueType = DiscreteValueType.ONEBITVALUE

TWOBITVALUE: DiscreteValueType = DiscreteValueType.TWOBITVALUE

FOURBITVALUE: DiscreteValueType = DiscreteValueType.FOURBITVALUE

EIGHTBITVALUE: DiscreteValueType = DiscreteValueType.EIGHTBITVALUE

SIXTEENBITVALUE: DiscreteValueType = DiscreteValueType.SIXTEENBITVALUE

class DiscreteValueVect:
    """
    A container class for storing unsigned integer
    values within a particular range.

    The length of the vector and type of its elements (determines the maximum value
    that can be stored) are both set at construction time.

    As you would expect, _DiscreteValueVects_ support a set of binary operations
    so you can do things like:
      dvv3 = dvv1 & dvv2  the result contains the smallest value in each entry
      dvv3 = dvv1 | dvv2  the result contains the largest value in each entry
      dvv1 += dvv2     values are truncated when necessary
      dvv3 = dvv1 + dvv2    values are truncated when necessary
      dvv1 -= dvv3    would-be negative values are set to zero
      dvv3 = dvv1 - dvv2    would-be negative values are set to zero

    Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, valType: DiscreteValueType, length: int) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    def __len__(self) -> int:
        """Get the number of entries in the vector"""

    def __setitem__(self, i: int, val: int) -> None:
        """Set the value at a specified location"""

    def __getitem__(self, i: int) -> int:
        """Get the value at a specified location"""

    def __and__(self, arg: DiscreteValueVect, /) -> DiscreteValueVect: ...

    def __or__(self, arg: DiscreteValueVect, /) -> DiscreteValueVect: ...

    def __sub__(self, arg: DiscreteValueVect, /) -> DiscreteValueVect: ...

    def __isub__(self, arg: DiscreteValueVect, /) -> DiscreteValueVect: ...

    def __add__(self, arg: DiscreteValueVect, /) -> DiscreteValueVect: ...

    def __iadd__(self, arg: DiscreteValueVect, /) -> DiscreteValueVect: ...

    def GetValueType(self) -> DiscreteValueType:
        """Get the type of value stored in the vector"""

    def GetTotalVal(self) -> int:
        """Get the sum of the values in the vector, basically L1 norm"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

@overload
def ComputeL1Norm(v1: DiscreteValueVect, v2: DiscreteValueVect) -> int:
    """Compute the distance between two discrete vector values"""

@overload
def ComputeL1Norm(arg0: RealValueVect, arg1: RealValueVect, /) -> float:
    """Compute the distance between two real vector values"""

class RealValueVect:
    """
    A container class for storing real
    values.

    The length of the vector is set at construction time.

    As you would expect, _RealValueVects_ support a set of binary operations
    so you can do things like:
      rvv3 = rvv1 & rvv2  the result contains the smallest value in each entry
      rvv3 = rvv1 | rvv2  the result contains the largest value in each entry
      rvv1 += rvv2     
      rvv3 = rvv1 + rvv2    
      rvv1 -= rvv3    
      rvv3 = rvv1 - rvv2    

    Elements can be set and read using indexing (i.e. bv[i] = 4 or val=bv[i])
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, length: int) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    def __len__(self) -> int:
        """Get the number of entries in the vector"""

    def __setitem__(self, arg0: int, arg1: float, /) -> None:
        """Set the value at a specified location"""

    def __getitem__(self, arg: int, /) -> float:
        """Get the value at a specified location"""

    def __and__(self, arg: RealValueVect, /) -> RealValueVect: ...

    def __or__(self, arg: RealValueVect, /) -> RealValueVect: ...

    def __sub__(self, arg: RealValueVect, /) -> RealValueVect: ...

    def __isub__(self, arg: RealValueVect, /) -> RealValueVect: ...

    def __add__(self, arg: RealValueVect, /) -> RealValueVect: ...

    def __iadd__(self, arg: RealValueVect, /) -> RealValueVect: ...

    def GetTotalVal(self) -> float:
        """Get the sum of the values in the vector, basically L1 norm"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class IntSparseIntVect:
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    @overload
    def __init__(self, arg: int, /) -> None:
        """Constructor"""

    def __setitem__(self, arg0: int, arg1: int, /) -> None:
        """Set the value at a specified location"""

    def __getitem__(self, arg: int, /) -> int:
        """Get the value at a specified location"""

    def __and__(self, arg: IntSparseIntVect, /) -> IntSparseIntVect: ...

    def __or__(self, arg: IntSparseIntVect, /) -> IntSparseIntVect: ...

    def __sub__(self, arg: IntSparseIntVect, /) -> IntSparseIntVect: ...

    @overload
    def __isub__(self, arg: IntSparseIntVect, /) -> IntSparseIntVect: ...

    @overload
    def __isub__(self, arg: int, /) -> IntSparseIntVect: ...

    def __add__(self, arg: IntSparseIntVect, /) -> IntSparseIntVect: ...

    @overload
    def __iadd__(self, arg: IntSparseIntVect, /) -> IntSparseIntVect: ...

    @overload
    def __iadd__(self, arg: int, /) -> IntSparseIntVect: ...

    def __eq__(self, arg: IntSparseIntVect, /) -> bool: ...

    def __ne__(self, arg: IntSparseIntVect, /) -> bool: ...

    def __itruediv__(self, arg: int, /) -> IntSparseIntVect: ...

    def __imul__(self, arg: int, /) -> IntSparseIntVect: ...

    def GetTotalVal(self, useAbs: bool = False) -> int:
        """Get the sum of the values in the vector, basically L1 norm"""

    def GetLength(self) -> int:
        """Returns the length of the vector"""

    def ToBinary(self) -> bytes:
        """returns a binary (pickle) representation of the vector"""

    def UpdateFromSequence(self, seq: Iterable[int]) -> None:
        """update the vector based on the values in the list or tuple"""

    def GetNonzeroElements(self) -> dict[int, int]:
        """returns a dictionary of the nonzero elements"""

    def ToList(self) -> list[int]:
        """Return the SparseIntVect as a python list"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class LongSparseIntVect:
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    @overload
    def __init__(self, arg: int, /) -> None:
        """Constructor"""

    def __setitem__(self, arg0: int, arg1: int, /) -> None:
        """Set the value at a specified location"""

    def __getitem__(self, arg: int, /) -> int:
        """Get the value at a specified location"""

    def __and__(self, arg: LongSparseIntVect, /) -> LongSparseIntVect: ...

    def __or__(self, arg: LongSparseIntVect, /) -> LongSparseIntVect: ...

    def __sub__(self, arg: LongSparseIntVect, /) -> LongSparseIntVect: ...

    @overload
    def __isub__(self, arg: LongSparseIntVect, /) -> LongSparseIntVect: ...

    @overload
    def __isub__(self, arg: int, /) -> LongSparseIntVect: ...

    def __add__(self, arg: LongSparseIntVect, /) -> LongSparseIntVect: ...

    @overload
    def __iadd__(self, arg: LongSparseIntVect, /) -> LongSparseIntVect: ...

    @overload
    def __iadd__(self, arg: int, /) -> LongSparseIntVect: ...

    def __eq__(self, arg: LongSparseIntVect, /) -> bool: ...

    def __ne__(self, arg: LongSparseIntVect, /) -> bool: ...

    def __itruediv__(self, arg: int, /) -> LongSparseIntVect: ...

    def __imul__(self, arg: int, /) -> LongSparseIntVect: ...

    def GetTotalVal(self, useAbs: bool = False) -> int:
        """Get the sum of the values in the vector, basically L1 norm"""

    def GetLength(self) -> int:
        """Returns the length of the vector"""

    def ToBinary(self) -> bytes:
        """returns a binary (pickle) representation of the vector"""

    def UpdateFromSequence(self, seq: Iterable[int]) -> None:
        """update the vector based on the values in the list or tuple"""

    def GetNonzeroElements(self) -> dict[int, int]:
        """returns a dictionary of the nonzero elements"""

    def ToList(self) -> list[int]:
        """Return the SparseIntVect as a python list"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class UIntSparseIntVect:
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    @overload
    def __init__(self, arg: int, /) -> None:
        """Constructor"""

    def __setitem__(self, arg0: int, arg1: int, /) -> None:
        """Set the value at a specified location"""

    def __getitem__(self, arg: int, /) -> int:
        """Get the value at a specified location"""

    def __and__(self, arg: UIntSparseIntVect, /) -> UIntSparseIntVect: ...

    def __or__(self, arg: UIntSparseIntVect, /) -> UIntSparseIntVect: ...

    def __sub__(self, arg: UIntSparseIntVect, /) -> UIntSparseIntVect: ...

    @overload
    def __isub__(self, arg: UIntSparseIntVect, /) -> UIntSparseIntVect: ...

    @overload
    def __isub__(self, arg: int, /) -> UIntSparseIntVect: ...

    def __add__(self, arg: UIntSparseIntVect, /) -> UIntSparseIntVect: ...

    @overload
    def __iadd__(self, arg: UIntSparseIntVect, /) -> UIntSparseIntVect: ...

    @overload
    def __iadd__(self, arg: int, /) -> UIntSparseIntVect: ...

    def __eq__(self, arg: UIntSparseIntVect, /) -> bool: ...

    def __ne__(self, arg: UIntSparseIntVect, /) -> bool: ...

    def __itruediv__(self, arg: int, /) -> UIntSparseIntVect: ...

    def __imul__(self, arg: int, /) -> UIntSparseIntVect: ...

    def GetTotalVal(self, useAbs: bool = False) -> int:
        """Get the sum of the values in the vector, basically L1 norm"""

    def GetLength(self) -> int:
        """Returns the length of the vector"""

    def ToBinary(self) -> bytes:
        """returns a binary (pickle) representation of the vector"""

    def UpdateFromSequence(self, seq: Iterable[int]) -> None:
        """update the vector based on the values in the list or tuple"""

    def GetNonzeroElements(self) -> dict[int, int]:
        """returns a dictionary of the nonzero elements"""

    def ToList(self) -> list[int]:
        """Return the SparseIntVect as a python list"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class ULongSparseIntVect:
    """
    A container class for storing integer
    values within a particular range.

    The length of the vector is set at construction time.

    As you would expect, _SparseIntVects_ support a set of binary operations
    so you can do things like:
      Arithmetic:
      siv1 += siv2
      siv3 = siv1 + siv2
      siv1 -= siv3
      siv3 = siv1 - siv2
      "Fuzzy" binary operations:
      siv3 = siv1 & siv2  the result contains the smallest value in each entry
      siv3 = siv1 | siv2  the result contains the largest value in each entry

    Elements can be set and read using indexing (i.e. siv[i] = 4 or val=siv[i])
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pkl: bytes) -> None: ...

    @overload
    def __init__(self, pkl: str) -> None: ...

    @overload
    def __init__(self, arg: int, /) -> None:
        """Constructor"""

    def __setitem__(self, arg0: int, arg1: int, /) -> None:
        """Set the value at a specified location"""

    def __getitem__(self, arg: int, /) -> int:
        """Get the value at a specified location"""

    def __and__(self, arg: ULongSparseIntVect, /) -> ULongSparseIntVect: ...

    def __or__(self, arg: ULongSparseIntVect, /) -> ULongSparseIntVect: ...

    def __sub__(self, arg: ULongSparseIntVect, /) -> ULongSparseIntVect: ...

    @overload
    def __isub__(self, arg: ULongSparseIntVect, /) -> ULongSparseIntVect: ...

    @overload
    def __isub__(self, arg: int, /) -> ULongSparseIntVect: ...

    def __add__(self, arg: ULongSparseIntVect, /) -> ULongSparseIntVect: ...

    @overload
    def __iadd__(self, arg: ULongSparseIntVect, /) -> ULongSparseIntVect: ...

    @overload
    def __iadd__(self, arg: int, /) -> ULongSparseIntVect: ...

    def __eq__(self, arg: ULongSparseIntVect, /) -> bool: ...

    def __ne__(self, arg: ULongSparseIntVect, /) -> bool: ...

    def __itruediv__(self, arg: int, /) -> ULongSparseIntVect: ...

    def __imul__(self, arg: int, /) -> ULongSparseIntVect: ...

    def GetTotalVal(self, useAbs: bool = False) -> int:
        """Get the sum of the values in the vector, basically L1 norm"""

    def GetLength(self) -> int:
        """Returns the length of the vector"""

    def ToBinary(self) -> bytes:
        """returns a binary (pickle) representation of the vector"""

    def UpdateFromSequence(self, seq: Iterable[int]) -> None:
        """update the vector based on the values in the list or tuple"""

    def GetNonzeroElements(self) -> dict[int, int]:
        """returns a dictionary of the nonzero elements"""

    def ToList(self) -> list[int]:
        """Return the SparseIntVect as a python list"""

    def __getstate__(self) -> tuple[bytes, dict]: ...

    def __setstate__(self, arg: object, /) -> None: ...

class FPBReader:
    """
    A class for reading and searching FPB files from Andrew Dalke's chemfp.
        Note that this functionality is still experimental and the API may
        change in future releases.
    """

    def __init__(self, filename: str, lazy: bool = False) -> None:
        """docstring"""

    def Init(self) -> None:
        """Read the fingerprints from the file. This can take a while."""

    def __len__(self) -> int: ...

    def __getitem__(self, which: int) -> tuple[ExplicitBitVect, str]: ...

    def GetNumBits(self) -> int:
        """returns the number of bits in a fingerprint"""

    def GetFP(self, idx: int) -> ExplicitBitVect:
        """returns a particular fingerprint as an ExplicitBitVect"""

    def GetBytes(self, which: int) -> bytes:
        """returns a particular fingerprint as bytes"""

    def GetId(self, idx: int) -> str:
        """returns the id of a particular fingerprint"""

    def GetTanimoto(self, which: int, bytes: bytes) -> float:
        """
        return the tanimoto similarity of a particular fingerprint to the bytes provided
        """

    def GetTanimotoNeighbors(self, bv: bytes, threshold: float = 0.7) -> tuple[tuple[float, int], ...]:
        """
        returns tanimoto similarities to and indices of all neighbors above the specified threshold
        """

    def GetTversky(self, which: int, bytes: bytes, ca: float, cb: float) -> float:
        """
        return the Tverksy similarity of a particular fingerprint to the bytes provided
        """

    def GetTverskyNeighbors(self, bv: bytes, ca: float, cb: float, threshold: float = 0.7) -> tuple[tuple[float, int], ...]:
        """
        returns Tversky similarities to and indices of all neighbors above the specified threshold
        """

    def GetContainingNeighbors(self, bv: bytes) -> tuple[int, ...]:
        """
        returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)
        """

class MultiFPBReader:
    """
    A class for reading and searching multiple FPB files from Andrew Dalke's chemfp.
        Note that this functionality is still experimental and the API may
        change in future releases.
    """

    def __init__(self, initOnSearch: bool = False) -> None:
        """docstring"""

    def Init(self) -> None:
        """Call Init() on each of our children. This can take a while."""

    def __len__(self) -> int: ...

    def GetNumBits(self) -> int:
        """returns the number of bits in a fingerprint"""

    def AddReader(self, rdr: FPBReader) -> int:
        """adds an FPBReader to our set of readers"""

    def GetReader(self, which: int) -> FPBReader:
        """returns one of our readers"""

    def GetTanimotoNeighbors(self, bv: bytes, threshold: float = 0.7, numThreads: int = 1) -> tuple[tuple[float, int, int], ...]:
        """
        returns tanimoto similarities to and indices of all neighbors above the specified threshold
        """

    def GetTverskyNeighbors(self, bv: bytes, ca: float, cb: float, threshold: float = 0.7, numThreads: int = 1) -> tuple[tuple[float, int, int], ...]:
        """
        returns Tversky similarities to and indices of all neighbors above the specified threshold
        """

    def GetContainingNeighbors(self, bv: bytes, numThreads: int = 1) -> tuple[tuple[int, int], ...]:
        """
        returns indices of neighbors that contain this fingerprint (where all bits from this fingerprint are also set)
        """

@overload
def ConvertToNumpyArray(bv: DiscreteValueVect, destArray: object) -> None: ...

@overload
def ConvertToNumpyArray(rvv: RealValueVect, destArray: object) -> None: ...

@overload
def ConvertToNumpyArray(bv: ExplicitBitVect, destArray: object) -> None: ...

@overload
def ConvertToNumpyArray(bv: IntSparseIntVect, destArray: object) -> None: ...

@overload
def ConvertToNumpyArray(bv: LongSparseIntVect, destArray: object) -> None: ...

@overload
def ConvertToNumpyArray(bv: UIntSparseIntVect, destArray: object) -> None: ...

@overload
def ConvertToNumpyArray(bv: ULongSparseIntVect, destArray: object) -> None: ...
