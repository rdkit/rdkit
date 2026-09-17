"""Module containing geometry objects like points, grids, etc."""

from collections.abc import Sequence
from typing import overload

import rdkit.DataStructs.cDataStructs


class Point2D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, x: float, y: float) -> None: ...

    @overload
    def __init__(self, arg: Point3D, /) -> None: ...

    @property
    def x(self) -> float: ...

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float: ...

    @y.setter
    def y(self, arg: float, /) -> None: ...

    def __getitem__(self, arg: int, /) -> float: ...

    def __len__(self) -> int: ...

    def __sub__(self, arg: Point2D, /) -> Point2D: ...

    def __isub__(self, arg: Point2D, /) -> Point2D: ...

    def __add__(self, arg: Point2D, /) -> Point2D: ...

    def __iadd__(self, arg: Point2D, /) -> Point2D: ...

    def __mul__(self, arg: float, /) -> Point2D: ...

    def __truediv__(self, arg: float, /) -> Point2D: ...

    def __imul__(self, arg: float, /) -> Point2D: ...

    def __itruediv__(self, arg: float, /) -> Point2D: ...

    def Normalize(self) -> None:
        """Normalize the vector (using L2 norm)"""

    def Length(self) -> float:
        """Length of the vector"""

    def LengthSq(self) -> float:
        """Square of the length"""

    def DotProduct(self, other: Point2D) -> float:
        """Dot product with another point"""

    def AngleTo(self, other: Point2D) -> float:
        """determines the angle between a vector to this point (between 0 and PI)"""

    def SignedAngleTo(self, other: Point2D) -> float:
        """
        determines the signed angle between a vector to this point (between 0 and 2*PI)
        """

    def DirectionVector(self, other: Point2D) -> Point2D:
        """return a normalized direction vector from this point to another"""

    def __getstate__(self) -> tuple[float, float]: ...

    def __setstate__(self, arg: tuple[float, float], /) -> None: ...

class Point3D:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, x: float, y: float, z: float) -> None: ...

    @overload
    def __init__(self, arg: Point3D) -> None: ...

    @property
    def x(self) -> float: ...

    @x.setter
    def x(self, arg: float, /) -> None: ...

    @property
    def y(self) -> float: ...

    @y.setter
    def y(self, arg: float, /) -> None: ...

    @property
    def z(self) -> float: ...

    @z.setter
    def z(self, arg: float, /) -> None: ...

    def __getitem__(self, arg: int, /) -> float: ...

    def __len__(self) -> int: ...

    def __sub__(self, arg: Point3D, /) -> Point3D: ...

    def __isub__(self, arg: Point3D, /) -> Point3D: ...

    def __add__(self, arg: Point3D, /) -> Point3D: ...

    def __iadd__(self, arg: Point3D, /) -> Point3D: ...

    def __mul__(self, arg: float, /) -> Point3D: ...

    def __truediv__(self, arg: float, /) -> Point3D: ...

    def __imul__(self, arg: float, /) -> Point3D: ...

    def __itruediv__(self, arg: float, /) -> Point3D: ...

    def Normalize(self) -> None:
        """Normalize the vector (using L2 norm)"""

    def Length(self) -> float:
        """Length of the vector"""

    def LengthSq(self) -> float:
        """Square of the length"""

    def DotProduct(self, other: Point3D) -> float:
        """Dot product with another point"""

    def CrossProduct(self, other: Point3D) -> Point3D:
        """Get the cross product between two points"""

    def AngleTo(self, other: Point3D) -> float:
        """determines the angle between a vector to this point (between 0 and PI)"""

    def SignedAngleTo(self, other: Point3D) -> float:
        """
        determines the signed angle between a vector to this point (between 0 and 2*PI)
        """

    def DirectionVector(self, other: Point3D) -> Point3D:
        """return a normalized direction vector from this point to another"""

    def Distance(self, pt2: Point3D) -> float:
        """Distance from this point to another point"""

    def __getstate__(self) -> tuple[float, float, float]: ...

    def __setstate__(self, arg: tuple[float, float, float], /) -> None: ...

class PointND:
    def __init__(self, dim: int) -> None: ...

    def __getitem__(self, arg: int, /) -> float: ...

    def __setitem__(self, arg0: int, arg1: float, /) -> float: ...

    def __len__(self) -> int: ...

    def __sub__(self, arg: PointND, /) -> PointND: ...

    def __isub__(self, arg: PointND, /) -> PointND: ...

    def __add__(self, arg: PointND, /) -> PointND: ...

    def __iadd__(self, arg: PointND, /) -> PointND: ...

    def __mul__(self, arg: float, /) -> PointND: ...

    def __truediv__(self, arg: float, /) -> PointND: ...

    def __imul__(self, arg: float, /) -> PointND: ...

    def __itruediv__(self, arg: float, /) -> PointND: ...

    def Normalize(self) -> None:
        """Normalize the vector (using L2 norm)"""

    def Length(self) -> float:
        """Length of the vector"""

    def LengthSq(self) -> float:
        """Square of the length"""

    def DotProduct(self, other: PointND) -> float:
        """Dot product with another point"""

    def AngleTo(self, other: PointND) -> float:
        """determines the angle between a vector to this point (between 0 and PI)"""

    def DirectionVector(self, other: PointND) -> PointND:
        """return a normalized direction vector from this point to another"""

    def __getstate__(self) -> list[float]: ...

    def __setstate__(self, arg: Sequence[float], /) -> None: ...

def ComputeDihedralAngle(pt1: Point3D, pt2: Point3D, pt3: Point3D, pt4: Point3D) -> float:
    """calculates the dihedral angle determined by four Point3D objects"""

def ComputeSignedDihedralAngle(pt1: Point3D, pt2: Point3D, pt3: Point3D, pt4: Point3D) -> float:
    """
    calculates the signed dihedral angle determined by four Point3D objects
    """

class UniformGrid3D:
    """
    Class to represent a uniform three-dimensional
        cubic grid. Each grid point can store a positive integer value. For the sake
        of efficiency these value can either be binary or fit in 2, 4, 8 or 16 bits
    """

    @overload
    def __init__(self, dimX: float, dimY: float, dimZ: float, spacing: float = 0.5, valType: rdkit.DataStructs.cDataStructs.DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, offSet: Point3D | None = None) -> None:
        """Constructor for a UniformGrid3D object"""

    @overload
    def __init__(self, pkl: str) -> None:
        """pickle constructor"""

    @overload
    def __init__(self, arg: bytes, /) -> None: ...

    def GetGridPointIndex(self, point: Point3D) -> int:
        """Get the index to the grid point closest to the specified point"""

    def GetGridIndex(self, xi: int, yi: int, zi: int) -> int:
        """
        Get the index to the grid point with the three integer indices provided
        """

    def GetGridIndices(self, idx: int) -> tuple[int, int, int]:
        """Returns the integer indices of the grid index provided."""

    def GetValPoint(self, pt: Point3D) -> int:
        """Get the value at the closest grid point"""

    def GetVal(self, id: int) -> int:
        """Get the value at the specified grid index"""

    def SetVal(self, id: int, val: int) -> None:
        """Set the value at the specified grid index"""

    def SetValPoint(self, pt: Point3D, val: int) -> None:
        """Set the value at grid point closest to the specified point"""

    def GetGridPointLoc(self, pointId: int) -> Point3D:
        """Get the location of the specified grid point"""

    def GetSize(self) -> int:
        """Get the size of the grid (number of grid points)"""

    def GetNumX(self) -> int:
        """Get the number of grid points along x-axis"""

    def GetNumY(self) -> int:
        """Get the number of grid points along y-axis"""

    def GetNumZ(self) -> int:
        """Get the number of grid points along z-axis"""

    def GetOffset(self) -> Point3D:
        """Get the location of the center of the grid"""

    def GetSpacing(self) -> float:
        """Get the grid spacing"""

    def GetOccupancyVect(self) -> rdkit.DataStructs.cDataStructs.DiscreteValueVect:
        """Get the occupancy vector for the grid"""

    def CompareParams(self, other: UniformGrid3D) -> bool:
        """Compare the parameters between two grid object"""

    def SetSphereOccupancy(self, center: Point3D, radius: float, stepSize: float, maxLayers: int = -1, ignoreOutOfBound: bool = True) -> None:
        """
        Set the occupancy on the grid for a sphere or specified radius
         and multiple layers around this sphere, with decreasing values of 
        occupancy
        """

    def __iadd__(self, arg: UniformGrid3D, /) -> UniformGrid3D: ...

    def __iand__(self, arg: UniformGrid3D, /) -> UniformGrid3D: ...

    def __ior__(self, arg: UniformGrid3D, /) -> UniformGrid3D: ...

    def __isub__(self, arg: UniformGrid3D, /) -> UniformGrid3D: ...

    def __getstate__(self) -> tuple[bytes]: ...

    def __setstate__(self, arg: tuple[bytes], /) -> None: ...

def WriteGridToFile(grid: UniformGrid3D, filename: str) -> None:
    """Write the grid to a grid file"""

def TverskyIndex(grid1: UniformGrid3D, grid2: UniformGrid3D, alpha: float, beta: float) -> float:
    """Compute the tversky index between two grid objects"""

def TanimotoDistance(grid1: UniformGrid3D, grid2: UniformGrid3D) -> float:
    """Compute the tanimoto distance between two grid objects"""

def ProtrudeDistance(grid1: UniformGrid3D, grid2: UniformGrid3D) -> float:
    """Compute the protrude distance between two grid objects"""

def ComputeGridCentroid(grid: UniformGrid3D, pt: Point3D, windowRadius: float) -> tuple[float, Point3D]:
    """Compute the grid point at the center of sphere around a Point3D"""

def FindGridTerminalPoints(grid: UniformGrid3D, windowRadius: float, inclusionFraction: float) -> list[Point3D]:
    """Find a grid's terminal points (defined in the subshape algorithm)."""

class UniformRealValueGrid3D:
    """
    Class to represent a uniform three-dimensional
        cubic grid. Each grid point can store a floating point value.
    """

    @overload
    def __init__(self) -> None:
        """Default constructor"""

    @overload
    def __init__(self, arg: UniformRealValueGrid3D) -> None:
        """Copy constructor"""

    @overload
    def __init__(self, dimX: float, dimY: float, dimZ: float, spacing: float = 0.5, offSet: Point3D | None = None) -> None:
        """Constructor"""

    @overload
    def __init__(self, arg: bytes, /) -> None: ...

    def GetGridPointIndex(self, arg: Point3D, /) -> int:
        """Get the index to the grid point closest to the specified point"""

    def GetGridIndex(self, arg0: int, arg1: int, arg2: int, /) -> int:
        """
        Get the index to the grid point with the three integer indices provided
        """

    def GetGridIndices(self, idx: int) -> tuple[int, int, int]:
        """Returns the integer indices of the grid index provided."""

    def GetValPoint(self, pt: Point3D) -> float:
        """Get the value at the closest grid point"""

    def GetVal(self, id: int) -> float:
        """Get the value at the specified grid index"""

    def SetVal(self, id: int, val: float) -> None:
        """Set the value at the specified grid index"""

    def SetValPoint(self, pt: Point3D, val: float) -> None:
        """Set the value at grid point closest to the specified point"""

    def GetGridPointLoc(self, pointId: int) -> Point3D:
        """Get the location of the specified grid point"""

    def GetSize(self) -> int:
        """Get the size of the grid (number of grid points)"""

    def GetNumX(self) -> int:
        """Get the number of grid points along x-axis"""

    def GetNumY(self) -> int:
        """Get the number of grid points along y-axis"""

    def GetNumZ(self) -> int:
        """Get the number of grid points along z-axis"""

    def GetOffset(self) -> Point3D:
        """Get the location of the center of the grid"""

    def GetSpacing(self) -> float:
        """Get the grid spacing"""

    def GetOccupancyVect(self) -> rdkit.DataStructs.cDataStructs.RealValueVect:
        """Get the occupancy vector for the grid"""

    def CompareVectors(self, arg: UniformRealValueGrid3D, /) -> bool:
        """Compare the vector values between two grid objects."""

    def CompareParams(self, arg: UniformRealValueGrid3D, /) -> bool:
        """Compare the parameters between two grid object."""

    def CompareGrids(self, arg: UniformRealValueGrid3D, /) -> bool:
        """Compare the parameters and values between two grid objects."""

    def __iand__(self, arg: UniformRealValueGrid3D, /) -> UniformRealValueGrid3D: ...

    def __ior__(self, arg: UniformRealValueGrid3D, /) -> UniformRealValueGrid3D: ...

    def __iadd__(self, arg: UniformRealValueGrid3D, /) -> UniformRealValueGrid3D: ...

    def __isub__(self, arg: UniformRealValueGrid3D, /) -> UniformRealValueGrid3D: ...

    def __getstate__(self) -> tuple[bytes]: ...

    def __setstate__(self, arg: tuple[bytes], /) -> None: ...
