"""
Module containing functions to encode and compare the shapes of molecules
"""

from typing import Annotated

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs
import rdkit.Geometry.rdGeometry


def EncodeShape(mol: rdkit.Chem.rdchem.Mol, grid: rdkit.Geometry.rdGeometry.UniformGrid3D, confId: int = -1, trans: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)] | None = None, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> None:
    """
    Encode the shape of a molecule (one of its conformer) onto a grid

    ARGUMENTS:

       - mol : the molecule of interest
       - grid : grid onto which the encoding is written
       - confId : id of the conformation of interest on mol (defaults to the first one)
       - trans : any transformation that needs to be used to encode onto the grid (note the molecule remains unchanged)
       - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
                    used in the encoding - grid points inside this sphere carry the maximum occupancy
       - setpSize : thickness of the layers outside the base radius, the occupancy value is decreased
                    from layer to layer from the maximum value
       - maxLayers : the maximum number of layers - defaults to the number of bits
                     used per grid point - e.g. two bits per grid point will allow 3 layers
       - ignoreHs : when set, the contribution of Hs to the shape will be ignored
    """

def ShapeTverskyIndex(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol, alpha: float, beta: float, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: rdkit.DataStructs.cDataStructs.DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> float:
    """
    Compute the shape tversky index between two molecule based on a predefined alignment

    ARGUMENTS:
       - mol1 : The first molecule of interest
       - mol2 : The second molecule of interest
       - alpha : first parameter of the Tversky index
       - beta : second parameter of the Tversky index
       - confId1 : Conformer in the first molecule (defaults to first conformer)
       - confId2 : Conformer in the second molecule (defaults to first conformer)
       - gridSpacing : resolution of the grid used to encode the molecular shapes
       - bitsPerPoint : number of bits used to encode the occupancy at each grid point
                             defaults to two bits per grid point
       - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
                   used in the encoding - grid points inside this sphere carry the maximum occupancy
       - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
                    from layer to layer from the maximum value
       - maxLayers : the maximum number of layers - defaults to the number of bits
                     used per grid point - e.g. two bits per grid point will allow 3 layers
       - ignoreHs : when set, the contribution of Hs to the shape will be ignored
    """

def ShapeTanimotoDist(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: rdkit.DataStructs.cDataStructs.DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True) -> float:
    """
    Compute the shape tanimoto distance between two molecule based on a predefined alignment

    ARGUMENTS:
       - mol1 : The first molecule of interest
       - mol2 : The second molecule of interest
       - confId1 : Conformer in the first molecule (defaults to first conformer)
       - confId2 : Conformer in the second molecule (defaults to first conformer)
       - gridSpacing : resolution of the grid used to encode the molecular shapes
       - bitsPerPoint : number of bits used to encode the occupancy at each grid point
                             defaults to two bits per grid point
       - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
                   used in the encoding - grid points inside this sphere carry the maximum occupancy
       - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
                    from layer to layer from the maximum value
       - maxLayers : the maximum number of layers - defaults to the number of bits
                     used per grid point - e.g. two bits per grid point will allow 3 layers
       - ignoreHs : when set, the contribution of Hs to the shape will be ignored
    """

def ShapeProtrudeDist(mol1: rdkit.Chem.rdchem.Mol, mol2: rdkit.Chem.rdchem.Mol, confId1: int = -1, confId2: int = -1, gridSpacing: float = 0.5, bitsPerPoint: rdkit.DataStructs.cDataStructs.DiscreteValueType = rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE, vdwScale: float = 0.8, stepSize: float = 0.25, maxLayers: int = -1, ignoreHs: bool = True, allowReordering: bool = True) -> float:
    """
    Compute the shape protrude distance between two molecule based on a predefined alignment

    ARGUMENTS:
       - mol1 : The first molecule of interest
       - mol2 : The second molecule of interest
       - confId1 : Conformer in the first molecule (defaults to first conformer)
       - confId2 : Conformer in the second molecule (defaults to first conformer)
       - gridSpacing : resolution of the grid used to encode the molecular shapes
       - bitsPerPoint : number of bit used to encode the occupancy at each grid point
                             defaults to two bits per grid point
       - vdwScale : Scaling factor for the radius of the atoms to determine the base radius
                   used in the encoding - grid points inside this sphere carry the maximum occupancy
       - stepSize : thickness of the each layer outside the base radius, the occupancy value is decreased
                    from layer to layer from the maximum value
       - maxLayers : the maximum number of layers - defaults to the number of bits
                     used per grid point - e.g. two bits per grid point will allow 3 layers
       - ignoreHs : when set, the contribution of Hs to the shape will be ignored
       - allowReordering : when set, the order will be automatically updated so that the value calculated
                           is the protrusion of the smaller shape from the larger one.
    """

def ComputeConfDimsAndOffset(conf: rdkit.Chem.rdchem.Conformer, trans: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)] | None = None, padding: float = 2.0) -> tuple:
    """
    Compute the size of the box that can fit the conformations, and offset
    of the box from the origin
    """

def ComputeConfBox(conf: rdkit.Chem.rdchem.Conformer, trans: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)] | None = None, padding: float = 2.0) -> tuple:
    """
    Compute the lower and upper corners of a cuboid that will fit the conformer
    """

def ComputeUnionBox(box1: tuple[rdkit.Geometry.rdGeometry.Point3D, rdkit.Geometry.rdGeometry.Point3D], box2: tuple[rdkit.Geometry.rdGeometry.Point3D, rdkit.Geometry.rdGeometry.Point3D]) -> tuple:
    """
    Compute the union of two boxes, so that all the points in both boxes are
    contained in the new box
    """
