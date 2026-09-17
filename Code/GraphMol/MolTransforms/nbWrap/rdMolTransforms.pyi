"""
Module containing functions to perform 3D operations like rotate and
translate conformations
"""

from collections.abc import Sequence
from typing import Annotated

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


def ComputeCentroid(conf: rdkit.Chem.rdchem.Conformer, ignoreHs: bool = True, weights: Sequence[float] | None = None) -> rdkit.Geometry.rdGeometry.Point3D:
    """
    Compute the centroid of the conformation - hydrogens are ignored and no attention
    is paid to the difference in sizes of the heavy atoms; however,
    an optional vector of weights can be passed.
    """

def ComputeCanonicalTransform(conf: rdkit.Chem.rdchem.Conformer, center: rdkit.Geometry.rdGeometry.Point3D | None = None, normalizeCovar: bool = False, ignoreHs: bool = True) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
    """
    Compute the transformation required to align a conformer so that
    the principal axes align up with the x,y,z axes.
    The conformer itself is left unchanged.
    ARGUMENTS:
      - conf : the conformer of interest
      - center : optional center point to compute the principal axes around (defaults to the centroid)
      - normalizeCovar : optionally normalize the covariance matrix by the number of atoms
    """

def ComputePrincipalAxesAndMoments(conf: rdkit.Chem.rdchem.Conformer, ignoreHs: bool = True, weights: Sequence[float] | None = None) -> object:
    """
    Compute principal axes and moments of inertia for a conformer.
    These values are calculated from the inertia tensor:
      Iij = - sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
      Iii = sum_{s=1..N} sum_{j!=i} (w_s * r_{sj} * r_{sj})
    where the coordinates are relative to the center of mass.

    ARGUMENTS:
      - conf : the conformer of interest
      - ignoreHs : if True, ignore hydrogen atoms
      - weights : if present, used to weight the atomic coordinates

    Returns a (principal axes, principal moments) tuple
    """

def ComputePrincipalAxesAndMomentsFromGyrationMatrix(conf: rdkit.Chem.rdchem.Conformer, ignoreHs: bool = True, weights: Sequence[float] | None = None) -> object:
    """
    Compute principal axes and moments from the gyration matrix of a conformer.
    These values are calculated from the gyration matrix/tensor:
      Iij = sum_{s=1..N}(w_s * r_{si} * r_{sj}) i != j
      Iii = sum_{s=1..N} sum_{t!=s}(w_s * r_{si} * r_{ti})
    where the coordinates are relative to the center of mass.

    ARGUMENTS:
      - conf : the conformer of interest
      - ignoreHs : if True, ignore hydrogen atoms
      - weights : if present, used to weight the atomic coordinates

    Returns a (principal axes, principal moments) tuple
    """

def TransformConformer(conf: rdkit.Chem.rdchem.Conformer, trans: Annotated[NDArray[numpy.float64], dict(shape=(None, None), order='C', writable=False)]) -> None:
    """Transform the coordinates of a conformer"""

def CanonicalizeConformer(conf: rdkit.Chem.rdchem.Conformer, center: rdkit.Geometry.rdGeometry.Point3D | None = None, normalizeCovar: bool = False, ignoreHs: bool = True) -> None:
    """
    Canonicalize the orientation of a conformer so that its principal axes
    around the specified center point coincide with the x, y, z axes.

    ARGUMENTS:
      - conf : conformer of interest
      - center : optionally center point about which the principal axes are computed;
    if not specified the centroid of the conformer will be used
      - normalizeCovar : Optionally normalize the covariance matrix by the number of atoms
    """

def CanonicalizeMol(mol: rdkit.Chem.rdchem.Mol, normalizeCovar: bool = False, ignoreHs: bool = True) -> None:
    """
    Loop over the conformers in a molecule and canonicalize their orientation
    """

def GetBondLength(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int) -> float:
    """Returns the bond length in angstrom between atoms i, j"""

def SetBondLength(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, value: float) -> None:
    """
    Sets the bond length in angstrom between atoms i, j; all atoms bonded to atom j are moved
    """

def GetAngleRad(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int) -> float:
    """Returns the angle in radians between atoms i, j, k"""

def GetAngleDeg(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int) -> float:
    """Returns the angle in degrees between atoms i, j, k"""

def SetAngleRad(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float) -> None:
    """
    Sets the angle in radians between atoms i, j, k; all atoms bonded to atom k are moved
    """

def SetAngleDeg(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int, value: float) -> None:
    """
    Sets the angle in degrees between atoms i, j, k; all atoms bonded to atom k are moved
    """

def GetDihedralRad(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int) -> float:
    """Returns the dihedral angle in radians between atoms i, j, k, l"""

def GetDihedralDeg(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int) -> float:
    """Returns the dihedral angle in degrees between atoms i, j, k, l"""

def SetDihedralRad(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int, value: float) -> None:
    """
    Sets the dihedral angle in radians between atoms i, j, k, l; all atoms bonded to atom l are moved
    """

def SetDihedralDeg(conf: rdkit.Chem.rdchem.Conformer, iAtomId: int, jAtomId: int, kAtomId: int, lAtomId: int, value: float) -> None:
    """
    Sets the dihedral angle in degrees between atoms i, j, k, l; all atoms bonded to atom l are moved
    """
