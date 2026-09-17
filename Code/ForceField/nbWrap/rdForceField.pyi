"""Exposes the ForceField class"""

from collections.abc import Sequence
import enum

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


class Snapshot:
    """A snapshot of atomic coordinates from a minimization trajectory"""

    def GetPoint2D(self, pointNum: int) -> rdkit.Geometry.rdGeometry.Point2D:
        """
        Returns the coordinates at pointNum as a Point2D object; requires the Trajectory dimension to be == 2
        """

    def GetPoint3D(self, pointNum: int) -> rdkit.Geometry.rdGeometry.Point3D:
        """
        Returns the coordinates at pointNum as a Point3D object; requires the Trajectory dimension to be >= 2
        """

    def GetEnergy(self) -> float:
        """Returns the energy for this Snapshot"""

    def SetEnergy(self, energy: float) -> None:
        """Sets the energy for this Snapshot"""

class ForceField:
    """A force field"""

    def CalcEnergy(self, pos: Sequence[float] | None = None) -> float:
        """
        Returns the energy (in kcal/mol) of the current arrangement
        or of the supplied coordinate list (if non-empty)
        """

    def CalcGrad(self, pos: Sequence[float] | None = None) -> tuple[float, ...]:
        """
        Returns a tuple filled with the per-coordinate gradients
        of the current arrangement or of the supplied coordinate list (if non-empty)
        """

    def Positions(self) -> tuple[float, ...]:
        """
        Returns a tuple filled with the coordinates of the
        points the ForceField is handling
        """

    def Dimension(self) -> int:
        """Returns the dimension of the ForceField"""

    def NumPoints(self) -> int:
        """Returns the number of points the ForceField is handling"""

    def Minimize(self, maxIts: int = 200, forceTol: float = 0.0001, energyTol: float = 1e-06) -> int:
        """
        Runs some minimization iterations.

          Returns 0 if the minimization succeeded.
        """

    def MinimizeTrajectory(self, snapshotFreq: int, maxIts: int = 200, forceTol: float = 0.0001, energyTol: float = 1e-06) -> tuple[int, list[Snapshot]]:
        """
        Runs some minimization iterations, recording the minimization
        trajectory every snapshotFreq steps.

        Returns a (int, []) tuple; the int is 0 if the minimization succeeded,
        while the list contains Snapshot objects.
        """

    def AddDistanceConstraint(self, idx1: int, idx2: int, minLen: float, maxLen: float, forceConstant: float) -> None:
        """
        Adds a distance constraint to the UFF force field (deprecated, use UFFAddDistanceConstraint instead).
        """

    def AddFixedPoint(self, idx: int) -> None:
        """Adds a fixed point to the force field."""

    def UFFAddDistanceConstraint(self, idx1: int, idx2: int, relative: bool, minLen: float, maxLen: float, forceConstant: float) -> None:
        """
        Adds a distance constraint to the UFF force field; if relative == True, then minLen and maxLen are intended as relative to the current distance.
        """

    def UFFAddAngleConstraint(self, idx1: int, idx2: int, idx3: int, relative: bool, minAngleDeg: float, maxAngleDeg: float, forceConstant: float) -> None:
        """
        Adds an angle constraint to the UFF force field; if relative == True, then minAngleDeg and maxAngleDeg are intended as relative to the current angle.
        """

    def UFFAddTorsionConstraint(self, idx1: int, idx2: int, idx3: int, idx4: int, relative: bool, minDihedralDeg: float, maxDihedralDeg: float, forceConstant: float) -> None:
        """
        Adds a dihedral angle constraint to the UFF force field; if relative == True, then minDihedralDeg and maxDihedralDeg are intended as relative to the current dihedral angle.
        """

    def UFFAddPositionConstraint(self, idx: int, maxDispl: float, forceConstant: float) -> None:
        """Adds a position constraint to the UFF force field."""

    def MMFFAddDistanceConstraint(self, idx1: int, idx2: int, relative: bool, minLen: float, maxLen: float, forceConstant: float) -> None:
        """
        Adds a distance constraint to the MMFF force field; if relative == True, then minLen and maxLen are intended as relative to the current distance.
        """

    def MMFFAddAngleConstraint(self, idx1: int, idx2: int, idx3: int, relative: bool, minAngleDeg: float, maxAngleDeg: float, forceConstant: float) -> None:
        """
        Adds an angle constraint to the MMFF force field; if relative == True, then minAngleDeg and maxAngleDeg are intended as relative to the current angle.
        """

    def MMFFAddTorsionConstraint(self, idx1: int, idx2: int, idx3: int, idx4: int, relative: bool, minDihedralDeg: float, maxDihedralDeg: float, forceConstant: float) -> None:
        """
        Adds a dihedral angle constraint to the MMFF force field; if relative == True, then minDihedralDeg and maxDihedralDeg are intended as relative to the current dihedral angle.
        """

    def MMFFAddPositionConstraint(self, idx: int, maxDispl: float, forceConstant: float) -> None:
        """Adds a position constraint to the MMFF force field."""

    def Initialize(self) -> None:
        """initializes the force field (call this before minimizing)"""

    def AddExtraPoint(self, x: float, y: float, z: float, fixed: bool = True) -> int:
        """Adds an extra point, this can be useful for adding constraints."""

    def GetExtraPointPos(self, idx: int) -> tuple[float, float, float]:
        """returns the location of an extra point as a tuple"""

class MMFFVerbosity(enum.Enum):
    MMFF_VERBOSITY_NONE = 0

    MMFF_VERBOSITY_LOW = 1

    MMFF_VERBOSITY_HIGH = 2

class MMFFMolProperties:
    """MMFF molecular properties"""

    def __init__(self, mol: rdkit.Chem.rdchem.Mol, mmffVariant: str = 'MMFF94', verbosity: int = MMFFVerbosity.MMFF_VERBOSITY_NONE) -> None: ...

    def GetMMFFAtomType(self, idx: int) -> int:
        """Retrieves MMFF atom type for atom with index idx"""

    def GetMMFFFormalCharge(self, idx: int) -> float:
        """Retrieves MMFF formal charge for atom with index idx"""

    def GetMMFFPartialCharge(self, idx: int) -> float:
        """Retrieves MMFF partial charge for atom with index idx"""

    def GetMMFFBondStretchParams(self, mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int) -> object:
        """
        Retrieves MMFF bond stretch parameters for atoms with indexes idx1, idx2 as a (bondType, kb, r0) tuple, or None if no parameters could be found
        """

    def GetMMFFAngleBendParams(self, mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int, idx3: int) -> object:
        """
        Retrieves MMFF angle bend parameters for atoms with indexes idx1, idx2, idx3 as a (angleType, ka, theta0) tuple, or None if no parameters could be found
        """

    def GetMMFFStretchBendParams(self, mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int, idx3: int) -> object:
        """
        Retrieves MMFF stretch-bend parameters for atoms with indexes idx1, idx2, idx3 as a (stretchBendType, kbaIJK, kbaKJI) tuple, or None if no parameters could be found
        """

    def GetMMFFTorsionParams(self, mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object:
        """
        Retrieves MMFF torsion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a (torsionType, V1, V2, V3) tuple, or None if no parameters could be found
        """

    def GetMMFFOopBendParams(self, mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object:
        """
        Retrieves MMFF out-of-plane bending force constant for atoms with indexes idx1, idx2, idx3, idx4 as a koop float value
        """

    def GetMMFFVdWParams(self, idx1: int, idx2: int) -> object:
        """
        Retrieves MMFF van der Waals parameters for atoms with indexes idx1, idx2 as a (R_ij_starUnscaled, epsilonUnscaled, R_ij_star, epsilon) tuple, or None if no parameters could be found
        """

    def SetMMFFDielectricModel(self, dielModel: int = 1) -> None:
        """
        Sets the DielModel MMFF property (1: constant; 2: distance-dependent; defaults to constant)
        """

    def GetMMFFDielectricModel(self) -> int:
        """
        Returns the currently configured MMFF dielectric model (1: constant; 2: distance-dependent).
        """

    def SetMMFFDielectricConstant(self, dielConst: float = 1.0) -> None:
        """Sets the DielConst MMFF property (defaults to 1.0)"""

    def GetMMFFDielectricConstant(self) -> float:
        """Returns the currently configured MMFF dielectric constant."""

    def SetMMFFBondTerm(self, state: bool = True) -> None:
        """
        Sets the bond term to be included in the MMFF equation (defaults to True)
        """

    def GetMMFFBondTerm(self) -> bool:
        """Returns whether the bond term is included in the MMFF equation."""

    def SetMMFFAngleTerm(self, state: bool = True) -> None:
        """
        Sets the angle term to be included in the MMFF equation (defaults to True)
        """

    def GetMMFFAngleTerm(self) -> bool:
        """Returns whether the angle term is included in the MMFF equation."""

    def SetMMFFStretchBendTerm(self, state: bool = True) -> None:
        """
        Sets the stretch-bend term to be included in the MMFF equation (defaults to True)
        """

    def GetMMFFStretchBendTerm(self) -> bool:
        """
        Returns whether the stretch-bend term is included in the MMFF equation.
        """

    def SetMMFFOopTerm(self, state: bool = True) -> None:
        """
        Sets the out-of-plane bend term to be included in the MMFF equation (defaults to True)
        """

    def GetMMFFOopTerm(self) -> bool:
        """
        Returns whether the out-of-plane bend term is included in the MMFF equation.
        """

    def SetMMFFTorsionTerm(self, state: bool = True) -> None:
        """
        Sets the torsional term to be included in the MMFF equation (defaults to True)
        """

    def GetMMFFTorsionTerm(self) -> bool:
        """Returns whether the torsional term is included in the MMFF equation."""

    def SetMMFFVdWTerm(self, state: bool = True) -> None:
        """
        Sets the Van der Waals term to be included in the MMFF equation (defaults to True)
        """

    def GetMMFFVdWTerm(self) -> bool:
        """
        Returns whether the Van der Waals term is included in the MMFF equation.
        """

    def SetMMFFEleTerm(self, state: bool = True) -> None:
        """
        Sets the electrostatic term to be included in the MMFF equation (defaults to True)
        """

    def GetMMFFEleTerm(self) -> bool:
        """
        Returns whether the electrostatic term is included in the MMFF equation.
        """

    def SetMMFFVariant(self, mmffVariant: str = 'MMFF94') -> None:
        """
        Sets the MMFF variant to be used ("MMFF94" or "MMFF94s"; defaults to "MMFF94")
        """

    def GetMMFFVariant(self) -> str:
        """Returns the currently configured MMFF variant ("MMFF94" or "MMFF94s")."""

    def SetMMFFVerbosity(self, verbosity: int = 0) -> None:
        """Sets the MMFF verbosity (0: none; 1: low; 2: high; defaults to 0)"""
