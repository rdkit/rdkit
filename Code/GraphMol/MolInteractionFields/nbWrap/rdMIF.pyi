"""
Module containing functions for calculating molecular interaction fields (MIFs)
NOTE: This functionality is experimental and the API and/or results may change in future releases.
"""

import os
from typing import overload

import rdkit.Chem.rdchem
import rdkit.Geometry.rdGeometry


class Coulomb:
    """
    Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in
            vacuum (no dielectric).
    """

    @overload
    def __init__(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1, probeCharge: float = 1.0, absVal: bool = False, chargeKey: str = '_GasteigerCharge', softcoreParam: float = 0.0, cutoff: float = 1.0) -> None:
        """
        Constructor for Coulomb class.

                ARGUMENTS:
                - mol:           the molecule of interest
                - confId:        the ID of the conformer to be used (defaults to -1)
                - probeCharge    charge of probe [e] (defaults to 1.0 e)
                - absVal:        if True, absolute values of interactions are calculated (defaults to False)
                - chargeKey      property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')
                - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
                - cutoff         minimum cutoff distance [A] (defaults to 1.0)
        """

    @overload
    def __init__(self, charges: object, positions: object, probeCharge: float = 1.0, absVal: bool = False, softcoreParam: float = 0.0, cutoff: float = 1.0) -> None:
        """
        Alternative constructor for Coulomb class.

                ARGUMENTS:
                - charges:       array of partial charges of a molecule's atoms
                - positions:     array of positions of a molecule's atoms
                - probeCharge    charge of probe [e] (defaults to 1.0 e)
                - absVal:        if True, absolute values of interactions are calculated (defaults to False)
                - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
                - cutoff         minimum cutoff distance [A] (defaults to 1.0)
        """

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
        Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in
                vacuum (no dielectric).

                ARGUMENTS:
                - x, y, z:   coordinates of probe position for energy calculation
                - threshold: maximal distance until which interactions are calculated
                RETURNS:
                - electrostatic potential in [kJ mol^-1]
        """

class CoulombDielectric:
    """
    Class for calculation of electrostatic interaction (Coulomb energy) between probe and molecule in
            by taking a distance-dependent dielectric into account.
            Same energy term as used in GRID MIFs.
            References:
            - J. Med. Chem. 1985, 28, 849.
            - J. Comp. Chem. 1983, 4, 187.
    """

    @overload
    def __init__(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1, probeCharge: float = 1.0, absVal: bool = False, chargeKey: str = '_GasteigerCharge', softcoreParam: float = 0.0, cutoff: float = 1.0, epsilon: float = 80.0, xi: float = 4.0) -> None:
        """
        Constructor for CoulombDielectric class.

                ARGUMENTS:
                - mol:           the molecule of interest
                - confId:        the ID of the conformer to be used (defaults to -1)
                - probeCharge    charge of probe [e] (defaults to 1.0 e)
                - absVal:        if True, absolute values of interactions are calculated (defaults to False)
                - chargeKey       property key for retrieving partial charges of atoms from molecule (defaults to '_GasteigerCharge')
                - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
                - cutoff         minimum cutoff distance [A] (defaults to 1.0)
                - epsilon        relative permittivity of solvent (defaults to 80.0)
                - xi             relative permittivity of solute (defaults to 4.0)
        """

    @overload
    def __init__(self, charges: object, positions: object, probeCharge: float = 1.0, absVal: bool = False, softcoreParam: float = 0.0, cutoff: float = 1.0, epsilon: float = 80.0, xi: float = 4.0) -> None:
        """
        Alternative constructor for CoulombDielectric class.

              - charges:       array of partial charges of a molecule's atoms
              - positions:     array of positions of a molecule's atoms
              - probeCharge    charge of probe [e] (defaults to 1.0 e)
              - absVal:        if True, absolute values of interactions are calculated (defaults to False)
              - softcoreParam  softcore interaction parameter [A^2], if zero, a minimum cutoff distance is used (defaults to 0.0)
              - cutoff         minimum cutoff distance [A] (defaults to 1.0)
              - epsilon        relative permittivity of solvent (defaults to 80.0)
              - xi             relative permittivity of solute (defaults to 4.0)
        """

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
        Calculates the electrostatic interaction (Coulomb energy) between probe and molecule in
                by taking a distance-dependent dielectric into account.

                ARGUMENTS:
                - x, y, z:   coordinates of probe position for energy calculation
                - threshold: maximal distance until which interactions are calculated
                RETURNS:
                - electrostatic potential in [kJ mol^-1]
        """

class MMFFVdWaals:
    """
    Class for calculating van der Waals interactions between molecule and a probe at a gridpoint        based on the MMFF forcefield.
    """

    def __init__(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1, probeAtomType: int = 6, scaling: bool = False, cutoff: float = 1.0) -> None:
        """
        ARGUMENTS:
                - mol           molecule object
                - confId        conformation id which is used to get positions of atoms (default=-1)
                - probeAtomType MMFF94 atom type for the probe atom (default=6, sp3 oxygen)
                - cutoff        minimum cutoff distance [A] (default:1.0)
                - scaling       scaling of VdW parameters to take hydrogen bonds into account (default=False)
        """

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
        Calculates the van der Waals interaction between molecule and a probe at a gridpoint.

                ARGUMENTS:
                - x, y, z:   coordinates of probe position for energy calculation
                - threshold: maximal distance until which interactions are calculated
                RETURNS:
                - van der Waals potential in [kJ mol^-1]
        """

class UFFVdWaals:
    """
    Class for calculating van der Waals interactions between molecule and a probe at a gridpoint        based on the UFF forcefield.
    """

    def __init__(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1, probeAtomType: str = 'O_3', cutoff: float = 1.0) -> None:
        """
        ARGUMENTS:
                - mol           molecule object
                - confId        conformation id which is used to get positions of atoms (default=-1)
                - probeAtomType UFF atom type for the probe atom (default='O_3', sp3 oxygen)
                - cutoff        minimum cutoff distance [A] (default:1.0)
        """

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
        Calculates the van der Waals interaction between molecule and a probe at a gridpoint.

                ARGUMENTS:
                - x, y, z:   coordinates of probe position for energy calculation
                - threshold: maximal distance until which interactions are calculated
                RETURNS:
                - van der Waals potential in [kJ mol^-1]
        """

class HBond:
    """
    Class for calculation of hydrogen bonding energy between a probe and a molecule.

            Similar to GRID hydrogen bonding descriptors.
            References:
            - J.Med.Chem. 1989, 32, 1083.
            - J.Med.Chem. 1993, 36, 140.
            - J.Med.Chem. 1993, 36, 148.
    """

    def __init__(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1, probeAtomType: str = 'OH', fixed: bool = True, cutoff: float = 1.0) -> None:
        """
        Constructor for HBond class.

                ARGUMENTS:
                - mol:           the molecule of interest
                - confId:        the ID of the conformer to be used (defaults to -1)
                - probeAtomType: atom type for the probe atom (either 'OH', 'O', 'NH' or 'N') (defaults to 'OH')
                - fixed:         for some groups, two different angle dependencies are defined:
                                 one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)
                                 into account and one for strictly fixed conformations
                                 if True, strictly fixed conformations (defaults to True)
                - cutoff         minimum cutoff distance [A] (defaults to 1.0)
        """

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
        Calculates the hydrogen bonding energy between probe and molecule in

                ARGUMENTS:
                - x, y, z:   coordinates of probe position for energy calculation
                - threshold: maximal distance until which interactions are calculated
                RETURNS:
                hydrogen bonding energy in [kJ mol^-1]
        """

class Hydrophilic:
    """
    Class for calculation of a hydrophilic potential of a molecule at a point.

            The interaction energy of hydrogen and oxygen of water is calculated at each point as a
            hydrogen bond interaction (either OH or O probe). The favored interaction is returned.
    """

    def __init__(self, mol: rdkit.Chem.rdchem.Mol, confId: int = -1, fixed: bool = True, cutoff: float = 1.0) -> None:
        """
        Constructor for Hydrophilic class.

                ARGUMENTS:
                - mol:         the molecule of interest
                - confId:      the ID of the conformer to be used (defaults to -1)
                - fixed:       for some groups, two different angle dependencies are defined:
                               one which takes some flexibility of groups (rotation/swapping of lone pairs and hydrogen)
                               into account and one for strictly fixed conformations
                               if True, strictly fixed conformations (defaults to True)
                - cutoff       minimum cutoff distance [A] (default:1.0)
        """

    def __call__(self, x: float, y: float, z: float, threshold: float) -> float:
        """
        Calculates the hydrophilic field energy at a point.

                ARGUMENTS:
                - x, y, z:   coordinates of probe position for energy calculation
                - threshold: maximal distance until which interactions are calculated
                RETURNS:
                hydrophilic field energy in [kJ mol^-1]
        """

def ConstructGrid(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, margin: float = 5.0, spacing: float = 0.5) -> rdkit.Geometry.rdGeometry.UniformRealValueGrid3D:
    """
    Constructs a UniformRealValueGrid3D (3D grid with real values at gridpoints) fitting to a molecule.

            ARGUMENTS:
            - mol:     molecule of interest
            - confId:  the ID of the conformer to be used (defaults to -1)
            - margin:  minimum distance of molecule to surface of grid [A] (defaults to 5.0 A)
            - spacing: grid spacing [A] (defaults to 0.5 A)
    """

@overload
def CalculateDescriptors(grid: rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, descriptor: Coulomb, threshold: float = -1.0) -> None: ...

@overload
def CalculateDescriptors(grid: rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, descriptor: CoulombDielectric, threshold: float = -1.0) -> None: ...

@overload
def CalculateDescriptors(grid: rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, descriptor: MMFFVdWaals, threshold: float = -1.0) -> None: ...

@overload
def CalculateDescriptors(grid: rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, descriptor: UFFVdWaals, threshold: float = -1.0) -> None: ...

@overload
def CalculateDescriptors(grid: rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, descriptor: HBond, threshold: float = -1.0) -> None: ...

@overload
def CalculateDescriptors(grid: rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, descriptor: Hydrophilic, threshold: float = -1.0) -> None:
    """
    Calculates descriptors (to be specified as parameter) of a molecule at every gridpoint of a grid.

            ARGUMENTS:
            - grid:      UniformRealValueGrid3D which get the MIF values
            - descriptor:  Descriptor class which is used to calculate values
    """

def WriteToCubeFile(grid: rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, filename: str | os.PathLike, mol: rdkit.Chem.rdchem.Mol | None = None, confId: int = -1) -> None:
    """
    Writes Grid to a file in Gaussian CUBE format.

            ARGUMENTS:
            - grid:      UniformRealValueGrid3D to be stored
            - filename:  filename of file to be written
            - mol:       associated molecule (defaults to None)
            - confId:    the ID of the conformer to be used (defaults to -1)
    """

def ReadFromCubeFile(filename: str | os.PathLike) -> tuple[rdkit.Geometry.rdGeometry.UniformRealValueGrid3D, rdkit.Chem.rdchem.Mol]:
    """
    Reads Grid from a file in Gaussian CUBE format.

            ARGUMENTS:
            - filename:  filename of file to be read
            RETURNS:
            a tuple where the first element is the grid and
            the second element is the molecule object associated to the grid
            (only atoms and coordinates, no bonds;
            None if no molecule was associated to the grid)
    """
