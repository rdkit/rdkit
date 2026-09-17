"""
Module containing interface to the YAeHMOP extended Hueckel library.
Please note that this interface should still be considered experimental and may
change from one release to the next.
"""

from typing import Annotated

import numpy
from numpy.typing import NDArray

import rdkit.Chem.rdchem


class EHTResults:
    @property
    def numOrbitals(self) -> int: ...

    @property
    def numElectrons(self) -> int: ...

    @property
    def fermiEnergy(self) -> float: ...

    @property
    def totalEnergy(self) -> float: ...

    def GetReducedChargeMatrix(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
        """returns the reduced charge matrix"""

    def GetReducedOverlapPopulationMatrix(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
        """returns the reduced overlap population matrix"""

    def GetAtomicCharges(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
        """returns the calculated atomic charges"""

    def GetHamiltonian(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
        """returns the symmetric Hamiltonian matrix"""

    def GetOverlapMatrix(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None, None))]:
        """returns the symmetric overlap matrix"""

    def GetOrbitalEnergies(self) -> Annotated[NDArray[numpy.float64], dict(shape=(None,))]:
        """returns the energies of the molecular orbitals as a vector"""

def RunMol(mol: rdkit.Chem.rdchem.Mol, confId: int = -1, keepOverlapAndHamiltonianMatrices: bool = False) -> tuple:
    """
    Runs an extended Hueckel calculation for a molecule.
    The molecule should have at least one conformation

    ARGUMENTS:
       - mol: molecule to use
       - confId: (optional) conformation to use
       - keepOverlapAndHamiltonianMatrices: (optional) triggers storing the overlap
         and hamiltonian matrices in the EHTResults object

    RETURNS: a 2-tuple:
       - a boolean indicating whether or not the calculation succeeded
       - an EHTResults object with the results
    """
