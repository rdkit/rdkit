"""
Module containing a C++ implementation of the xyz2mol algorithm. This is based on xyz2mol: https://github.com/jensengroup/xyz2mol
"""

import rdkit.Chem.rdchem


class MaxFindBondOrdersItersExceeded(RuntimeError):
    pass

def DetermineConnectivity(mol: rdkit.Chem.rdchem.Mol, useHueckel: bool = False, charge: int = 0, covFactor: float = 1.3, useVdw: bool = False) -> None:
    """
    Assigns atomic connectivity to a molecule using atomic coordinates,
    disregarding pre-existing bonds

    Args:
       mol : the molecule of interest; it must have a 3D conformer
       useHueckel : (optional) if this is  \\c true, extended Hueckel theory
           will be used to determine connectivity rather than the van der Waals
           or connect-the-dots methods
       charge : (optional) the charge of the molecule; it must be provided if
           the Hueckel method is used and charge is non-zero
       covFactor : (optional) the factor with which to multiply each covalent
           radius if the van der Waals method is used
       useVdw: (optional) if this is false, the connect-the-dots method
           will be used instead of the van der Waals method
    """

def DetermineBondOrders(mol: rdkit.Chem.rdchem.Mol, charge: int = 0, allowChargedFragments: bool = True, embedChiral: bool = True, useAtomMap: bool = False, maxIterations: int = 0) -> None:
    """
    Assigns atomic connectivity to a molecule using atomic coordinates,
    disregarding pre-existing bonds

    Args:
       mol : the molecule of interest; it must have a 3D conformer
       charge : (optional) the charge of the molecule; it must be provided if
           the Hueckel method is used and charge is non-zero
       allowChargedFragments : (optional) if this is true, formal charges
           will be placed on atoms according to their valency; otherwise, radical
           electrons will be placed on the atoms
       embedChiral : (optional) if this is true,
          chirality information will be embedded into the molecule; the function calls
          sanitizeMol() when this is true
       useAtomMap : (optional) if this is true, an atom map will be created for the
          molecule
       maxIterations: (optional) maximum number of iterations to run in the bond order
       determination algorithm, after which a MaxFindBondOrdersItersExceeded
       exception will be thrown. Defaults to 0 (no limit)
    """

def DetermineBonds(mol: rdkit.Chem.rdchem.Mol, useHueckel: bool = False, charge: int = 0, covFactor: float = 1.3, allowChargedFragments: bool = True, embedChiral: bool = True, useAtomMap: bool = False, useVdw: bool = False, maxIterations: int = 0) -> None:
    """
    Assigns atomic connectivity to a molecule using atomic coordinates,
    disregarding pre-existing bonds

    Args:
       mol : the molecule of interest; it must have a 3D conformer
       useHueckel : (optional) if this is true, extended Hueckel theory
           will be used to determine connectivity rather than the van der Waals
           or connect-the-dots methods
       charge : (optional) the charge of the molecule; it must be provided if
           the Hueckel method is used and charge is non-zero
       covFactor : (optional) the factor with which to multiply each covalent
           radius if the van der Waals method is used
       allowChargedFragments : (optional) if this is true, formal charges
           will be placed on atoms according to their valency; otherwise, radical
           electrons will be placed on the atoms
       embedChiral : (optional) if this is true,
          chirality information will be embedded into the molecule; the function calls
          sanitizeMol() when this is true
       useAtomMap : (optional) if this is true, an atom map will be created for the
          molecule
       useVdw: (optional) if this is false, the connect-the-dots method
           will be used instead of the van der Waals method
       maxIterations: (optional) maximum number of iterations to run in the bond order
       determination algorithm, after which a MaxFindBondOrdersItersExceeded
       exception will be thrown. Defaults to 0 (no limit)
    """

def hueckelEnabled() -> bool:
    """whether or not the RDKit was compiled with YAeHMOP support"""
