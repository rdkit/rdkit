"""Module containing functions to handle force fields"""

import rdkit.Chem.rdchem
import rdkit.ForceField.rdForceField


def UFFOptimizeMolecule(mol: rdkit.Chem.rdchem.Mol, maxIters: int = 200, vdwThresh: float = 10.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> int:
    """
    uses UFF to optimize a molecule's structure

    ARGUMENTS:

     - mol : the molecule of interest
     - maxIters : the maximum number of iterations (defaults to 200)
     - vdwThresh : used to exclude long-range van der Waals interactions
                   (defaults to 10.0)
     - confId : indicates which conformer to optimize
     - ignoreInterfragInteractions : if true, nonbonded terms between
                   fragments will not be added to the forcefield.

    RETURNS: 0 if the optimization converged, 1 if more iterations are required.
    """

def UFFOptimizeMoleculeConfs(mol: rdkit.Chem.rdchem.Mol, numThreads: int = 1, maxIters: int = 200, vdwThresh: float = 10.0, ignoreInterfragInteractions: bool = True) -> list[tuple[int, float]]:
    """
    uses UFF to optimize all of a molecule's conformations

    ARGUMENTS:

     - mol : the molecule of interest
     - numThreads : the number of threads to use, only has an effect if the RDKit
                    was built with thread support (defaults to 1)
                    If set to zero, the max supported by the system will be used.
     - maxIters : the maximum number of iterations (defaults to 200)
     - vdwThresh : used to exclude long-range van der Waals interactions
                   (defaults to 10.0)
     - ignoreInterfragInteractions : if true, nonbonded terms between
                   fragments will not be added to the forcefield.

    RETURNS: a list of (not_converged, energy) 2-tuples.
        If not_converged is 0 the optimization converged for that conformer.
    """

def UFFGetMoleculeForceField(mol: rdkit.Chem.rdchem.Mol, vdwThresh: float = 10.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> rdkit.ForceField.rdForceField.ForceField:
    """
    returns a UFF force field for a molecule

    ARGUMENTS:

     - mol : the molecule of interest
     - vdwThresh : used to exclude long-range van der Waals interactions
                   (defaults to 10.0)
     - confId : indicates which conformer to optimize
     - ignoreInterfragInteractions : if true, nonbonded terms between
                   fragments will not be added to the forcefield.
    """

def UFFHasAllMoleculeParams(mol: rdkit.Chem.rdchem.Mol) -> bool:
    """
    checks if UFF parameters are available for all of a molecule's atoms

    ARGUMENTS:

     - mol : the molecule of interest.
    """

def MMFFOptimizeMolecule(mol: rdkit.Chem.rdchem.Mol, mmffVariant: str = 'MMFF94', maxIters: int = 200, nonBondedThresh: float = 100.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> int:
    """
    uses MMFF to optimize a molecule's structure

    ARGUMENTS:

     - mol : the molecule of interest
     - mmffVariant : "MMFF94" or "MMFF94s"
     - maxIters : the maximum number of iterations (defaults to 200)
     - nonBondedThresh : used to exclude long-range non-bonded
                    interactions (defaults to 100.0)
     - confId : indicates which conformer to optimize
     - ignoreInterfragInteractions : if true, nonbonded terms between
                    fragments will not be added to the forcefield

    RETURNS: 0 if the optimization converged, -1 if the forcefield could
             not be set up, 1 if more iterations are required.
    """

def MMFFSanitizeMolecule(mol: rdkit.Chem.rdchem.Mol) -> int:
    """
    sanitizes a molecule according to MMFF requirements.

     - mol : the molecule of interest.
    """

def MMFFGetMoleculeProperties(mol: rdkit.Chem.rdchem.Mol, mmffVariant: str = 'MMFF94', mmffVerbosity: int = 0) -> rdkit.ForceField.rdForceField.MMFFMolProperties:
    """
    returns an MMFFMolProperties object for a
    molecule, which is required by MMFFGetMoleculeForceField()
    and can be used to get/set MMFF properties

    ARGUMENTS:

     - mol : the molecule of interest
     - mmffVariant : "MMFF94" or "MMFF94s"
                   (defaults to "MMFF94")
     - mmffVerbosity : 0: none; 1: low; 2: high (defaults to 0).
    """

def MMFFGetMoleculeForceField(mol: rdkit.Chem.rdchem.Mol, MMFFMolProperties: rdkit.ForceField.rdForceField.MMFFMolProperties | None, nonBondedThresh: float = 100.0, confId: int = -1, ignoreInterfragInteractions: bool = True) -> rdkit.ForceField.rdForceField.ForceField:
    """
    returns a MMFF force field for a molecule

    ARGUMENTS:

     - mol : the molecule of interest
     - MMFFMolProperties : MMFFMolProperties object as returned
                   by MMFFGetMoleculeProperties()
     - nonBondedThresh : used to exclude long-range non-bonded
                   interactions (defaults to 100.0)
     - confId : indicates which conformer to optimize
     - ignoreInterfragInteractions : if true, nonbonded terms between
                   fragments will not be added to the forcefield
    """

def CreateEmptyForceFieldForMol(mol: rdkit.Chem.rdchem.Mol, confId: int = -1) -> rdkit.ForceField.rdForceField.ForceField:
    """
    Get An empty Force Field, with only the positions of the atoms but no Contributions.

    ARGUMENTS :

        - mol : the molecule of interest
        - confId: the conformer which positions should be added to the force field.
    """

def MMFFHasAllMoleculeParams(mol: rdkit.Chem.rdchem.Mol) -> bool:
    """
    checks if MMFF parameters are available for all of a molecule's atoms

    ARGUMENTS:

     - mol : the molecule of interest
    """

def MMFFOptimizeMoleculeConfs(mol: rdkit.Chem.rdchem.Mol, numThreads: int = 1, maxIters: int = 200, mmffVariant: str = 'MMFF94', nonBondedThresh: float = 100.0, ignoreInterfragInteractions: bool = True) -> list[tuple[int, float]]:
    """
    uses MMFF to optimize all of a molecule's conformations

    ARGUMENTS:

     - mol : the molecule of interest
     - numThreads : the number of threads to use, only has an effect if the RDKit
                    was built with thread support (defaults to 1)
                    If set to zero, the max supported by the system will be used.
     - maxIters : the maximum number of iterations (defaults to 200)
     - mmffVariant : "MMFF94" or "MMFF94s"
     - nonBondedThresh : used to exclude long-range non-bonded
                   interactions (defaults to 100.0)
     - ignoreInterfragInteractions : if true, nonbonded terms between
                   fragments will not be added to the forcefield.

    RETURNS: a list of (not_converged, energy) 2-tuples.
        If not_converged is 0 the optimization converged for that conformer.
    """

def GetUFFBondStretchParams(mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int) -> object:
    """
    Retrieves UFF bond stretch parameters for atoms with indexes idx1, idx2 as a (kb, r0) tuple, or None if no parameters could be found
    """

def GetUFFAngleBendParams(mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int, idx3: int) -> object:
    """
    Retrieves UFF angle bend parameters for atoms with indexes idx1, idx2, idx3 as a (ka, theta0) tuple, or None if no parameters could be found
    """

def GetUFFTorsionParams(mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object:
    """
    Retrieves UFF torsion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a V float value, or None if no parameters could be found
    """

def GetUFFInversionParams(mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int, idx3: int, idx4: int) -> object:
    """
    Retrieves UFF inversion parameters for atoms with indexes idx1, idx2, idx3, idx4 as a K float value, or None if no parameters could be found
    """

def GetUFFVdWParams(mol: rdkit.Chem.rdchem.Mol, idx1: int, idx2: int) -> object:
    """
    Retrieves UFF van der Waals parameters for atoms with indexes idx1, idx2 as a (x_ij, D_ij) tuple, or None if no parameters could be found
    """

def OptimizeMolecule(ff: rdkit.ForceField.rdForceField.ForceField, maxIters: int = 200) -> int:
    """
    uses the supplied force field to optimize a molecule's structure

    ARGUMENTS:

     - ff : the force field
     - maxIters : the maximum number of iterations (defaults to 200)

    RETURNS: 0 if the optimization converged, 1 if more iterations are required.
    """

def OptimizeMoleculeConfs(mol: rdkit.Chem.rdchem.Mol, ff: rdkit.ForceField.rdForceField.ForceField, numThreads: int = 1, maxIters: int = 200) -> list[tuple[int, float]]:
    """
    uses the supplied force field to optimize all of a molecule's conformations

    ARGUMENTS:

     - mol : the molecule of interest
     - ff : the force field
     - numThreads : the number of threads to use, only has an effect if the RDKit
                    was built with thread support (defaults to 1)
                    If set to zero, the max supported by the system will be used.
     - maxIters : the maximum number of iterations (defaults to 200)

    RETURNS: a list of (not_converged, energy) 2-tuples.
        If not_converged is 0 the optimization converged for that conformer.
    """
