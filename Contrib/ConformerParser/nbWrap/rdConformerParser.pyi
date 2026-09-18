"""
Module containing functions to read conformations of a molecule from MD trajectories
"""

import rdkit.Chem.rdchem


def AddConformersFromAmberTrajectory(mol: rdkit.Chem.rdchem.Mol, traj: str, numConfs: int = -1, clearConfs: bool = True) -> list[int]:
    """
    Read conformations of a molecule from
    an Amber trajectory

    ARGUMENTS:

       - mol : the molecule of interest
       - traj : the filename of the trajectory
       - numConfs : number of conformations to read
                   The default (-1) reads all.
       - clearConfs : clear all existing conformations on the molecule
                      The default is true.

    RETURNS:

       IDs of the new conformations added to the molecule
    """
