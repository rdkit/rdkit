from collections.abc import Iterable, Sequence

import rdkit.Chem.rdchem
import rdkit.DataStructs.cDataStructs


class MHFPEncoder:
    def __init__(self, n_permutations: int = 2048, seed: int = 42) -> None: ...

    def FromStringArray(self, vec: Iterable[str]) -> list[int]:
        """Creates a MHFP vector from a list of arbitrary strings."""

    def FromArray(self, vec: Iterable[int]) -> list[int]:
        """Creates a MHFP vector from a list of unsigned integers."""

    def CreateShinglingFromSmiles(self, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1) -> list[str]:
        """
        Creates a shingling (a list of circular n-grams / substructures) from a SMILES string.
        """

    def CreateShinglingFromMol(self, mol: rdkit.Chem.rdchem.Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1) -> list[str]:
        """
        Creates a shingling (a list of circular n-grams / substructures) from a RDKit Mol instance.
        """

    def EncodeSmiles(self, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1) -> list[int]:
        """Creates a MHFP vector from a SMILES string."""

    def EncodeMol(self, mol: rdkit.Chem.rdchem.Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1) -> list[int]:
        """Creates a MHFP vector from an RDKit Mol instance."""

    def EncodeSmilesBulk(self, smiles: Iterable[str], radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1) -> tuple:
        """Creates a MHFP vector from a list of SMILES strings."""

    def EncodeMolsBulk(self, mols: Iterable[rdkit.Chem.rdchem.Mol], radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1) -> tuple:
        """Creates a MHFP vector from a list of RDKit Mol instances."""

    def EncodeSECFPSmiles(self, smiles: str, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1, length: int = 2048) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
        """Creates a SECFP binary vector from a SMILES string."""

    def EncodeSECFPMol(self, mol: rdkit.Chem.rdchem.Mol, radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1, length: int = 2048) -> rdkit.DataStructs.cDataStructs.ExplicitBitVect:
        """Creates a SECFP binary vector from an RDKit Mol instance."""

    def EncodeSECFPSmilesBulk(self, smiles: Iterable[str], radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1, length: int = 2048) -> tuple:
        """Creates a SECFP binary vector from a list of SMILES strings."""

    def EncodeSECFPMolsBulk(self, mols: Iterable[rdkit.Chem.rdchem.Mol], radius: int = 3, rings: bool = True, isomeric: bool = False, kekulize: bool = True, min_radius: int = 1, length: int = 2048) -> tuple:
        """Creates a SECFP binary vector from a list of RDKit Mol instances."""

    @staticmethod
    def Distance(a: Sequence[int], b: Sequence[int]) -> float: ...
