"""Module containing functions to generate hashes for molecules"""

import enum

import rdkit.Chem.rdchem


class HashFunction(enum.Enum):
    AnonymousGraph = 1

    ElementGraph = 2

    CanonicalSmiles = 3

    MurckoScaffold = 4

    ExtendedMurcko = 5

    MolFormula = 6

    AtomBondCounts = 7

    DegreeVector = 8

    Mesomer = 9

    HetAtomTautomer = 10

    HetAtomProtomer = 11

    RedoxPair = 12

    Regioisomer = 13

    NetCharge = 14

    SmallWorldIndexBR = 15

    SmallWorldIndexBRL = 16

    ArthorSubstructureOrder = 17

    HetAtomTautomerv2 = 18

    HetAtomProtomerv2 = 19

AnonymousGraph: HashFunction = HashFunction.AnonymousGraph

ElementGraph: HashFunction = HashFunction.ElementGraph

CanonicalSmiles: HashFunction = HashFunction.CanonicalSmiles

MurckoScaffold: HashFunction = HashFunction.MurckoScaffold

ExtendedMurcko: HashFunction = HashFunction.ExtendedMurcko

MolFormula: HashFunction = HashFunction.MolFormula

AtomBondCounts: HashFunction = HashFunction.AtomBondCounts

DegreeVector: HashFunction = HashFunction.DegreeVector

Mesomer: HashFunction = HashFunction.Mesomer

HetAtomTautomer: HashFunction = HashFunction.HetAtomTautomer

HetAtomProtomer: HashFunction = HashFunction.HetAtomProtomer

RedoxPair: HashFunction = HashFunction.RedoxPair

Regioisomer: HashFunction = HashFunction.Regioisomer

NetCharge: HashFunction = HashFunction.NetCharge

SmallWorldIndexBR: HashFunction = HashFunction.SmallWorldIndexBR

SmallWorldIndexBRL: HashFunction = HashFunction.SmallWorldIndexBRL

ArthorSubstructureOrder: HashFunction = HashFunction.ArthorSubstructureOrder

HetAtomTautomerv2: HashFunction = HashFunction.HetAtomTautomerv2

HetAtomProtomerv2: HashFunction = HashFunction.HetAtomProtomerv2

def MolHash(mol: rdkit.Chem.rdchem.Mol, func: HashFunction, useCxSmiles: bool = False, cxFlagsToSkip: int = 0) -> str:
    """
    Generate a hash for a molecule. The func argument determines which hash is generated.
    """
