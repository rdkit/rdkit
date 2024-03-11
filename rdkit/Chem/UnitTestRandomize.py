import pytest
from contextlib import contextmanager
from rdkit import Chem
from rdkit.Chem import Randomize
from rdkit import RDRandom as random


def _get_bond_indices(mol):
    return [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]


@contextmanager
def set_random_seed():
    random.seed(42)
    yield


@set_random_seed()
@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)[O-]",
        "[NH3+]CC(=O)[O-]",
        "[O-][O-]",
        "[O-]C[O-]",
    ],
)
def test_RandomizeMol(smiles):
    for _ in range(5):
        mol = Chem.MolFromSmiles(smiles)
        mol_randomized = Randomize.RandomizeMol(mol)

        assert _get_bond_indices(mol) != _get_bond_indices(mol_randomized)


@set_random_seed()
@pytest.mark.parametrize(
    "smiles",
    ["CON", "c1ccccn1", "C/C=C/F"],
)
def test_smiles_canonicalization(smiles):
    mol = Chem.MolFromSmiles(smiles)
    canonical_reference_smiles = Chem.MolToSmiles(mol, False)
    for _ in range(5):
        mol_randomized = Randomize.RandomizeMol(mol)
        canonical_smiles = Chem.MolToSmiles(mol_randomized, False)

        assert canonical_reference_smiles == canonical_smiles
