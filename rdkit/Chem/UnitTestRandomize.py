import unittest
import random
from rdkit import Chem
from rdkit.Chem import Randomize


def _get_bond_indices(mol):
    return [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) for bond in mol.GetBonds()]


class TestCase(unittest.TestCase):
    random.seed(42)

    def test_RandomizeMol(self):
        smiles = [
            "CC(=O)[O-]",
            "[NH3+]CC(=O)[O-]",
            "[O-][O-]",
            "[O-]C[O-]",
        ]
        for _ in range(5):
            for smi in smiles:
                mol = Chem.MolFromSmiles(smi)
                mol_randomized = Randomize.RandomizeMol(mol)

                self.assertNotEqual(
                    _get_bond_indices(mol),
                    _get_bond_indices(mol_randomized),
                    msg=f"Failed to randomize charged mol {smi}",
                )

    def test_smiles_canonicalization(self):
        smiles = ["CON", "c1ccccn1", "C/C=C/F"]
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            canonical_reference_smiles = Chem.MolToSmiles(mol, False)
            for _ in range(5):
                mol_randomized = Randomize.RandomizeMol(mol)
                canonical_smiles = Chem.MolToSmiles(mol_randomized, False)
                self.assertEqual(
                    canonical_reference_smiles,
                    canonical_smiles,
                    msg=f"Canonicalization of {smi} resulted in {canonical_smiles} and {canonical_reference_smiles}",
                )


if __name__ == "__main__":  # pragma: nocover
    unittest.main()
