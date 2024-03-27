import unittest
import random
import textwrap
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

    def test_RandomizeMolBlock(self):
        molblock = textwrap.dedent(
            """
                 RDKit          2D

              6  6  0  0  0  0  0  0  0  0999 V2000
                1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
               -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
                0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
              1  2  2  0
              2  3  1  0
              3  4  2  0
              4  5  1  0
              5  6  2  0
              6  1  1  0
            M  END
            """
        )
        molblock_randomized = Randomize.RandomizeMolBlock(molblock)
        mol = Chem.MolFromMolBlock(molblock)
        mol_randomized = Chem.MolFromMolBlock(molblock_randomized)

        self.assertNotEqual(
            _get_bond_indices(mol),
            _get_bond_indices(mol_randomized),
            msg="Failed to randomize molblock",
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
