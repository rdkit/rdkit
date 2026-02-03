from rdkit import Chem
import sys

print(f"Python version: {sys.version}")
print(f"rdkit Version: {Chem.rdBase.rdkitVersion}")

# Test ring molecule
m_ring = Chem.MolFromSmiles("C1CC1")
print(f"Ring molecule SMILES: '{Chem.MolToSmiles(m_ring)}'")
print(f"Number of atoms: {m_ring.GetNumAtoms()}")
print(f"Number of bonds: {m_ring.GetNumBonds()}")

# Get ring information
ring_info = m_ring.GetRingInfo()
print(f"Ring info: {ring_info}")
print(f"Number of rings: {ring_info.NumRings()}")

# Check each atom's ring bond count
for i, atom in enumerate(m_ring.GetAtoms()):
    ring_bonds = 0
    for bond in atom.GetBonds():
        if bond.IsInRing():
            ring_bonds += 1
    print(f"Atom {i} (atomic num 6): {ring_bonds} ring bonds")