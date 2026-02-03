from rdkit import Chem
import sys

print(f"Python version: {sys.version}")
print(f"rdkit Version: {Chem.rdBase.rdkitVersion}")

# Test molecule that should match
m = Chem.MolFromSmiles("C1CC1C")
print(f"Non-ring molecule SMILES: '{Chem.MolToSmiles(m)}'")
print(f"Number of atoms: {m.GetNumAtoms()}")
print(f"Number of bonds: {m.GetNumBonds()}")

# Get ring information
ring_info = m.GetRingInfo()
print(f"Ring info: {ring_info}")
print(f"Number of rings: {ring_info.NumRings()}")

# Check each atom's ring bond count
for i, atom in enumerate(m.GetAtoms()):
    ring_bonds = 0
    for bond in atom.GetBonds():
        if bond.IsInRing():
            ring_bonds += 1
    print(f"Atom {i} (atomic num {atom.GetAtomicNum()}): {ring_bonds} ring bonds")