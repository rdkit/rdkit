from rdkit import Chem
import sys

print(f"Python version: {sys.version}")
print(f"rdkit Version: {Chem.rdBase.rdkitVersion}")

qs = "C!@c"
q = Chem.MolFromSmarts(qs)
print(f"Query SMARTS: '{qs}'")

# Create a simple non-ring molecule
m = Chem.MolFromSmiles("C1CC1C")
print(f"Non-ring molecule SMILES: '{Chem.MolToSmiles(m)}'")
print(f"Should NOT match: {m.GetSubstructMatches(q)}")

# Create a ring molecule
m_ring = Chem.MolFromSmiles("C1CC1")
print(f"Ring molecule SMILES: '{Chem.MolToSmiles(m_ring)}'")
print(f"Should match: {m_ring.GetSubstructMatches(q)}")