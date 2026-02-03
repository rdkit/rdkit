from rdkit import Chem
import sys

print(f"Python version: {sys.version}")
print(f"rdkit Version: {Chem.rdBase.rdkitVersion}")

qs = "C!@c"
q = Chem.MolFromSmarts(qs)
print(f"Query SMARTS: '{qs}'")

# Test molecule that should match
m = Chem.MolFromSmiles("C1CC1C")
print(f"SMILES: '{Chem.MolToSmiles(m)}', matches: {m.GetSubstructMatches(q)}")