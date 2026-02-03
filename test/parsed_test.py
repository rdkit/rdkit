from rdkit import Chem
import sys

print(f"Python version: {sys.version}")
print(f"rdkit Version: {Chem.rdBase.rdkitVersion}")

qs = "C!@c"
q = Chem.MolFromSmarts(qs)
print(f"Query SMARTS: '{qs}'")
print(f"Parsed query: {q}")
print(f"Query atoms: {q.GetNumAtoms()}")
print(f"Query bonds: {q.GetNumBonds()}")

# Test with a simple ring
m = Chem.MolFromSmiles("C1CC1")
print(f"Ring molecule: '{Chem.MolToSmiles(m)}'")
print(f"Matches: {m.GetSubstructMatches(q)}")