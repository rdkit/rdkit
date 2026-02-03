from rdkit import Chem
import sys

print(f"Python version: {sys.version}")
print(f"rdkit Version: {Chem.rdBase.rdkitVersion}")

qs = "C!@c"
q = Chem.MolFromSmarts(qs)
print(f"Query SMARTS: '{qs}'")

# Test molecule that should match
m1 = Chem.MolFromSmiles("C1CC1C")
print(f"SMILES: '{Chem.MolToSmiles(m1)}', matches: {m1.GetSubstructMatches(q)}")

# Test molecule from issue (simplified)
m2 = Chem.MolFromMolBlock("""C1C2CC3CCC(C3)c3ccccc3C12

 0 0 0 0 0 999 V3000
M V30 BEGIN CTAB
M V30 COUNTS 16 20 0 0 0
M V30 BEGIN ATOM
M V30 1 C 14.0076 -9.07 0 0
M V30 2 C 14.6193 -8.0101 0 0
M V30 3 C 15.7208 -7.2455 0 0
M V30 4 C 15.7476 -6.0691 0 0
M V30 5 C 15.1326 -5.1066 0 0
M V30 6 C 13.3948 -5.0531 0 0
M V30 7 C 12.5296 -5.9194 0 0
M V30 8 C 11.4687 -5.3077 0 0
M V30 9 C 10.4089 -5.9194 0 0
M V30 10 C 10.4089 -7.1439 0 0
M V30 11 C 11.4687 -7.7567 0 0
M V30 12 C 12.5296 -7.1439 0 0
M V30 13 C 13.3948 -8.0101 0 0
M V30 14 C 13.4215 -6.176 0 0
M V30 15 C 14.3573 -7.1385 0 0
M V30 16 C 15.1326 -6.283 0 0
M V30 END ATOM
M V30 BEGIN BOND
M V30 1 1 1 2
M V30 2 1 2 3
M V30 3 1 3 4
M V30 4 1 4 5
M V30 5 1 5 6
M V30 6 1 6 7
M V30 7 1 7 8
M V30 8 2 8 9
M V30 9 1 9 10
M V30 10 2 10 11
M V30 11 1 11 12
M V30 12 2 7 12
M V30 13 1 12 13
M V30 14 1 2 13
M V30 15 1 1 13
M V30 16 1 6 14
M V30 17 1 14 15
M V30 18 1 3 15
M V30 19 1 15 16
M V30 20 1 5 16
M V30 END BOND
M V30 END CTAB
M END
""")
print(f"SMILES: '{Chem.MolToSmiles(m2)}', matches: {m2.GetSubstructMatches(q)}")