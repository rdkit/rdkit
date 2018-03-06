from rdkit.Chem import AllChem


mol = AllChem.MolFromSmiles("CCCC(N)(O)")
info ={}
fp1 = AllChem.GetMorganFingerprint(mol, 3,
          bitInfo=info, includeRedundantEnvironments=True)
fp = fp.GetNonzeroElements()  # convert to a dict
print fp
print info