from rdkit import Chem
from rdkit.Chem import rdRascalMCES

mol1 = Chem.MolFromSmiles("c1ccccc1Cl")
mol2 = Chem.MolFromSmiles("c1ccccc1F")
opts = rdRascalMCES.RascalOptions()

results = rdRascalMCES.FindMCES(mol1, mol2, opts)

for res in results:
    print(res.smartsString, res.similarity, res.numFrags, res.timedOut)
    print(res.bondMatches(), res.atomMatches())
    
ad1 = Chem.MolFromSmiles("CN(C)c1ccc(CC(=O)NCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL153934")
ad2 = Chem.MolFromSmiles("N(C)c1ccc(CC(=O)NCCCCCCCCCCCCNC23CC4CC(C2)CC(C3)C4)cc1 CHEMBL157336")

results = rdRascalMCES.FindMCES(ad1, ad2, opts)

for res in results:
    print(res.smartsString, res.similarity, res.numFrags, res.timedOut)
    print(res.bondMatches(), res.atomMatches())
    res.largestFragmentOnly()
    print(res.smartsString, res.similarity, res.numFrags, res.timedOut)
    print(res.bondMatches(), res.atomMatches())
