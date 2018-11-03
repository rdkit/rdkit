from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem
import numpy as np
import time 
import scipy.io as sio

# convert smile to mol
def mols(smile):
	m = Chem.MolFromSmiles(smile) #,sanitize=False)
	#Chem.SanitizeMol(m,sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL^Chem.SanitizeFlags.SANITIZE_KEKULIZE^Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
	ps = AllChem.ETKDG()
	ps.randomSeed = 0xf00d
	AllChem.EmbedMolecule(m,ps)
	AllChem.MMFFOptimizeMolecule(m)
	return m

# generate from rdkit coulomb matrix from smiles
def generateX(m, smi, MaxBags):
    r = rdMD.CalcBagOfBondVector(m, alpha = 2, MaxBags = MaxBags)
    return r

def getmols(smi):
    Mols = []
    for smile in smi:
        Mols.append(mols(smile))
    return Mols

# read the file and get the data
def getXmols(mols, MaxBags,smi):
    # get BoB maaping
    X=[]
    for m in mols:
            x = generateX(m,smi,MaxBags)
            X.append(x)
    return X

# read the file and get the data only first 10 values!
def loadsmi(f):
    smi=[]
    fo = open(f)
    i=0
    for line in fo:
        i+=1
        if i>50001:
            break
        if i>1:
            rl=line.split(',')
            smi.append(rl[0])
    return smi

if __name__ == "__main__":
    f='/Users/tgg/Github/rdkit/Code/GraphMol/Descriptors/test_data/qm7smiles.csv'
    startTime = time.time()
    smi=loadsmi(f)
    print(len(smi))
    print ('Smiles reading took {0} second !'.format(time.time() - startTime))
    startTime = time.time()
    mols = getmols(smi)
    print ('Smiles to 3d mols script took {0} second !'.format(time.time() - startTime))
    startTime = time.time()
    MaxBags = rdMD.CalcBagOfBondsMap(smi)
    print(MaxBags)
    print ('smiles to BoB map script took {0} second !'.format(time.time() - startTime))
    startTime = time.time()
    X =getXmols(mols, MaxBags,smi)
    print ('mols to BoB matrix script took {0} second !'.format(time.time() - startTime))
    print(len(X))
    print(X[0:2])
    print ("---------")
    Mat = np.concatenate(X,axis=0)
    sio.savemat('XBoB.mat', {'X': Mat})
