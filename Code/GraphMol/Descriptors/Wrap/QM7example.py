
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem
import numpy as np
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
def generateX(m):
	vs = rdMD.CalcCoulombMat(m,confId = -1, nbmats= 1, seed = 1, padding=23, local = False, sorted=True, eigenval=False) # true, false
	for v in vs:
		natoms = len(v)
		arr = np.reshape(v,(1,natoms))
	return arr

# read the file and get the data
def loaddata(f):
	X=[]
	Y=[]
	fo = open(f)
	i=0
	for line in fo:
		i+=1
		if i>1:
			rl=line.split(',')
			m = mols(rl[0])
			x= generateX(m)
			X.append(x)
			Y.append(float(rl[1].replace('\n','')))
		if i%100==0:
			print i
		if i>5500:
			break
	return X,Y

if __name__ == "__main__":
	f='/pathto/qm7smiles.csv'
	S,Y=loaddata(f)
	print len(Y)
	print "---------"
	Mat = np.concatenate(S,axis=0)
	sio.savemat('Xoptsortedsautoanitized.mat', {'X':Mat,'Y':Y})


