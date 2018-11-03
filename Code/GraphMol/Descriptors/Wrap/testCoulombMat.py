from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem
import numpy as np
import time 

startTime = time.time()

def isqrt(x):
    if x < 0:
        raise ValueError('square root not defined for negative numbers')
    n = int(x)
    if n == 0:
        return 0
    a, b = divmod(n.bit_length(), 2)
    x = 2**(a+b)
    while True:
        y = (x + n//x)//2
        if y >= x:
            return x
        x = y

m = Chem.MolFromSmiles('NCCCCCO')
m = Chem.AddHs(m)
ps = AllChem.ETKDG()
ps.randomSeed = 0xf00d
AllChem.EmbedMolecule(m,ps)
#vs = rdMD.CalcEEMcharges(m)
#print vs
print ('Initialization script took {0} second !'.format(time.time() - startTime))

startTime = time.time()


#CalcCoulombMat(mol, confId, nbmats, seed, rcut, local = false, decaying = false, reduced = false)
for i in range(1):
	# randomized sorted global matrix
	vs = rdMD.CalcCoulombMat(m,confId = -1, nbmats= 10, seed = 0xf00d, local = False)
	for v in vs:
		#natoms = isqrt(len(v))
		natoms = len(v)
		arr = np.reshape(v,(1,natoms))
		#print(arr)
#print ('The mat script took {0} second !'.format(time.time() - startTime))

startTime = time.time()

for i in range(1):
	# local matrix 
	vs = rdMD.CalcCoulombMat(m,-1,10,1, rcut = 2.0, local = True)
	for v in vs:
		#natoms = isqrt(len(v))
		natoms = len(v)
		arr = np.reshape(v,(1,natoms))
		print("------------------")
		#print(arr)
#print ('The mat script took {0} second !'.format(time.time() - startTime))

startTime = time.time()

for i in range(1):
	vs = rdMD.CalcCoulombMat(m,-1,10,1, rcut = 2.0, local = True, decaying =  True)
	for v in vs:
		#natoms = isqrt(len(v))
		natoms = len(v)
		arr = np.reshape(v,(1,natoms))
		print(arr)
print ('The mat script took {0} second !'.format(time.time() - startTime))

startTime = time.time()

for i in range(1):
	vs = rdMD.CalcCoulombMat(m,-1,10,1, rcut = 6.5, local = True, decaying =  True, reduced = False, alpha = 6)
    # check if the limit 10 is realy apply yes 
	for v in vs:
		#natoms = isqrt(len(v))
		natoms = len(v)
		arr = np.reshape(v,(1,natoms))
		#print(arr)
#print ('The mat script took {0} second !'.format(time.time() - startTime))


for i in range(1):
	vs = rdMD.CalcCoulombMat(m,-1,10,1, rcut = 6.5, local = True, decaying =  True, reduced = True, alpha = 6) 
	# check if alpha is also apply to diag values... strange 
	for v in vs:
		natoms = len(v)
		arr = np.reshape(v,(1,natoms))
		#print(arr)
#print ('The mat script took {0} second !'.format(time.time() - startTime))