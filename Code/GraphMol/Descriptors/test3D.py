from rdkit import Chem
from rdkit import rdBase

from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem
import time
print rdBase.rdkitVersion
print rdBase.boostVersion


def get3D(m,is3d):
	if not is3d:
		m = Chem.AddHs(m)
		AllChem.EmbedMolecule(m)
		AllChem.MMFFOptimizeMolecule(m)
	r= rdMD.CalcAUTOCORR3D(m)+rdMD.CalcRDF(m)+rdMD.CalcMORSE(m)+rdMD.CalcWHIM(m)+rdMD.CalcGETAWAY(m)
	return r


m = Chem.MolFromSmiles('Cc1ccccc1')
thefile = open('test.txt', 'w')
filename="/Users/mbp/Github/rdkit_mine/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf"
suppl = Chem.SDMolSupplier(filename,removeHs=False)
mols = [x for x in suppl]
start = time.time()
for m in mols:
	r= get3D(m,True)
	for item in r:
  		thefile.write("%.3f," % item)
  	thefile.write("\n")

end = time.time()
print end - start
