from rdkit import Chem
from rdkit import rdBase

from rdkit.Chem import rdMolDescriptors as rdMD
from rdkit.Chem import AllChem

print rdBase.rdkitVersion
print rdBase.boostVersion


def get3D(m):
	m = Chem.AddHs(m)
	AllChem.EmbedMolecule(m)
	AllChem.MMFFOptimizeMolecule(m)
	r= rdMD.CalcAUTOCORR3D(m)+rdMD.CalcRDF(m)+rdMD.CalcMORSE(m)+rdMD.CalcWHIM(m)+rdMD.CalcGETAWAY(m)
	return r


m = Chem.MolFromSmiles('Cc1ccccc1')


r= get3D(m)
print len(r)


# from rdkit import Chem
# from rdkit.Chem import Descriptors
# #from rdkit import rdBase
# #print rdBase.rdkitVersion
# #print rdBase.boostVersion

# m = Chem.MolFromSmiles('c1ccccc1C(=O)O')

# #print dir(Descriptors)

# #print Descriptors.CalcRadiusOfGyration(m)


# from rdkit.Chem import rdMolDescriptors

# from rdkit.Chem import AllChem

# m=Chem.AddHs(m)
# AllChem.EmbedMolecule(m) #,useRandomCoords=True)
# AllChem.MMFFOptimizeMolecule(m)

# # print rdMolDescriptors.CalcRadiusOfGyration(m)
# # print rdMolDescriptors.CalcTPSA(m)
# # print rdMolDescriptors.CalcLabuteASA(m)
# # print rdMolDescriptors.CalcAsphericity(m)
# # print rdMolDescriptors.CalcInertialShapeFactor(m)
# # print rdMolDescriptors.CalcNPR1(m)
# # print rdMolDescriptors.CalcNPR2(m)
# # print rdMolDescriptors.CalcPBF(m)
# # print rdMolDescriptors.CalcPMI1(m)
# # print rdMolDescriptors.CalcPMI2(m)
# # print rdMolDescriptors.CalcPMI3(m)
# # print rdMolDescriptors.CalcEccentricity(m)
# # print rdMolDescriptors.CalcSpherocityIndex(m)
# print "------"
# print rdMolDescriptors.CalcRDF(m)
# print "------"

# print rdMolDescriptors.CalcMORSE(m)