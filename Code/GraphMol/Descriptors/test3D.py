import os
import time

from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors as rdMD


def get3D(m, is3d):
  if not is3d:
    m = Chem.AddHs(m)
    # define the new code from RDKit Molecule 3D ETKDG.
    ps = AllChem.ETKDG()
    ps.randomSeed = 0xf00d
    AllChem.EmbedMolecule(m, ps)
  r = rdMD.CalcAUTOCORR3D(m) + rdMD.CalcRDF(m) + rdMD.CalcMORSE(m) + rdMD.CalcWHIM(
    m) + rdMD.CalcGETAWAY(m, precision=0.001)
  return r


def generateAll():
  filename = '/Users/GVALMTGG/Github/rdkit_mine/Code/GraphMol/Descriptors/test_data/PBF_egfr.sdf'
  suppl = Chem.SDMolSupplier(filename, removeHs=False)
  mols = [x for x in suppl]
  start = time.time()
  for m in mols:
    r = get3D(m, True)
    print(r)
  end = time.time()
  print(end - start)


def simple_case():
  start = time.time()
  smi = 'CCC(C)COCCCC'
  m = Chem.MolFromSmiles(smi)
  T = get3D(m, False)
  print(T)
  end = time.time()
  print(end - start)


if (__name__ == '__main__'):
  # FIX: We need to actually add some tests here, but this doees not need to
  # to be done until the C++ code and tests are straightened out.
  generateAll()

  start = time.time()
  smi = 'CCC(C)COCCCC'
  m = Chem.MolFromSmiles(smi)
  T = get3D(m, False)
  print(T)
  end = time.time()
  print(end - start)
