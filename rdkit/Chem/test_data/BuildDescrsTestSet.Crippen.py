
from rdkit import RDConfig
import os.path
import pickle
from rdkit import Chem
from rdkit.Chem import Descriptors

descrs = ['SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6',
          'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11',
          'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7',
          'SlogP_VSA8', 'SlogP_VSA9']


def runIt(inFileName, outFileName, smiCol=0, maxMols=-1, delim=','):
  outF = open(outFileName, 'w+')
  outF.write('#' + ','.join(['SMILES'] + descrs))
  outF.write('\n')
  mols = []
  nDone = 0
  for line in inD:
    if line[0] != '#':
      splitL = line.strip().split(delim)
      if not splitL:
        continue
      smi = splitL[smiCol].strip()
      mol = Chem.MolFromSmiles(smi)
      print(smi)
      if mol:
        vals = []
        for descr in descrs:
          fn = getattr(Descriptors, descr)
          try:
            v = fn(mol)
          except Exception:
            v = 666
          vals.append(v)
        outF.write(','.join([smi] + ['%.4f' % x for x in vals]))
        outF.write('\n')
      nDone += 1
      if maxMols > 0 and nDone >= maxMols:
        break
  outF.close()


if __name__ == '__main__':
  if 1:
    inD = file(
      os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'PP_descrs_regress.VSA.csv'),
      'r').readlines()
    outFileName = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'PP_descrs_regress.VSA.csv')
    runIt(inD, outFileName, smiCol=0, delim=',')
