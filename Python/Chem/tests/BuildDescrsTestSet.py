import RDConfig
import os.path,cPickle
import Chem
from Chem import AvailDescriptors

descrs=['BertzCT','Chi0','HallKierAlpha','Ipc','Kappa3','LabuteASA','SMR_VSA1',
        'SMR_VSA10','SMR_VSA2','SMR_VSA3','SMR_VSA4','SMR_VSA5','SMR_VSA6','SMR_VSA7',
        'SMR_VSA8','SMR_VSA9','SlogP_VSA1','SlogP_VSA10','SlogP_VSA11','SlogP_VSA12',
        'SlogP_VSA2','SlogP_VSA3','SlogP_VSA5','SlogP_VSA6','SlogP_VSA7','SlogP_VSA8',
        'SlogP_VSA9','TPSA','NumHAcceptors','NumHDonors','NumHeteroatoms',
        'NumRotatableBonds','MolLogP']

def runIt(inFileName,outFileName,smiCol=0,maxMols=-1,delim=','):
  inF = open(inFileName,'r')
  outF = open(outFileName,'w+')
  outF.write(','.join(['SMILES']+descrs))
  outF.write('\n')
  mols = []
  nDone = 0
  for line in inF.readlines():
    if line[0] != '#':
      splitL = line.strip().split(delim)
      smi = splitL[smiCol].strip()
      print smi
      mol = Chem.MolFromSmiles(smi)
      if mol:
        vals = []
        for descr in descrs:
          fn = AvailDescriptors.descDict[descr]
          try:
            v = fn(mol)
          except:
            v = 666
          vals.append(v)  
        outF.write(','.join([smi]+['%.4f'%x for x in vals]))
        outF.write('\n')
      nDone += 1
      if maxMols>0 and nDone>=maxMols:
        break
  outF.close()    
    
if __name__ == '__main__':
  if 1:
    # 100 PhysProp compounds
    inFileName = os.path.join(RDConfig.RDDataDir,'PhysProp',
                              'Limits-pyridine100.smi')
    outFileName = os.path.join(RDConfig.RDCodeDir,'Chem','tests',
                               'PP_descrs_regress.csv')
    runIt(inFileName,outFileName,smiCol=1,delim='\t')

    
