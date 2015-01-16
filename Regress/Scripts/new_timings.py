from __future__ import print_function
import time,gzip,random,os,sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Recap
from rdkit.RDLogger import logger
logger = logger()

tests=[1]*1001
if len(sys.argv)>1:
    tests=[0]*1001
    tests[1]=1
    for x in sys.argv[1:]:
        x = int(x)
        tests[x] = 1
ts = []
mols = []
lines = gzip.open('../Data/znp.50k.smi.gz','rb').readlines()
logger.info('mols from smiles')
nMols=0
nBad=0
t1=time.time()
for line in lines:
    line = line.strip().split(' ')
    m = Chem.MolFromSmiles(line[0])
    if m:
        nMols+=1
        mols.append(m)
    else:
        nBad += 1
        
t2=time.time()
logger.info('Results1: %.2f seconds, %d passed, %d failed'%(t2-t1,nMols,nBad))
ts.append(t2-t1)

if tests[1]:
    logger.info('Writing: Canonical SMILES')
    t1=time.time()
    for mol in mols:
        smi = Chem.MolToSmiles(mol,True)
    t2 = time.time()
    logger.info('Results2: %.2f seconds'%(t2-t1))
    ts.append(t2-t1)

print('times: ',' || '.join(['%.1f'%x for x in ts]))
