from rdkit import Chem, DataStructs, RDConfig, __version__, rdBase

import time
import os

ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Wrap', 'test_data',
                         'chiral_test_mols.tsv')
# with open('./test_data/chiral_test_mols/chiral_test_mols.tsv', 'r') as f:
with open(ofile, 'r') as f:
  lines = f.readlines()
lines = [x.strip().split('\t') for x in lines]
lines.pop(0)

taccum = 0
cx_taccum = 0
res = []
cx_res = []
failed_res = []
n_fail = 0

smiList = ""
cxsmiList = ""
failedList = ""
for line in lines:
  print("line:" +str(line))
  cid, smiles = line
  mol = Chem.MolFromSmiles(smiles)
  t0 = time.time()
  smi = Chem.MolToSmiles(mol)
  t1 = time.time()
  taccum += t1 - t0
  res.append((cid, smiles, smi))
  smiList += cid + ':' + smi + '\n'
  if 0:
    t0 = time.time()
    cxsmi = Chem.MolToCXSmiles(mol)
    t1 = time.time()
  else:
    t0 = time.time()
    try:
      mol = Chem.CanonicalizeStereoGroups(mol)
      cxsmi = Chem.MolToCXSmiles(mol)
    # except ValueError:
    except:
      cxsmi = ''
      n_fail += 1
      failed_res.append((cid, smiles))
      failedList += cid + ':'  + smiles + '\n'
    t1 = time.time()

  cx_taccum += t1 - t0
  cx_res.append((cid, smiles, smi))
  cxsmiList += cid + ':' + cxsmi + '\n'

print('Time for MolToSmiles:', taccum)
print('Time for MolToCXSmiles:', cx_taccum)
print('Number of failures:', n_fail)
with open('./run.pkl', 'wb') as f:
  import pickle
  pickle.dump((taccum, res, cx_taccum, cx_res), f)

with open('./chiral_test_mols_smiOut.txt', 'w') as f:
  f.write(smiList)
  
with open('./chiral_test_mols_cxsmiOut.txt', 'w') as f:
  f.write(cxsmiList)
  
with open('./chiral_test_mols_failedOut.txt', 'w') as f:
  f.write(failedList)
  