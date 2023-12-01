import gzip
import pickle
import sys
import time

from rdkit import Chem, RDConfig
from rdkit.RDLogger import logger

logger = logger()

logger.info('reading smarts')
qs = []
smas = []
for line in file(RDConfig.RDDataDir + '/SmartsLib/RLewis_smarts.txt', 'r').readlines():
  if line[0] == '#':
    continue
  line = line.split(' ')
  p = Chem.MolFromSmarts(line[0])
  if not p:
    print(line[0], file=sys.stderr)
    continue
  smas.append(line[0])
  qs.append(p)

logger.info('reading target counts')
refFps = pickle.loads(gzip.open('fps.1000.counts.pkl.gz', 'rb').read())

fps = []
logger.info('reading mols:')
ms = pickle.loads(gzip.open('mols.1000.pkl.gz', 'rb').read())
t1 = time.time()
nFail = 0
for i, m in enumerate(ms):
  fp = [0] * len(qs)
  for j, q in enumerate(qs):
    o = m.GetSubstructMatches(q)
    if len(o) != refFps[i][j]:
      print('  >', i, j, o, refFps[i][j], Chem.MolToSmiles(m), smas[j])
      nFail += 1
      if nFail == 10:
        raise ValueError
    fp[j] = len(o)
  fps.append(fp)
  if not i % 50:
    logger.info('Done %d' % i)
t2 = time.time()
print('%.2f' % (t2 - t1))

#pickle.dump(fps,file('fps.1000.counts.pkl','wb+'))
nFail = 0
for i, fp in enumerate(fps):
  if fp != refFps[i]:
    nFail += 1
print('%d mismatches' % nFail)
