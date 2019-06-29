#
#  Copyright (C) 2019 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import time
import gzip
import rdkit
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.RDLogger import logger
import argparse

logger = logger()

parser = argparse.ArgumentParser(
    description="benchmark and test fingerprint screenout and substructure searching")
parser.add_argument("--validate", dest='validateResults', default=False, action='store_true',
                    help="validate that the screenout isn't missing anything")
parser.add_argument("--short", dest='doShort', default=False, action='store_true',
                    help="run a small subset of the molecules")
args = parser.parse_args()

ts = []

logger.info('mols from smiles')
mols = []
t1 = time.time()
# find this file here: https://raw.githubusercontent.com/greglandrum/rdkit_blog/master/data/chembl21_25K.pairs.txt.gz
with gzip.open('../Data/chembl21_25K.pairs.txt.gz', 'rb') as inf:
    for line in inf:
        line = line.decode().strip().split()
        smi1 = line[1]
        smi2 = line[3]
        mols.append(Chem.MolFromSmiles(smi1))
        mols.append(Chem.MolFromSmiles(smi2))
        if args.doShort and len(mols) >= 1000:
            break
t2 = time.time()
ts.append(t2 - t1)
logger.info(f'Results{len(ts)}: {t2-t1 : .2f} seconds, {len(mols)} mols')


logger.info('queries from smiles')
t1 = time.time()
# find this file here: https://raw.githubusercontent.com/greglandrum/rdkit_blog/master/data/zinc.frags.500.q.smi
frags = [Chem.MolFromSmiles(x.split()[0]) for x in open('../Data/zinc.frags.500.q.smi', 'r')]
# find this file here: https://raw.githubusercontent.com/greglandrum/rdkit_blog/master/data/zinc.leads.500.q.smi
leads = [Chem.MolFromSmiles(x.split()[0]) for x in open('../Data/zinc.leads.500.q.smi', 'r')]
# find this file here: https://raw.githubusercontent.com/greglandrum/rdkit_blog/master/data/fragqueries.q.txt
pieces = [Chem.MolFromSmiles(x) for x in open('../Data/fragqueries.q.txt', 'r')]
t2 = time.time()
ts.append(t2 - t1)
logger.info(f'Results{len(ts)}: {t2-t1 : .2f} seconds')


logger.info('generating pattern fingerprints for mols')
t1 = time.time()
mfps = [Chem.PatternFingerprint(m) for m in mols]
t2 = time.time()
ts.append(t2 - t1)
logger.info(f'Results{len(ts)}: {t2-t1 : .2f} seconds')


logger.info('generating pattern fingerprints for queries')
t1 = time.time()
fragsfps = [Chem.PatternFingerprint(m, 2048) for m in frags]
leadsfps = [Chem.PatternFingerprint(m, 2048) for m in leads]
piecesfps = [Chem.PatternFingerprint(m, 2048) for m in pieces]
t2 = time.time()
ts.append(t2 - t1)
logger.info(f'Results{len(ts)}: {t2-t1 : .2f} seconds')

for nm, qs, qfps in [('frags', frags, fragsfps), ('leads', leads, leadsfps), ('pieces', pieces, piecesfps)]:
    logger.info(f'testing {nm} queries')
    t1 = time.time()
    nPossible = 0
    nTested = 0
    nFound = 0
    nErrors = 0
    for i, fragfp in enumerate(qfps):
        for j, mfp in enumerate(mfps):
            nPossible += 1
            if args.validateResults:
                matched = mols[j].HasSubstructMatch(qs[i])
                fpMatch = DataStructs.AllProbeBitsMatch(fragfp, mfp)
                if fpMatch:
                    nTested += 1
                if matched:
                    nFound += 1
                    if not fpMatch:
                        nErrors += 1
                        logger.error(f"ERROR: mol {j} query {i}")
            else:
                if DataStructs.AllProbeBitsMatch(fragfp, mfp):
                    nTested += 1
                    if mols[j].HasSubstructMatch(qs[i]):
                        nFound += 1
    t2 = time.time()
    ts.append(t2 - t1)
    logger.info(
        f'Results{len(ts)}: {t2-t1 : .2f} seconds. {nTested} tested ({nTested/nPossible :.4f} of total), {nFound} found, {nFound/nTested : .2f} accuracy. {nErrors} errors.')


print(f"| {rdkit.__version__} | {' | '.join(['%.1f' % x for x in ts])} |")
