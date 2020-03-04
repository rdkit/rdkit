
import time
import gzip
import random
import os
import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.RDLogger import logger

dname = os.path.dirname(__file__)


def data(fname):
    return os.path.join(dname, '..', 'Data', fname)


logger = logger()

tests = [1] * 1001
if len(sys.argv) > 1:
    tests = [0] * 1001
    for x in sys.argv[1:]:
        x = int(x)
        tests[x] = 1
ts = []
mols = []

if tests[0]:
    lines = gzip.open(data('znp.50k.smi.gz'), 'rt').readlines()
    logger.info('mols from smiles')
    nMols = 0
    nBad = 0
    t1 = time.time()
    for line in lines:
        line = line.strip().split(' ')
        m = Chem.MolFromSmiles(line[0])
        if m:
            nMols += 1
            mols.append(m)
        else:
            nBad += 1

    t2 = time.time()
    logger.info('Results1: %.2f seconds, %d passed, %d failed' % (t2 - t1, nMols, nBad))
    ts.append(t2 - t1)

if tests[1]:
    logger.info('Writing: Canonical SMILES')
    t1 = time.time()
    for mol in mols:
        smi = Chem.MolToSmiles(mol, True)
    t2 = time.time()
    logger.info('Results2: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[2]:
    sdData = gzip.open(data('mols.1000.sdf.gz'), 'rb').read()
    logger.info('mols from sdf')
    suppl = Chem.SDMolSupplier()
    suppl.SetData(sdData)
    nMols = 0
    nBad = 0
    t1 = time.time()
    for i in range(10):
        for m in suppl:
            if m:
                nMols += 1
                # mols.append(m)
            else:
                nBad += 1
    t2 = time.time()
    logger.info('Results1: %.2f seconds, %d passed, %d failed' % (t2 - t1, nMols, nBad))
    ts.append(t2 - t1)

if tests[3] or tests[4] or tests[5]:
    pattData = gzip.open(data('queries.txt.gz'), 'rt').readlines()
    pattData = [x.strip().replace('[H]', '').replace('()', '') for x in pattData]
    logger.info('patterns from smiles')
    patts = []
    nMols = 0
    t1 = time.time()
    for line in pattData:
        m = Chem.MolFromSmarts(line)
        if m:
            nMols += 1
            patts.append(m)
        else:
            nBad += 1
    t2 = time.time()
    logger.info('Results3: %.2f seconds, %d passed, %d failed' % (t2 - t1, nMols, nBad))
    ts.append(t2 - t1)
    random.seed(23)
    random.shuffle(patts)
    patts = patts[:100]

if tests[4]:
    logger.info('Matching1: HasSubstructMatch')
    t1 = time.time()
    for mol in mols:
        for patt in patts:
            mol.HasSubstructMatch(patt)
    t2 = time.time()
    logger.info('Results4: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[5]:
    logger.info('Matching2: GetSubstructMatches')
    t1 = time.time()
    for mol in mols:
        for patt in patts:
            mol.GetSubstructMatches(patt)
    t2 = time.time()
    logger.info('Results5: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[6] or tests[7] or tests[8]:
    logger.info('reading SMARTS')
    patts = []
    t1 = time.time()
    for line in open(data('RLewis_smarts.txt')):
        line = line.strip()
        if line == '' or line[0] == '#':
            continue
        splitL = line.split(' ')
        sma = splitL[0]
        m = Chem.MolFromSmarts(sma)
        if m:
            patts.append(m)
    t2 = time.time()
    logger.info('Results6: %.2f seconds for %d patterns' % (t2 - t1, len(patts)))
    ts.append(t2 - t1)

if tests[7]:
    logger.info('Matching3: HasSubstructMatch')
    t1 = time.time()
    for mol in mols:
        for patt in patts:
            mol.HasSubstructMatch(patt)
    t2 = time.time()
    logger.info('Results7: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[8]:
    logger.info('Matching4: GetSubstructMatches')
    t1 = time.time()
    for mol in mols:
        for patt in patts:
            mol.GetSubstructMatches(patt)
    t2 = time.time()
    logger.info('Results8: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[9]:
    logger.info('Writing: Mol blocks')
    t1 = time.time()
    for mol in mols:
        mb = Chem.MolToMolBlock(mol)
    t2 = time.time()
    logger.info('Results10: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[10]:
    from rdkit.Chem import BRICS
    logger.info('BRICS decomposition')
    t1 = time.time()
    for mol in mols:
        d = BRICS.BreakBRICSBonds(mol)
    t2 = time.time()
    logger.info('Results11: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[11]:
    logger.info('Generate 2D coords')
    t1 = time.time()
    for mol in mols:
        AllChem.Compute2DCoords(mol)
    t2 = time.time()
    logger.info('Results12: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[12]:
    logger.info('Generate topological fingerprints')
    t1 = time.time()
    for mol in mols:
        Chem.RDKFingerprint(mol)
    t2 = time.time()
    logger.info('Results16: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[13]:
    logger.info('Generate morgan fingerprints')
    t1 = time.time()
    for mol in mols:
        AllChem.GetMorganFingerprint(mol, radius=2)
    t2 = time.time()
    logger.info('Results16: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

print(f"| {rdkit.__version__} | {' | '.join(['%.1f' % x for x in ts])} |")
