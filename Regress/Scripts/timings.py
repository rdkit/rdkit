
import time
import gzip
import random
import os
import sys
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Recap
from rdkit.RDLogger import logger

logger = logger()

tests = [1] * 1001
if len(sys.argv) > 1:
    tests = [0] * 1001
    tests[1] = 1
    for x in sys.argv[1:]:
        x = int(x)
        tests[x] = 1
ts = []

sdData = gzip.open('../Data/mols.1000.sdf.gz').read()
logger.info('mols from sdf')
suppl = Chem.SDMolSupplier()
suppl.SetData(sdData)
mols = []
nMols = 0
nBad = 0
t1 = time.time()
for m in suppl:
    if m:
        nMols += 1
        mols.append(m)
    else:
        nBad += 1
t2 = time.time()
logger.info('Results1: %.2f seconds, %d passed, %d failed' % (t2 - t1, nMols, nBad))
ts.append(t2 - t1)

if tests[2]:
    lines = gzip.open('../Data/mols.1000.txt.gz').readlines()
    logger.info('mols from smiles')
    nMols = 0
    nBad = 0
    t1 = time.time()
    for line in lines:
        line = line.decode().strip().split(' ')
        m = Chem.MolFromSmiles(line[1])
        if m:
            nMols += 1
        else:
            nBad += 1
    t2 = time.time()
    logger.info('Results2: %.2f seconds, %d passed, %d failed' % (t2 - t1, nMols, nBad))
    ts.append(t2 - t1)

if tests[3] or tests[4] or tests[5]:
    pattData = gzip.open('../Data/queries.txt.gz').readlines()
    pattData = [x.decode().strip().replace('[H]', '').replace('()', '') for x in pattData]
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
    for line in open('../Data/RLewis_smarts.txt'):
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
    logger.info('Writing: Canonical SMILES')
    t1 = time.time()
    for mol in mols:
        smi = Chem.MolToSmiles(mol)
    t2 = time.time()
    logger.info('Results9: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[10]:
    logger.info('Generate 2D coords')
    t1 = time.time()
    for mol in mols:
        AllChem.Compute2DCoords(mol)
    t2 = time.time()
    logger.info('Results10: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[11]:
    logger.info('Writing: Mol blocks')
    t1 = time.time()
    for mol in mols:
        mb = Chem.MolToMolBlock(mol)
    t2 = time.time()
    logger.info('Results11: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[12]:
    logger.info('RECAP decomposition')
    t1 = time.time()
    for mol in mols:
        d = Recap.RecapDecompose(mol)
    t2 = time.time()
    logger.info('Results12: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[13]:
    logger.info('Generate 3D coords for 50 molecules with ETKDG')
    mols3d = mols[200:250]
    t1 = time.time()
    nBad = 0
    for mol in mols3d:
        cid = AllChem.EmbedMolecule(mol, randomSeed=0xF00D,
                                    useExpTorsionAnglePrefs=True, useBasicKnowledge=True)
        if cid < 0:
            nBad += 1
    t2 = time.time()
    logger.info('Results13: %.2f seconds %d failures' % (t2 - t1, nBad))
    ts.append(t2 - t1)

if tests[14]:
    logger.info('UFF optimizing those:')
    t1 = time.time()
    for mol in mols3d:
        if not mol.GetNumConformers():
            continue
        mol = Chem.Mol(mol)
        needMore = 1
        while needMore:
            needMore = AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    t2 = time.time()
    logger.info('Results14: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[15]:
    logger.info('MMFF optimizing the molecules:')
    t1 = time.time()
    for i, mol in enumerate(mols3d):
        mol = Chem.Mol(mol)
        if not mol.GetNumConformers():
            continue
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            continue
        needMore = 1
        while needMore:
            try:
                needMore = AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except ValueError:
                logger.warning('Problems with MMFF and mol %d' % i)
                break
    t2 = time.time()
    logger.info('Results15: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[16]:
    logger.info('Find unique subgraphs')
    t1 = time.time()
    for mol in mols:
        Chem.FindUniqueSubgraphsOfLengthN(mol, 6)
    t2 = time.time()
    logger.info('Results16: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[17]:
    logger.info('Generate topological fingerprints')
    t1 = time.time()
    for mol in mols:
        Chem.RDKFingerprint(mol)
    t2 = time.time()
    logger.info('Results17: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

if tests[18]:
    logger.info('Generate morgan fingerprints')
    t1 = time.time()
    for mol in mols:
        AllChem.GetMorganFingerprint(mol, radius=2)
    t2 = time.time()
    logger.info('Results18: %.2f seconds' % (t2 - t1))
    ts.append(t2 - t1)

print(f"| {rdkit.__version__} | {' | '.join(['%.1f' % x for x in ts])} |")
