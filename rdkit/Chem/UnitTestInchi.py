# $Id$

# Copyright (C) 2011-2011 Novartis Institutes for BioMedical Research, Inc
#
#  @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

from rdkit.Chem import *
from rdkit import RDConfig
from rdkit.Chem import rdDepictor
import unittest
import gzip
import tempfile
import os
import gzip
from subprocess import Popen, PIPE
import re
from rdkit.rdBase import DisableLog, EnableLog
from pickle import load

curdir = os.path.dirname(os.path.abspath(__file__))

esc = chr(27)
red = esc + '[31m' 
green = esc + '[32m' 
reset = esc + '[0m'

class NoReentrySDMolSupplier(object):
    def __init__(self, filepath, sanitize=False):
        self.ms = SDMolSupplier()
        self.valid = True
        self.filepath = filepath
        self.f = self.fileopen()
        self.sanitize = sanitize

    def fileopen(self):
        return file(self.filepath, 'r')

    def fileclose(self):
        if self.f:
            self.f.close()

    def reset(self):
        self.fileclose()
        self.f = gzip.open(self.filepath, 'r')
        self.valid = True
        
    def next(self):
        buf = ''
        for line in self.f:
            buf += line
            if line.strip() == '$$$$':
                break
        if buf:
            self.sdf = buf
            self.ms.SetData(buf, sanitize=self.sanitize, removeHs=False)
            return self.ms[0]
        else:
            self.fileclose()
            self.f = None
            raise StopIteration

    def dumpCurMol(self):
        f = file(self.ms[0].GetProp('PUBCHEM_COMPOUND_CID') + '.sdf', 'w')
        f.write(self.sdf)
        f.close()
        return f.name

    def __iter__(self):
        if not self.valid:
            raise ValueError("No reentry allowed")
        self.valid = False
        return self

class GzippedSDMolSupplier(NoReentrySDMolSupplier):
    def fileopen(self):
        return gzip.open(self.filepath, 'r')


def inchiDiffPrefix(inchi1, inchi2):
    inchi1 = inchi1.split('/')
    inchi2 = inchi2.split('/')
    for i in range(len(inchi1) + 1):
        if i == len(inchi1):
            break
        if i == len(inchi2) or inchi1[i] != inchi2[i]:
            break
    if len(inchi1) >= i:
        return inchi1[i][0]
    else:
        return inchi2[i][0]

def inchiDiff(inchi1, inchi2):
    inchi1 = inchi1.split('/')
    inchi2 = inchi2.split('/')
    for i in range(len(inchi1) + 1):
        if i == len(inchi1):
            break
        if i == len(inchi2) or inchi1[i] != inchi2[i]:
            break
    return '/'.join(inchi1[:i]) + red + '/' + '/'.join(inchi1[i:]) + \
        '\n' + reset + \
        '/'.join(inchi2[:i]) + red + '/' + '/'.join(inchi2[i:]) + \
        reset 

class TestCase(unittest.TestCase):
    def setUp(self):
        self.dataset = dict()
        self.dataset_inchi = dict()
        self.dataset['problematic'] = GzippedSDMolSupplier(
                os.path.join(RDConfig.RDCodeDir, 'Chem/test_data',
                    'pubchem-hard-set.sdf.gz'))
        _ = file(os.path.join(RDConfig.RDCodeDir, 'Chem/test_data',
            'pubchem-hard-set.inchi'))
        self.dataset_inchi['problematic'] = load(_)
        _.close()
        # disable logging
        DisableLog('rdApp.warning')

    def tearDown(self):
        EnableLog('rdApp.warning')

    def test0InchiWritePubChem(self):
        for fp, f in self.dataset.items():
            inchi_db = self.dataset_inchi[fp]
            same, diff, reasonable = 0, 0, 0
            for m in f:
                if m is None:
                    continue
                ref_inchi = inchi_db[m.GetProp('PUBCHEM_COMPOUND_CID')]
                x, y = MolToInchi(m), ref_inchi
                if x != y:
                    if re.search(r'.[1-9]?ClO4', x) is not None:
                        reasonable += 1
                        continue
                    SanitizeMol(m)
                    if filter(lambda i: i >= 8, [len(r) for r in m.GetRingInfo().AtomRings()]):
                        reasonable += 1
                        continue
                    # if it is because RDKit does not think the bond is stereo
                    z = MolToInchi(MolFromMolBlock(MolToMolBlock(m)))
                    if y != z and inchiDiffPrefix(y, z) == 'b':
                        reasonable += 1
                        continue
                    # some warning
                    try:
                        MolToInchi(m, treatWarningAsError=True)
                    except InchiReadWriteError as inst:
                        inchi, error = inst
                        if 'Metal' in error:
                            reasonable += 1
                            continue

                    diff += 1
                    print 'InChI mismatach for PubChem Compound ' + \
                            m.GetProp('PUBCHEM_COMPOUND_CID') + \
                            '\n' + MolToSmiles(m) + '\n' + inchiDiff(x, y)
                    print

                else:
                    same += 1

            print green + "InChI write Summary: %d identical, %d suffix variance, %d reasonable" % (same, diff, reasonable) + reset
            self.assertEqual(same, 1174)
            self.assertEqual(diff, 0)
            self.assertEqual(reasonable, 7)
            

    def test1InchiReadPubChem(self):
        for fp, f in self.dataset.items():
            inchi_db = self.dataset_inchi[fp]
            same, diff, reasonable = 0, 0, 0
            for m in f:
                if m is None:
                    continue
                ref_inchi = inchi_db[m.GetProp('PUBCHEM_COMPOUND_CID')]
                cid = m.GetProp('PUBCHEM_COMPOUND_CID')
                x = MolToInchi(m)
                try:
                    y = MolToInchi(
                            MolFromSmiles(
                                MolToSmiles(
                                    MolFromInchi(x), isomericSmiles=True
                                    )
                                )
                            )
                except:
                    y = ''
                if y == '':
                    # metal involved?
                    try:
                        MolToInchi(m, treatWarningAsError=True)
                    except InchiReadWriteError as inst:
                        _, error = inst
                        if 'Metal' in error or \
                                'Charges were rearranged' in error:
                            reasonable += 1
                            continue
                    # RDKit does not like the SMILES? use MolBlock instead
                    inchiMol = MolFromInchi(x)
                    if inchiMol:
                        rdDepictor.Compute2DCoords(inchiMol)
                        z = MolToInchi(MolFromMolBlock(MolToMolBlock(inchiMol)))
                        if x == z:
                            reasonable += 1
                            continue
                    # InChI messed up the radical?
                    unsanitizedInchiMol = MolFromInchi(x, sanitize=False)
                    if sum([a.GetNumRadicalElectrons() * a.GetAtomicNum()
                            for a in m.GetAtoms()
                            if a.GetNumRadicalElectrons() != 0]
                            ) != sum([a.GetNumRadicalElectrons() * a.GetAtomicNum()
                            for a in unsanitizedInchiMol.GetAtoms()
                            if a.GetNumRadicalElectrons() != 0]):
                        reasonable += 1
                        continue

                    diff += 1
                    print green + 'Empty mol for PubChem Compound ' + \
                        cid + '\n' + reset
                    continue
                if x != y:
                    # if there was warning in the first place, then this is
                    # tolerable
                    try:
                        MolToInchi(m, treatWarningAsError=True)
                        MolFromInchi(x, treatWarningAsError=True)
                    except InchiReadWriteError as inst:
                        reasonable += 1
                        continue
                    # or if there are big rings
                    SanitizeMol(m)
                    if filter(lambda i: i >= 8, [len(r) for r in m.GetRingInfo().AtomRings()]):
                        reasonable += 1
                        continue
                    # or if RDKit loses bond stereo
                    s = MolToSmiles(m, True)
                    if MolToSmiles(MolFromSmiles(s), True) != s:
                        reasonable += 1
                        continue
                    # or if it is RDKit SMILES writer unhappy about the mol
                    inchiMol = MolFromInchi(x)
                    rdDepictor.Compute2DCoords(inchiMol)
                    z = MolToInchi(MolFromMolBlock(MolToMolBlock(inchiMol)))
                    if x == z:
                        reasonable += 1
                        continue

                    diff += 1
                    print green + 'Molecule mismatch for PubChem Compound ' + \
                            cid + '\n' + reset + \
                            inchiDiff(x, y) + reset
                    print
                else:
                    same += 1
            print green + "InChI Read Summary: %d identical, %d  variance, %d reasonable variance" % (same, diff, reasonable) + reset
            self.assertEqual(same, 428)
            self.assertEqual(diff, 1)
            self.assertEqual(reasonable, 752)

    def test2InchiOptions(self):
        m = MolFromSmiles("CC=C(N)C")
        inchi1 = MolToInchi(m).split('/', 1)[1]
        inchi2 = MolToInchi(m, "/SUU").split('/', 1)[1]
        self.assertEqual(inchi1 + '/b4-3?', inchi2)

    def test3InchiKey(self):
        inchi = 'InChI=1S/C9H12/c1-2-6-9-7-4-3-5-8-9/h3-5,7-8H,2,6H2,1H3'
        self.assertEqual(InchiToInchiKey(inchi), 'ODLMAHJVESYWTB-UHFFFAOYSA-N')

if __name__ == '__main__':
    unittest.main()

