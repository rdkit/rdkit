# $Id$
#
#  Copyright (c) 2011, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: 
#
#     * Redistributions of source code must retain the above copyright 
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following 
#       disclaimer in the documentation and/or other materials provided 
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc. 
#       nor the names of its contributors may be used to endorse or promote 
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
from __future__ import print_function
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
from rdkit.six.moves.cPickle import load

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
            line = str(line)
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
    def __next__(self):
        return self.next()

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

class RegressionTest(unittest.TestCase):
    def testPrechloricAcid(self):
        examples = (
                ('OCl(=O)(=O)=O', 'InChI=1S/ClHO4/c2-1(3,4)5/h(H,2,3,4,5)'),
                ('CC1=CC2=NCC(CN2C=C1)C(=O)c3ccc4cc(C)ccc4c3.OCl(=O)(=O)=O', 
                    'InChI=1S/C21H20N2O.ClHO4/c1-14-3-4-17-11-18(6-5-16(17)9-14)21(24)19-12-22-20-10-15(2)7-8-23(20)13-19;2-1(3,4)5/h3-11,19H,12-13H2,1-2H3;(H,2,3,4,5)'),
                ('CNc1ccc2nc3ccccc3[n+](C)c2c1.[O-]Cl(=O)(=O)=O',
                    'InChI=1S/C14H13N3.ClHO4/c1-15-10-7-8-12-14(9-10)17(2)13-6-4-3-5-11(13)16-12;2-1(3,4)5/h3-9H,1-2H3;(H,2,3,4,5)'),
                )
        for smiles, expected in examples:
            m = MolFromSmiles(smiles)
            inchi = MolToInchi(m)
            self.assertEqual(inchi, expected)


class TestCase(unittest.TestCase):
    def setUp(self):
        self.dataset = dict()
        self.dataset_inchi = dict()
        inf = gzip.open(os.path.join(RDConfig.RDCodeDir, 'Chem/test_data',
                                     'pubchem-hard-set.sdf.gz'),'r')
        self.dataset['problematic'] = ForwardSDMolSupplier(inf,sanitize=False,removeHs=False)
        with open(os.path.join(RDConfig.RDCodeDir, 'Chem/test_data',
                               'pubchem-hard-set.inchi'),'rb') as inF:
            self.dataset_inchi['problematic'] = load(inF)
        # disable logging
        DisableLog('rdApp.warning')

    def tearDown(self):
        EnableLog('rdApp.warning')

    def test0InchiWritePubChem(self):
        for fp, f in self.dataset.items():
            inchi_db = self.dataset_inchi[fp]
            same, diff, reasonable = 0, 0, 0
            for i_,m in enumerate(f):
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
                        inchi, error = inst.args
                        if 'Metal' in error:
                            reasonable += 1
                            continue

                    diff += 1
                    print('InChI mismatch for PubChem Compound ' + 
                          m.GetProp('PUBCHEM_COMPOUND_CID') + 
                          '\n' + MolToSmiles(m,True) + '\n' + inchiDiff(x, y))
                    print()

                else:
                    same += 1

            print(green + "InChI write Summary: %d identical, %d suffix variance, %d reasonable" % (same, diff, reasonable) + reset)
            self.assertEqual(same, 1164)
            self.assertEqual(diff, 0)
            self.assertEqual(reasonable, 17)
            

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
                        _, error = inst.args
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
                    print(green + 'Empty mol for PubChem Compound ' + 
                          cid + '\n' + reset)
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
                    print(green + 'Molecule mismatch for PubChem Compound ' + 
                          cid + '\n' + reset + 
                          inchiDiff(x, y) + reset)
                    print()
                else:
                    same += 1
            print(green + "InChI Read Summary: %d identical, %d  variance, %d reasonable variance" % (same, diff, reasonable) + reset)
            self.assertEqual(same, 550)
            self.assertEqual(diff, 0)
            self.assertEqual(reasonable, 631)

    def test2InchiOptions(self):
        m = MolFromSmiles("CC=C(N)C")
        inchi1 = MolToInchi(m).split('/', 1)[1]
        inchi2 = MolToInchi(m, "/SUU").split('/', 1)[1]
        self.assertEqual(inchi1 + '/b4-3?', inchi2)

    def test3InchiKey(self):
        inchi = 'InChI=1S/C9H12/c1-2-6-9-7-4-3-5-8-9/h3-5,7-8H,2,6H2,1H3'
        self.assertEqual(InchiToInchiKey(inchi), 'ODLMAHJVESYWTB-UHFFFAOYSA-N')

if __name__ == '__main__':
    # only run the test if InChI is available
    if inchi.INCHI_AVAILABLE:
        unittest.main()

