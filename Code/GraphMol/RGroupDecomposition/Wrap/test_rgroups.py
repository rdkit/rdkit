#  Copyright (c) 2018, Novartis Institutes for BioMedical Research Inc.
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
#       products derived from this software without specific prior written
#       permission.
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

import unittest
import os,sys, copy

from rdkit.six.moves import cPickle

from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem.rdRGroupDecomposition import RGroupDecompose, RGroupDecomposition, RGroupDecompositionParameters
from collections import OrderedDict

class TestCase(unittest.TestCase) :
    def test_multicores(self):
        cores_smi_easy = OrderedDict()
        cores_smi_hard = OrderedDict()

        #cores_smi_easy['cephem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])CS2')
        cores_smi_easy['cephem'] = Chem.MolFromSmarts('O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CS2')
        cores_smi_hard['cephem'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])CS2')

        #cores_smi_easy['carbacephem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])CC2')
        cores_smi_easy['carbacephem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])CC2')
        cores_smi_hard['carbacephem'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])CC2')

        #cores_smi_easy['oxacephem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])CO2')
        cores_smi_easy['oxacephem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])CO2')
        cores_smi_hard['oxacephem'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])CO2')

        #cores_smi_easy['carbapenem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])C2')
        cores_smi_easy['carbapenem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])C2')
        cores_smi_hard['carbapenem'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])C2')

        #cores_smi_easy['carbapenam'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])C2')
        cores_smi_easy['carbapenam'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)C([3*])([4*])C2')
        cores_smi_hard['carbapenam'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])C2')

        #cores_smi_easy['penem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])S2')
        cores_smi_easy['penem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])S2')
        cores_smi_hard['penem'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])S2')

        #cores_smi_easy['penam'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])S2')
        cores_smi_easy['penam'] = Chem.MolFromSmarts('O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])S2')
        cores_smi_hard['penam'] = Chem.MolFromSmarts('O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])S2')

        #cores_smi_easy['oxapenam'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])O2')
        cores_smi_easy['oxapenam'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)C([3*])([4*])O2')
        cores_smi_hard['oxapenam'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])O2')

        cores_smi_easy['monobactam'] = Chem.MolFromSmarts('O=C1C([1*])C([5*])N1')
        cores_smi_hard['monobactam'] = Chem.MolFromSmarts('O=C1C([2*])([1*])C([6*])([5*])N1')
        rg_easy = RGroupDecomposition(cores_smi_easy.values())
        rg_stereo = RGroupDecomposition(cores_smi_hard.values())

    def test_stereo(self):
        smiles = """C1CCO[C@@H](N)1
C1CCO[C@H](N)1
C1CCO[C@@](N)(O)1
C1CCO[C@@](N)(P)1
C1CCO[C@@](N)(S)1
C1CCO[C@@H](O)1
C1CCO[C@H](O)1
C1CCO[C@@](O)(N)1
C1CCO[C@@](O)(P)1
C1CCO[C@@](O)(S)1
C1CCO[C@@H](P)1
C1CCO[C@H](P)1
C1CCO[C@@](P)(N)1
C1CCO[C@@](P)(O)1
C1CCO[C@@](P)(S)1
C1CCO[C@@H](S)1
C1CCO[C@H](S)1
C1CCO[C@@](S)(N)1
C1CCO[C@@](S)(O)1
C1CCO[C@@](S)(P)1
"""
        mols = []
        for smi in smiles.split():
            m = Chem.MolFromSmiles(smi)
            assert m, smi
            mols.append(m)
        core = Chem.MolFromSmarts("C1CCOC1")
        rgroups = RGroupDecomposition(core)
        for m in mols:
            rgroups.Add(m)
        rgroups.Process()
        columns = rgroups.GetRGroupsAsColumns()
        data = {}
        for k,v in columns.items():
            data[k] = [Chem.MolToSmiles(m,True) for m in v]

        rgroups2,unmatched = RGroupDecompose([core], mols)
        columns2,unmatched = RGroupDecompose([core], mols, asRows=False)
        data2 = {}
        for k,v in columns2.items():
            data2[k] = [Chem.MolToSmiles(m,True) for m in v]

        self.assertEqual(data, data2)
        columns3, unmatched = RGroupDecompose([core], mols, asRows=False, asSmiles=True)
        self.assertEqual(data, columns3)
        

    def test_h_options(self):
        core = Chem.MolFromSmiles("O=c1oc2ccccc2cc1")
        smiles = ("O=c1cc(Cn2ccnc2)c2ccc(Oc3ccccc3)cc2o1",
                  "O=c1oc2ccccc2c(Cn2ccnc2)c1-c1ccccc1",
                  "COc1ccc2c(Cn3cncn3)cc(=O)oc2c1")
        params = RGroupDecompositionParameters()
        rgd = RGroupDecomposition(core,params)
        for smi in smiles:
            m = Chem.MolFromSmiles(smi)
            rgd.Add(m)
        rgd.Process()
        columns = rgd.GetRGroupsAsColumns()
        self.assertEqual(columns['R2'][0].GetNumAtoms(),12)

        params.removeHydrogensPostMatch = True
        rgd = RGroupDecomposition(core,params)
        for smi in smiles:
            m = Chem.MolFromSmiles(smi)
            rgd.Add(m)
        rgd.Process()
        columns = rgd.GetRGroupsAsColumns()
        self.assertEqual(columns['R2'][0].GetNumAtoms(),7)

    def test_unmatched(self):
        cores = [Chem.MolFromSmiles("N")]
        mols = [Chem.MolFromSmiles("CC"),
                Chem.MolFromSmiles("CC"),
                Chem.MolFromSmiles("CC"),
                Chem.MolFromSmiles("N"),
                Chem.MolFromSmiles("CC")]

        res, unmatched = RGroupDecompose(cores, mols)
        self.assertEqual(len(res), 1)
        self.assertEqual(unmatched, [0,1,2,4])

    def test_userlabels(self):
        smis = ["C(Cl)N(N)O(O)"]
        mols = [Chem.MolFromSmiles(smi) for smi in smis]
        smarts = 'C([*:1])N([*:5])O([*:6])'
        core = Chem.MolFromSmarts(smarts)
        rg = RGroupDecomposition(core)
        for m in mols:
            rg.Add(m)
        rg.Process()
        self.assertEqual(rg.GetRGroupsAsColumns(asSmiles=True),
                         {'Core': ['C(N(O[*:6])[*:5])[*:1]'],
                          'R1': ['Cl[*:1]'],
                          'R5': ['[H]N([H])[*:5]'],
                          'R6': ['[H]O[*:6]']})

        smarts = 'C([*:4])N([*:5])O([*:6])'

        core = Chem.MolFromSmarts(smarts)
        rg = RGroupDecomposition(core)
        for m in mols:
            rg.Add(m)
        rg.Process()
        self.assertEqual(rg.GetRGroupsAsColumns(asSmiles=True),
                         {'Core': ['C(N(O[*:6])[*:5])[*:4]'],
                          'R4': ['Cl[*:4]'],
                          'R5': ['[H]N([H])[*:5]'],
                          'R6': ['[H]O[*:6]']})

    def test_match_only_at_rgroups(self):
        smiles = ['c1ccccc1']#, 'c1(Cl)ccccc1', 'c1(Cl)cc(Br)ccc1']
        mols = [Chem.MolFromSmiles(smi) for smi in smiles]

        core1 = Chem.MolFromSmiles("c1([*:5])cc([*:6])ccc1")
        params = RGroupDecompositionParameters()
        params.onlyMatchAtRGroups=True
        rg = RGroupDecomposition(core1, params)
        for smi,m in zip(smiles,mols):
            self.assertTrue(rg.Add(m)!=-1, smi)



if __name__ == '__main__':
  unittest.main()
