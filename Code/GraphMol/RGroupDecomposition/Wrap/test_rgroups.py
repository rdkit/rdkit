#  Copyright (c) 2017, Novartis Institutes for BioMedical Research Inc.
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
from rdkit.Chem.rdRGroupDecomposition import RGroupDecomposition
from collections import OrderedDict

class TestCase(unittest.TestCase) :
    def atest_multicores(self):
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
            print (m, Chem.MolToSmiles(m))
            rgroups.Add(m)
        rgroups.Process()
        columns = rgroups.GetRGroupsAsColumns()
        for k,v in columns.items():
            print("%s:%s"%(k, " ".join([Chem.MolToSmiles(m,True) for m in v])))
            
if __name__ == '__main__':
  unittest.main()
