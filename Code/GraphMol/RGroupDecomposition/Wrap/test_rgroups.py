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

import unittest
import os, sys, copy

import pickle

from rdkit import rdBase
from rdkit import Chem
from rdkit.Chem.rdRGroupDecomposition import RGroupDecompose, RGroupDecomposition, RGroupDecompositionParameters
from collections import OrderedDict

# the RGD code can generate a lot of warnings. disable them
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.warning")


class TestCase(unittest.TestCase):

  def test_multicores(self):
    cores_smi_easy = OrderedDict()
    cores_smi_hard = OrderedDict()

    #cores_smi_easy['cephem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])CS2')
    cores_smi_easy['cephem'] = Chem.MolFromSmarts('O=C1C([*:1])C2N1C(C(O)=O)=C([*:3])CS2')
    cores_smi_hard['cephem'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])CS2')

    #cores_smi_easy['carbacephem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])CC2')
    cores_smi_easy['carbacephem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])CC2')
    cores_smi_hard['carbacephem'] = Chem.MolFromSmarts(
      'O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])CC2')

    #cores_smi_easy['oxacephem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])CO2')
    cores_smi_easy['oxacephem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])CO2')
    cores_smi_hard['oxacephem'] = Chem.MolFromSmarts(
      'O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])CO2')

    #cores_smi_easy['carbapenem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])C2')
    cores_smi_easy['carbapenem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])C2')
    cores_smi_hard['carbapenem'] = Chem.MolFromSmarts(
      'O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])C2')

    #cores_smi_easy['carbapenam'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])C2')
    cores_smi_easy['carbapenam'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)C([3*])([4*])C2')
    cores_smi_hard['carbapenam'] = Chem.MolFromSmarts(
      'O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])C2')

    #cores_smi_easy['penem'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)=C([3*])S2')
    cores_smi_easy['penem'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)=C([3*])S2')
    cores_smi_hard['penem'] = Chem.MolFromSmarts('O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)=C([3*])S2')

    #cores_smi_easy['penam'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])S2')
    cores_smi_easy['penam'] = Chem.MolFromSmarts('O=C1C([*:1])C2N1C(C(O)=O)C([*:3])([*:4])S2')
    cores_smi_hard['penam'] = Chem.MolFromSmarts(
      'O=C1C([*:2])([*:1])[C@@H]2N1C(C(O)=O)C([*:3])([*:4])S2')

    #cores_smi_easy['oxapenam'] = Chem.MolFromSmiles('O=C1C([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])O2')
    cores_smi_easy['oxapenam'] = Chem.MolFromSmarts('O=C1C([1*])C2N1C(C(O)=O)C([3*])([4*])O2')
    cores_smi_hard['oxapenam'] = Chem.MolFromSmarts(
      'O=C1C([2*])([1*])[C@@H]2N1C(C(O)=O)C([3*])([4*])O2')

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
    for k, v in columns.items():
      data[k] = [Chem.MolToSmiles(m, True) for m in v]

    rgroups2, unmatched = RGroupDecompose([core], mols)
    columns2, unmatched = RGroupDecompose([core], mols, asRows=False)
    data2 = {}
    for k, v in columns2.items():
      data2[k] = [Chem.MolToSmiles(m, True) for m in v]

    self.assertEqual(data, data2)
    columns3, unmatched = RGroupDecompose([core], mols, asRows=False, asSmiles=True)
    self.assertEqual(data, columns3)

  def test_h_options(self):
    core = Chem.MolFromSmiles("O=c1oc2ccccc2cc1")
    smiles = ("O=c1cc(Cn2ccnc2)c2ccc(Oc3ccccc3)cc2o1", "O=c1oc2ccccc2c(Cn2ccnc2)c1-c1ccccc1",
              "COc1ccc2c(Cn3cncn3)cc(=O)oc2c1")
    params = RGroupDecompositionParameters()
    rgd = RGroupDecomposition(core, params)
    for smi in smiles:
      m = Chem.MolFromSmiles(smi)
      rgd.Add(m)
    rgd.Process()
    columns = rgd.GetRGroupsAsColumns()
    self.assertEqual(columns['R2'][0].GetNumAtoms(), 7)

    params.removeHydrogensPostMatch = False
    rgd = RGroupDecomposition(core, params)
    for smi in smiles:
      m = Chem.MolFromSmiles(smi)
      rgd.Add(m)
    rgd.Process()
    columns = rgd.GetRGroupsAsColumns()
    self.assertEqual(columns['R2'][0].GetNumAtoms(), 12)

  def test_unmatched(self):
    cores = [Chem.MolFromSmiles("N")]
    mols = [
      Chem.MolFromSmiles("CC"),
      Chem.MolFromSmiles("CC"),
      Chem.MolFromSmiles("CC"),
      Chem.MolFromSmiles("N"),
      Chem.MolFromSmiles("CC")
    ]

    res, unmatched = RGroupDecompose(cores, mols)
    self.assertEqual(len(res), 1)
    self.assertEqual(unmatched, [0, 1, 2, 4])

  def test_userlabels(self):
    smis = ["C(Cl)N(N)O(O)"]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    smarts = 'C([*:1])N([*:5])O([*:6])'
    core = Chem.MolFromSmarts(smarts)
    rg = RGroupDecomposition(core)
    for m in mols:
      rg.Add(m)
    rg.Process()
    self.assertEqual(rg.GetRGroupsAsColumns(asSmiles=True), {
      'Core': ['C(N(O[*:6])[*:5])[*:1]'],
      'R1': ['Cl[*:1]'],
      'R5': ['N[*:5]'],
      'R6': ['O[*:6]']
    })

    smarts = 'C([*:4])N([*:5])O([*:6])'

    core = Chem.MolFromSmarts(smarts)
    rg = RGroupDecomposition(core)
    for m in mols:
      rg.Add(m)
    rg.Process()
    self.assertEqual(rg.GetRGroupsAsColumns(asSmiles=True), {
      'Core': ['C(N(O[*:6])[*:5])[*:4]'],
      'R4': ['Cl[*:4]'],
      'R5': ['N[*:5]'],
      'R6': ['O[*:6]']
    })

  def test_match_only_at_rgroups(self):
    smiles = ['c1ccccc1']  #, 'c1(Cl)ccccc1', 'c1(Cl)cc(Br)ccc1']
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]

    core1 = Chem.MolFromSmiles("c1([*:5])cc([*:6])ccc1")
    params = RGroupDecompositionParameters()
    params.onlyMatchAtRGroups = True
    rg = RGroupDecomposition(core1, params)
    for smi, m in zip(smiles, mols):
      self.assertTrue(rg.Add(m) != -1, smi)

  def test_incorrect_multiple_rlabels(self):
    mols = [Chem.MolFromSmiles(smi) for smi in (
      "C1NC(Cl)CC1",
      "C1OC(Cl)CC1",
      "C1(Cl)OCCC1",
    )]

    scaffolds = [Chem.MolFromSmarts(x) for x in ("C1NCCC1", )]
    groups, unmatched = RGroupDecompose(scaffolds, mols, asSmiles=True, asRows=False)

    self.assertEqual(groups, {'Core': ['C1CNC([*:1])C1'], 'R1': ['Cl[*:1].[H][*:1]']})

    scaffolds = [Chem.MolFromSmarts(x) for x in ("C1OCCC1", )]
    groups, unmatched = RGroupDecompose(scaffolds, mols, asSmiles=True, asRows=False)
    self.assertEqual(groups, {
      'Core': ['C1COC([*:1])C1', 'C1COC([*:1])C1'],
      'R1': ['Cl[*:1].[H][*:1]', 'Cl[*:1].[H][*:1]']
    })
    scaffolds = [Chem.MolFromSmarts(x) for x in ("C1NCCC1", "C1OCCC1")]
    groups, unmatched = RGroupDecompose(scaffolds, mols, asSmiles=True, asRows=False)
    self.assertEqual(
      groups, {
        'Core': ['C1CNC([*:1])C1', 'C1COC([*:1])C1', 'C1COC([*:1])C1'],
        'R1': ['Cl[*:1].[H][*:1]', 'Cl[*:1].[H][*:1]', 'Cl[*:1].[H][*:1]']
      })

  def test_aligned_cores(self):
    scaffolds = [Chem.MolFromSmarts(x) for x in ("C1NC1OC", "C1NC1NC")]
    mols = [Chem.MolFromSmiles(smi) for smi in ("C1NC1OCCC", "C1NC1NCCC")]
    groups, unmatched = RGroupDecompose(scaffolds, mols, asSmiles=True, asRows=True)
    #print("test_aligned_Cores")
    #print("groups:", groups)
    self.assertEqual(groups, [{
      'Core': 'C1NC1OC[*:1]',
      'R1': 'CC[*:1].[H][*:1].[H][*:1]'
    }, {
      'Core': 'C1NC1NC[*:2]',
      'R2': 'CC[*:2].[H][*:2].[H][*:2]'
    }])

  def test_aligned_cores2(self):
    scaffolds = [Chem.MolFromSmarts(x) for x in ("C1NCC1", "C1SCC1")]
    mols = [Chem.MolFromSmiles(smi) for smi in ("C1N(P)CC1", "C1S(P)CC1")]
    #print("test_aligned_Cores2")
    groups, unmatched = RGroupDecompose(scaffolds, mols, asSmiles=True, asRows=True)
    #print("groups: ", groups)
    self.assertEqual(groups, [{
      'Core': 'C1CN([*:1])C1',
      'R1': 'P[*:1]'
    }, {
      'Core': 'C1C[SH]([*:2])C1',
      'R2': 'P[*:2].[H][*:2]'
    }])

  def test_getrgrouplabels(self):
    smis = ["C(Cl)N(N)O(O)"]
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    smarts = 'C([*:1])N([*:5])O([*:6])'
    core = Chem.MolFromSmarts(smarts)
    rg = RGroupDecomposition(core)
    for m in mols:
      rg.Add(m)
    rg.Process()
    self.assertEqual(set(rg.GetRGroupLabels()), set(rg.GetRGroupsAsColumns()))

  def test_timeout(self):
    smis = '''CN(C)Cc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CNc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
Cc1cc2cc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NC4CCN(C)CC4)c([N+](=O)[O-])c3)ccc2[nH]1
Cc1cc2cc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)ccc2[nH]1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccccc1Cl
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(Cl)c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(Cl)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc([N+](=O)[O-])c1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(CO)c1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2ccccc2Cl)cc1[N+](=O)[O-]
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(Cl)c1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(Cl)cc1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc(Cl)c2)cc1[N+](=O)[O-]
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2ccc(Cl)cc2)cc1[N+](=O)[O-]
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc3c2ccn3C)cc1[N+](=O)[O-]
CC(=O)Nc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
Nc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
Nc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
COc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CN(C)c1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
N#Cc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
Cc1nc2ccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCC4CCOCC4)c([N+](=O)[O-])c3)cc2s1
Cc1nc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)ccc2s1
Cc1nc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN(C)C)c([N+](=O)[O-])c3)ccc2s1
CN(C)C(=O)CCc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)C(=O)Cc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)CCCc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)CCc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)C(=O)c1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)Cc1ccccc1Oc1cc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)ccc1C(=O)NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc(N3CCOCC3)c2)cc1[N+](=O)[O-]
Cc1nc(C)c(-c2cccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)c2)s1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cc(Cl)cc(Cl)c3)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(Cl)c3)cc2[N+](=O)[O-])CC1
CN(C)CCOc1ccc(-c2ccc(Cl)cc2)c(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NC5CCN(C)CC5)c([N+](=O)[O-])c4)c(Oc4cccc(Cl)c4)c3)CC2)c1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)OC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN(C)CCOc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccc(N)c(Cl)c3)cc2[N+](=O)[O-])CC1
CC(C)N1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Br)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CCCC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
Cc1n[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)c12
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(F)c3F)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(Br)c3)cc2[N+](=O)[O-])CC1
CCN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN1C(C)(C)CC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1(C)C
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NC5CCN(C6CCOCC6)CC5)c([N+](=O)[O-])c4)c(Oc4cccc(F)c4F)c3)CC2)=C(c2ccc(Cl)cc2)C1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cc(F)c4[nH]ccc4c3)cc2[N+](=O)[O-])CC1
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NCCCN5CCOCC5)c([N+](=O)[O-])c4)c(Oc4cccc(F)c4F)c3)CC2)=C(c2ccc(Cl)cc2)C1
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NCCCN5CCOCC5)c([N+](=O)[O-])c4)c(Oc4ccc(N)c(Cl)c4)c3)CC2)=C(c2ccc(Cl)cc2)C1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3cc(CCN4CCCC4)ccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(Cl)c1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(Cl)c3Cl)cc2[N+](=O)[O-])CC1
Cc1n[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NC4CCN(C)CC4)c([N+](=O)[O-])c3)c12
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CCCCC5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(C(F)(F)F)c3)cc2[N+](=O)[O-])CC1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc3c2CCC(=O)N3)cc1[N+](=O)[O-]
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccccc3Cl)cc2S(=O)(=O)C(F)(F)F)CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cc(Cl)ccc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3ccc(F)cc3Cl)cc2[N+](=O)[O-])CC1
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)C5)CC4)cc3Oc3ccccc3Cl)cc2[N+](=O)[O-])CC1
Cc1c[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)c12
CN1CCC(Nc2ccc(S(=O)(=O)NC(=O)c3ccc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)cc3Oc3cccc(C(F)(F)F)c3Cl)cc2[N+](=O)[O-])CC1
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NC5CCN(C6CC6)CC5)c([N+](=O)[O-])c4)c(Oc4ccccc4Cl)c3)CC2)=C(c2ccc(Cl)cc2)C1
Cc1c[nH]c2cccc(Oc3cc(N4CCN(CC5=C(c6ccc(Cl)cc6)CC(C)(C)CC5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NC4CCN(C)CC4)c([N+](=O)[O-])c3)c12
CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NCCCN5CCOCC5)c([N+](=O)[O-])c4)c(Oc4cc(Cl)ccc4Cl)c3)CC2)=C(c2ccc(Cl)cc2)C1
Cn1ccc2c(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)cccc21
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(N2CCOCC2)c1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2ccc3[nH]cc(CCC(=O)N4CCOCC4)c3c2)cc1[N+](=O)[O-]
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(OCc2ccccc2)c1
N#Cc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc2[nH]cc(CCC(=O)N3CCOCC3)c2c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc2[nH]cc(CCCN3CCOCC3)c2c1
CN(C)Cc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(-n2ccnc2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)cc1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc([N+](=O)[O-])c1
CCN(Cc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1)C(=O)OC(C)(C)C
CCN(Cc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1)C(=O)OC(C)(C)C
CCNCc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
CCNCc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
CC(=O)Nc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
CC(C)(C)OC(=O)Nc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccccc1-c1ccccc1
CC(C)(C)OC(=O)Nc1cccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(-c2ccccc2)c1
CN(C)CCc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(OCc2ccccc2)cc1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(N2CCOCC2)c1
Cc1nc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCC4CCOCC4)c([N+](=O)[O-])c3)ccc2s1
CC(C)(C)OC(=O)N1CCN(c2cccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCC4CCOCC4)c([N+](=O)[O-])c3)c2)CC1
CN(C)CCCNc1ccc(S(=O)(=O)NC(=O)c2ccc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)cc2Oc2cccc(OCc3ccccc3)c2)cc1[N+](=O)[O-]
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(OCc2ccccc2)c1
O=C(NS(=O)(=O)c1ccc(NCC2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(OCCN2CCOCC2)cc1
O=C1CCc2c(cccc2Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)N1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(OCc2ccccc2)cc1
CC(C)(C)OC(=O)N1CCN(c2ccc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)cc2)CC1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1cccc(-c2ccncc2)c1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(-c2ccncc2)cc1
O=C(NS(=O)(=O)c1ccc(NCCCN2CCOCC2)c([N+](=O)[O-])c1)c1ccc(N2CCN(Cc3ccccc3-c3ccc(Cl)cc3)CC2)cc1Oc1ccc(-c2cccnc2)cc1
CN(C)C(=O)COc1ccc(Oc2cc(N3CCN(Cc4ccccc4-c4ccc(Cl)cc4)CC3)ccc2C(=O)NS(=O)(=O)c2ccc(NCC3CCOCC3)c([N+](=O)[O-])c2)cc1
Cn1cnc2cc(Oc3cc(N4CCN(Cc5ccccc5-c5ccc(Cl)cc5)CC4)ccc3C(=O)NS(=O)(=O)c3ccc(NCCCN4CCOCC4)c([N+](=O)[O-])c3)ccc21'''
    mols = [Chem.MolFromSmiles(x) for x in smis.split('\n')]

    core = Chem.MolFromSmarts('O=C(NS(=O)(=O)c1ccccc1)c1ccccc1Oc1ccccc1')
    ps = RGroupDecompositionParameters()
    ps.timeout = 1.0
    with self.assertRaises(RuntimeError):
      rg = RGroupDecomposition(core, ps)
      for m in mols:
        rg.Add(m)

    with self.assertRaises(RuntimeError):
      columns2, unmatched = RGroupDecompose([core], mols, asRows=False, options=ps)


if __name__ == '__main__':
  unittest.main()
