#
#  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
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
# Created by Nadine Schneider, July 2016

import copy
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import AllChem


def transferAgentsToReactants(rxn):
  for a in range(rxn.GetNumAgentTemplates()):
    agent = rxn.GetAgentTemplate(a)
    rxn.AddReactantTemplate(agent)


def removeAgentsAndTransferToReactants(rxn):
  tmp = []
  rxn.RemoveAgentTemplates(tmp)
  for a in tmp:
    rxn.AddReactantTemplate(a)


def getNumPositiveCounts(fp):
  count = 0
  for k, v in fp.GetNonzeroElements().items():
    if v > 0:
      count += v
  return count


def getNumNegativeCounts(fp):
  count = 0
  for k, v in fp.GetNonzeroElements().items():
    if v < 0:
      count += abs(v)
  return count


def getNumPositiveBitCountsOfRadius0(fp, bitinfo):
  count = 0
  bitsUnmappedAtoms = []
  for k in bitinfo:
    if bitinfo[k][0][1] == 0:
      v = fp[k]
      if v > 0:
        count += 1
        bitsUnmappedAtoms.append((k, v))
  return count, bitsUnmappedAtoms


def getSumFps(fps):
  summedFP = copy.deepcopy(fps[0])
  for fp in fps[1:]:
    summedFP += fp
  return summedFP


def uniqueMolecules(mols):
  smiles = [Chem.MolToSmiles(mol) for mol in mols]
  uniqueMolecules = defaultdict(int)
  for n, smi in enumerate(smiles):
    uniqueMolecules[n] = smiles.index(smi)
  return uniqueMolecules, smiles
