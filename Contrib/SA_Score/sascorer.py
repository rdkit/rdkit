#
# calculation of synthetic accessibility score as described in:
#
# Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
# Peter Ertl and Ansgar Schuffenhauer
# Journal of Cheminformatics 1:8 (2009)
# http://www.jcheminf.com/content/1/1/8
#
# several small modifications to the original paper are included
# particularly slightly different formula for marocyclic penalty
# and taking into account also molecule symmetry (fingerprint density)
#
# for a set of 10k diverse molecules the agreement between the original method
# as implemented in PipelinePilot and this implementation is r2 = 0.97
#
# peter ertl & greg landrum, september 2013
#
# The scoring itself now lives in C++ as RDKit::Descriptors::calcSAScore. This
# module remains the supported entry point; the C++ code is reached through a
# private binding and is not exposed as an RDKit descriptor.
#

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

import gzip
import pickle

import os.path as op

from rdkit import RDConfig

_fragmentTable = None

_DEFAULT_SCORES = "fpscores.pkl.gz"


def readFragmentScores(name=_DEFAULT_SCORES):
  """ reads the fragment contributions used by calculateScore()

    The default reads the binary table shipped with RDKit. A .bin path is read
    directly; anything else is treated as a pickle from build_sascore_db.py,
    which is slower to load than converting it once with
    $RDBASE/Data/SA_Score/build_fpscores_bin.py.

  """
  global _fragmentTable

  if name == _DEFAULT_SCORES:
    _fragmentTable = rdMolDescriptors._SAScoreFragmentTable(
      op.join(RDConfig.RDDataDir, 'SA_Score', 'fpscores.bin'))
    return

  if str(name).endswith('.bin'):
    _fragmentTable = rdMolDescriptors._SAScoreFragmentTable(str(name))
    return

  bitIds = []
  contributions = []
  for row in pickle.load(gzip.open(name)):
    bitIds.extend(row[1:])
    contributions.extend([float(row[0])] * (len(row) - 1))
  _fragmentTable = rdMolDescriptors._SAScoreFragmentTable(bitIds, contributions)


def numBridgeheadsAndSpiro(mol, ri=None):
  nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
  nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
  return nBridgehead, nSpiro


def calculateScore(m):
  """ returns a score between 1 (easy to make) and 10 (hard to make), or None
      for a molecule with no atoms

  """
  if not m.GetNumAtoms():
    return None

  if _fragmentTable is None:
    readFragmentScores()

  return rdMolDescriptors._CalcSAScore(m, _fragmentTable)


def processMols(mols):
  print('smiles\tName\tsa_score')
  for i, m in enumerate(mols):
    if m is None:
      continue

    s = calculateScore(m)

    smiles = Chem.MolToSmiles(m)
    if s is None:
      print(f"{smiles}\t{m.GetProp('_Name')}\t{s}")
    else:
      print(f"{smiles}\t{m.GetProp('_Name')}\t{s:3f}")


if __name__ == '__main__':
  import sys
  import time

  t1 = time.time()
  if len(sys.argv) == 2:
    readFragmentScores()
  else:
    readFragmentScores(sys.argv[2])
  t2 = time.time()

  molFile = sys.argv[1]
  if molFile.endswith("smi"):
    suppl = Chem.SmilesMolSupplier(molFile)
  elif molFile.endswith("sdf"):
    suppl = Chem.SDMolSupplier(molFile)
  else:
    print(f"Unrecognized file extension for {molFile}")
    sys.exit()

  t3 = time.time()
  processMols(suppl)
  t4 = time.time()

  print('Reading took %.2f seconds. Calculating took %.2f seconds' % ((t2 - t1), (t4 - t3)),
        file=sys.stderr)

#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
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
