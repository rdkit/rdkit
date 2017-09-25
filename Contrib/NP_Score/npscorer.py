#
# calculation of natural product-likeness as described in:
#
# Natural Product-likeness Score and Its Application for Prioritization of Compound Libraries
# Peter Ertl, Silvio Roggo, and Ansgar Schuffenhauer
# Journal of Chemical Information and Modeling, 48, 68-74 (2008)
# http://pubs.acs.org/doi/abs/10.1021/ci700286x
#
# for the training of this model only openly available data have been used
# ~50,000 natural products collected from various open databases
# ~1 million drug-like molecules from ZINC as a "non-NP background"
#
# peter ertl, august 2015
#

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import sys, math, gzip, pickle
import os.path


def readNPModel(filename=os.path.join(os.path.dirname(__file__), 'publicnp.model.gz')):
  sys.stderr.write("reading NP model ...\n")
  fscore = pickle.load(gzip.open(filename))
  sys.stderr.write("model in\n")
  return fscore


def scoreMol(mol, fscore):
  if mol is None:
    raise ValueError('invalid molecule')
  fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
  bits = fp.GetNonzeroElements()

  # calculating the score
  score = 0.
  for bit in bits:
    score += fscore.get(bit, 0)
  score /= float(mol.GetNumAtoms())

  # preventing score explosion for exotic molecules
  if score > 4:
    score = 4. + math.log10(score - 4. + 1.)
  if score < -4:
    score = -4. - math.log10(-4. - score + 1.)
  return score


def processMols(fscore, suppl):
  sys.stderr.write("calculating ...\n")
  count = {}
  n = 0
  for i, m in enumerate(suppl):
    if m is None:
      continue

    n += 1
    score = "%.3f" % scoreMol(m, fscore)

    smiles = Chem.MolToSmiles(m, True)
    name = m.GetProp('_Name')
    print(smiles + "\t" + name + "\t" + score)

  sys.stderr.write("finished, " + str(n) + " molecules processed\n")


if __name__ == '__main__':

  fscore = readNPModel()  # fills fscore

  suppl = Chem.SmilesMolSupplier(sys.argv[1], smilesColumn=0, nameColumn=1, titleLine=False)
  processMols(fscore, suppl)

#
# Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
# All rights reserved.
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
