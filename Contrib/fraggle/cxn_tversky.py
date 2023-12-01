#! /usr/bin/env jython

# Copyright (c) 2013, GlaxoSmithKline Research & Development Ltd.
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
#     * Neither the name of GlaxoSmithKline Research & Development Ltd.
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
# Created by Jameed Hussain, July 2013

import array
import re
import sys

from chemaxon.descriptors import (CFParameters, ChemicalFingerprint,
                                  SimilarityCalculatorFactory)
from chemaxon.struc import Molecule
from chemaxon.util import MolHandler


def desalt(mol):
  parmol = mol
  smi = mol.toFormat("smiles")
  parcount = 0
  msmi = smi.split('.')
  for smi in msmi:
    mol = MolHandler(smi).getMolecule()
    count = mol.getAtomCount()
    if count > parcount:
      parcount = count
      parmol = mol
  return parmol


cfp = CFParameters(
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?><ChemicalFingerprintConfiguration Version =\"0.3\" schemaLocation=\"cfp.xsd\">    <Parameters Length=\"1024\" BondCount=\"7\" BitCount=\"4\"/>    <StandardizerConfiguration Version =\"0.1\"><Actions><Action ID=\"aromatize\" Act=\"aromatize\"/> </Actions> </StandardizerConfiguration><ScreeningConfiguration><ParametrizedMetrics><ParametrizedMetric Name=\"Tversky\" ActiveFamily=\"Generic\" Metric=\"Tversky\" Threshold=\"0.5\" TverskyAlpha=\"0.1\" TverskyBeta=\"0.9\"/></ParametrizedMetrics></ScreeningConfiguration></ChemicalFingerprintConfiguration>"
)
cfp.setLength(1024)
cfp.setBondCount(7)
cfp.setBitCount(4)

#output needs to look like this:
#qSubs,qSmi,qID,inSmi,id,tversky

#first read in queries
q_split_input = open("frag_q_split_out", 'r')

queries = []
for line in q_split_input:
  info = line.rstrip().split(",")
  #print info

  #generate fp for query
  #print info[2]
  mol = MolHandler(info[2]).getMolecule()
  mol.aromatize(Molecule.AROM_GENERAL)

  qfp = ChemicalFingerprint(cfp)
  qfp.generate(mol)
  qintfp = array.array('i', list(map(int, qfp.toFloatArray())))

  queries.append((qintfp, info[0], info[1], info[2]))

#print queries

for line in sys.stdin:

  line_fields = re.split(r'\s|,', line)
  dbsmi = line_fields[0]
  dbid = line_fields[1]

  mol = MolHandler(dbsmi).getMolecule()
  mol_desalted = desalt(mol)
  mol_desalted.aromatize(Molecule.AROM_GENERAL)
  #print mol_desalted.toFormat("smiles")

  fp = ChemicalFingerprint(cfp)
  fp.generate(mol)
  intfp = array.array('i', list(map(int, fp.toFloatArray())))

  for q in queries:
    qsmi = q[1]
    qid = q[2]
    qsub = q[3]

    sc = SimilarityCalculatorFactory.create("Tversky;0.95;0.05")
    sc.setQueryFingerprint(q[0])

    tversky = sc.getSimilarity(intfp)

    if (tversky >= 0.9):
      print("%s,%s,%s,%s,%s,%s" % (qsub, qsmi, qid, dbsmi, dbid, tversky))
