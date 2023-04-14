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
# Created by Jameed Hussain, October 2013

import re
import sys
from optparse import OptionParser

from rdkit import Chem, DataStructs

parser = OptionParser(
  description="Program to Tversky search results as part of Fraggle",
  epilog="Format of input file: whole mol smiles,ID,fraggle split smiles\t\t"
  "Output: query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim")
parser.add_option('-f', '--frags', dest='f_file', type='string',
                  help="File containing the query fragmentations from Fraggle")
parser.add_option(
  '-c', '--cutoff', dest='cutoff', type='float', default=0.8, help=
  "Cutoff for Tversy similarity. Only Tversky results with similarity greater than the cutoff will be output. DEFAULT = 0.8"
)

#parse the command line options
(options, args) = parser.parse_args()

if (options.f_file is None):
  print("Please specify the file containing the Fraggle fragmentations")
  sys.exit(1)

if ((options.cutoff < 0) or (options.cutoff > 1)):
  print("Please specify a Tversky cut-off between 0-1")
  sys.exit(1)

#sys.exit(1)
queries = []
query_info = []

q_split_input = open(options.f_file, "r")

for line in q_split_input:
  info = line.rstrip().split(",")

  qmol = Chem.MolFromSmiles(info[2])
  #print info[2]
  if qmol is None:
    sys.stderr.write("Can't generate mol for: %s\n" % (info[2]))
    continue

  #generate fp of query_substructs
  qfp = Chem.RDKFingerprint(qmol, maxPath=5, fpSize=1024, nBitsPerHash=2)

  queries.append(qfp)
  query_info.append((info[0], info[1], info[2]))

fragments = len(query_info)

for line in sys.stdin:

  line = line.rstrip()
  smi, id = re.split(r'\s|,', line)
  #print smi,id

  mol = Chem.MolFromSmiles(smi)
  if mol is None:
    sys.stderr.write("Can't generate mol for: %s\n" % (smi))
    continue

  mfp = Chem.RDKFingerprint(mol, maxPath=5, fpSize=1024, nBitsPerHash=2)
  #print smi

  res = DataStructs.BulkTverskySimilarity(mfp, queries, 0, 1, False)

  #query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim
  for i in range(fragments):
    if (res[i] >= options.cutoff):
      print("%s,%s,%s,%s,%s,%s" %
            (query_info[i][2], query_info[i][0], query_info[i][1], smi, id, res[i]))
