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
# Created by Jameed Hussain, May 2013

from rdkit import Chem
from rdkit.Chem.Fraggle import FraggleSim

if __name__ == '__main__':
  import re
  import sys
  if (len(sys.argv) >= 2):
    print(
      "Program to run the first part of Fraggle. Program splits the molecule\nready for the search\n"
    )
    print("USAGE: ./fraggle.py <file_of_smiles")
    print("Format of smiles file: SMILES ID (space or comma separated)")
    print("Output: whole mol smiles,ID,fraggle split smiles\n")
    sys.exit(1)

  #read the STDIN
  for line in sys.stdin:
    line = line.rstrip()
    smi, id_ = re.split(r'\s|,', line)
    #print smi,id_

    mol = Chem.MolFromSmiles(smi)

    if mol is None:
      sys.stderr.write("Can't generate mol for: %s\n" % (smi))
      continue

    out_fragments = FraggleSim.generate_fraggle_fragmentation(mol)
    #print out the unique fragments
    for x in out_fragments:
      #cansmi
      temp = Chem.MolFromSmiles(x)

      print("%s,%s,%s" % (smi, id_, Chem.MolToSmiles(temp)))
