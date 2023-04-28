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

import sys
from collections import defaultdict
from optparse import OptionParser

from rdkit import Chem, DataStructs
from rdkit.Chem.Fraggle import FraggleSim

#input format
#query_substructs,query_smiles,SMILES,ID,Tversky_sim

#algorithm
#read in query_substructs and smiles
#feed to atomcontrib function to return generalised_SMILES
#use Tanimoto to compare generalised_SMILES with query smiles to give fraggle similarity


parser = OptionParser(
  description="Program to post-process Tversky search results as part of Fraggle", epilog=
  "Format of input file: query_frag_smiles,query_smiles,query_id,retrieved_smi,retrieved_id,tversky_sim\t"
  "Output: SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,RDK5_Similarity")
parser.add_option(
  '-c', '--cutoff', action='store', dest='cutoff', type='float', default=0.7, help=
  "Cutoff for fraggle similarity. Only results with similarity greater than the cutoff will be output. DEFAULT = 0.7"
)
parser.add_option('-p', '--pfp', action='store', dest='pfp', type='float', default=0.8,
                  help="Cutoff for partial fp similarity. DEFAULT = 0.8")

if __name__ == '__main__':
  #parse the command line options
  (options, args) = parser.parse_args()

  if ((options.cutoff >= 0) and (options.cutoff <= 1)):
    fraggle_cutoff = options.cutoff
  else:
    print("Fraggle cutoff must be in range 0-1")
    sys.exit(1)

  print("SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,RDK5_Similarity")

  #create some data structure to store results
  id_to_smi = {}
  day_sim = {}
  frag_sim = {}
  query_size = {}
  query_mols = {}

  #generate dummy mol object which generates empty fp
  emptyMol = Chem.MolFromSmiles('*')

  #read the STDIN
  for line in sys.stdin:
    line = line.rstrip()
    qSubs, qSmi, qID, inSmi, id_, tversky = line.split(",")

    #add query to id_to_smi
    id_to_smi[qID] = qSmi
    id_to_smi[id_] = inSmi

    #add query to data structures
    frag_sim.setdefault(qID, defaultdict(float))
    day_sim.setdefault(qID, {})

    if (qID not in query_size):
      qMol = Chem.MolFromSmiles(qSmi)
      if qMol is None:
        sys.stderr.write("Can't generate mol for: %s\n" % (qSmi))
        continue
      query_mols[qID] = qMol
      query_size[qID] = qMol.GetNumAtoms()

    iMol = Chem.MolFromSmiles(inSmi)

    if iMol is None:
      sys.stderr.write("Can't generate mol for: %s\n" % (inSmi))
      continue

    #discard based on atom size
    if (iMol.GetNumAtoms() < query_size[qID] - 3):
      #sys.stderr.write("Too small: %s\n" % (inSmi) )
      continue

    if (iMol.GetNumAtoms() > query_size[qID] + 4):
      #sys.stderr.write("Too large: %s\n" % (inSmi) )
      continue

    #print '>>>',id_
    rdkit_sim, fraggle_sim = FraggleSim.compute_fraggle_similarity_for_subs(
      iMol, query_mols[qID], qSmi, qSubs, options.pfp)
    day_sim[qID][id_] = rdkit_sim
    frag_sim[qID][id_] = max(frag_sim[qID][id_], fraggle_sim)

    #check if you have the fp for the modified query
    #and generate if need to

  #right, print out the results for the query
  #Format: SMILES,ID,QuerySMI,QueryID,Fraggle_Similarity,Daylight_Similarity
  for qID in frag_sim:
    for id_ in frag_sim[qID]:
      if (frag_sim[qID][id_] >= fraggle_cutoff):
        print("%s,%s,%s,%s,%s,%s" %
              (id_to_smi[id_], id_, id_to_smi[qID], qID, frag_sim[qID][id_], day_sim[qID][id_]))
