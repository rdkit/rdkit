#
#  Copyright (c) 2009, Novartis Institutes for BioMedical Research Inc.
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
# Created by Greg Landrum and Anna Vulpetti, March 2009

import pickle
import sys

from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

# sims is the list of similarity thresholds used to generate clusters
sims = [.9, .8, .7, .6]
smis = []
uniq = []
uFps = []

for fileN in sys.argv[1:]:
  inF = file(sys.argv[1], 'r')
  cols = pickle.load(inF)
  fps = pickle.load(inF)

  for row in fps:
    nm, smi, fp = row[:3]
    if smi not in smis:
      try:
        fpIdx = uFps.index(fp)
      except ValueError:
        fpIdx = len(uFps)
        uFps.append(fp)
      uniq.append([fp, nm, smi, 'FP_%d' % fpIdx] + row[3:])
      smis.append(smi)


def distFunc(a, b):
  return 1. - DataStructs.DiceSimilarity(a[0], b[0])


for sim in sims:
  clusters = Butina.ClusterData(uniq, len(uniq), 1. - sim, False, distFunc)
  print('Sim: %.2f, nClusters: %d' % (sim, len(clusters)), file=sys.stderr)
  for i, cluster in enumerate(clusters):
    for pt in cluster:
      uniq[pt].append(str(i + 1))
  cols.append('cluster_thresh_%d' % (int(100 * sim)))
print(' '.join(cols))
for row in uniq:
  print(' '.join(row[1:]))
