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
from CreateFps import GetMolFingerprint

from rdkit import Chem, DataStructs
from rdkit.RDLogger import logger

logger = logger()
import sys

# maxPathLength is the maximum path length in atoms
# maxPathLength=6 corresponds to F-FP-5
# maxPathLength=7 corresponds to F-FP-6
# maxPathLength=8 corresponds to F-FP-7
maxPathLength = 8

# nameField is the name of the property (from the SD file) that has molecule
# names... If the molecules have names in the first row of the file, use "_Name"
nameField = 'Compound_orig'
#nameField = '_Name'

# propField is the name of the property (from the SD file) you want to use
# as the "activity"
propField = 'chemical_shift_1'

# similarity threshold for a pair to be considered interesting.
# (i.e. pairs with a similarity below this value will not be
# added to the output.
similarityThreshold = 0.5

if __name__ == '__main__':
  suppl = Chem.SDMolSupplier(sys.argv[1])
  outF = file(sys.argv[2], 'w+')

  data = []
  logger.info('reading molecules and generating fingeprints')
  for i, mol in enumerate(suppl):
    if not mol:
      continue
    smi = Chem.MolToSmiles(mol, True)
    nm = mol.GetProp(nameField)
    property = float(mol.GetProp(propField))
    fp = GetMolFingerprint(mol, maxPathLength)
    data.append((nm, smi, property, fp))

  logger.info('  got %d molecules' % len(data))

  logger.info('calculating pairs')
  pairs = []
  for i in range(len(data)):
    for j in range(i + 1, len(data)):
      if DataStructs.DiceSimilarity(data[i][-1], data[j][-1]) > similarityThreshold:
        pairs.append((i, j))
    if not (i + 1) % 100:
      logger.info('Done %d molecules' % (i + 1))

  logger.info('  got %d reasonable pairs' % len(pairs))

  logger.info('creating output file')
  print >> outF, 'nameA|nameB|nameAB|smilesA|smilesB|smilesAB|actA|actB|dAct|dist|disparity'
  for i, j in pairs:
    if data[i][2] < data[j][2]:
      i, j = j, i
    nmi, smii, propi, fpi = data[i]
    nmj, smij, propj, fpj = data[j]
    dAct = propi - propj
    dist = 1. - DataStructs.DiceSimilarity(fpi, fpj)
    if dist != 0:
      disparity = dAct / dist
    else:
      disparity = 1000
    print >> outF, '%s|%s|%s_%s|%s|%s|%s.%s|%f|%f|%f|%f|%f' % (
      nmi, nmj, nmi, nmj, smii, smij, smii, smij, propi, propj, dAct, dist, disparity)
