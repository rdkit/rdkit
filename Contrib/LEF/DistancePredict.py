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
from rdkit.ML.KNN.KNNRegressionModel import KNNRegressionModel
from rdkit.RDLogger import logger

logger = logger()
import sys

# nameField is the name of the property (from the SD file) that has molecule
# names...If the molecules have names in the first row of the file, use "_Name"
nameField = 'Compound_orig'
#nameField = '_Name'

# propField is the name of the property (from the SD file) you want to generate
# predictions for
propField = 'chemical_shift_1'

weightedAverage = True
import copy
import types
from optparse import Option, OptionParser, OptionValueError


def check_floatlist(option, opt, value):
  try:
    v = eval(value)
    if type(v) not in (types.ListType, types.TupleType):
      raise ValueError
    v = [float(x) for x in v]
  except ValueError:
    raise OptionValueError("option %s : invalid float list value: %r" % (opt, value))
  return v


class MyOption(Option):
  TYPES = Option.TYPES + ("floatlist", )
  TYPE_CHECKER = copy.copy(Option.TYPE_CHECKER)
  TYPE_CHECKER["floatlist"] = check_floatlist


parser = OptionParser("distance predict", version='%prog', option_class=MyOption)
parser.add_option('--maxPathLength', '--max', default=8, type=int,
                  help='maximum length path for the fingerprint')
parser.add_option('--similarityThreshold', '--sim', default=[0.9], type='floatlist',
                  help='threshold for similarity')
parser.add_option('--numNeighbors', '--num', '-n', '-k', default=50, type=int,
                  help='number of neighbors to consider')
parser.add_option('--neighborsFile', '--nbrs', default='',
                  help='name of an output file to hold the neighbor lists')
parser.add_option('--scan', default=False, action="store_true")

if __name__ == '__main__':
  options, args = parser.parse_args()
  outF = file(args[-1], 'w+')

  logger.info('reading training molecules and generating fingerprints')
  suppl = Chem.SDMolSupplier(args[0])
  train = []
  for i, mol in enumerate(suppl):
    if not mol:
      continue
    smi = Chem.MolToSmiles(mol, True)
    nm = mol.GetProp(nameField)
    property = float(mol.GetProp(propField))
    fp = GetMolFingerprint(mol, options.maxPathLength)
    train.append((nm, smi, fp, property))
  logger.info('  got %d molecules' % len(train))

  if len(args) > 2:
    suppl = Chem.SDMolSupplier(args[1])
    haveTest = True
    logger.info('reading testing molecules and generating fingerprints')
    test = []
    for i, mol in enumerate(suppl):
      if not mol:
        continue
      smi = Chem.MolToSmiles(mol, True)
      nm = mol.GetProp(nameField)
      if mol.HasProp(propField):
        property = float(mol.GetProp(propField))
      else:
        property = 0
      fp = GetMolFingerprint(mol, options.maxPathLength)
      test.append((nm, smi, fp, property))
    logger.info('  got %d molecules' % len(test))
  else:
    haveTest = False
    test = train

  results = [None] * len(test)
  for i in range(len(test)):
    results[i] = [None] * len(options.similarityThreshold)
  if options.neighborsFile:
    nbrFile = file(options.neighborsFile, 'w+')
    print >> nbrFile, 'ID|CompoundName|CompoundSmiles|NeighborName|NeighborSmiles|NeighborShift|Similarity'
    id = 1
  else:
    nbrFile = None
  for j, thresh in enumerate(options.similarityThreshold):
    if not haveTest:
      logger.info('Doing cross validation with threshold %.2f' % thresh)
    else:
      logger.info('Doing prediction with threshold %.2f' % thresh)
    for i in range(len(test)):
      if not haveTest:
        localTrain = [train[x] for x in range(len(train)) if x != i]
      else:
        localTrain = train
      localTest = test[i]
      mdl = KNNRegressionModel(options.numNeighbors, [],
                               lambda x, y, *args: 1 - DataStructs.DiceSimilarity(x[-2], y[-2]),
                               radius=1. - thresh)
      mdl.SetTrainingExamples(localTrain)
      nbrs = []
      pred = mdl.PredictExample(localTest, weightedAverage=weightedAverage, neighborList=nbrs)
      nm, smi, fp, prop = test[i]

      if nbrFile:
        for dist, data in nbrs:
          if data is None:
            continue
          nnm, nsmi, nfp, nproperty = data
          outRow = [str(id), nm, smi, nnm, nsmi, str(nproperty), str(dist - 1.)]
          id += 1
          print >> nbrFile, '|'.join(outRow)
      nbrs = [x for x in nbrs if x[1] is not None]
      results[i][j] = (nm, smi, prop, pred, len(nbrs))
      if not (i + 1) % 100:
        logger.info('Done %d molecules' % (i + 1))
    logger.info('  done')
  numNeighbors = options.numNeighbors
  maxPathLength = options.maxPathLength - 1
  logger.info('creating output file')
  headers = ['name', 'smiles', 'shift']
  for thresh in options.similarityThreshold:
    headers.append('predShift_%(maxPathLength)d_%(numNeighbors)d_%(thresh).2f' % locals())
    headers.append('dPred_%(maxPathLength)d_%(numNeighbors)d_%(thresh).2f' % locals())
    headers.append('nbrs_%(maxPathLength)d_%(numNeighbors)d_%(thresh).2f' % locals())
  print >> outF, '|'.join(headers)
  for i in range(len(test)):
    nm = results[i][0][0]
    smi = results[i][0][1]
    prop = results[i][0][2]
    row = [nm, smi, str(prop)]
    for j in range(len(options.similarityThreshold)):
      nbrs = results[i][j][4]
      pred = results[i][j][3]
      row.append(str(pred))
      row.append(str(abs(prop - pred)))
      row.append(str(nbrs))
    print >> outF, '|'.join(row)
