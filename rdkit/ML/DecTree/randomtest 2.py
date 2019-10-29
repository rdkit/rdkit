import random

import numpy

from rdkit.ML.DecTree import ID3


def GenRandomExamples(nVars=10, randScale=0.3, bitProb=0.5, nExamples=500, seed=(0, 0),
                      addResults=1):
  random.seed(seed[0])
  varWeights = numpy.array([random.random() for _ in range(nVars)]) * randScale
  examples = [None] * nExamples

  for i in range(nExamples):
    varVals = [random.random() > bitProb for _ in range(nVars)]
    temp = numpy.array(varVals) * varWeights
    res = sum(temp)
    if addResults:
      varVals.append(res >= 1.)
    examples[i] = varVals

  nPossibleVals = [2] * (nExamples + 1)
  attrs = list(range(nVars))

  return (examples, attrs, nPossibleVals)


if __name__ == '__main__':  # pragma: nocover
  import pickle
  examples, attrs, nPossibleVals = GenRandomExamples()
  outF = open('random.dat.pkl', 'wb+')
  pickle.dump(examples, outF)
  pickle.dump(attrs, outF)
  pickle.dump(nPossibleVals, outF)

  tree = ID3.ID3Boot(examples, attrs, nPossibleVals)
  tree.Pickle('save.pkl')
