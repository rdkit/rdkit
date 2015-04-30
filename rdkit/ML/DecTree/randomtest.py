import numpy
import random
from rdkit.ML.DecTree import ID3
from rdkit.ML.DecTree import CrossValidate
from rdkit.six.moves import xrange

def GenRandomExamples(nVars=10,randScale=0.3,bitProb=0.5,nExamples=500,seed=(0,0),
                      addResults=1):
    random.seed(seed[0])
    varWeights = numpy.array([random.random() for x in range(nVars)])*randScale
    examples = [None]*nExamples

    for i in xrange(nExamples):
        varVals=[random.random()>bitProb for x in range(nVars)]
        temp = numpy.array(varVals) * varWeights
        res = sum(temp)
        if addResults:
            varVals.append(res>=1.)
        examples[i] = varVals

    nPossibleVals = [2]*(nExamples+1)
    attrs = list(range(nVars))

    return (examples,attrs,nPossibleVals)
    
if __name__ == '__main__':
    from rdkit.six.moves import cPickle
    examples,attrs,nPossibleVals = GenRandomExamples()
    outF = open('random.dat.pkl','wb+')
    cPickle.dump(examples,outF)
    cPickle.dump(attrs,outF)
    cPickle.dump(nPossibleVals,outF)
    
    tree = ID3.ID3Boot(examples,attrs,nPossibleVals)
    tree.Pickle('save.pkl')



