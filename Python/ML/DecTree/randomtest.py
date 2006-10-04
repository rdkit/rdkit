## Automatically adapted for numpy.oldnumeric Sep 23, 2006 by alter_code1.py
import numpy
import numpy.oldnumeric.random_array as RandomArray

def GenRandomExamples(nVars=10,randScale=0.3,bitProb=0.5,nExamples=500,seed=(0,0),
                      addResults=1):
    apply(RandomArray.seed,seed)
    varWeights = RandomArray.random(nVars)*randScale
    examples = [None]*nExamples

    for i in xrange(nExamples):
        varVals = map(lambda x,y=bitProb:int(x>y),RandomArray.random(nVars))
        temp = numpy.array(varVals) * varWeights
        res = sum(temp)
        if addResults:
            varVals.append(int(res>=1.))
        examples[i] = varVals

    nPossibleVals = [2]*(nExamples+1)
    attrs = range(nVars)

    return (examples,attrs,nPossibleVals)


