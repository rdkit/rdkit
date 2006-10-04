## Automatically adapted for numpy.oldnumeric Sep 23, 2006 by alter_code1.py

from Numeric import *
import numpy.oldnumeric.random_array as RandomArray
from ML.DecTree import ID3,CrossValidate,Forest
from ML.DecTree import randomtest
    
if __name__ == '__main__':
    examples,attrs,nPossibleVals = randomtest.GenRandomExamples(nExamples=100)
    forest = Forest.Forest()
    forest.Grow(examples,attrs,nPossibleVals,nTries=10)
    print str(forest)
    forest.Pickle('regress/RandomForest.pkl')
    



