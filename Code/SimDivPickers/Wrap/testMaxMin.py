from __future__ import print_function
from rdkit.SimDivFilters import rdSimDivPickers as rdsimdiv
import numpy
from rdkit import RDRandom
RDRandom.seed(23)


pkr = rdsimdiv.MaxMinPicker()

n = 1000
m = 80
dataPts = []
for i in range(n) :
    pt = numpy.zeros(2, 'd')
    pt[0] = 10.*RDRandom.random()
    pt[1] = 10.*RDRandom.random()
    dataPts.append(pt)

# compute the distance matrix
distMat = numpy.zeros(n*(n-1)/2, 'd')
for i in range(n-1) :
    itab = n*i - ((i+1)*(i+2))/2
    pt1 = dataPts[i]
    for j in range(i+1, n) :
        id = itab + j
        pt2 = dataPts[j]
        diff = pt2 - pt1
        
        dist = numpy.sqrt(numpy.dot(diff, diff))
        distMat[id] = dist
        


# now do the picking
res = pkr.Pick(distMat, n, m)

print("Results:")
for k in res :
    print(dataPts[k][0], dataPts[k][1])


