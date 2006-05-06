from SimDivFilters import rdSimDivPickers as rdsimdiv
from Numeric import *
from random import random

pkr = rdsimdiv.MaxMinPicker()

n = 1000
m = 80
dataPts = []
for i in range(n) :
    pt = zeros(2, 'd')
    pt[0] = 10.*random()
    pt[1] = 10.*random()
    print pt[0], pt[1]
    dataPts.append(pt)

# compute the distance matrix
distMat = zeros(n*(n-1)/2, 'd')
for i in range(n-1) :
    itab = n*i - ((i+1)*(i+2))/2
    pt1 = dataPts[i]
    for j in range(i+1, n) :
        id = itab + j
        pt2 = dataPts[j]
        diff = pt2 - pt1
        
        dist = sqrt(dot(diff, diff))
        #print pt1, pt2, diff, dist
        distMat[id] = dist
        


# now do the picking
res = pkr.Pick(distMat, n, m)

print "Results:"
for k in res :
    print dataPts[k][0], dataPts[k][1]


